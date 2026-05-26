# Bootstrap and noisy-data stability comparison for seven network methods.
#
# This workflow fits a fixed single-model specification for PhyMapNet. It is
# separate from the reliability-consensus analysis in 2-sensitivity_analysis.
#
# The script defaults to a short smoke test. Run the full publication workflow
# only with the explicit environment settings documented in the root README.

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages(library(dplyr))

base_dir <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
source(file.path(base_dir, "3-bootstrap", "bootstrap_phymapnet_helpers.R"))

parse_values <- function(value) {
  values <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  values[nzchar(values)]
}

is_true <- function(value) {
  tolower(value) %in% c("true", "t", "yes", "y", "1")
}

smoke_test <- is_true(Sys.getenv("PHYMAPNET_BOOTSTRAP_SMOKE_TEST", unset = "true"))

all_datasets <- c("smoking", "caffeine", "hmp_stool")
dataset_default <- if (smoke_test) "caffeine" else paste(all_datasets, collapse = ",")
datasets <- parse_values(Sys.getenv("PHYMAPNET_BOOTSTRAP_DATASETS", unset = dataset_default))
invalid_datasets <- setdiff(datasets, all_datasets)
if (length(invalid_datasets)) {
  stop("Unknown datasets: ", paste(invalid_datasets, collapse = ", "), call. = FALSE)
}

max_reps_default <- if (smoke_test) "1" else "100"
max_reps <- as.integer(Sys.getenv("PHYMAPNET_BOOTSTRAP_MAX_REPS", unset = max_reps_default))
if (is.na(max_reps) || max_reps < 1L) {
  stop("PHYMAPNET_BOOTSTRAP_MAX_REPS must be a positive integer.", call. = FALSE)
}

all_method_keys <- c(
  "phymap_phylo", "phymap_uniform", "spieceasi_mb", "spieceasi_glasso",
  "sparcc", "cclasso", "c_mi"
)
method_default <- if (smoke_test) "phymap_phylo,phymap_uniform" else {
  paste(all_method_keys, collapse = ",")
}
method_keys <- parse_values(Sys.getenv("PHYMAPNET_BOOTSTRAP_METHODS", unset = method_default))
invalid_methods <- setdiff(method_keys, all_method_keys)
if (length(invalid_methods)) {
  stop("Unknown method keys: ", paste(invalid_methods, collapse = ", "), call. = FALSE)
}

if (any(method_keys %in% c("spieceasi_mb", "spieceasi_glasso"))) {
  suppressPackageStartupMessages(library(SpiecEasi))
}
if (any(method_keys %in% c("sparcc", "cclasso", "c_mi"))) {
  suppressPackageStartupMessages(library(CMiNet))
}

dir_filter <- file.path(base_dir, "data", "filter")
dir_boot <- file.path(base_dir, "data", "generate", "bootstrap")
dir_noisy <- file.path(base_dir, "data", "generate", "noisy")
output_default <- if (smoke_test) {
  file.path(tempdir(), "phymapnet_bootstrap_smoke")
} else {
  file.path(base_dir, "result", "bootstrap")
}
dir_out <- Sys.getenv("PHYMAPNET_BOOTSTRAP_OUTPUT_DIR", unset = output_default)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# Fixed settings used for the reported stability comparison.
PHYMAP_PARAMS <- list(
  alpha = 0.1,
  kernel = "laplacian",
  nue_k = 10,
  epsilon1 = 1,
  epsilon2 = 0.4,
  th = 0.99
)

SE_MB_PARAMS <- list(
  method = "mb",
  lambda.min.ratio = 1e-2,
  nlambda = 15,
  pulsar.params = list(rep.num = 20, ncores = 4)
)

SE_GL_PARAMS <- list(
  method = "glasso",
  lambda.min.ratio = 1e-2,
  nlambda = 15,
  pulsar.params = list(rep.num = 50, ncores = 4)
)

SPARCC_PARAMS <- list(imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4)
CCLASSO_PARAMS <- list(
  counts = FALSE, pseudo = 0.5, k_cv = 3, lam_int = c(1e-4, 1),
  k_max = 20, n_boot = 20
)
CMI_PARAMS <- list(quantitative = TRUE, q1 = 0.9, q2 = 0.95)

count_edges <- function(B) {
  sum(B[upper.tri(B)])
}

clean_B <- function(B) {
  B <- as.matrix(B)
  diag(B) <- 0L
  B[lower.tri(B)] <- 0L
  (B != 0) * 1L
}

get_f1 <- function(ev) {
  metrics <- ev[[2]]
  hit <- which(tolower(names(metrics)) %in% c("f_score", "f.score", "f1", "fscore"))
  if (length(hit)) {
    return(as.numeric(metrics[[hit[1]]]))
  }
  as.numeric(metrics[[ncol(metrics)]])
}

make_row <- function(method, dataset, split, replicate, n_ref, n_est, f1) {
  data.frame(
    method = method, dataset = dataset, split = split, replicate = replicate,
    n_edges_ref = n_ref, n_edges_est = n_est, F1 = f1,
    stringsAsFactors = FALSE
  )
}

read_filtered_inputs <- function(dataset) {
  otu_file <- file.path(dir_filter, paste0(dataset, "_otu_filtered.csv"))
  dist_file <- file.path(dir_filter, paste0(dataset, "_dist_matrix.csv"))
  if (!file.exists(otu_file) || !file.exists(dist_file)) {
    stop("Missing filtered inputs for dataset: ", dataset, call. = FALSE)
  }
  otu <- as.matrix(read.csv(otu_file, row.names = 1, check.names = FALSE))
  dist <- as.matrix(read.csv(dist_file, row.names = 1, check.names = FALSE))
  if (!identical(sort(colnames(otu)), sort(rownames(dist))) ||
      !identical(sort(colnames(otu)), sort(colnames(dist)))) {
    stop("OTU and distance-matrix taxa do not match for dataset: ", dataset, call. = FALSE)
  }
  if (!isTRUE(all.equal(dist, t(dist), tolerance = 1e-10, check.attributes = FALSE))) {
    stop("Distance matrix is not symmetric for dataset: ", dataset, call. = FALSE)
  }
  taxa <- colnames(otu)
  list(otu = otu, dist = dist[taxa, taxa, drop = FALSE], taxa = taxa)
}

load_all_replicates <- function(dataset) {
  flatten_rds <- function(files) {
    if (!length(files)) {
      return(NULL)
    }
    out <- list()
    for (file in files) {
      object <- readRDS(file)
      if (is.list(object)) {
        out <- c(out, object)
      } else {
        out[[length(out) + 1L]] <- object
      }
    }
    out
  }
  list(
    boot = flatten_rds(Sys.glob(file.path(dir_boot, paste0(dataset, "_bootstrap_*.rds")))),
    noisy = flatten_rds(Sys.glob(file.path(dir_noisy, paste0(dataset, "_noisy_*.rds"))))
  )
}

select_subset <- function(reps, n = max_reps, seed = 42) {
  set.seed(seed)
  select_one <- function(values) {
    if (is.null(values) || !length(values)) {
      return(list(values = NULL, indices = integer()))
    }
    indices <- sort(sample(seq_along(values), min(n, length(values)), replace = FALSE))
    list(values = values[indices], indices = indices)
  }
  boot <- select_one(reps$boot)
  noisy <- select_one(reps$noisy)
  list(
    boot = boot$values,
    noisy = noisy$values,
    boot_indices = boot$indices,
    noisy_indices = noisy$indices
  )
}

extract_cminet_B <- function(res, TT = 0.95) {
  mat <- NULL
  for (slot in c("score", "weighted_adj", "adj", "network", "W", "result")) {
    if (is.matrix(res[[slot]]) && nrow(res[[slot]]) > 1) {
      mat <- res[[slot]]
      break
    }
  }
  if (is.null(mat)) {
    for (name in names(res)) {
      if (is.matrix(res[[name]]) && nrow(res[[name]]) > 1) {
        mat <- res[[name]]
        break
      }
    }
  }
  if (is.null(mat)) {
    stop("Cannot extract a matrix from CMiNet result.", call. = FALSE)
  }
  upper_positive <- mat[upper.tri(mat)]
  upper_positive <- upper_positive[upper_positive > 0]
  if (!length(upper_positive)) {
    return(clean_B(mat * 0))
  }
  threshold <- quantile(upper_positive, TT, na.rm = TRUE)
  clean_B((mat >= threshold) * 1L)
}

run_cminet_single <- function(X, method_name, params) {
  methods <- list(
    pearson = list(enabled = FALSE),
    spearman = list(enabled = FALSE),
    bicor = list(enabled = FALSE),
    sparcc = list(enabled = FALSE),
    spiecEasi_mb = list(enabled = FALSE),
    spiecEasi_glasso = list(enabled = FALSE),
    spring = list(enabled = FALSE),
    gcoda = list(enabled = FALSE),
    c_MI = list(enabled = FALSE),
    cclasso = list(enabled = FALSE)
  )
  methods[[method_name]] <- list(enabled = TRUE, params = params)
  do.call(CMiNet, c(list(data = X, quantitative = TRUE, TT = 0.95), methods))
}

# Preserved normalization used by this existing fixed-model stability analysis.
norm_log <- function(X) {
  X[X == 0] <- 0.5
  sweep(log(X), 1, rowMeans(log(X)), "-")
}

fit_phymapnet_phylo <- function(X, D) {
  p <- ncol(X)
  C <- exp(-D / PHYMAP_PARAMS$alpha)
  y <- norm_log(X)
  Iso <- phymapnet(
    y, PHYMAP_PARAMS$nue_k * p, PHYMAP_PARAMS$epsilon1,
    PHYMAP_PARAMS$epsilon2, C
  )
  clean_B(sparse_quantile(y, Iso, PHYMAP_PARAMS$th)$theta)
}

fit_phymapnet_uniform <- function(X) {
  p <- ncol(X)
  taxa <- colnames(X)
  D_uniform <- matrix(1 / p, p, p, dimnames = list(taxa, taxa))
  diag(D_uniform) <- 0
  C <- exp(-D_uniform / PHYMAP_PARAMS$alpha)
  y <- norm_log(X)
  Iso <- phymapnet(
    y, PHYMAP_PARAMS$nue_k * p, PHYMAP_PARAMS$epsilon1,
    PHYMAP_PARAMS$epsilon2, C
  )
  clean_B(sparse_quantile(y, Iso, PHYMAP_PARAMS$th)$theta)
}

fit_se_mb <- function(X) {
  fit <- do.call(spiec.easi, c(list(data = as.matrix(X)), SE_MB_PARAMS))
  clean_B(as.matrix(fit$refit$stars))
}

fit_se_glasso <- function(X) {
  fit <- do.call(spiec.easi, c(list(data = as.matrix(X)), SE_GL_PARAMS))
  clean_B(as.matrix(fit$refit$stars))
}

fit_sparcc <- function(X) {
  extract_cminet_B(run_cminet_single(X, "sparcc", SPARCC_PARAMS), TT = 0.95)
}

fit_cclasso <- function(X) {
  extract_cminet_B(run_cminet_single(X, "cclasso", CCLASSO_PARAMS), TT = 0.95)
}

fit_cmi <- function(X) {
  extract_cminet_B(run_cminet_single(X, "c_MI", CMI_PARAMS), TT = 0.95)
}

run_method <- function(method_name, dataset, X_orig, taxa, fit_fn, replicates) {
  cat("    Fitting reference...\n")
  ref_B <- tryCatch(
    fit_fn(X_orig),
    error = function(e) {
      message("    REF ERROR: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ref_B)) {
    message("    Skipping ", method_name, ": reference failed.")
    return(NULL)
  }

  n_ref <- count_edges(ref_B)
  ev0 <- evaluate(ncol(ref_B), ref_B, ref_B)
  rows <- list(make_row(method_name, dataset, "original", 0, n_ref, n_ref, get_f1(ev0)))

  run_reps <- function(values, label) {
    if (is.null(values) || !length(values)) {
      return()
    }
    cat(sprintf("    %s: %d replicates\n", label, length(values)))
    for (i in seq_along(values)) {
      matrix_i <- as.matrix(values[[i]])
      if (!all(taxa %in% colnames(matrix_i))) {
        stop("Replicate is missing taxa for dataset: ", dataset, call. = FALSE)
      }
      Xi <- matrix_i[, taxa, drop = FALSE]
      Bi <- tryCatch(
        fit_fn(Xi),
        error = function(e) {
          message("    replicate ", i, " error: ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(Bi)) {
        next
      }
      ev <- evaluate(ncol(Bi), Bi, ref_B)
      rows[[length(rows) + 1L]] <<- make_row(
        method_name, dataset, label, i, n_ref, count_edges(Bi), get_f1(ev)
      )
    }
  }

  run_reps(replicates$boot, "bootstrap")
  run_reps(replicates$noisy, "noisy")
  bind_rows(rows)
}

cat("Bootstrap stability workflow\n")
cat("  Mode:", if (smoke_test) "smoke test" else "full configured run", "\n")
cat("  Datasets:", paste(datasets, collapse = ", "), "\n")
cat("  Methods:", paste(method_keys, collapse = ", "), "\n")
cat("  Maximum replicates per split:", max_reps, "\n")
cat("  Output:", dir_out, "\n")

all_results <- list()
for (dataset in datasets) {
  cat("\nDataset:", dataset, "\n")
  inputs <- read_filtered_inputs(dataset)
  X_orig <- inputs$otu
  D_orig <- inputs$dist
  taxa <- inputs$taxa
  cat("  Taxa:", length(taxa), "\n")

  all_reps <- load_all_replicates(dataset)
  cat("  Available bootstrap:", length(all_reps$boot),
      "| noisy:", length(all_reps$noisy), "\n")
  reps <- select_subset(all_reps, n = max_reps, seed = 42)
  cat("  Selected bootstrap:", length(reps$boot),
      "| noisy:", length(reps$noisy), " (seed = 42)\n")

  saveRDS(
    list(
      bootstrap_indices = reps$boot_indices,
      noisy_indices = reps$noisy_indices,
      max_reps = max_reps,
      seed = 42
    ),
    file.path(dir_out, paste0(dataset, "_selected_replicate_indices.rds"))
  )

  method_specs <- list(
    phymap_phylo = list(
      name = "PhyMapNet (phylo)",
      fn = function(X) fit_phymapnet_phylo(X, D_orig)
    ),
    phymap_uniform = list(name = "PhyMapNet (uniform)", fn = fit_phymapnet_uniform),
    spieceasi_mb = list(name = "SPIEC-EASI mb", fn = fit_se_mb),
    spieceasi_glasso = list(name = "SPIEC-EASI glasso", fn = fit_se_glasso),
    sparcc = list(name = "SparCC", fn = fit_sparcc),
    cclasso = list(name = "ccLasso", fn = fit_cclasso),
    c_mi = list(name = "c_MI", fn = fit_cmi)
  )

  dataset_results <- list()
  for (method_key in method_keys) {
    method <- method_specs[[method_key]]
    cat("\n  Method:", method$name, "\n")
    result <- run_method(method$name, dataset, X_orig, taxa, method$fn, reps)
    if (!is.null(result)) {
      dataset_results[[length(dataset_results) + 1L]] <- result
    }
  }

  if (length(dataset_results)) {
    dataset_result <- bind_rows(dataset_results)
    write.csv(
      dataset_result,
      file.path(dir_out, paste0(dataset, "_all_methods_results.csv")),
      row.names = FALSE
    )
    all_results[[dataset]] <- dataset_result
  }
}

combined <- bind_rows(all_results)
write.csv(
  combined,
  file.path(dir_out, "all_datasets_all_methods_results.csv"),
  row.names = FALSE
)
cat("\nCompleted. Results written to:", dir_out, "\n")
