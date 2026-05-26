# =============================================================================
# Runtime comparison across seven network methods
#
# Methods:
#   1. PhyMapNet (phylo)    — phylogenetic kernel, best params
#   2. PhyMapNet (uniform)  — uniform 1/p kernel, ablation
#   3. SPIEC-EASI mb        — neighbourhood selection
#   4. SPIEC-EASI glasso    — graphical lasso
#   5. SparCC               — compositional correlation
#   6. ccLasso              — compositional graphical model
#   7. c_MI                 — mutual information
#
# Each method: fit reference on original → fit on all replicates → F1 vs ref
# =============================================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages(library(dplyr))

BASE_DIR <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
source(file.path(BASE_DIR, "3-bootstrap", "bootstrap_phymapnet_helpers.R"))

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(CMiNet)
})

# =============================================================================
# USER CONFIG
# =============================================================================
smoke_test <- tolower(Sys.getenv("PHYMAPNET_RUNTIME_SMOKE_TEST", unset = "true")) %in%
  c("true", "t", "yes", "y", "1")
DATASETS <- trimws(strsplit(
  Sys.getenv(
    "PHYMAPNET_RUNTIME_DATASETS",
    unset = if (smoke_test) "caffeine" else "smoking,caffeine,hmp_stool"
  ),
  ",",
  fixed = TRUE
)[[1]])
MAX_REPS <- as.integer(Sys.getenv("PHYMAPNET_RUNTIME_MAX_REPS", unset = "1"))
if (is.na(MAX_REPS) || MAX_REPS < 1) {
  stop("PHYMAPNET_RUNTIME_MAX_REPS must be a positive integer.", call. = FALSE)
}
METHOD_KEYS <- trimws(strsplit(
  Sys.getenv(
    "PHYMAPNET_RUNTIME_METHODS",
    unset = if (smoke_test) "phymap_phylo" else
      "phymap_phylo,phymap_uniform,spieceasi_mb,spieceasi_glasso,sparcc,cclasso,c_mi"
  ),
  ",",
  fixed = TRUE
)[[1]])

PHYMAP_PARAMS <- list(
  alpha    = 0.1,
  kernel   = "laplacian",
  nue_k    = 10,
  epsilon1 = 1,
  epsilon2 = 0.4,
  th       = 0.99
)

SE_MB_PARAMS <- list(
  method           = "mb",
  lambda.min.ratio = 1e-2,
  nlambda          = 15,
  pulsar.params    = list(rep.num = 20, ncores = 4)
)

SE_GL_PARAMS <- list(
  method           = "glasso",
  lambda.min.ratio = 1e-2,
  nlambda          = 15,
  pulsar.params    = list(rep.num = 50, ncores = 4)
)

SPARCC_PARAMS <- list(
  imax  = 20,
  kmax  = 10,
  alpha = 0.1,
  Vmin  = 1e-4
)

CCLASSO_PARAMS <- list(
  counts  = FALSE,
  pseudo  = 0.5,
  k_cv    = 3,
  lam_int = c(1e-4, 1),
  k_max   = 20,
  n_boot  = 20
)

CMI_PARAMS <- list(
  quantitative = TRUE,
  q1           = 0.9,
  q2           = 0.95
)
# =============================================================================

dir_filter <- file.path(BASE_DIR, "data", "filter")
dir_boot   <- file.path(BASE_DIR, "data", "generate", "bootstrap")
dir_noisy  <- file.path(BASE_DIR, "data", "generate", "noisy")
dir_out_default <- if (smoke_test) {
  file.path(tempdir(), "phymapnet_method_runtime_smoke")
} else {
  file.path(BASE_DIR, "result", "running_time", "method_comparison")
}
dir_out <- Sys.getenv("PHYMAPNET_RUNTIME_OUTPUT_DIR", unset = dir_out_default)
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPERS
# =============================================================================
count_edges <- function(B) sum(B[upper.tri(B)])

clean_B <- function(B) {
  B <- as.matrix(B)
  diag(B) <- 0L
  B[lower.tri(B)] <- 0L
  (B != 0) * 1L
}

get_f1 <- function(ev) {
  m   <- ev[[2]]
  hit <- which(tolower(names(m)) %in% c("f_score","f.score","f1","fscore"))
  if (length(hit)) return(as.numeric(m[[hit[1]]]))
  as.numeric(m[[ncol(m)]])
}

make_row <- function(method, dataset, split, rep, n_ref, n_est, f1) {
  data.frame(method=method, dataset=dataset, split=split,
             replicate=rep, n_edges_ref=n_ref,
             n_edges_est=n_est, F1=f1, stringsAsFactors=FALSE)
}

# Load ALL replicates from all RDS files for a dataset
load_all_replicates <- function(ds) {
  flatten_rds <- function(files) {
    if (!length(files)) return(NULL)
    out <- list()
    for (f in files) {
      obj <- readRDS(f)
      if (is.list(obj)) out <- c(out, obj) else out[[length(out)+1]] <- obj
    }
    out
  }
  list(
    boot  = flatten_rds(Sys.glob(file.path(dir_boot,
                                           paste0(ds,"_bootstrap_*.rds")))),
    noisy = flatten_rds(Sys.glob(file.path(dir_noisy,
                                           paste0(ds,"_noisy_*.rds"))))
  )
}

# Select a fixed random subset of MAX_REPS replicates
# Called ONCE per dataset — same indices used for ALL methods
select_subset <- function(reps, n = MAX_REPS, seed = 42) {
  set.seed(seed)
  subset_list <- function(lst) {
    if (is.null(lst) || !length(lst)) return(list(data = NULL, indices = integer()))
    idx <- sort(sample(seq_along(lst), min(n, length(lst)), replace=FALSE))
    list(data = lst[idx], indices = idx)
  }
  boot <- subset_list(reps$boot)
  noisy <- subset_list(reps$noisy)
  list(
    boot = boot$data,
    noisy = noisy$data,
    boot_indices = boot$indices,
    noisy_indices = noisy$indices
  )
}

# Extract binary adjacency from CMiNet result
extract_cminet_B <- function(res, TT = 0.95) {
  mat <- NULL
  for (sl in c("score","weighted_adj","adj","network","W","result")) {
    if (is.matrix(res[[sl]]) && nrow(res[[sl]]) > 1) { mat <- res[[sl]]; break }
  }
  if (is.null(mat)) {
    for (nm in names(res))
      if (is.matrix(res[[nm]]) && nrow(res[[nm]]) > 1) { mat <- res[[nm]]; break }
  }
  if (is.null(mat))
    stop("Cannot extract matrix from CMiNet result. Slots: ",
         paste(names(res), collapse=", "))
  ut_pos <- mat[upper.tri(mat)]
  ut_pos <- ut_pos[ut_pos > 0]
  if (!length(ut_pos)) return(clean_B(mat * 0))
  thr <- quantile(ut_pos, TT, na.rm=TRUE)
  clean_B((mat >= thr) * 1L)
}

# Shared CMiNet call — only the requested method enabled
run_cminet_single <- function(X, method_name, params) {
  all_off <- list(
    pearson          = list(enabled=FALSE),
    spearman         = list(enabled=FALSE),
    bicor            = list(enabled=FALSE),
    sparcc           = list(enabled=FALSE),
    spiecEasi_mb     = list(enabled=FALSE),
    spiecEasi_glasso = list(enabled=FALSE),
    spring           = list(enabled=FALSE),
    gcoda            = list(enabled=FALSE),
    c_MI             = list(enabled=FALSE),
    cclasso          = list(enabled=FALSE)
  )
  all_off[[method_name]] <- list(enabled=TRUE, params=params)

  do.call(CMiNet, c(
    list(data=X, quantitative=TRUE, TT=0.95),
    all_off
  ))
}

# =============================================================================
# NORMALISATION
# =============================================================================
norm_log <- function(X) {
  X[X == 0] <- 0.5
  sweep(log(X), 1, rowMeans(log(X)), "-")
}

# =============================================================================
# FIT FUNCTIONS
# =============================================================================
fit_phymapnet_phylo <- function(X, D) {
  p   <- ncol(X)
  nue <- PHYMAP_PARAMS$nue_k * p
  C   <- exp(-D / PHYMAP_PARAMS$alpha)
  y   <- norm_log(X)
  Iso <- phymapnet(y, nue, PHYMAP_PARAMS$epsilon1,
                   PHYMAP_PARAMS$epsilon2, C)
  clean_B(sparse_quantile(y, Iso, PHYMAP_PARAMS$th)$theta)
}

fit_phymapnet_uniform <- function(X) {
  p     <- ncol(X)
  taxa  <- colnames(X)
  D_uni <- matrix(1/p, p, p, dimnames=list(taxa,taxa))
  diag(D_uni) <- 0
  C   <- exp(-D_uni / PHYMAP_PARAMS$alpha)
  y   <- norm_log(X)
  Iso <- phymapnet(y, PHYMAP_PARAMS$nue_k*p,
                   PHYMAP_PARAMS$epsilon1, PHYMAP_PARAMS$epsilon2, C)
  clean_B(sparse_quantile(y, Iso, PHYMAP_PARAMS$th)$theta)
}

fit_se_mb <- function(X) {
  fit <- do.call(spiec.easi, c(list(data=as.matrix(X)), SE_MB_PARAMS))
  clean_B(as.matrix(fit$refit$stars))
}

fit_se_glasso <- function(X) {
  fit <- do.call(spiec.easi, c(list(data=as.matrix(X)), SE_GL_PARAMS))
  clean_B(as.matrix(fit$refit$stars))
}

fit_sparcc <- function(X) {
  res <- run_cminet_single(X, "sparcc", SPARCC_PARAMS)
  extract_cminet_B(res, TT=0.95)
}

fit_cclasso <- function(X) {
  res <- run_cminet_single(X, "cclasso", CCLASSO_PARAMS)
  extract_cminet_B(res, TT=0.95)
}

fit_cmi <- function(X) {
  res <- run_cminet_single(X, "c_MI", CMI_PARAMS)
  extract_cminet_B(res, TT=0.95)
}

# =============================================================================
# GENERIC RUNNER
# =============================================================================
run_method <- function(method_name, dataset, taxa, fit_fn, replicates) {

  X_file <- file.path(dir_filter, paste0(dataset, "_otu_filtered.csv"))
  X_orig <- as.matrix(read.csv(X_file, row.names=1, check.names=FALSE))
  X_orig <- X_orig[, taxa, drop=FALSE]

  cat("    Fitting reference...\n")
  ref_B <- tryCatch(
    fit_fn(X_orig),
    error = function(e) { message("    REF ERROR: ", conditionMessage(e)); NULL }
  )
  if (is.null(ref_B)) {
    message("    Skipping ", method_name, " — reference failed.")
    return(NULL)
  }

  n_ref <- count_edges(ref_B)
  ev0   <- evaluate(ncol(ref_B), ref_B, ref_B)
  rows  <- list(make_row(method_name, dataset, "original", 0,
                         n_ref, n_ref, get_f1(ev0)))

  run_reps <- function(lst, label) {
    if (is.null(lst) || !length(lst)) return()
    cat(sprintf("    %s: %d replicates\n", label, length(lst)))
    for (i in seq_along(lst)) {
      Xi <- tryCatch(
        { m <- lst[[i]]; m[, taxa, drop=FALSE] },
        error = function(e) NULL
      )
      if (is.null(Xi)) next
      Bi <- tryCatch(
        fit_fn(Xi),
        error = function(e) {
          message("    rep ", i, " error: ", conditionMessage(e)); NULL
        }
      )
      if (is.null(Bi)) next
      ev <- evaluate(ncol(Bi), Bi, ref_B)
      rows[[length(rows)+1]] <<-
        make_row(method_name, dataset, label, i,
                 n_ref, count_edges(Bi), get_f1(ev))
    }
  }

  run_reps(replicates$boot,  "bootstrap")
  run_reps(replicates$noisy, "noisy")
  bind_rows(rows)
}

# =============================================================================
# MAIN LOOP
# =============================================================================
all_results <- list()
all_timing  <- list()

cat("Method runtime comparison workflow\n")
cat("  Mode:", if (smoke_test) "smoke test" else "full configured run", "\n")
cat("  Datasets:", paste(DATASETS, collapse = ", "), "\n")
cat("  Methods:", paste(METHOD_KEYS, collapse = ", "), "\n")
cat("  Output:", dir_out, "\n")

for (ds in DATASETS) {
  cat("\n", strrep("=", 60), "\n")
  cat(" Dataset:", ds, "\n")
  cat(strrep("=", 60), "\n")

  X_file <- file.path(dir_filter, paste0(ds, "_otu_filtered.csv"))
  D_file <- file.path(dir_filter, paste0(ds, "_dist_matrix.csv"))
  X_orig <- as.matrix(read.csv(X_file, row.names=1, check.names=FALSE))
  D_orig <- as.matrix(read.csv(D_file, row.names=1, check.names=FALSE))
  if (!identical(sort(colnames(X_orig)), sort(rownames(D_orig))) ||
      !identical(sort(colnames(X_orig)), sort(colnames(D_orig)))) {
    stop("OTU and distance-matrix taxa do not match for dataset: ", ds, call. = FALSE)
  }
  if (!isTRUE(all.equal(D_orig, t(D_orig), tolerance = 1e-10, check.attributes = FALSE))) {
    stop("Distance matrix is not symmetric for dataset: ", ds, call. = FALSE)
  }
  taxa   <- colnames(X_orig)
  X_orig <- X_orig[, taxa, drop=FALSE]
  D_orig <- D_orig[taxa, taxa, drop=FALSE]
  cat("  taxa:", length(taxa), "\n")

  all_reps <- load_all_replicates(ds)
  cat("  Available — bootstrap:", length(all_reps$boot),
      "| noisy:", length(all_reps$noisy), "\n")

  reps <- select_subset(all_reps, n=MAX_REPS, seed=42)
  cat("  Selected  — bootstrap:", length(reps$boot),
      "| noisy:", length(reps$noisy),
      " (seed=42, same for all methods)\n")

  saveRDS(list(
    bootstrap_indices = reps$boot_indices,
    noisy_indices = reps$noisy_indices,
    max_reps = MAX_REPS,
    seed = 42
  ), file.path(dir_out, paste0(ds, "_selected_replicate_indices.rds")))

  methods <- list(
    phymap_phylo = list(name = "PhyMapNet (phylo)",
         fn   = function(X) fit_phymapnet_phylo(X, D_orig)),
    phymap_uniform = list(name = "PhyMapNet (uniform)",
         fn   = fit_phymapnet_uniform),
    spieceasi_mb = list(name = "SPIEC-EASI mb",
         fn   = fit_se_mb),
    spieceasi_glasso = list(name = "SPIEC-EASI glasso",
         fn   = fit_se_glasso),
    sparcc = list(name = "SparCC",
         fn   = fit_sparcc),
    cclasso = list(name = "ccLasso",
         fn   = fit_cclasso),
    c_mi = list(name = "c_MI",
         fn   = fit_cmi)
  )
  invalid_methods <- setdiff(METHOD_KEYS, names(methods))
  if (length(invalid_methods)) {
    stop("Unknown method keys: ", paste(invalid_methods, collapse = ", "), call. = FALSE)
  }
  methods <- methods[METHOD_KEYS]

  ds_rows   <- list()
  ds_timing <- list()

  for (m_idx in seq_along(methods)) {
    m <- methods[[m_idx]]
    cat(sprintf("\n  [%d/%d] %s\n", m_idx, length(methods), m$name))

    t0      <- proc.time()
    res     <- run_method(m$name, ds, taxa, m$fn, reps)
    elapsed <- (proc.time() - t0)[["elapsed"]]

    cat(sprintf("    Time: %.1f s (%.2f min)\n", elapsed, elapsed / 60))

    ds_timing[[length(ds_timing) + 1]] <- data.frame(
      dataset     = ds,
      method      = m$name,
      n_taxa      = length(taxa),
      n_boot_reps = length(reps$boot),
      n_noisy_reps = length(reps$noisy),
      elapsed_sec = round(elapsed, 1),
      elapsed_min = round(elapsed / 60, 2),
      stringsAsFactors = FALSE
    )

    if (!is.null(res)) ds_rows[[length(ds_rows) + 1]] <- res
  }

  if (length(ds_rows)) {
    ds_result <- bind_rows(ds_rows)
    out_csv   <- file.path(dir_out, paste0(ds, "_all_methods_results.csv"))
    write.csv(ds_result, out_csv, row.names=FALSE)
    cat("\n  Saved:", basename(out_csv), "\n")
    all_results[[ds]] <- ds_result
  }

  timing_ds <- bind_rows(ds_timing)
  write.csv(timing_ds,
            file.path(dir_out, paste0(ds, "_method_timing.csv")),
            row.names = FALSE)
  cat("  Timing saved:", paste0(ds, "_method_timing.csv"), "\n")
  all_timing[[ds]] <- timing_ds
}

# =============================================================================
# COMBINE + SAVE
# =============================================================================
combined <- bind_rows(all_results)
write.csv(combined,
          file.path(dir_out, "all_datasets_all_methods_results.csv"),
          row.names = FALSE)
cat("\nCombined results saved.\n")

timing_combined <- bind_rows(all_timing)
write.csv(timing_combined,
          file.path(dir_out, "all_datasets_method_timing.csv"),
          row.names = FALSE)
cat("Combined timing saved to all_datasets_method_timing.csv\n")

# Print summary table to console
cat("\n--- Runtime Summary ---\n")
print(
  timing_combined %>%
    select(dataset, method, n_taxa, elapsed_sec, elapsed_min) %>%
    arrange(dataset, elapsed_sec),
  row.names = FALSE
)
