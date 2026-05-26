# ============================================================
# Run PhyMapNet runtime evaluation on different numbers of HMP taxa
# Save reliability matrices and edge tables
# Includes safer numerical grid for large taxa sets
# ============================================================

rm(list = ls())

library(phymapnet)
library(ape)
library(dplyr)

# -----------------------------
# Run mode and paths
# -----------------------------

base_dir <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
smoke_test <- tolower(Sys.getenv("PHYMAPNET_RUNTIME_SMOKE_TEST", unset = "true")) %in%
  c("true", "t", "yes", "y", "1")
output_default <- if (smoke_test) {
  file.path(tempdir(), "phymapnet_hmp_runtime_smoke")
} else {
  file.path(base_dir, "result", "hmp_taxa_sensitivity")
}
outdir <- Sys.getenv("PHYMAPNET_RUNTIME_OUTPUT_DIR", unset = output_default)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load original HMP OTU table and tree
# -----------------------------

load(file.path(base_dir, "data", "original_data", "data_hmp.RData"))

otu <- as.matrix(data_hmp$otu.tab)
tree <- data_hmp$tree

colnames(otu) <- gsub("^X", "", colnames(otu))

common <- intersect(colnames(otu), tree$tip.label)

otu <- otu[, common, drop = FALSE]
tree <- ape::keep.tip(tree, common)
otu <- otu[, tree$tip.label, drop = FALSE]

stopifnot(identical(colnames(otu), tree$tip.label))

cat("Original aligned HMP data:", nrow(otu), "samples x", ncol(otu), "taxa\n")

# Remove all-zero taxa
otu <- otu[, colSums(otu) > 0, drop = FALSE]

common <- intersect(colnames(otu), tree$tip.label)

otu <- otu[, common, drop = FALSE]
tree <- ape::keep.tip(tree, common)
otu <- otu[, tree$tip.label, drop = FALSE]

stopifnot(identical(colnames(otu), tree$tip.label))

cat("After removing all-zero taxa:", nrow(otu), "samples x", ncol(otu), "taxa\n")

# -----------------------------
# Settings
# -----------------------------

sizes <- if (smoke_test) {
  as.integer(Sys.getenv("PHYMAPNET_RUNTIME_TAXA_SIZES", unset = "20"))
} else {
  as.integer(strsplit(
    Sys.getenv("PHYMAPNET_RUNTIME_TAXA_SIZES", unset = "100,300,500,1000"),
    ",",
    fixed = TRUE
  )[[1]])
}
if (any(is.na(sizes)) || any(sizes < 2)) {
  stop("PHYMAPNET_RUNTIME_TAXA_SIZES must contain integers greater than 1.", call. = FALSE)
}

norm_use <- "log"
kernel_use <- "laplacian"

# Original grid
alpha_range_default <- seq(0.01, 0.12, 0.01)
epsilon_default <- seq(0, 1, 0.1)

# Safer grid for larger taxa sets
alpha_range_large <- seq(0.05, 0.12, 0.01)
epsilon_large <- seq(0.1, 1, 0.1)

# Fallback grid if singular matrix still happens
alpha_range_fallback <- seq(0.08, 0.20, 0.02)
epsilon_fallback <- seq(0.2, 1, 0.2)

k_range_default <- 2:10
k_range_fallback <- 2:8

th_fixed <- 0.95
consensus_cut <- 0.65

if (smoke_test) {
  alpha_range_default <- 0.1
  epsilon_default <- 0.1
  k_range_default <- 2
}

# -----------------------------
# Helper: extract reliability matrix from result
# -----------------------------

extract_reliability_matrix <- function(res) {
  candidate_names <- c(
    "reliability_matrix",
    "reliability",
    "R",
    "rel_mat",
    "consensus_matrix",
    "consensus",
    "score_matrix",
    "matrix"
  )

  if (is.matrix(res) || is.data.frame(res)) {
    return(as.matrix(res))
  }

  if (is.list(res)) {
    for (nm in candidate_names) {
      if (nm %in% names(res)) {
        x <- res[[nm]]
        if (is.matrix(x) || is.data.frame(x)) {
          message("Using result element: ", nm)
          return(as.matrix(x))
        }
      }
    }

    matrix_candidates <- names(res)[
      sapply(res, function(x) is.matrix(x) || is.data.frame(x))
    ]

    if (length(matrix_candidates) > 0) {
      message("Using matrix-like result element: ", matrix_candidates[1])
      return(as.matrix(res[[matrix_candidates[1]]]))
    }
  }

  stop(
    "Could not find a reliability matrix inside phymapnet_reliability output. ",
    "Inspect the object with names(res) and str(res)."
  )
}

# -----------------------------
# Helper: matrix to edge table
# -----------------------------

matrix_to_edge_table <- function(mat) {
  mat <- as.matrix(mat)
  diag(mat) <- NA

  idx <- which(
    upper.tri(mat) & !is.na(mat),
    arr.ind = TRUE
  )

  data.frame(
    taxon1 = rownames(mat)[idx[, 1]],
    taxon2 = colnames(mat)[idx[, 2]],
    reliability = mat[idx],
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Rank taxa by prevalence and abundance
# -----------------------------

prevalence <- colSums(otu > 0)
abundance <- colSums(otu)

taxa_rank <- names(
  sort(
    prevalence + abundance / max(abundance),
    decreasing = TRUE
  )
)

taxa_rank_table <- data.frame(
  OTU_ID = taxa_rank,
  prevalence = prevalence[taxa_rank],
  total_abundance = abundance[taxa_rank],
  stringsAsFactors = FALSE
)

write.csv(
  taxa_rank_table,
  file.path(outdir, "hmp_taxa_rank_by_prevalence_abundance.csv"),
  row.names = FALSE
)

# -----------------------------
# Run sensitivity analysis
# -----------------------------

timing <- list()

cat("PhyMapNet HMP taxa runtime workflow\n")
cat("  Mode:", if (smoke_test) "smoke test" else "full configured run", "\n")
cat("  Taxa sizes:", paste(sizes, collapse = ", "), "\n")
cat("  Output:", outdir, "\n")

for (p in sizes) {
  if (p > ncol(otu)) {
    warning(
      "Skipping p = ", p,
      " because OTU table has only ", ncol(otu), " taxa."
    )
    next
  }

  message("\n==============================")
  message("Running p = ", p)
  message("==============================")

  taxa_use <- taxa_rank[1:p]

  otu_sub <- otu[, taxa_use, drop = FALSE]
  tree_sub <- ape::keep.tip(tree, taxa_use)
  otu_sub <- otu_sub[, tree_sub$tip.label, drop = FALSE]

  stopifnot(identical(colnames(otu_sub), tree_sub$tip.label))

  prefix <- paste0(
    "hmp_p", p,
    "_", norm_use,
    "_", kernel_use
  )

  write.csv(
    data.frame(OTU_ID = colnames(otu_sub)),
    file.path(outdir, paste0(prefix, "_taxa_used.csv")),
    row.names = FALSE
  )

  # Choose grid
  if (p >= 1000) {
    alpha_range_use <- alpha_range_large
    epsilon_use <- epsilon_large
    k_range_use <- k_range_default
    grid_used <- "large"
  } else {
    alpha_range_use <- alpha_range_default
    epsilon_use <- epsilon_default
    k_range_use <- k_range_default
    grid_used <- "default"
  }

  t0 <- proc.time()

  fallback_used <- FALSE

  res <- tryCatch(
    {
      phymapnet_reliability(
        otu = otu_sub,
        tree = tree_sub,
        th_fixed = th_fixed,
        alpha_range = alpha_range_use,
        k_range = k_range_use,
        epsilon1_range = epsilon_use,
        epsilon2_range = epsilon_use,
        kernels = kernel_use,
        normalizations = norm_use,
        consensus_cut = consensus_cut,
        prune_tree = TRUE,
        progress = FALSE
      )
    },
    error = function(e) {
      message("First run failed for p = ", p)
      message("Error: ", e$message)
      message("Retrying with safer fallback numerical grid...")

      fallback_used <<- TRUE
      grid_used <<- "fallback"

      phymapnet_reliability(
        otu = otu_sub,
        tree = tree_sub,
        th_fixed = th_fixed,
        alpha_range = alpha_range_fallback,
        k_range = k_range_fallback,
        epsilon1_range = epsilon_fallback,
        epsilon2_range = epsilon_fallback,
        kernels = kernel_use,
        normalizations = norm_use,
        consensus_cut = consensus_cut,
        prune_tree = TRUE,
        progress = FALSE
      )
    }
  )

  elapsed <- (proc.time() - t0)[3]

  rel_mat <- extract_reliability_matrix(res)

  # Make sure matrix names are present and aligned
  if (is.null(rownames(rel_mat))) {
    rownames(rel_mat) <- colnames(otu_sub)
  }

  if (is.null(colnames(rel_mat))) {
    colnames(rel_mat) <- colnames(otu_sub)
  }

  if (
    nrow(rel_mat) == ncol(otu_sub) &&
    ncol(rel_mat) == ncol(otu_sub)
  ) {
    rownames(rel_mat) <- colnames(otu_sub)
    colnames(rel_mat) <- colnames(otu_sub)
  }

  edge_table <- matrix_to_edge_table(rel_mat)

  saveRDS(
    res,
    file.path(outdir, paste0(prefix, "_phymapnet_result.rds"))
  )

  saveRDS(
    rel_mat,
    file.path(outdir, paste0(prefix, "_reliability_matrix.rds"))
  )

  write.csv(
    rel_mat,
    file.path(outdir, paste0(prefix, "_reliability_matrix.csv")),
    row.names = TRUE
  )

  write.csv(
    edge_table,
    file.path(outdir, paste0(prefix, "_all_reliability_edges.csv")),
    row.names = FALSE
  )

  timing[[length(timing) + 1]] <- data.frame(
    taxa = p,
    norm = norm_use,
    kernel = kernel_use,
    th_fixed = th_fixed,
    consensus_cut = consensus_cut,
    grid_used = grid_used,
    fallback_used = fallback_used,
    elapsed_sec = elapsed,
    n_edges_total = nrow(edge_table),
    stringsAsFactors = FALSE
  )

  message(sprintf(
    "size=%d norm=%s kernel=%s grid=%s fallback=%s time=%.1fs",
    p,
    norm_use,
    kernel_use,
    grid_used,
    fallback_used,
    elapsed
  ))
}

# -----------------------------
# Save timing summary
# -----------------------------

timing_df <- dplyr::bind_rows(timing)

write.csv(
  timing_df,
  file.path(outdir, "hmp_taxa_sensitivity_runtime.csv"),
  row.names = FALSE
)

cat("\nDone. Files saved to:\n")
cat(outdir, "\n")
