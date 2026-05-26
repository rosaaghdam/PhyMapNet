# =============================================================================
# PhyMapNet - generate reliability masters for CMiNet overlap analysis
#
# Automatically loops over:
#   Datasets   : smoking, caffeine, hmp_stool
#   Norms      : log, tss, clr     (run separately — not robust)
#   Kernels    : gaussian, laplacian (run separately — not robust)
#
# For each combination (dataset × norm × kernel) saves ONE master CSV:
#   taxa1 | taxa2 | reliability | cminet_096 | cminet_097 | cminet_098
#
# Output files named: <dataset>_<norm>_<kernel>_master.csv
# Example: smoking_log_gaussian_master.csv
#
# Total runs: 3 datasets × 3 norms × 2 kernels = 18 files
#
# This is a time-consuming network-inference script. Do not rerun it when
# completed master files already exist under result/reliability_master/phymap.
# =============================================================================

rm(list = ls())
library(phymapnet)
library(ape)
library(dplyr)

# =============================================================================
# USER CONFIG — edit only this block
# =============================================================================
BASE_DIR <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)


DATASETS <- c("smoking", "caffeine", "hmp_stool")
NORMS    <- c("log", "tss", "clr")
KERNELS  <- c("gaussian", "laplacian")

# Robust hyperparameter ranges (confirmed from sensitivity analysis)
ALPHA_RANGE <- seq(0.05, 0.10, 0.01)
K_RANGE     <- 2:10
EPS1_RANGE  <- seq(0, 1, 0.1)
EPS2_RANGE  <- seq(0, 1, 0.1)
TH_SPARSITY <- 0.95

# CMiNet threshold subfolders
CMINET_THS <- c("th=0.96", "th=0.97", "th=0.98")
# =============================================================================

data_dir <- file.path(BASE_DIR, "data")
out_dir  <- file.path(BASE_DIR, "result", "reliability_master", "phymap")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# HELPERS
# =============================================================================
make_edge_id <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  ifelse(a < b, paste(a, b, sep="__"), paste(b, a, sep="__"))
}

load_dataset <- function(ds) {
  if (ds == "smoking") {
    otu_df <- read.csv(file.path(data_dir,"filter","smoking_otu_filtered.csv"),
                       row.names=1)
    load(file.path(data_dir,"original_data","data_smoking.RData"))
    colnames(otu_df) <- gsub("^X","", colnames(otu_df))
    otu <- as.matrix(otu_df); tree <- data_smoking$tree

  } else if (ds == "caffeine") {
    otu_df <- read.csv(file.path(data_dir,"filter","caffeine_otu_filtered.csv"),
                       row.names=1)
    load(file.path(data_dir,"original_data","data_caff.RData"))
    colnames(otu_df) <- gsub("^X","", colnames(otu_df))
    otu <- as.matrix(otu_df); tree <- data_caff$tree

  } else if (ds == "hmp_stool") {
    otu_df <- read.csv(file.path(data_dir,"filter","hmp_stool_otu_filtered.csv"),
                       row.names=1, check.names=FALSE)
    load(file.path(data_dir,"original_data","data_hmp.RData"))
    colnames(otu_df) <- gsub("^X","", colnames(otu_df))
    otu <- as.matrix(otu_df); tree <- data_hmp$tree

  } else stop("Unknown dataset: ", ds)

  # Prune tree to OTU columns — applied to ALL datasets
  common <- intersect(colnames(otu), tree$tip.label)
  if (!length(common)) stop("No common taxa for: ", ds)
  tree <- ape::keep.tip(tree, common)
  otu  <- otu[, tree$tip.label, drop=FALSE]
  stopifnot(identical(colnames(otu), tree$tip.label))
  list(otu=otu, tree=tree)
}

read_cminet <- function(file) {
  if (!file.exists(file)) { warning("Not found: ", file); return(NULL) }
  df <- read.csv(file, check.names=FALSE, stringsAsFactors=FALSE)
  cn <- tolower(colnames(df))
  fc <- which(cn %in% c("from","taxa1","node1","source"))
  tc <- which(cn %in% c("to","taxa2","node2","target"))
  wc <- which(cn %in% c("weight","value","score","methods"))
  if (!length(fc)||!length(tc)||!length(wc)) {
    df <- df[,1:3]; colnames(df) <- c("taxa1","taxa2","cminet_value")
  } else {
    df <- df[,c(fc[1],tc[1],wc[1])]; colnames(df) <- c("taxa1","taxa2","cminet_value")
  }
  df %>%
    mutate(taxa1        = as.character(taxa1),
           taxa2        = as.character(taxa2),
           cminet_value = suppressWarnings(as.numeric(cminet_value)),
           edge         = make_edge_id(taxa1, taxa2)) %>%
    filter(!is.na(cminet_value)) %>%
    distinct(edge, .keep_all=TRUE) %>%
    select(edge, cminet_value)
}

# =============================================================================
# MAIN LOOP: dataset × norm × kernel
# =============================================================================
# Track progress
total_runs   <- length(DATASETS) * length(NORMS) * length(KERNELS)
current_run  <- 0
run_log      <- list()

for (ds in DATASETS) {

  # Load dataset ONCE per dataset (reused across all norm × kernel combos)
  cat("\n", strrep("#",60), "\n")
  cat(" Loading dataset:", ds, "\n")
  cat(strrep("#",60), "\n")
  dat    <- load_dataset(ds)
  otu    <- dat$otu
  tree   <- dat$tree
  cat("  samples:", nrow(otu), "| taxa:", ncol(otu), "\n")

  # Load CMiNet files ONCE per dataset (reused across norm × kernel combos)
  cmi_data <- list()
  for (th in CMINET_THS) {
    cmi_file <- file.path(BASE_DIR,"results","cminet",th,
                          paste0(ds,"_edge_list.csv"))
    cmi_data[[th]] <- read_cminet(cmi_file)
    if (!is.null(cmi_data[[th]])) {
      cat(sprintf("  CMiNet %s: %d edges loaded\n", th, nrow(cmi_data[[th]])))
    }
  }

  for (norm in NORMS) {
    for (kernel in KERNELS) {

      current_run <- current_run + 1
      run_label   <- paste0(ds, "_", norm, "_", kernel)
      out_file    <- file.path(out_dir, paste0(run_label, "_master.csv"))

      cat("\n", strrep("=",55), "\n")
      cat(sprintf(" [%d/%d] %s\n", current_run, total_runs, run_label))
      cat(strrep("=",55), "\n")

      # Skip if already done
      if (file.exists(out_file)) {
        cat("  Already exists - skipping.\n")
        run_log[[run_label]] <- "skipped"
        next
      }

      if (any(vapply(cmi_data, is.null, logical(1)))) {
        stop(
          "Cannot generate a missing master result: raw CMiNet edge-list input files ",
          "are missing for dataset ", ds, ". Do not replace existing master results ",
          "with zero-filled CMiNet columns.",
          call. = FALSE
        )
      }

      # Run reliability
      res <- tryCatch({
        phymapnet_reliability(
          otu            = otu,
          tree           = tree,
          th_fixed       = TH_SPARSITY,
          alpha_range    = ALPHA_RANGE,
          k_range        = K_RANGE,
          epsilon1_range = EPS1_RANGE,
          epsilon2_range = EPS2_RANGE,
          kernels        = kernel,
          normalizations = norm,
          consensus_cut  = 0.01,
          prune_tree     = TRUE,
          progress       = TRUE
        )
      }, error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        NULL
      })

      if (is.null(res)) {
        run_log[[run_label]] <- "error"
        next
      }

      # Build edge list
      edge_df <- as.data.frame(res$edge_list)
      colnames(edge_df)[1:3] <- c("taxa1","taxa2","reliability")
      edge_df <- edge_df %>%
        mutate(
          taxa1       = as.character(taxa1),
          taxa2       = as.character(taxa2),
          reliability = as.numeric(reliability),
          edge        = make_edge_id(taxa1, taxa2)
        ) %>%
        filter(!is.na(reliability)) %>%
        distinct(edge, .keep_all=TRUE) %>%
        arrange(desc(reliability)) %>%
        select(taxa1, taxa2, edge, reliability)

      cat("  Reliability edges:", nrow(edge_df), "\n")

      # Attach CMiNet scores from each threshold
      for (th in CMINET_THS) {
        col_name <- paste0("cminet_", gsub("[^0-9]","", th))  # cminet_096/097/098
        cmi      <- cmi_data[[th]]

        if (!is.null(cmi)) {
          edge_df <- edge_df %>%
            left_join(cmi %>% rename(!!col_name := cminet_value), by="edge") %>%
            mutate(!!col_name := as.integer(
              replace(.data[[col_name]], is.na(.data[[col_name]]), 0L)
            ))
          cat(sprintf("  CMiNet %s: %d edges with score > 0\n",
                      th, sum(edge_df[[col_name]] > 0)))
        } else {
          edge_df[[col_name]] <- 0L
        }
      }

      # Save master CSV: 6 columns as specified
      master <- edge_df %>%
        select(taxa1, taxa2, reliability,
               cminet_096, cminet_097, cminet_098)

      write.csv(master, out_file, row.names=FALSE)
      cat("  Saved:", basename(out_file), "\n")
      run_log[[run_label]] <- "done"
    }
  }
}

# =============================================================================
# RUN SUMMARY
# =============================================================================
cat("\n", strrep("=",60), "\n")
cat("ALL RUNS COMPLETE\n")
cat(strrep("=",60), "\n")
cat(sprintf("Total: %d | Done: %d | Skipped: %d | Errors: %d\n",
            total_runs,
            sum(unlist(run_log) == "done"),
            sum(unlist(run_log) == "skipped"),
            sum(unlist(run_log) == "error")))

cat("\nFiles saved to:", out_dir, "\n\n")
cat("File naming: <dataset>_<norm>_<kernel>_master.csv\n")
cat("Example files:\n")
for (ds in DATASETS)
  for (norm in NORMS)
    for (kernel in KERNELS)
      cat(sprintf("  %s_%s_%s_master.csv\n", ds, norm, kernel))
