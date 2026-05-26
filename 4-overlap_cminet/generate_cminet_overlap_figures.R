# =============================================================================
# PhyMapNet - CMiNet overlap figures and statistics from existing master results
#
# Layout:
#   Column header  = "Dataset\n(CMiNet n=X, top-N% n=Y)"
#   Above each bar = p-value only (hypergeometric)
#   Legend         = Score 0–9 with simple labels (no method-count text)
# =============================================================================

rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# =============================================================================
# USER CONFIG
# =============================================================================
BASE_DIR <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)

DATASETS  <- c("smoking", "caffeine", "hmp_stool")
DS_LABELS <- c(smoking = "Smoking", caffeine = "Caffeine",
               hmp_stool = "HMP stool")
NORMS     <- c("log", "tss", "clr")
KERNELS   <- c("gaussian", "laplacian")

CMI_COL    <- "cminet_096"
CMI_LABEL  <- "th=0.96"
CMI_STRONG <- 6

TOP_PERCENTS <- c(0.5, 1.0, 1.5, 2.0, 3.0)
# =============================================================================

master_dir <- file.path(BASE_DIR, "result", "reliability_master", "phymap")
out_dir    <- file.path(BASE_DIR, "result", "cminet_comparison", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DS_COL <- c(Smoking = "#D32F2F", Caffeine = "#1565C0", `HMP stool` = "#E65100")

SCORE_COLS <- c(
  "0" = "#F5F5F5", "1" = "#E3F2FD", "2" = "#BBDEFB", "3" = "#90CAF9",
  "4" = "#64B5F6", "5" = "#42A5F5", "6" = "#1E88E5", "7" = "#1565C0",
  "8" = "#0D47A1", "9" = "#01237A"
)

NORM_LABELS   <- c(log = "Log", tss = "TSS", clr = "CLR")
KERNEL_LABELS <- c(gaussian = "Gaussian", laplacian = "Laplacian")

fmt_p <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "p < 0.0001",
    p < 0.001 ~ "p < 0.001",
    p < 0.01  ~ sprintf("p = %.3f", p),
    p < 0.05  ~ sprintf("p = %.3f", p),
    TRUE      ~ sprintf("p = %.2f",  p)
  )
}

# =============================================================================
# LOAD MASTER FILES
# =============================================================================
cat("Loading master files...\n")
master_all <- list()

for (ds in DATASETS) {
  for (norm in NORMS) {
    for (kernel in KERNELS) {
      key  <- paste(ds, norm, kernel, sep = "_")
      file <- file.path(master_dir,
                        paste0(ds, "_", norm, "_", kernel, "_master.csv"))
      if (!file.exists(file)) { cat("  Missing:", basename(file), "\n"); next }
      df <- read.csv(file, stringsAsFactors = FALSE)
      required_cols <- c("taxa1", "taxa2", "reliability", CMI_COL)
      if (!all(required_cols %in% colnames(df))) {
        stop("Required columns missing from: ", file, call. = FALSE)
      }
      df <- df %>%
        mutate(dataset = DS_LABELS[[ds]], norm = norm, kernel = kernel)
      master_all[[key]] <- df
      cat("  Loaded:", basename(file), "| edges:", nrow(df), "\n")
    }
  }
}
expected_file_count <- length(DATASETS) * length(NORMS) * length(KERNELS)
if (length(master_all) != expected_file_count) {
  stop(
    "Expected ", expected_file_count, " master files but found ", length(master_all),
    " in: ", master_dir,
    call. = FALSE
  )
}

# =============================================================================
# COMPUTE STATS
# =============================================================================
dist_rows <- list()
pval_rows <- list()

for (key in names(master_all)) {
  df     <- master_all[[key]]
  n_tot  <- nrow(df)
  ds_lab <- df$dataset[1]
  norm   <- df$norm[1]
  kernel <- df$kernel[1]

  p_taxa  <- length(unique(c(df$taxa1, df$taxa2)))
  E_total <- p_taxa * (p_taxa - 1) / 2
  n_B     <- sum(df[[CMI_COL]] >= CMI_STRONG)

  for (top_pct in TOP_PERCENTS) {
    n_A     <- ceiling(n_tot * top_pct / 100)
    net_top <- df %>% arrange(desc(reliability)) %>% slice(1:n_A)
    scores  <- net_top[[CMI_COL]]

    # Score distribution
    sc_tab <- table(factor(scores, levels = 0:9))
    for (sc in 0:9) {
      dist_rows[[length(dist_rows) + 1]] <- data.frame(
        dataset = ds_lab, norm = norm, kernel = kernel,
        top_pct = top_pct, n_A = n_A, n_B = n_B,
        score = sc,
        count = as.integer(sc_tab[[as.character(sc)]]),
        stringsAsFactors = FALSE
      )
    }

    # Hypergeometric p-value
    TP      <- sum(scores >= CMI_STRONG)
    p_hyper <- phyper(TP - 1, m = n_B, n = E_total - n_B,
                      k = n_A, lower.tail = FALSE)
    jaccard <- if ((n_A + n_B - TP) > 0)
      round(TP / (n_A + n_B - TP), 4) else 0

    pval_rows[[length(pval_rows) + 1]] <- data.frame(
      dataset = ds_lab, norm = norm, kernel = kernel,
      top_pct = top_pct, n_A = n_A, n_B = n_B,
      TP = TP, E_total = E_total, jaccard = jaccard,
      p_hyper = p_hyper, p_fmt = fmt_p(p_hyper),
      stringsAsFactors = FALSE
    )
  }
}

dist_df <- bind_rows(dist_rows) %>%
  group_by(dataset, norm, kernel, top_pct) %>%
  mutate(pct = count / n_A * 100) %>%
  ungroup() %>%
  mutate(
    dataset      = factor(dataset,        levels = DS_LABELS),
    norm_label   = factor(NORM_LABELS[norm],     levels = NORM_LABELS),
    kernel_label = factor(KERNEL_LABELS[kernel], levels = KERNEL_LABELS),
    score        = factor(score, levels = 0:9)
  )

pval_df <- bind_rows(pval_rows) %>%
  mutate(
    dataset      = factor(dataset,        levels = DS_LABELS),
    norm_label   = factor(NORM_LABELS[norm],     levels = NORM_LABELS),
    kernel_label = factor(KERNEL_LABELS[kernel], levels = KERNEL_LABELS)
  )

write.csv(pval_df, file.path(out_dir, "overlap_statistics.csv"), row.names = FALSE)

# =============================================================================
# PRODUCE FIGURES
# =============================================================================
for (top_pct in TOP_PERCENTS) {
  pct_str  <- gsub("\\.", "p", as.character(top_pct))

  dist_sub <- dist_df %>% filter(top_pct == !!top_pct)
  pval_sub <- pval_df %>% filter(top_pct == !!top_pct)

  # Build column header: "Dataset\n(CMiNet n=X, top-N% n=Y)"
  # n_B is same across norm×kernel for a given dataset×top_pct
  # n_A varies by norm×kernel — use the value per norm×kernel combination
  # For the header we show n_B (CMiNet) which is dataset-level
  # n_A per bar is shown in p-value label below — so header just has n_B
  ds_nb <- pval_sub %>%
    distinct(dataset, n_B, n_A) %>%
    # n_A can differ across norm×kernel; take representative (e.g. median)
    group_by(dataset, n_B) %>%
    summarise(n_A_rep = as.integer(median(n_A)), .groups = "drop") %>%
    mutate(
      ds_facet = paste0(
        as.character(dataset), "\n",
        "(CMiNet n = ", n_B,
        ",  top-", top_pct, "% reliable n = ", n_A_rep, ")"
      )
    ) %>%
    select(dataset, ds_facet)

  # Factor levels preserve dataset order
  ds_nb <- ds_nb %>%
    mutate(ds_facet = factor(ds_facet, levels = ds_nb$ds_facet))

  dist_sub <- dist_sub %>% left_join(ds_nb, by = "dataset") %>%
    mutate(ds_facet = factor(ds_facet, levels = levels(ds_nb$ds_facet)))
  pval_sub <- pval_sub %>% left_join(ds_nb, by = "dataset") %>%
    mutate(ds_facet = factor(ds_facet, levels = levels(ds_nb$ds_facet)))

  # Above each bar: p-value only
  y_ann <- 103

  fig <- ggplot() +

    # Stacked bars
    geom_col(
      data  = dist_sub,
      aes(x = norm_label, y = pct, fill = score),
      width = 0.80, colour = "white", linewidth = 0.25
    ) +

    # p-value above each bar (no n_A — it is in the header)
    geom_text(
      data = pval_sub,
      aes(x = norm_label, y = y_ann, label = p_fmt),
      size = 2.6, vjust = 0, colour = "grey15"
    ) +

    facet_grid(kernel_label ~ ds_facet) +

    # Legend: score 0–9, simple numeric labels only
    scale_fill_manual(
      values = SCORE_COLS,
      name   = "CMiNet\nscore",
      breaks = as.character(9:0),
      labels = as.character(9:0)
    ) +
    scale_y_continuous(
      limits = c(0, 116),
      breaks = seq(0, 100, 25),
      labels = paste0(seq(0, 100, 25), "%"),
      expand = c(0, 0)
    ) +
    labs(
      # title = paste0(
      #   "CMiNet score distribution of top-", top_pct,
      #   "% reliability edges  |  CMiNet: ", CMI_LABEL
      # ),
      # subtitle = paste0(
      #   "Fill colour: CMiNet score 0\u20139  ",
      #   "(0 = not confirmed by any method;  ",
      #   "9 = confirmed by all 9 methods).\n",
      #   "Column header: number of CMiNet edges with score \u2265 ", CMI_STRONG,
      #   " (n_B) and number of top-", top_pct, "% reliable edges (n_A).\n",
      #   "Above each bar: hypergeometric p-value  P(X \u2265 TP),  ",
      #   "E\u1d57\u1d52\u1d57\u1d43\u02e1 = p(p\u22121)/2."
      # ),
      x = "Normalization",
      y = paste0("% of top-", top_pct, "% reliable edges")
    ) +
    theme_bw(11) +
    theme(
      plot.title       = element_text(face = "bold", size = 11,
                                      margin = margin(b = 3)),
      plot.subtitle    = element_text(size = 7.8, colour = "grey25",
                                      lineheight = 1.35,
                                      margin = margin(b = 8)),
      strip.text.x     = element_text(face = "bold", size = 8.5,
                                      lineheight = 1.15),
      strip.text.y     = element_text(face = "bold", size = 9.5),
      strip.background = element_rect(fill = "grey94", colour = NA),
      panel.grid       = element_blank(),
      panel.spacing.x  = unit(0.6, "lines"),
      panel.spacing.y  = unit(0.8, "lines"),
      axis.text.x      = element_text(size = 10, face = "bold"),
      axis.text.y      = element_text(size = 9),
      axis.ticks.x     = element_blank(),
      legend.position  = "right",
      legend.key.size  = unit(0.40, "cm"),
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8.5),
      plot.margin      = margin(8, 8, 8, 8)
    )

  fig_base <- file.path(out_dir, paste0("Fig_score_dist_", pct_str, "pct"))
  ggsave(paste0(fig_base, ".pdf"), fig,
         width = 14, height = 7, device = cairo_pdf)
  ggsave(paste0(fig_base, ".png"), fig,
         width = 14, height = 7, dpi = 300)
  cat("Saved:", basename(fig_base), "\n")
}

cat("\n=== DONE ===\n")
cat("Figures saved to:", out_dir, "\n")
cat("Statistics table: overlap_statistics.csv\n")
