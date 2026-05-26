# Generate the bootstrap/noisy-data stability comparison figure from completed
# result files. This script does not rerun network inference.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

base_dir <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
dir_out <- Sys.getenv(
  "PHYMAPNET_BOOTSTRAP_OUTPUT_DIR",
  unset = file.path(base_dir, "result", "bootstrap")
)
results_file <- file.path(dir_out, "all_datasets_all_methods_results.csv")
if (!file.exists(results_file)) {
  stop(
    "Bootstrap result file not found: ", results_file,
    "\nThis figure requires completed F1-score results from the bootstrap analysis.",
    "\nThe repository currently contains replicate input data, not this completed summary.",
    "\nFor a short code check, run the smoke-test commands in the root README.",
    "\nFor the paper figure, run the full analysis or add your existing completed summary file.",
    call. = FALSE
  )
}
combined <- read.csv(results_file, stringsAsFactors = FALSE)

dataset_labels <- c(smoking = "Smoking", caffeine = "Caffeine", hmp_stool = "HMP stool")
method_order <- c(
  "PhyMapNet (phylo)", "PhyMapNet (uniform)", "SPIEC-EASI mb",
  "SPIEC-EASI glasso", "SparCC", "ccLasso", "c_MI"
)
method_colors <- c(
  "PhyMapNet (phylo)" = "#1B5E20",
  "PhyMapNet (uniform)" = "#81C784",
  "SPIEC-EASI mb" = "#1565C0",
  "SPIEC-EASI glasso" = "#42A5F5",
  "SparCC" = "#6A1B9A",
  "ccLasso" = "#E65100",
  "c_MI" = "#D32F2F"
)

plot_df <- combined %>%
  filter(split %in% c("bootstrap", "noisy"), !is.na(F1)) %>%
  mutate(
    dataset = factor(dataset_labels[dataset], levels = dataset_labels),
    method = factor(method, levels = method_order),
    split = factor(split, levels = c("bootstrap", "noisy"), labels = c("Bootstrap", "Noisy"))
  )
if (!nrow(plot_df)) {
  stop("No bootstrap or noisy replicate results available for plotting.", call. = FALSE)
}

figure <- ggplot(plot_df, aes(x = method, y = F1, fill = method)) +
  geom_boxplot(
    outlier.shape = 16, outlier.size = 1.0, outlier.alpha = 0.4,
    width = 0.72, linewidth = 0.35, colour = "grey20"
  ) +
  facet_grid(dataset ~ split) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_discrete(labels = c(
    "PhyMapNet (phylo)" = "PhyMapNet",
    "PhyMapNet (uniform)" = "PhyMapNet (uniform)",
    "SPIEC-EASI mb" = "SPIEC-EASI_mb",
    "SPIEC-EASI glasso" = "SPIEC-EASI_glasso",
    "SparCC" = "SparCC",
    "ccLasso" = "CClasso",
    "c_MI" = "CMIMN"
  )) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = c("0", "0.25", "0.50", "0.75", "1.0")
  ) +
  labs(x = NULL, y = "F1-score", title = "", subtitle = "") +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 3)),
    plot.subtitle = element_text(
      size = 7.5, colour = "grey30", lineheight = 1.3, margin = margin(b = 6)
    ),
    strip.background = element_rect(fill = "grey94", colour = NA),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 9, angle = 0),
    axis.text.x = element_text(size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(dir_out, "Fig_bootstrap_boxplot.pdf"), figure, width = 13, height = 10)
ggsave(file.path(dir_out, "Fig_bootstrap_boxplot.png"), figure, width = 13, height = 10, dpi = 300)
cat("Saved bootstrap stability figure files to:", dir_out, "\n")
