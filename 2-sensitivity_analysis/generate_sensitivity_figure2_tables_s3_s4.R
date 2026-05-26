################################################################################
# Sensitivity analysis outputs for the revised paper
#
# Purpose:
#   Generate the sensitivity-analysis outputs used in the revised manuscript:
#     - Main Figure 2: pairwise within-parameter-group weighted Jaccard
#       similarity
#     - Supplementary Table S3: Kruskal-Wallis tests comparing weighted Jaccard
#       similarity across parameter groups
#     - Supplementary Table S4: post hoc pairwise Mann-Whitney tests comparing
#       weighted Jaccard similarity between parameter groups
#
# Input:
#   Existing precomputed reliability edge lists in result/reliable_score_all/.
#   This script does not rerun the time-consuming network generation step.
#
# Notes:
#   This combines the relevant logic from steps 2-7 without changing the
#   scientific calculations. It evaluates sensitivity to normalization, kernel,
#   alpha, k, epsilon1, and epsilon2 using the existing parameter-specific
#   reliability networks.
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ctx <- rstudioapi::getActiveDocumentContext()
    if (!is.null(ctx$path) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}

script_dir <- get_script_dir()
if (file.exists(file.path(getwd(), "result", "reliable_score_all"))) {
  repo_root <- normalizePath(getwd())
} else {
  repo_root <- normalizePath(file.path(script_dir, ".."))
}

input_dir <- file.path(repo_root, "result", "reliable_score_all")
output_dir <- file.path(repo_root, "result", "sensitivity_analysis")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figure")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- c("caffeine", "hmp_stool", "smoking")
expected_group_counts <- c(
  normalization = 3,
  kernel = 2,
  alpha = 10,
  epsilon1 = 10,
  k = 9,
  epsilon2 = 10
)

dataset_labels <- c(
  caffeine = "Caffeine",
  hmp_stool = "HMP stool",
  smoking = "Smoking"
)

group_order <- c("normalization", "kernel", "alpha", "epsilon1", "k", "epsilon2")
group_labels <- c(
  normalization = "norm",
  kernel = "kernel",
  alpha = "alpha",
  epsilon1 = "epsilon1",
  k = "k",
  epsilon2 = "epsilon2"
)

figure_group_labels <- c(
  normalization = "norm",
  kernel = "kernel",
  alpha = "alpha",
  epsilon1 = "eps1",
  k = "k",
  epsilon2 = "eps2"
)

sensitivity_palette <- c(
  "Very Robust" = "#2E7D32",
  "Robust" = "#1E88E5",
  "Moderate" = "#FBC02D",
  "Sensitive" = "#E53935"
)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

extract_info <- function(fname) {
  if (grepl("^norm_", fname, ignore.case = TRUE)) {
    list(group = "normalization", parameter = sub("^norm_", "", fname, ignore.case = TRUE))
  } else if (grepl("^alpha_", fname, ignore.case = TRUE)) {
    list(group = "alpha", parameter = sub("^alpha_", "", fname, ignore.case = TRUE))
  } else if (grepl("^k_", fname, ignore.case = TRUE)) {
    list(group = "k", parameter = sub("^k_", "", fname, ignore.case = TRUE))
  } else if (grepl("^eps1_|^epsilon1_", fname, ignore.case = TRUE)) {
    list(group = "epsilon1", parameter = sub("^eps1_|^epsilon1_", "", fname, ignore.case = TRUE))
  } else if (grepl("^eps2_|^epsilon2_", fname, ignore.case = TRUE)) {
    list(group = "epsilon2", parameter = sub("^eps2_|^epsilon2_", "", fname, ignore.case = TRUE))
  } else if (grepl("^kernel_", fname, ignore.case = TRUE)) {
    list(group = "kernel", parameter = sub("^kernel_", "", fname, ignore.case = TRUE))
  } else if (grepl("^baseline", fname, ignore.case = TRUE)) {
    list(group = "baseline", parameter = "baseline")
  } else {
    list(group = "other", parameter = fname)
  }
}

make_edge_id <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  ifelse(a < b, paste(a, b, sep = "__"), paste(b, a, sep = "__"))
}

read_network <- function(file) {
  df <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1:3] <- c("taxa1", "taxa2", "weight")
  df %>%
    mutate(
      taxa1 = as.character(taxa1),
      taxa2 = as.character(taxa2),
      weight = suppressWarnings(as.numeric(weight)),
      edge = make_edge_id(taxa1, taxa2)
    ) %>%
    filter(!is.na(weight)) %>%
    distinct(edge, .keep_all = TRUE) %>%
    select(taxa1, taxa2, edge, weight)
}

weighted_jaccard <- function(x, y) {
  all_edges <- union(names(x), names(y))
  x_full <- x[all_edges]
  y_full <- y[all_edges]
  x_full[is.na(x_full)] <- 0
  y_full[is.na(y_full)] <- 0

  denominator <- sum(pmax(x_full, y_full))
  if (denominator == 0) return(NA_real_)
  sum(pmin(x_full, y_full)) / denominator
}

sig_code <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

format_p <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}

# -------------------------------------------------------------------------
# Load networks and compute similarity tables
# -------------------------------------------------------------------------

all_baseline <- list()
all_pairwise <- list()

for (ds in datasets) {
  message("Processing dataset: ", ds)
  ds_dir <- file.path(input_dir, ds)
  files <- list.files(ds_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No CSV files found in ", ds_dir, call. = FALSE)

  networks <- list()
  network_info <- list()

  for (f in files) {
    fname <- tools::file_path_sans_ext(basename(f))
    info <- extract_info(fname)
    if (identical(info$group, "other")) next

    net <- read_network(f)
    networks[[fname]] <- setNames(net$weight, net$edge)
    network_info[[fname]] <- info
  }

  baseline_file <- file.path(ds_dir, "baseline.csv")
  if (!file.exists(baseline_file)) stop("Missing baseline.csv for ", ds, call. = FALSE)

  baseline_net <- read_network(baseline_file)
  baseline_weights <- setNames(baseline_net$weight, baseline_net$edge)

  param_names <- setdiff(names(networks), "baseline")
  baseline_rows <- list()

  for (fname in param_names) {
    info <- network_info[[fname]]
    baseline_rows[[fname]] <- data.frame(
      dataset = ds,
      group = info$group,
      parameter = info$parameter,
      weighted_jaccard = weighted_jaccard(networks[[fname]], baseline_weights),
      n_edges = length(networks[[fname]]),
      stringsAsFactors = FALSE
    )
  }

  baseline_df <- bind_rows(baseline_rows)
  observed_group_counts <- table(baseline_df$group)
  observed_group_counts <- observed_group_counts[names(expected_group_counts)]
  observed_group_counts[is.na(observed_group_counts)] <- 0L
  if (!identical(as.integer(observed_group_counts), as.integer(expected_group_counts))) {
    stop(
      "Incomplete sensitivity network set for ", ds, ". Expected group counts: ",
      paste(names(expected_group_counts), expected_group_counts, sep = "=", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  all_baseline[[ds]] <- baseline_df

  pairwise_rows <- list()
  for (grp in unique(baseline_df$group)) {
    grp_names <- param_names[vapply(param_names, function(nm) {
      identical(network_info[[nm]]$group, grp)
    }, logical(1))]

    if (length(grp_names) < 2) next

    for (i in seq_len(length(grp_names) - 1)) {
      for (j in (i + 1):length(grp_names)) {
        name1 <- grp_names[i]
        name2 <- grp_names[j]
        pairwise_rows[[paste(ds, grp, name1, name2, sep = "__")]] <- data.frame(
          dataset = ds,
          group = grp,
          parameter1 = network_info[[name1]]$parameter,
          parameter2 = network_info[[name2]]$parameter,
          pairwise_wj = weighted_jaccard(networks[[name1]], networks[[name2]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  all_pairwise[[ds]] <- bind_rows(pairwise_rows)
}

combined_baseline <- bind_rows(all_baseline)
combined_pairwise <- bind_rows(all_pairwise)

write.csv(combined_baseline,
          file.path(table_dir, "combined_parameter_similarity.csv"),
          row.names = FALSE)
write.csv(combined_pairwise,
          file.path(table_dir, "combined_pairwise_similarity.csv"),
          row.names = FALSE)

# -------------------------------------------------------------------------
# Table S3: Kruskal-Wallis tests
# -------------------------------------------------------------------------

kw_results <- list()

for (ds in datasets) {
  ds_data <- combined_baseline %>% filter(dataset == ds, group %in% group_order)
  groups <- unique(ds_data$group)
  group_data <- lapply(groups, function(g) {
    ds_data %>% filter(group == g) %>% pull(weighted_jaccard)
  })

  kw_test <- kruskal.test(group_data)

  kw_results[[ds]] <- data.frame(
    dataset = ds,
    n = nrow(ds_data),
    n_groups = length(groups),
    H_statistic = as.numeric(kw_test$statistic),
    p_value = kw_test$p.value,
    significance = sig_code(kw_test$p.value),
    stringsAsFactors = FALSE
  )
}

table_s3 <- bind_rows(kw_results) %>%
  mutate(
    Dataset = dataset_labels[dataset],
    `Kruskal-Wallis H` = sprintf("%.2f", H_statistic),
    `p-value` = format_p(p_value),
    Significance = significance
  ) %>%
  select(
    Dataset,
    n,
    `No. groups` = n_groups,
    `Kruskal-Wallis H`,
    `p-value`,
    Significance
  )

write.csv(table_s3, file.path(table_dir, "Table_S3_kruskal_wallis.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Table S4: post hoc pairwise Mann-Whitney tests
# -------------------------------------------------------------------------

pairwise_results <- list()

for (ds in datasets) {
  ds_data <- combined_baseline %>% filter(dataset == ds, group %in% group_order)
  groups <- sort(unique(ds_data$group))

  for (i in seq_len(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]
      g2 <- groups[j]
      data_g1 <- ds_data %>% filter(group == g1) %>% pull(weighted_jaccard)
      data_g2 <- ds_data %>% filter(group == g2) %>% pull(weighted_jaccard)
      mw_test <- wilcox.test(data_g1, data_g2)

      pairwise_results[[length(pairwise_results) + 1]] <- data.frame(
        dataset = ds,
        group1 = g1,
        group2 = g2,
        median_wj_g1 = median(data_g1, na.rm = TRUE),
        median_wj_g2 = median(data_g2, na.rm = TRUE),
        U_statistic = as.numeric(mw_test$statistic),
        p_value = mw_test$p.value,
        n_g1 = length(data_g1),
        n_g2 = length(data_g2),
        stringsAsFactors = FALSE
      )
    }
  }
}

table_s4_raw <- bind_rows(pairwise_results) %>%
  mutate(
    p_value_BH = p.adjust(p_value, method = "BH"),
    significance = sig_code(p_value_BH)
  )

comparison_label <- function(g1, g2) {
  paste(group_labels[g1], "vs.", group_labels[g2])
}

table_s4 <- table_s4_raw %>%
  mutate(
    Dataset = dataset_labels[dataset],
    Comparison = comparison_label(group1, group2),
    `Median WJ 1` = sprintf("%.2f", median_wj_g1),
    `Median WJ 2` = sprintf("%.2f", median_wj_g2),
    `U statistic` = U_statistic,
    `BH-adjusted p-value` = format_p(p_value_BH),
    Significance = significance
  ) %>%
  select(
    Dataset,
    Comparison,
    `Median WJ 1`,
    `Median WJ 2`,
    `U statistic`,
    `BH-adjusted p-value`,
    Significance
  )

write.csv(table_s4_raw, file.path(table_dir, "Table_S4_pairwise_mann_whitney_raw.csv"), row.names = FALSE)
write.csv(table_s4, file.path(table_dir, "Table_S4_pairwise_mann_whitney.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Figure 2: pairwise within-parameter-group sensitivity
# -------------------------------------------------------------------------

plot_pairwise <- combined_pairwise %>%
  filter(group %in% group_order) %>%
  mutate(
    group = factor(group, levels = group_order, labels = figure_group_labels[group_order]),
    dataset = factor(dataset, levels = names(dataset_labels), labels = dataset_labels),
    sensitivity_class = case_when(
      pairwise_wj >= 0.85 ~ "Very Robust",
      pairwise_wj >= 0.70 ~ "Robust",
      pairwise_wj >= 0.55 ~ "Moderate",
      TRUE ~ "Sensitive"
    ),
    sensitivity_class = factor(
      sensitivity_class,
      levels = c("Very Robust", "Robust", "Moderate", "Sensitive")
    )
  )

figure2 <- ggplot(plot_pairwise, aes(x = group, y = pairwise_wj)) +
  geom_vline(
    xintercept = seq(1.5, length(group_order) - 0.5, by = 1),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_boxplot(
    fill = "white",
    alpha = 0.6,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = sensitivity_class),
    width = 0.20,
    size = 2.2,
    alpha = 0.9,
    shape = 19
  ) +
  facet_wrap(~ dataset, ncol = 3, scales = "fixed") +
  scale_color_manual(values = sensitivity_palette) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(
    x = "Parameter Group",
    y = "Pairwise Weighted Jaccard",
    color = "Sensitivity Class"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(file.path(figure_dir, "Figure2_pairwise_within_group.pdf"),
       figure2, width = 14, height = 10, dpi = 600)
ggsave(file.path(figure_dir, "Figure2_pairwise_within_group.png"),
       figure2, width = 14, height = 10, dpi = 300)

message("Done.")
message("Tables:  ", table_dir)
message("Figures: ", figure_dir)
