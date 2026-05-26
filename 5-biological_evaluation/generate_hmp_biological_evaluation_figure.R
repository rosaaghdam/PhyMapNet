# ============================================================
# HMP biological-evaluation network figure from completed master results
# Paper-ready HMP filtered 312-OTU network visualization
# Isolated nodes omitted from visualization: degree > 0
# Gray edges: reliability > 0.70
# Bold black edges: reliability > 0.9999999 and cminet_096 > 8
# Node color: FAMILY or GENUS
# Node shape: Hub / IVI / Hub + IVI
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(igraph)
  library(ggplot2)
  library(influential)
})

# ============================================================
# Parameters
# ============================================================

expected_n_taxa <- 312

taxonomic_rank_to_color <- "GENUS"   # Supported choices: PHYLUM, CLASS, ORDER, FAMILY, GENUS
top_taxonomic_groups_to_color <- 20

reliability_background_cut <- 0.70
reliability_strong_cut <- 0.9999999
cminet_strong_cut <- 8
cminet_column <- "cminet_096"

top_n <- 20

# Core expansion parameters
expand_core <- TRUE
core_top_degree_n <- 20
expand_radius <- 1.8
expand_strength <- 3

base_dir <- normalizePath(
  Sys.getenv("PHYMAPNET_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
master_file <- file.path(
  base_dir, "result", "reliability_master", "phymap",
  "hmp_stool_log_laplacian_master.csv"
)
otu_file <- file.path(base_dir, "data", "filter", "hmp_stool_otu_filtered.csv")
taxonomy_file <- file.path(base_dir, "data", "filter", "hmp_stool_nonisolated_taxonomy.csv")

outdir <- file.path(base_dir, "result", "biological_evaluation_filtered312")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: Read filtered OTU table
# ============================================================

otu_df <- read.csv(
  otu_file,
  row.names = 1,
  check.names = FALSE
)

colnames(otu_df) <- gsub("^X", "", colnames(otu_df))
otu <- as.matrix(otu_df)
all_taxa <- colnames(otu)

cat("Filtered OTU table:", nrow(otu), "samples x", length(all_taxa), "taxa\n")

if (length(all_taxa) != expected_n_taxa) {
  stop(
    "Expected ", expected_n_taxa, " filtered taxa, but found ",
    length(all_taxa),
    ". Please check the filtered OTU file."
  )
}

# ============================================================
# STEP 2: Read master edge file
# ============================================================

master <- read.csv(
  master_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("taxa1", "taxa2", "reliability", cminet_column)
if (!all(required_columns %in% colnames(master))) {
  stop(
    "Required columns missing from master result: ",
    paste(setdiff(required_columns, colnames(master)), collapse = ", "),
    call. = FALSE
  )
}

edge_all <- master[, required_columns]
colnames(edge_all) <- c("taxon1", "taxon2", "reliability", "cminet")

edge_all$reliability <- as.numeric(edge_all$reliability)
edge_all$cminet <- as.numeric(edge_all$cminet)

edge_background <- edge_all %>%
  dplyr::filter(
    reliability > reliability_background_cut,
    taxon1 %in% all_taxa,
    taxon2 %in% all_taxa
  )

edge_strong <- edge_all %>%
  dplyr::filter(
    reliability > reliability_strong_cut,
    cminet > cminet_strong_cut,
    taxon1 %in% all_taxa,
    taxon2 %in% all_taxa
  )

cat("Background edges reliability >", reliability_background_cut, ":", nrow(edge_background), "\n")
cat(
  "Strong edges reliability >", reliability_strong_cut,
  "and", cminet_column, ">", cminet_strong_cut, ":", nrow(edge_strong), "\n"
)

if (nrow(edge_background) == 0) {
  stop("No background edges passed reliability threshold.")
}

# ============================================================
# STEP 3: Build graphs
# ============================================================

g_all <- igraph::graph_from_data_frame(
  d = edge_background,
  directed = FALSE,
  vertices = data.frame(name = all_taxa)
)

E(g_all)$weight <- E(g_all)$reliability
E(g_all)$cminet <- E(g_all)$cminet

# Keep only non-isolated nodes for visualization and centrality
nonisolated_taxa <- V(g_all)$name[igraph::degree(g_all) > 0]

g <- igraph::induced_subgraph(
  g_all,
  vids = nonisolated_taxa
)

edge_background_plot <- edge_background %>%
  dplyr::filter(
    taxon1 %in% nonisolated_taxa,
    taxon2 %in% nonisolated_taxa
  )

edge_strong_plot <- edge_strong %>%
  dplyr::filter(
    taxon1 %in% nonisolated_taxa,
    taxon2 %in% nonisolated_taxa
  )

g_strong <- igraph::graph_from_data_frame(
  d = edge_strong_plot,
  directed = FALSE,
  vertices = data.frame(name = nonisolated_taxa)
)

E(g_strong)$weight <- E(g_strong)$reliability
E(g_strong)$cminet <- E(g_strong)$cminet

cat("All graph nodes:", igraph::vcount(g_all), "\n")
cat("Non-isolated plotted nodes:", igraph::vcount(g), "\n")
cat("Background plotted edges:", igraph::ecount(g), "\n")
cat("Strong plotted edges:", igraph::ecount(g_strong), "\n")

# ============================================================
# STEP 4: HMP taxonomy annotations used in the final saved analysis
# ============================================================

if (!file.exists(taxonomy_file)) {
  stop("Missing local HMP taxonomy annotation file: ", taxonomy_file, call. = FALSE)
}
tax_tab <- read.csv(taxonomy_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!setequal(tax_tab$OTU_ID, nonisolated_taxa)) {
  stop(
    "Local HMP taxonomy annotations do not match non-isolated network taxa.",
    call. = FALSE
  )
}

tax_tab <- tax_tab %>%
  dplyr::mutate(
    SUPERKINGDOM = ifelse(is.na(SUPERKINGDOM) | SUPERKINGDOM == "", "Unclassified", SUPERKINGDOM),
    PHYLUM = ifelse(is.na(PHYLUM) | PHYLUM == "", "Unclassified", PHYLUM),
    CLASS = ifelse(is.na(CLASS) | CLASS == "", "Unclassified", CLASS),
    ORDER = ifelse(is.na(ORDER) | ORDER == "", "Unclassified", ORDER),
    FAMILY = ifelse(is.na(FAMILY) | FAMILY == "", "Unclassified", FAMILY),
    GENUS = ifelse(is.na(GENUS) | GENUS == "", "Unclassified", GENUS),
    biological_label = dplyr::case_when(
      GENUS != "Unclassified" & FAMILY != "Unclassified" ~ paste0(FAMILY, " / ", GENUS),
      GENUS != "Unclassified" ~ GENUS,
      FAMILY != "Unclassified" ~ FAMILY,
      ORDER != "Unclassified" ~ ORDER,
      CLASS != "Unclassified" ~ CLASS,
      PHYLUM != "Unclassified" ~ PHYLUM,
      TRUE ~ OTU_ID
    )
  )

if (!taxonomic_rank_to_color %in% colnames(tax_tab)) {
  stop("taxonomic_rank_to_color must be one of: PHYLUM, CLASS, ORDER, FAMILY, GENUS")
}

# ============================================================
# STEP 5: Centrality and IVI using strong plotted network
# ============================================================

centrality_df <- tibble::tibble(
  OTU_ID = V(g_strong)$name,
  degree = igraph::degree(g_strong),
  strength = igraph::strength(g_strong, weights = E(g_strong)$weight),
  betweenness = igraph::betweenness(g_strong, directed = FALSE, weights = NULL),
  closeness = igraph::closeness(g_strong, normalized = TRUE),
  eigenvector = igraph::eigen_centrality(
    g_strong,
    directed = FALSE,
    weights = E(g_strong)$weight
  )$vector,
  pagerank = igraph::page_rank(
    g_strong,
    directed = FALSE,
    weights = E(g_strong)$weight
  )$vector
)

adj_strong <- igraph::as_adjacency_matrix(
  g_strong,
  attr = "weight",
  sparse = FALSE
)

g_ivi <- igraph::graph_from_adjacency_matrix(
  adj_strong,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

ivi_values <- influential::ivi(
  graph = g_ivi,
  weights = NULL,
  directed = FALSE,
  mode = "all",
  loops = FALSE,
  d = 3,
  scale = "range",
  ncores = 1
)

ivi_df <- tibble::tibble(
  OTU_ID = names(ivi_values),
  IVI = as.numeric(ivi_values)
)

abundance_df <- tibble::tibble(
  OTU_ID = colnames(otu),
  prevalence_n = colSums(otu > 0),
  prevalence_fraction = colSums(otu > 0) / nrow(otu),
  total_abundance = colSums(otu),
  mean_abundance = colMeans(otu),
  median_abundance = apply(otu, 2, median)
)

node_table <- centrality_df %>%
  dplyr::left_join(ivi_df, by = "OTU_ID") %>%
  dplyr::left_join(abundance_df, by = "OTU_ID") %>%
  dplyr::left_join(
    tax_tab %>%
      dplyr::select(
        OTU_ID,
        CONSENSUS_LINEAGE,
        SUPERKINGDOM,
        PHYLUM,
        CLASS,
        ORDER,
        FAMILY,
        GENUS,
        biological_label
      ),
    by = "OTU_ID"
  )

top_hub_ids <- node_table %>%
  dplyr::arrange(dplyr::desc(degree), dplyr::desc(betweenness), dplyr::desc(eigenvector)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(OTU_ID)

top_ivi_ids <- node_table %>%
  dplyr::arrange(dplyr::desc(IVI)) %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::pull(OTU_ID)

node_table <- node_table %>%
  dplyr::mutate(
    important_status = dplyr::case_when(
      OTU_ID %in% top_hub_ids & OTU_ID %in% top_ivi_ids ~ "Hub + IVI",
      OTU_ID %in% top_hub_ids ~ "Hub",
      OTU_ID %in% top_ivi_ids ~ "IVI",
      TRUE ~ "Other"
    )
  )

important_taxa <- node_table %>%
  dplyr::filter(important_status != "Other") %>%
  dplyr::arrange(dplyr::desc(degree), dplyr::desc(IVI))

# ============================================================
# STEP 6: Save CSV outputs
# ============================================================

write.csv(
  node_table,
  file.path(outdir, "hmp_filtered312_nonisolated_node_metrics_taxonomy.csv"),
  row.names = FALSE
)

write.csv(
  important_taxa,
  file.path(outdir, "hmp_filtered312_nonisolated_selected_important_otus.csv"),
  row.names = FALSE
)

write.csv(
  edge_background_plot,
  file.path(outdir, paste0("hmp_filtered312_nonisolated_background_edges_reliability_gt_", reliability_background_cut, ".csv")),
  row.names = FALSE
)

write.csv(
  edge_strong_plot,
  file.path(outdir, "hmp_filtered312_nonisolated_strong_edges_reliability_gt_0p9999999_cminet096_gt_8.csv"),
  row.names = FALSE
)

# ============================================================
# STEP 7: Prepare plotting data
# ============================================================

node_table$color_group <- node_table[[taxonomic_rank_to_color]]

node_table$color_group <- ifelse(
  is.na(node_table$color_group) | node_table$color_group == "",
  "Unclassified",
  node_table$color_group
)

group_counts <- sort(table(node_table$color_group), decreasing = TRUE)
top_groups <- names(group_counts)[seq_len(min(top_taxonomic_groups_to_color, length(group_counts)))]

node_table$color_group_plot <- ifelse(
  node_table$color_group %in% top_groups,
  node_table$color_group,
  "Other"
)

node_table$important_status <- factor(
  node_table$important_status,
  levels = c("Other", "Hub", "IVI", "Hub + IVI")
)

set.seed(123)

layout_mat <- igraph::layout_with_graphopt(
  g,
  niter = 3500,
  charge = 0.012,
  mass = 45,
  spring.length = 0.12,
  spring.constant = 1.2
)

rescale_robust <- function(x) {
  q <- quantile(x, probs = c(0.02, 0.98), na.rm = TRUE)
  x <- pmax(pmin(x, q[2]), q[1])
  as.numeric(scale(x))
}

layout_mat[, 1] <- rescale_robust(layout_mat[, 1])
layout_mat[, 2] <- rescale_robust(layout_mat[, 2])

node_df <- data.frame(
  OTU_ID = V(g)$name,
  x = layout_mat[, 1],
  y = layout_mat[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(node_table, by = "OTU_ID")

# ============================================================
# STEP 7B: Expand dense central region
# ============================================================

if (expand_core) {
  core_nodes <- node_df %>%
    dplyr::arrange(dplyr::desc(degree)) %>%
    dplyr::slice_head(n = core_top_degree_n)

  core_center_x <- mean(core_nodes$x, na.rm = TRUE)
  core_center_y <- mean(core_nodes$y, na.rm = TRUE)

  dx <- node_df$x - core_center_x
  dy <- node_df$y - core_center_y
  dist <- sqrt(dx^2 + dy^2)

  effect <- ifelse(
    dist < expand_radius,
    1 + expand_strength * (1 - dist / expand_radius)^2,
    1
  )

  node_df$x <- core_center_x + dx * effect
  node_df$y <- core_center_y + dy * effect
}

x_pad <- 0.03 * diff(range(node_df$x, na.rm = TRUE))
y_pad <- 0.03 * diff(range(node_df$y, na.rm = TRUE))

x_lim <- range(node_df$x, na.rm = TRUE) + c(-x_pad, x_pad)
y_lim <- range(node_df$y, na.rm = TRUE) + c(-y_pad, y_pad)

edge_bg_raw <- igraph::as_data_frame(g, what = "edges")

edge_bg_df <- edge_bg_raw %>%
  dplyr::left_join(
    node_df %>% dplyr::select(from = OTU_ID, x_from = x, y_from = y),
    by = "from"
  ) %>%
  dplyr::left_join(
    node_df %>% dplyr::select(to = OTU_ID, x_to = x, y_to = y),
    by = "to"
  )

edge_strong_df <- edge_strong_plot %>%
  dplyr::left_join(
    node_df %>% dplyr::select(taxon1 = OTU_ID, x_from = x, y_from = y),
    by = "taxon1"
  ) %>%
  dplyr::left_join(
    node_df %>% dplyr::select(taxon2 = OTU_ID, x_to = x, y_to = y),
    by = "taxon2"
  )

important_df <- node_df %>%
  dplyr::filter(important_status != "Other")

# ============================================================
# STEP 8: Plot network
# ============================================================

p_network <- ggplot() +

  geom_segment(
    data = edge_bg_df,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    color = "grey78",
    alpha = 0.30,
    linewidth = 0.18
  ) +

  geom_segment(
    data = edge_strong_df,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    color = "black",
    alpha = 0.90,
    linewidth = 1.05
  ) +

  geom_point(
    data = node_df,
    aes(x = x, y = y, color = color_group_plot),
    shape = 16,
    size = 1.5,
    alpha = 0.85
  ) +

  geom_point(
    data = important_df,
    aes(x = x, y = y, color = color_group_plot, shape = important_status),
    size = 2.5,
    stroke = 1.1,
    alpha = 1
  ) +

  scale_shape_manual(
    values = c(
      "Hub" = 18,
      "IVI" = 15,
      "Hub + IVI" = 17
    ),
    drop = TRUE
  ) +

  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.spacing.y = unit(0.9, "mm"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10.5, hjust = 0),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  labs(
    title = "",
    subtitle = "",
    color = taxonomic_rank_to_color,
    shape = "Important taxa"
  )

fig_base <- paste0(
  "fig_hmp_filtered312_nonisolated_compact_expanded_core_",
  taxonomic_rank_to_color,
  "_bold_reliability0p9999999_cminet096_gt8"
)

ggsave(
  file.path(outdir, paste0(fig_base, "_tight.pdf")),
  p_network,
  width = 7,
  height = 5,
  limitsize = FALSE
)

ggsave(
  file.path(outdir, paste0(fig_base, "_tight.png")),
  p_network,
  width = 7,
  height = 5,
  dpi = 400,
  limitsize = FALSE
)

cat("\nSaved files:\n")
cat(file.path(outdir, paste0(fig_base, "_tight.pdf")), "\n")
cat(file.path(outdir, paste0(fig_base, "_tight.png")), "\n")
cat(file.path(outdir, "hmp_filtered312_nonisolated_selected_important_otus.csv"), "\n")
