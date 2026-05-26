# 1. Install and load phymapnet version 0.1.3
install.packages("phymapnet_0.1.3.tar.gz", repos = NULL, type = "source")
library(phymapnet)
library(ape)

# 2. Read the OTU table
otu <- as.matrix(
  read.csv("data/filter/caffeine_otu_filtered.csv",
           row.names = 1, check.names = FALSE)
)
colnames(otu) <- gsub("^X", "", colnames(otu))

# 3. Read the phylogenetic tree
load("data/original_data/data_caff.RData")
tree <- data_caff$tree

# Align the tree and OTU table
common <- intersect(colnames(otu), tree$tip.label)
tree <- ape::keep.tip(tree, common)
otu <- otu[, tree$tip.label, drop = FALSE]

# 4. Find the phylogenetic distance matrix from the tree
DIS <- as.matrix(ape::cophenetic.phylo(tree))

# 5. Reliability inference using the tree and parameter ranges
res_tree <- phymapnet_reliability(
  otu = otu,
  tree = tree,
  th_fixed = 0.95,
  alpha_range = c(0.05, 0.10),
  k_range = c(2, 3),
  epsilon1_range = c(0.1, 0.2),
  epsilon2_range = c(0.1, 0.2),
  kernels = c("gaussian"),
  normalizations = c("log"),
  consensus_cut = 0.50,
  prune_tree = TRUE,
  progress = FALSE
)

# 6. Reliability inference using the distance matrix and parameter ranges
res_DIS <- phymapnet_reliability(
  otu = otu,
  tree = DIS,
  th_fixed = 0.95,
  alpha_range = c(0.05, 0.10),
  k_range = c(2, 3),
  epsilon1_range = c(0.1, 0.2),
  epsilon2_range = c(0.1, 0.2),
  kernels = c("gaussian"),
  normalizations = c("log"),
  consensus_cut = 0.50,
  progress = FALSE
)

# 7. Single-model inference using the tree and one parameter setting
fit_tree <- phymapnet_fit(
  otu = otu,
  tree = tree,
  alpha = 0.05,
  k = 3,
  epsilon1 = 0.1,
  epsilon2 = 0.1,
  kernel = "gaussian",
  th_sparsity = 0.95,
  normalization = "log",
  prune_tree = TRUE
)

# 8. Single-model inference using the distance matrix and one parameter setting
fit_DIS <- phymapnet_fit(
  otu = otu,
  tree = DIS,
  alpha = 0.05,
  k = 3,
  epsilon1 = 0.1,
  epsilon2 = 0.1,
  kernel = "gaussian",
  th_sparsity = 0.95,
  normalization = "log"
)

# Simple checks: tree and distance inputs should agree for the same settings
stopifnot(isTRUE(all.equal(res_tree$rel_mat, res_DIS$rel_mat)))
stopifnot(identical(fit_tree$adjacency, fit_DIS$adjacency))

cat("Reliability models fitted:", res_tree$N_models, "\n")
cat("Number of taxa:", ncol(otu), "\n")
cat("Tree and distance-input results agree.\n")
