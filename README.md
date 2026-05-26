# PhyMapNet

**PhyMapNet** is an R package for phylogeny-guided Bayesian microbial network
inference. This repository contains the package source and the reproducible
analysis workflow for the revised PhyMapNet study, including filtering,
sensitivity analyses, bootstrap/noisy-data stability analyses, CMiNet overlap
comparisons, HMP biological evaluation, and runtime scripts.

The reviewed package version in this repository is **phymapnet 0.1.3**.

## Main Features

- Inference of microbial association networks using a phylogenetic tree or a
  named phylogenetic distance matrix.
- Single-model network estimation with `phymapnet_fit()`.
- Ensemble edge reliability estimation with `phymapnet_reliability()`.
- Supported normalizations: `"log"`, `"clr"`, and `"tss"`.
- Supported phylogenetic kernels: `"gaussian"` and `"laplacian"`.
- Taxon alignment and optional tree pruning before kernel construction.

## Installation

### Install the reviewed local source package

From the repository root:

```r
install.packages("phymapnet_0.1.3.tar.gz", repos = NULL, type = "source")
library(phymapnet)
packageVersion("phymapnet")
```

### Install from CRAN

After version 0.1.3 is accepted on CRAN:

```r
install.packages("phymapnet")
library(phymapnet)
```

### Install from GitHub

Because the R package source is stored in the `phymapnet/` subdirectory:

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/PhyMapNet", subdir = "phymapnet")
```

Replace `YOUR_GITHUB_USERNAME` with the final repository owner after upload.

## Package Dependencies

The R package imports:

```r
stats
ape
compositions
```

Paper-analysis scripts additionally require selected packages including
`dplyr`, `tibble`, `tidyr`, `ggplot2`, `patchwork`, `igraph`, `influential`,
`SpiecEasi`, and `CMiNet`.

## Input Data

### OTU/ASV table

`otu` is a numeric matrix or data frame:

- rows: samples;
- columns: taxa or ASVs;
- column names: taxon identifiers.

### Phylogenetic tree input

The second input may be an `ape::phylo` tree:

```r
fit <- phymapnet_fit(otu = otu, tree = tree)
```

When `prune_tree = TRUE`, nonshared OTU/tree taxa are removed and the OTU
columns and tree tip labels are aligned before phylogenetic distances are
computed.

### Phylogenetic distance-matrix input

The second input may alternatively be a named, symmetric distance matrix:

```r
fit <- phymapnet_fit(otu = otu, tree = DIS)
```

`DIS` must have matching row and column taxon names, nonnegative values, a
zero diagonal, and taxa compatible with the OTU columns. The public argument
name remains `tree`, but the function recognizes either a tree or distance
matrix.

## Single-Model Inference

```r
library(phymapnet)

fit <- phymapnet_fit(
  otu = otu,
  tree = tree,
  alpha = 0.05,
  k = 3,
  epsilon1 = 0.1,
  epsilon2 = 0.1,
  kernel = "gaussian",
  th_sparsity = 0.90,
  normalization = "log",
  prune_tree = TRUE
)

fit$adjacency
fit$precision
fit$edge_list
```

The same call may use `tree = DIS` when a named phylogenetic distance matrix
is already available.

## Reliability Ensemble Inference

```r
res <- phymapnet_reliability(
  otu = otu,
  tree = tree,
  th_fixed = 0.90,
  alpha_range = c(0.03, 0.05),
  k_range = 2:10,
  epsilon1_range = c(0, 1),
  epsilon2_range = c(0, 1),
  kernels = c("gaussian"),
  normalizations = c("log"),
  consensus_cut = 0.50,
  prune_tree = TRUE,
  progress = FALSE
)

res$rel_mat
res$consensus_mat
res$edge_list
```

For distance-matrix input, replace `tree = tree` with `tree = DIS`.

## Main Parameters

| Parameter | Meaning |
|---|---|
| `alpha` / `alpha_range` | Kernel bandwidth controlling decay with phylogenetic distance. |
| `k` / `k_range` | Prior degree-of-freedom multiplier; internally used as `k * p`, where `p` is the number of taxa. |
| `epsilon1` / `epsilon1_range` | Regularization term applied during prior covariance construction. |
| `epsilon2` / `epsilon2_range` | Additional numerical regularization in precision estimation. |
| `kernel` / `kernels` | `"gaussian"` or `"laplacian"` phylogenetic kernel. |
| `normalization` / `normalizations` | `"log"`, `"clr"`, or `"tss"` transformation. |
| `th_sparsity` | Quantile threshold for sparsifying a single-model precision estimate. |
| `th_fixed` | Fixed quantile threshold applied to each model in an ensemble. |
| `consensus_cut` | Minimum edge-selection reliability required in the binary consensus network. |
| `prune_tree` | If `TRUE`, retains and aligns taxa shared by the OTU table and phylogenetic tree. |

## Revised Paper Analysis Design

GMPR normalization is not used in the revised package analysis workflow because
it produced problematic values for some datasets. The revised package interface
includes `"log"`, `"clr"`, and `"tss"` only.

For final reliability networks used downstream, normalization and kernel are
fixed first, and reliability is merged over:

```text
alpha, k, epsilon1, epsilon2
```

For example, the HMP biological evaluation reads the completed HMP
`log`/`laplacian` reliability master result.

The sensitivity analysis is different by design: it evaluates how the network
changes with normalization, kernel, `alpha`, `k`, `epsilon1`, and `epsilon2`.
It should not be interpreted as one final consensus network merged across all
normalizations and kernels.

## Reproduce Paper Outputs

The GitHub repository retains the two precomputed result sets needed for the
reported downstream analyses: `result/reliable_score_all/` and
`result/reliability_master/`. Generated figures and derived tables are not
tracked because they can be regenerated locally from these retained results.

From a terminal, install the plotting dependencies and run the reproducibility
script:

```bash
Rscript -e 'install.packages(c("ape", "dplyr", "tibble", "tidyr", "ggplot2", "patchwork", "igraph", "influential"))'
bash scripts/reproduce_paper_outputs.sh
```

This command regenerates filtered OTU/distance inputs, sensitivity Figure 2
and Tables S3/S4, CMiNet overlap outputs, and the HMP biological outputs from
the included precomputed result files.

The script intentionally does **not** rerun multi-day network generation, the
complete bootstrap analysis, or full runtime benchmarks. The bootstrap script
is supplied for reproducibility, but its completed summary results are not
included in the GitHub upload.

For a complete file-by-file guide and commands for expensive optional
regeneration, see [README_FILES.md](README_FILES.md).

## Quick Package Tests

```bash
Rscript test_local_phymapnet_0.1.3.R
Rscript test_local_phymapnet_0.1.3_real_data.R
Rscript example_run_phymapnet_tree_and_distance.R
```

The real-data test uses the three filtered study datasets and checks agreement
between tree-based and distance-matrix-based package calls with small test
grids.

## CRAN Source Package

The CRAN submission archive prepared for version 0.1.3 is:

```text
phymapnet_0.1.3.tar.gz
```

The editable package source is in `phymapnet/`. For GitHub, the source
directory should be uploaded. The `.tar.gz` archive may also be attached to a
GitHub Release for convenience; it is not required in the tracked repository
because it can be rebuilt from source and is excluded by `.gitignore`.

## Citation

Please cite the PhyMapNet manuscript when using this package or its study
workflow. 

` @article{shahdoust2026phymapnet,
  title={PhyMapNet: A Phylogeny-Guided Bayesian Framework for Reliable Microbiome Network Inference},
  author={Shahdoust, Maryam and Aghdam, Rosa and Taheri, Golnaz},
  journal={bioRxiv},
  pages={2026--02},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}`

## Contact

Rosa Aghdam
rosaaghdam@gmail.com

## License

The R package is distributed under the GPL-3 license.
