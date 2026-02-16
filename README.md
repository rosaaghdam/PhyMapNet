## PhyMapNet: A Phylogeny-Guided Bayesian Framework for Reliable Microbiome Network Inference
<img src="logo/logo.png" style="width:50%;" align=right>

[![GitHub license](https://img.shields.io/github/license/solislemuslab/CMiNet?color=yellow)](https://github.com/solislemuslab/CMiNet/blob/main/LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/solislemuslab/CMiNet)](https://github.com/solislemuslab/CMiNet/issues)
![Code Size](https://img.shields.io/github/languages/code-size/solislemuslab/CMiNet?color=white)
[![GitHub Releases](https://img.shields.io/github/v/release/solislemuslab/CMiNet?display_name=tag)](https://github.com/solislemuslab/CMiNet/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17459775.svg)](https://doi.org/10.5281/zenodo.17459775)

## Description
<div align="justify">
  
**PhyMapNet** is a Bayesian Gaussian Graphical Model framework for inferring microbiome interaction networks that explicitly integrates phylogenetic information. By incorporating evolutionary distances through kernel-based priors, PhyMapNet embeds biological structure directly into precision matrix estimation, enabling robust inference of conditional dependencies among microbial taxa. To reduce sensitivity to hyperparameter choices, PhyMapNet constructs reliability-driven consensus networks by aggregating results across a large hyperparameter ensemble, yielding stable and interpretable microbiome networks with controllable sparsity. The method is computationally efficient and supports multiple normalization strategies, making it suitable for real-world, high-dimensional microbiome datasets.
</div>

**Figure overview**. PhyMapNet takes normalized microbiome count data and phylogenetic relationships as input, embeds evolutionary structure into a Bayesian graphical model, and infers sparse microbial interaction networks. To ensure robustness, the procedure is repeated across multiple hyperparameter settings, and a stable consensus network is selected based on predefined criteria.

## Installation
```bash
# install devtools if needed
install.packages("devtools")

# install phymapnet from GitHub
devtools::install_github("rosaaghdam/PhyMapNet")
```
If required packages are missing, install them manually:
```bash
# Dependencies
install.packages(c("ape", "GMPR", "compositions"))
```

## Running `phymapnet`

The package provides two main analysis functions for phylogeny-aware microbial network inference.

---

### 1. `phymapnet_fit()` — Single-Model Inference

Fits a single PhyMapNet model using specified hyperparameters to estimate a sparse microbial association network.

#### Inputs

- `otu`: Samples × taxa abundance matrix  
  (rows = samples, columns = taxa)
- `tree`: Phylogenetic tree (`ape::phylo`)  
  with `tip.label` matching OTU column names

#### Outputs

- `precision_map`: Estimated precision matrix  
- `adjacency`: Binary network (0/1)  
- `threshold`: Sparsification threshold used  

---

### 2. `phymapnet_reliability()` — Ensemble Reliability Inference

Performs hyperparameter-ensemble inference to assess edge stability across multiple model configurations.

#### Inputs

- `otu`: Samples × taxa abundance matrix  
  (rows = samples, columns = taxa)
- `tree`: Phylogenetic tree (`ape::phylo`)  
  with `tip.label` matching OTU column names

#### Outputs

- `rel_mat`: Weighted reliability network (values in [0, 1])  
- `consensus_mat`: Binary consensus network based on `consensus_cut`  
- `edge_list`: Edges ranked by reliability  

---

### Model Parameters

Both functions allow customization of the following parameters:

- `alpha`: Bandwidth parameter controlling the decay of the phylogenetic kernel  
- `k`: Neighborhood scaling factor defining the effective prior sample size (`k × p`)  
- `epsilon1`: Diagonal regularization parameter added to the prior covariance matrix  
- `epsilon2`: Regularization parameter added during precision matrix inversion  
- `kernel`: Kernel type (`"gaussian"` or `"laplacian"`) defining how phylogenetic distances are transformed into similarity weights  
- `normalization`: Data transformation method (`"log"`, `"gmpr"`, `"clr"`, `"tss"`) applied before network inference  
---

### Conceptual Difference

- `phymapnet_fit()` constructs a network under a fixed set of hyperparameters.  
- `phymapnet_reliability()` evaluates network stability across multiple parameter configurations and identifies robust edges through ensemble consensus.


## Toy Example (Fully Reproducible)

This example generates a small synthetic OTU table and a matching phylogenetic tree, then runs both `phymapnet_fit()` and `phymapnet_reliability()`.

```{r, message=FALSE, warning=FALSE}
library(phymapnet)
library(ape)

set.seed(1)

# -------------------------------------------------
# 1. Generate a toy OTU table (10 samples × 6 taxa)
# -------------------------------------------------

otu <- matrix(rpois(10 * 6, lambda = 20), nrow = 10, ncol = 6)
rownames(otu) <- paste0("Sample", 1:10)
colnames(otu) <- paste0("Taxa", 1:6)

# -------------------------------------------------
# 2. Generate a matching phylogenetic tree
# -------------------------------------------------

tree <- read.tree(text = 
  "((Taxa1:0.1,Taxa2:0.1):0.2,
     (Taxa3:0.2,
        (Taxa4:0.1,
           (Taxa5:0.05,Taxa6:0.05):0.05
        ):0.1
     ):0.1
   );"
)
plot(tree, main = "Toy Phylogenetic Tree")
```

## Single-Model Inference
```bash
fit <- phymapnet_fit(
  otu,
  tree,
  alpha = 0.05,
  k = 3,
  epsilon1 = 0.1,
  epsilon2 = 0.1,
  kernel = "gaussian",
  th_sparsity = 0.90,
  normalization = "log"
)

# Binary adjacency matrix
fit$adjacency
```

## Ensemble Reliability Inference
```bash
res <- phymapnet_reliability(
  otu,
  tree,
  th_fixed = 0.90,
  alpha_range = c(0.03, 0.05),
  k_range = 2:3,
  epsilon1_range = c(0, 0.1),
  epsilon2_range = c(0, 0.1),
  kernels = c("gaussian"),
  normalizations = c("log"),
  consensus_cut = 0.50,
  progress = FALSE
)

# Weighted reliability matrix
res$rel_mat

# Binary consensus network
res$consensus_mat

# Top edges ranked by reliability
head(res$edge_list, 10)

```


## Reporting Issues and Asking Questions

If you encounter a bug, experience a failed function, or have a feature request, please open an issue in the GitHub [issue tracker](https://github.com/rosaaghdam/PhyMapNet/issues). 

## License

CMIMN is licensed under the [GNU General Public License v3.0 (GPL-3)](https://www.gnu.org/licenses/gpl-3.0.html). &copy; Solis-Lemus Lab (2024).


## Citation

If you use PhyMapNet in your work, we kindly ask that you cite the following paper:

```bibtex
@article{shahdoust2025simmapnet,
  title={SimMapNet: A Bayesian Framework for Gene Regulatory Network Inference Using Gene Ontology Similarities as External Hint},
  author={Shahdoust, Maryam and Aghdam, Rosa and Sadeghi, Mehdi},
  journal={bioRxiv},
  pages={2025--04},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}


```
