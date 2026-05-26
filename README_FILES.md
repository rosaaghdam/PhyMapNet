# PhyMapNet Repository Files and Reproducibility Workflow

This document describes the files to upload to GitHub and the order in which
the revised-paper outputs can be reproduced. It is intentionally separate from
the main front-page README so that package users can quickly find installation
and function documentation, while paper readers can find the analysis files.

## Repository Structure

```text
PhyMapNet/
|-- README.md                                # GitHub front page and package usage
|-- README_FILES.md                          # File guide and analysis workflow
|-- .gitignore                               # Excludes local R/macOS artifacts
|-- PhyMapNet_Git.Rproj                     # Optional RStudio project file
|-- phymapnet/                               # Editable R package source (v0.1.3)
|   |-- DESCRIPTION                          # Package metadata and dependencies
|   |-- NAMESPACE                            # Exported functions
|   |-- R/                                   # Package implementation
|   |-- man/                                 # Function documentation
|   `-- tests/testthat/                      # Package tests
|-- 1-filtering/
|   `-- 1_run_filtering.R                    # Creates filtered OTU and distance files
|-- 2-sensitivity_analysis/
|   |-- generate_sensitivity_reliability_networks.R
|   |                                        # Expensive upstream network generation
|   `-- generate_sensitivity_figure2_tables_s3_s4.R
|                                            # Figure 2 and Tables S3/S4 from results
|-- 3-bootstrap/
|   |-- bootstrap_phymapnet_helpers.R        # Fixed-model calculations used in stability run
|   |-- run_bootstrap_stability_analysis.R   # Expensive full run; smoke test by default
|   `-- generate_bootstrap_figure.R          # Figure from completed bootstrap summary
|-- 4-overlap_cminet/
|   |-- generate_cminet_reliability_master.R # Expensive upstream reliability generation
|   `-- generate_cminet_overlap_figures.R    # Overlap statistics and figures from masters
|-- 5-biological_evaluation/
|   `-- generate_hmp_biological_evaluation_figure.R
|                                            # HMP figure and biological tables
|-- 6-running_time/
|   |-- run_phymapnet_hmp_taxa_runtime.R     # HMP taxa-size runtime workflow
|   `-- run_method_runtime_comparison.R      # Multi-method runtime workflow
|-- scripts/
|   `-- reproduce_paper_outputs.sh           # One-command output reproduction
|-- data/
|   |-- original_data/                       # Original study R data objects
|   |-- filter/                              # Filtered OTU/distance inputs and HMP taxonomy
|   `-- generate/                            # Bootstrap and noisy replicate inputs
|-- result/
|   |-- reliable_score_all/                  # Precomputed sensitivity networks
|   |-- reliability_master/                  # Fixed-normalization/kernel masters
|   `-- selected_important_otus.csv
|                                            # Important genus-level taxa in the HMP biological network
|-- check_tree_distance_inputs.R             # Simple input-alignment helper
|-- example_run_phymapnet_tree_and_distance.R
|-- test_local_phymapnet_0.1.3.R
`-- test_local_phymapnet_0.1.3_real_data.R
```

## Terminal Reproduction Workflow

From the repository root:

```bash
bash scripts/reproduce_paper_outputs.sh
```

The script regenerates final outputs from retained precomputed results without
starting computationally intensive upstream inference. Figure files created
during a local run remain untracked through `.gitignore`.

### Outputs Regenerated

| Analysis | Script | Output location |
|---|---|---|
| Filtering | `1-filtering/1_run_filtering.R` | `data/filter/` |
| Sensitivity Figure 2 and Tables S3/S4 | `2-sensitivity_analysis/generate_sensitivity_figure2_tables_s3_s4.R` | `result/sensitivity_analysis/` |
| CMiNet overlap statistics and figures | `4-overlap_cminet/generate_cminet_overlap_figures.R` | Generated locally in `result/cminet_comparison/figures/` |
| HMP biological figure and tables | `5-biological_evaluation/generate_hmp_biological_evaluation_figure.R` | Generated locally in `result/biological_evaluation_filtered312/` |

### Retained Biological Result Table

`result/selected_important_otus.csv` reports the
important taxa identified in the HMP biological network at the genus level. It
includes network-importance measures and taxonomic annotations for the selected
non-isolated taxa used to interpret the biological network.

## Expensive Upstream Analyses

The following scripts are included for completeness but are intentionally not
called by the one-command reproduction workflow.

| Script | Reason not run automatically |
|---|---|
| `2-sensitivity_analysis/generate_sensitivity_reliability_networks.R` | Full sensitivity network generation takes approximately days; precomputed outputs are retained in `result/reliable_score_all/`. |
| `3-bootstrap/run_bootstrap_stability_analysis.R` with full settings | Complete bootstrap/noisy inference is very time-consuming; completed bootstrap summaries are not included in the GitHub upload. |
| `4-overlap_cminet/generate_cminet_reliability_master.R` | Completed master files are retained in `result/reliability_master/`; raw CMiNet edge-list inputs are not currently included. |
| `6-running_time/run_phymapnet_hmp_taxa_runtime.R` with full settings | Full runtime benchmarking is computationally expensive. |
| `6-running_time/run_method_runtime_comparison.R` with full settings | Full seven-method timing comparison is computationally expensive. |


## Package Archive Recommendation

Upload the editable `phymapnet/` source folder to GitHub. The CRAN submission
archive `phymapnet_0.1.3.tar.gz` is best attached as a GitHub Release asset or
deposited with a permanent study archive, rather than committed as tracked
source. The source folder remains the authoritative GitHub code.
