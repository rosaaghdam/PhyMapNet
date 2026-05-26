#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "PhyMapNet revised-paper reproducibility workflow"
echo "Repository: $ROOT"
echo

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript is not available on PATH." >&2
  exit 1
fi

Rscript -e '
required <- c("ape", "dplyr", "tibble", "tidyr", "ggplot2", "patchwork", "igraph", "influential")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Install required R packages before running reproduction: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}
cat("Required R packages are available.\n")
'

required_files=(
  "result/reliable_score_all/caffeine/baseline.csv"
  "result/reliability_master/phymap/hmp_stool_log_laplacian_master.csv"
  "data/filter/hmp_stool_nonisolated_taxonomy.csv"
)

for file in "${required_files[@]}"; do
  if [[ ! -f "$file" ]]; then
    echo "Error: required precomputed input is missing: $file" >&2
    exit 1
  fi
done

echo
echo "[1/4] Recreating filtered study inputs"
Rscript "1-filtering/1_run_filtering.R"

echo
echo "[2/4] Regenerating sensitivity Figure 2 and Supplementary Tables S3/S4"
Rscript "2-sensitivity_analysis/generate_sensitivity_figure2_tables_s3_s4.R"

echo
echo "[3/4] Regenerating CMiNet overlap statistics and figures"
Rscript "4-overlap_cminet/generate_cminet_overlap_figures.R"

echo
echo "[4/4] Regenerating HMP biological-evaluation figure and tables"
Rscript "5-biological_evaluation/generate_hmp_biological_evaluation_figure.R"

echo
echo "Completed reproducible output generation."
echo "Computationally expensive network-generation, bootstrap, and runtime workflows were not rerun."
echo "See README_FILES.md for full workflow details."
