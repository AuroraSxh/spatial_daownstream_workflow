# Spatial Transcriptomics Downstream Workflow

Code-only downstream spatial transcriptomics workflow for 10x Visium-style inputs using R and Seurat.

This repository is packaged for public sharing of workflow logic, analysis scripts, and output structure without bundling data, derived objects, or generated reports.

## Overview

The workflow covers:

- Visium input preparation
- spatial matrix loading
- QC metric computation and filtering
- SCTransform normalization
- PCA, graph construction, clustering, and UMAP
- marker extraction and heuristic annotation
- export of plots, tables, object checkpoints, and summary JSON

## Repository Layout

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
├── expected_results/
└── report_templates/
```

## Workflow Entry Points

- `scripts/install_r_packages.R`: install required R dependencies
- `scripts/run_spatial_analysis.R`: run the downstream workflow

## Quick Start

Install dependencies and run the analysis from the repository root:

```bash
Rscript scripts/install_r_packages.R
Rscript scripts/run_spatial_analysis.R
```

The workflow expects a local input directory with:

- one Visium feature-barcode H5 file
- one Visium spatial tarball

The repository does not ship the actual dataset.

## Repository Policy

Included:

- analysis scripts
- installation helper
- documentation
- report template
- expected results structure

Excluded:

- test data
- raw or processed spatial objects
- generated plots and tables
- local package libraries
- generated `results/` directories

## Expected Results

See:

- `expected_results/README.md`
- `report_templates/spatial_transcriptomics_analysis_report.template.md`

These files define the expected deliverable structure without publishing real outputs.
