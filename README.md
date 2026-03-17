# Spatial Downstream Workflow

This repository contains a downstream spatial transcriptomics workflow for 10x Visium-style inputs using R and Seurat.

It excludes:
- test data
- generated Seurat objects
- generated plots and tables

It includes:
- analysis scripts
- installation helper
- README
- expected results structure
- report template

## Repository Layout

```text
.
├── README.md
├── scripts/
├── expected_results/
└── report_templates/
```

## Included Scripts

- `scripts/install_r_packages.R`
  - install required R dependencies
- `scripts/run_spatial_analysis.R`
  - execute the spatial downstream workflow

## Workflow Summary

The downstream analysis performs:
- Visium input preparation
- spatial matrix loading
- QC metric computation
- filtering
- SCTransform normalization
- PCA / neighbors / clustering / UMAP
- marker extraction
- heuristic cluster annotation
- export of plots, tables, object checkpoints, and summary JSON

## Expected Input

The script expects a `test_data/`-style directory containing:
- one Visium feature-barcode H5 file
- one Visium spatial tarball

The actual dataset is not included in this repository.

## Expected Outputs

See:
- `expected_results/README.md`
- `report_templates/spatial_transcriptomics_analysis_report.template.md`
