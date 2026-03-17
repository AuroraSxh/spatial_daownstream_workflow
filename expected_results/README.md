# Expected Results: Spatial Workflow

Expected result layout:

```text
results/
├── prepared_visium_input/
├── plots/
├── tables/
└── objects/
```

## Directory Meaning

- `prepared_visium_input/`
  - normalized input layout prepared for Seurat loading

- `plots/`
  - QC violin and scatter plots
  - spatial QC overlays
  - UMAP cluster plots
  - spatial cluster plots
  - marker heatmaps
  - coarse and fine annotation maps
  - canonical marker spatial plots

- `tables/`
  - pre/post-filter QC tables
  - all marker tables
  - top marker summaries
  - cluster annotation tables

- `objects/`
  - saved Seurat object checkpoints

Additional top-level artifacts typically include:

```text
report/
└── spatial_transcriptomics_analysis_report.md

results/
└── analysis_summary.json
```

This repository includes only the expected structure, not real outputs.
