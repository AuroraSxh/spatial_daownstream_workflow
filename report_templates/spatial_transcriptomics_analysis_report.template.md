# Spatial Transcriptomics Analysis Report Template

## What This Workflow Does

- Prepare Visium-compatible input directories
- Load expression matrix and spatial image/coordinate assets
- Compute QC metrics and filter low-quality spots
- Run normalization, dimensionality reduction, clustering, and marker analysis
- Generate spatial and embedding-level annotation plots

## Input Summary

- Sample name:
- Input H5 matrix:
- Input spatial archive:

## QC Summary

- Spots before filtering:
- Spots after filtering:
- Genes detected:
- Cluster count:

## Key Outputs

- `results/plots/qc_violin_before_filtering.png`
- `results/plots/qc_scatter_before_filtering.png`
- `results/plots/spatial_qc_features.png`
- `results/plots/umap_clusters.png`
- `results/plots/spatial_clusters.png`
- `results/plots/marker_heatmap_top5.png`
- `results/tables/cluster_markers_all.csv`
- `results/tables/cluster_markers_top10.csv`
- `results/tables/cluster_cell_type_annotations.csv`
- `results/objects/<sample>_seurat.rds`

## Notes

- Cluster annotations are heuristic unless replaced by a stricter reference-based annotation stage.
- This template is included as a reporting scaffold, not as a real analysis output.
