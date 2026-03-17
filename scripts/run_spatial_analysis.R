options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(jsonlite)
})

set.seed(1234)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_root, "test_data")
results_dir <- file.path(project_root, "results")
prepared_dir <- file.path(results_dir, "prepared_visium_input")
plots_dir <- file.path(results_dir, "plots")
tables_dir <- file.path(results_dir, "tables")
objects_dir <- file.path(results_dir, "objects")
report_dir <- file.path(project_root, "report")
log_dir <- file.path(project_root, "logs")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(prepared_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

h5_name <- "Visium_Fresh_Frozen_Adult_Mouse_Brain_raw_feature_bc_matrix.h5"
h5_path <- file.path(input_dir, h5_name)
spatial_tar <- file.path(input_dir, "Visium_Fresh_Frozen_Adult_Mouse_Brain_spatial.tar.gz")
prepared_h5 <- file.path(prepared_dir, h5_name)
prepared_spatial <- file.path(prepared_dir, "spatial")

if (!file.exists(h5_path)) {
  stop("Missing H5 file: ", h5_path)
}
if (!file.exists(spatial_tar)) {
  stop("Missing spatial tarball: ", spatial_tar)
}

if (!file.exists(prepared_h5)) {
  ok <- file.symlink(from = h5_path, to = prepared_h5)
  if (!ok) {
    file.copy(h5_path, prepared_h5, overwrite = TRUE)
  }
}

if (!dir.exists(prepared_spatial) || length(list.files(prepared_spatial)) == 0) {
  utils::untar(spatial_tar, exdir = prepared_dir)
}

visium <- Load10X_Spatial(
  data.dir = prepared_dir,
  filename = h5_name,
  assay = "Spatial",
  slice = "adult_mouse_brain",
  filter.matrix = TRUE
)

visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "(?i)^mt-")

qc_table <- visium@meta.data %>%
  rownames_to_column("barcode") %>%
  as_tibble()

write_csv(qc_table, file.path(tables_dir, "spot_qc_metrics_before_filtering.csv"))

qc_violin <- VlnPlot(
  visium,
  features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
) + patchwork::plot_annotation(title = "QC metrics before filtering")

ggsave(
  filename = file.path(plots_dir, "qc_violin_before_filtering.png"),
  plot = qc_violin,
  width = 13,
  height = 4.5,
  dpi = 300
)

qc_scatter <- (
  FeatureScatter(visium, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
    ggtitle("UMI vs detected genes")
) | (
  FeatureScatter(visium, feature1 = "nCount_Spatial", feature2 = "percent.mt") +
    ggtitle("UMI vs mitochondrial fraction")
)

ggsave(
  filename = file.path(plots_dir, "qc_scatter_before_filtering.png"),
  plot = qc_scatter,
  width = 11,
  height = 4.5,
  dpi = 300
)

spatial_qc_plot <- wrap_plots(
  SpatialFeaturePlot(visium, features = "nCount_Spatial"),
  SpatialFeaturePlot(visium, features = "nFeature_Spatial"),
  SpatialFeaturePlot(visium, features = "percent.mt"),
  ncol = 3
) + patchwork::plot_annotation(title = "Spatial QC features")

ggsave(
  filename = file.path(plots_dir, "spatial_qc_features.png"),
  plot = spatial_qc_plot,
  width = 15,
  height = 5,
  dpi = 300
)

min_features <- max(100, unname(quantile(visium$nFeature_Spatial, 0.05)))
max_features <- unname(quantile(visium$nFeature_Spatial, 0.98))
max_counts <- unname(quantile(visium$nCount_Spatial, 0.98))
max_percent_mt <- min(20, unname(quantile(visium$percent.mt, 0.98)))

visium <- subset(
  visium,
  subset =
    nFeature_Spatial >= min_features &
    nFeature_Spatial <= max_features &
    nCount_Spatial <= max_counts &
    percent.mt <= max_percent_mt
)

write_csv(
  visium@meta.data %>%
    rownames_to_column("barcode") %>%
    as_tibble(),
  file.path(tables_dir, "spot_qc_metrics_after_filtering.csv")
)

visium <- SCTransform(visium, assay = "Spatial", verbose = FALSE)
visium <- RunPCA(visium, assay = "SCT", verbose = FALSE)
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30, verbose = FALSE)
visium <- FindClusters(visium, resolution = 0.5, verbose = FALSE)
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30, verbose = FALSE)

cluster_plot <- DimPlot(visium, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("UMAP clusters")
spatial_cluster_plot <- SpatialDimPlot(visium, label = TRUE, label.size = 3) +
  ggtitle("Spatial clusters")

ggsave(
  filename = file.path(plots_dir, "umap_clusters.png"),
  plot = cluster_plot,
  width = 7,
  height = 5.5,
  dpi = 300
)
ggsave(
  filename = file.path(plots_dir, "spatial_clusters.png"),
  plot = spatial_cluster_plot,
  width = 7,
  height = 6,
  dpi = 300
)

visium <- PrepSCTFindMarkers(visium)
markers <- FindAllMarkers(
  visium,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers <- markers %>%
  arrange(cluster, desc(avg_log2FC))

write_csv(markers, file.path(tables_dir, "cluster_markers_all.csv"))

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write_csv(top_markers, file.path(tables_dir, "cluster_markers_top10.csv"))

heatmap_features <- top_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  pull(gene) %>%
  unique()

heatmap_features <- intersect(heatmap_features, rownames(visium[["SCT"]]$scale.data))

if (length(heatmap_features) > 0) {
  marker_heatmap <- DoHeatmap(
    visium,
    features = heatmap_features,
    assay = "SCT",
    size = 3
  ) + NoLegend()

  ggsave(
    filename = file.path(plots_dir, "marker_heatmap_top5.png"),
    plot = marker_heatmap,
    width = 12,
    height = 9,
    dpi = 300
  )
}

reference_markers <- list(
  Excitatory_Neuron = c("Slc17a7", "Camk2a", "Satb2", "Tbr1"),
  Inhibitory_Neuron = c("Gad1", "Gad2", "Slc6a1", "Pvalb"),
  Astrocyte = c("Aqp4", "Gfap", "Aldh1l1", "Slc1a3"),
  Oligodendrocyte = c("Mbp", "Mobp", "Plp1", "Mag"),
  OPC = c("Pdgfra", "Cspg4", "Olig1", "Olig2"),
  Microglia = c("C1qa", "C1qb", "Ctss", "Tyrobp"),
  Endothelial = c("Cldn5", "Kdr", "Pecam1", "Emcn")
)

cluster_means <- AverageExpression(
  visium,
  assays = "SCT",
  slot = "data",
  group.by = "seurat_clusters",
  verbose = FALSE
)$SCT

available_genes <- rownames(cluster_means)
annotation_scores <- lapply(names(reference_markers), function(cell_type) {
  genes <- intersect(reference_markers[[cell_type]], available_genes)
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(cluster_means)))
  }
  colMeans(cluster_means[genes, , drop = FALSE])
})

annotation_scores <- do.call(rbind, annotation_scores)
rownames(annotation_scores) <- names(reference_markers)
predicted_cell_types <- apply(annotation_scores, 2, function(x) {
  if (all(is.na(x))) {
    return("Unassigned")
  }
  names(which.max(x))
})

annotation_table <- tibble(
  cluster = sub("^g", "", as.character(colnames(cluster_means))),
  predicted_cell_type = unname(predicted_cell_types)
)

score_table <- as_tibble(t(annotation_scores), rownames = "cluster") %>%
  mutate(cluster = sub("^g", "", cluster))
annotation_table <- left_join(annotation_table, score_table, by = "cluster")

write_csv(annotation_table, file.path(tables_dir, "cluster_cell_type_annotations.csv"))

cluster_ids <- as.character(visium$seurat_clusters)
predicted_labels <- annotation_table$predicted_cell_type[match(cluster_ids, annotation_table$cluster)]
predicted_labels[is.na(predicted_labels)] <- "Unassigned"
visium$predicted_cell_type <- factor(predicted_labels)

if (length(levels(visium$predicted_cell_type)) > 0) {
  cell_type_plot <- SpatialDimPlot(
    visium,
    group.by = "predicted_cell_type",
    label = TRUE,
    label.size = 3
  ) + ggtitle("Heuristic cell type labels")

  ggsave(
    filename = file.path(plots_dir, "spatial_predicted_cell_types.png"),
    plot = cell_type_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Finer cluster-level annotation derived from observed top markers per cluster.
fine_annotation_map <- c(
  "0" = "Neuronal_1_Gas7_Ncan",
  "1" = "Hypothalamic_Hcrt_Pmch",
  "2" = "Fezf2_Nr4a2_Neuron",
  "3" = "Penk_Gpr88_Neuron",
  "4" = "Mature_Oligodendrocyte",
  "5" = "Vip_Interneuron",
  "6" = "Prkcd_Vipr2_Neuron",
  "7" = "Choroid_Plexus",
  "8" = "Lefty1_Spink8_Ependymal",
  "9" = "VLMC_Fibroblast",
  "10" = "Avp_Gal_Neuroendocrine",
  "11" = "Cbln1_Ctxn3_Neuron",
  "12" = "Trh_Six3_Neuron"
)

fine_annotation_notes <- c(
  "0" = "Gas7/Ncan/Synpo-rich neuronal cluster",
  "1" = "Hcrt/Pmch/Lhx1os hypothalamic peptidergic neurons",
  "2" = "Fezf2/Nr4a2/Myl4 neuronal cluster",
  "3" = "Penk/Gpr88/C1ql3 neuronal cluster",
  "4" = "Plp1/Mobp/Enpp6 mature oligodendrocytes",
  "5" = "Vip-positive interneuron-like cluster",
  "6" = "Prkcd/Vipr2/Slitrk6 neuronal cluster",
  "7" = "Ttr/Folr1/Kl choroid plexus epithelial cells",
  "8" = "Lefty1/Spink8/Lct ependymal-like cells",
  "9" = "Ogn/Fmod/Aldh1a2 vascular leptomeningeal-like cells",
  "10" = "Avp/Gal/Otp neuroendocrine neurons",
  "11" = "Cbln1/Ctxn3/Fzd5 neuronal cluster",
  "12" = "Trh/Six3/Isl1 neuronal cluster"
)

fine_annotation_table <- annotation_table %>%
  distinct(cluster, .keep_all = TRUE) %>%
  mutate(
    fine_cell_type = unname(fine_annotation_map[cluster]),
    annotation_basis = unname(fine_annotation_notes[cluster])
  )

write_csv(
  fine_annotation_table,
  file.path(tables_dir, "cluster_fine_cell_type_annotations.csv")
)

fine_cluster_ids <- as.character(visium$seurat_clusters)
fine_labels <- fine_annotation_table$fine_cell_type[
  match(fine_cluster_ids, fine_annotation_table$cluster)
]
fine_labels[is.na(fine_labels)] <- "Unassigned"
visium$fine_cell_type <- factor(fine_labels)

fine_cell_type_plot <- SpatialDimPlot(
  visium,
  group.by = "fine_cell_type",
  label = TRUE,
  label.size = 2.8,
  repel = TRUE
) + ggtitle("Fine-grained spatial annotation")

ggsave(
  filename = file.path(plots_dir, "spatial_fine_cell_types.png"),
  plot = fine_cell_type_plot,
  width = 10,
  height = 8,
  dpi = 300
)

fine_umap_plot <- DimPlot(
  visium,
  reduction = "umap",
  group.by = "fine_cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("Fine-grained annotation on UMAP")

ggsave(
  filename = file.path(plots_dir, "umap_fine_cell_types.png"),
  plot = fine_umap_plot,
  width = 9,
  height = 7,
  dpi = 300
)

canonical_genes <- c("Snap25", "Slc17a7", "Gad1", "Mbp", "Aqp4", "Cldn5", "C1qa")
canonical_genes <- intersect(canonical_genes, rownames(visium))

if (length(canonical_genes) > 0) {
  feature_plot <- wrap_plots(
    lapply(canonical_genes, function(gene) SpatialFeaturePlot(visium, features = gene))
  ) + patchwork::plot_annotation(title = "Canonical marker spatial expression")

  ggsave(
    filename = file.path(plots_dir, "spatial_canonical_markers.png"),
    plot = feature_plot,
    width = 14,
    height = 8,
    dpi = 300
  )
}

saveRDS(visium, file.path(objects_dir, "adult_mouse_brain_visium_seurat.rds"))

summary_list <- list(
  sample_name = "Adult Mouse Brain Coronal Section (Fresh Frozen)",
  total_spots_before_filtering = nrow(qc_table),
  total_spots_after_filtering = ncol(visium),
  detected_genes = nrow(visium),
  clusters = length(unique(visium$seurat_clusters)),
  min_features_cutoff = unname(min_features),
  max_features_cutoff = unname(max_features),
  max_counts_cutoff = unname(max_counts),
  max_percent_mt_cutoff = unname(max_percent_mt),
  output_dirs = list(
    plots = plots_dir,
    tables = tables_dir,
    objects = objects_dir,
    report = report_dir
  )
)

write_json(
  summary_list,
  path = file.path(results_dir, "analysis_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

report_lines <- c(
  "# 空间转录组分析说明",
  "",
  "## 我完成了什么",
  "",
  "- 使用 `test_data` 中的 10x Visium 数据完成了标准 Seurat 空间转录组下游分析。",
  "- 自动整理输入目录，读取 `h5` 表达矩阵和 `spatial` 图像/坐标文件。",
  "- 计算 spot 级别质控指标并基于分位数阈值过滤低质量 spot。",
  "- 使用 `SCTransform` 归一化，随后进行了 PCA、邻居图构建、聚类和 UMAP 降维。",
  "- 输出了空间聚类图、UMAP 图、QC 图、marker 热图、代表性 marker 空间表达图。",
  "- 基于经典脑组织 marker 基因，对每个 cluster 给出了启发式细胞类型标签。",
  "",
  "## 分析步骤",
  "",
  "1. 准备 Visium 输入目录：解压 `spatial.tar.gz`，并把 `h5` 与 `spatial/` 组织成 Seurat 可直接读取的结构。",
  "2. 读取数据：`Load10X_Spatial()` 导入表达矩阵、组织切片图和 spot 坐标。",
  "3. 质控：计算 `nCount_Spatial`、`nFeature_Spatial`、`percent.mt`，生成 QC 图。",
  sprintf("4. 过滤：保留 `nFeature_Spatial >= %.0f`、`nFeature_Spatial <= %.0f`、`nCount_Spatial <= %.0f`、`percent.mt <= %.2f` 的 spot。", min_features, max_features, max_counts, max_percent_mt),
  "5. 归一化与降维：`SCTransform`、`RunPCA`、`FindNeighbors`、`FindClusters`、`RunUMAP`。",
  "6. marker 分析：`FindAllMarkers()` 识别每个 cluster 的上调基因，并绘制热图。",
  "7. 启发式注释：根据神经元、星形胶质、少突胶质、OPC、微胶质、内皮细胞等 marker 对 cluster 进行粗略标注。",
  "",
  "## 数据概况",
  "",
  sprintf("- 样本：%s", summary_list$sample_name),
  sprintf("- 过滤前 spot 数：%s", format(summary_list$total_spots_before_filtering, big.mark = ",")),
  sprintf("- 过滤后 spot 数：%s", format(summary_list$total_spots_after_filtering, big.mark = ",")),
  sprintf("- 检测到的基因数：%s", format(summary_list$detected_genes, big.mark = ",")),
  sprintf("- 聚类数：%s", summary_list$clusters),
  "",
  "## 输出位置",
  "",
  sprintf("- 图像文件：`%s`", plots_dir),
  sprintf("- 表格文件：`%s`", tables_dir),
  sprintf("- Seurat 对象：`%s`", objects_dir),
  sprintf("- 汇总 JSON：`%s`", file.path(results_dir, "analysis_summary.json")),
  sprintf("- 本说明文档：`%s`", file.path(report_dir, "spatial_transcriptomics_analysis_report.md")),
  "",
  "## 关键输出文件",
  "",
  "- `plots/qc_violin_before_filtering.png`",
  "- `plots/qc_scatter_before_filtering.png`",
  "- `plots/spatial_qc_features.png`",
  "- `plots/umap_clusters.png`",
  "- `plots/spatial_clusters.png`",
  "- `plots/marker_heatmap_top5.png`",
  "- `plots/spatial_predicted_cell_types.png`",
  "- `plots/spatial_fine_cell_types.png`",
  "- `plots/umap_fine_cell_types.png`",
  "- `plots/spatial_canonical_markers.png`（如果对应基因存在）",
  "- `tables/cluster_markers_all.csv`",
  "- `tables/cluster_markers_top10.csv`",
  "- `tables/cluster_cell_type_annotations.csv`",
  "- `tables/cluster_fine_cell_type_annotations.csv`",
  "- `objects/adult_mouse_brain_visium_seurat.rds`",
  "",
  "## 备注",
  "",
  "- 细胞类型注释属于基于 marker 的启发式结果，用于快速解释 cluster，不等同于严格监督注释。",
  "- `spatial_fine_cell_types.png` 是基于每个 cluster 的 top markers 进行的更细粒度手工注释，分辨率高于粗粒度大类注释。",
  "- 当前流程面向单张 Visium 切片的标准下游分析；若后续需要差异区域分析、空间变异基因检测、cell2location/SPOTlight 去卷积，可在此基础上继续扩展。"
)

writeLines(report_lines, con = file.path(report_dir, "spatial_transcriptomics_analysis_report.md"))
