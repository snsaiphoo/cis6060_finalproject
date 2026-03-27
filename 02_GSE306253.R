library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# 22240 features across 6045 samples within 1 assay

# Set seed for reproducibility
set.seed(42)


# Load Data
time_load <- system.time({
  data <- Read10X("GSM9195629_rep1/")
})

# Create Seurat Object
time_create <- system.time({
  obj <- CreateSeuratObject(
    counts = data,
    project = "pcl_injury",
    min.cells = 3,
    min.features = 200
  )
})

DefaultAssay(obj) <- "RNA"
options(Seurat.object.assay.version = "v5")
obj <- UpdateSeuratObject(obj)
print(obj)

# Mitochondrial Percentage
time_mt <- system.time({
  counts <- GetAssayData(obj, layer = "counts")
  mt_genes <- grep("^mt-", rownames(counts), value = TRUE)
  
  percent.mt <- colSums(counts[mt_genes, ]) / colSums(counts) * 100
  obj[["percent.mt"]] <- percent.mt
})

# QC Plots

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC Filtering
# QC Filtering:
# Remove low-quality cells (nFeature_RNA < 1000)
# Remove potential doublets (nFeature_RNA > 7000, nCount_RNA > 50000)
# Remove high mitochondrial cells (percent.mt > 10)
time_qc <- system.time({
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 1000 &
      nFeature_RNA < 7000 &
      nCount_RNA < 50000 &
      percent.mt < 10
  )
})

# Normalization
time_norm <- system.time({
  obj <- NormalizeData(obj)
})


# Variable Features
time_hvg <- system.time({
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000
  )
})

# Recalculate top features AFTER filtering
top10 <- head(VariableFeatures(obj), 10)

# Plot
plot <- VariableFeaturePlot(obj)

LabelPoints(
  plot = plot,
  points = top10,
  repel = TRUE
)

# Scaling
time_scale <- system.time({
  obj <- ScaleData(obj)
})

# PCA (IRLBA)
time_pca <- system.time({
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    approx = TRUE,
    seed.use = 42
  )
})

# PCA Visualization
ElbowPlot(obj, ndims = 10)

print(obj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(obj, dims = 1:2, reduction = "pca")

DimPlot(obj, reduction = "pca")

# PCA Evaluation Metrics

# Embeddings
time_embed <- system.time({
  pca_embeddings <- Embeddings(obj, "pca")[, 1:2]
})
write.csv(pca_embeddings, "seurat_02/seurat_pca_embeddings_02.csv")

# Variance explained
time_var <- system.time({
  stdev <- obj[["pca"]]@stdev
  var_explained <- (stdev^2) / sum(stdev^2)
})

variance_df <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained
)

write.csv(variance_df, "seurat_02/seurat_variance_explained_02.csv")

# Loadings
time_loadings <- system.time({
  loadings <- Loadings(obj, "pca")
})

pc1_vals <- sort(abs(loadings[,1]), decreasing = TRUE)[1:10]
pc2_vals <- sort(abs(loadings[,2]), decreasing = TRUE)[1:10]

pc1_df <- data.frame(Gene = names(pc1_vals), Loading = pc1_vals)
pc2_df <- data.frame(Gene = names(pc2_vals), Loading = pc2_vals)

write.csv(pc1_df, "seurat_02/seurat_pc1_genes_02.csv", row.names = FALSE)
write.csv(pc2_df, "seurat_02/seurat_pc2_genes_02.csv", row.names = FALSE)

# Clustering 
time_neighbors <- system.time({
  obj <- FindNeighbors(obj, dims = 1:10)
})

time_clusters <- system.time({
  obj <- FindClusters(obj)
})

# UMAP
time_umap <- system.time({
  obj <- RunUMAP(obj, dims = 1:10)
})

DimPlot(obj, reduction = "umap", label = TRUE)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "seurat_clusters")

# Marker Gene Analysis
time_markers <- system.time({
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
})

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top_markers, "seurat_02/seurat_top_markers_02.csv", row.names = FALSE)

top_markers %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  View()

library(dplyr)

cluster_summary <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  summarise(genes = paste(gene, collapse = ", "))

print(cluster_summary)
write.csv(cluster_summary, "seurat_02/seurat_cluster_summary_02.csv", row.names = FALSE)

new_labels <- c(
  "Fibroblasts (Remodeling)",          # 0
  "Inflammatory Myeloid Cells",        # 1
  "Macrophages (MMP12+)",              # 2
  "Activated Fibroblasts",             # 3
  "Mesenchymal Stromal Cells",         # 4
  "Epithelial-like Cells",             # 5
  "Inflammatory Monocytes",            # 6
  "Basement Membrane Fibroblasts",     # 7
  "Stress-Response Cells",             # 8
  "Neutrophils",                       # 9
  "Activated Endothelial Cells",       # 10
  "Activated Neutrophils",             # 11
  "Dendritic Cells",                   # 12
  "T Cells",                           # 13
  "Pericytes / Smooth Muscle Cells",   # 14
  "Fibroblasts (ECM)",                 # 15
  "Lymphatic Endothelial Cells",       # 16
  "Schwann / Neural Cells",            # 17
  "Skeletal Muscle Cells",             # 18
  "Endothelial Cells (Apln+)"          # 19
)

names(new_labels) <- levels(obj)
obj <- RenameIdents(obj, new_labels)

DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(obj, reduction = "pca", repel = TRUE) 

# Runtime Summary
runtime_summary <- data.frame(
  Step = c(
    "Load Data",
    "Create Object",
    "Mitochondrial %",
    "QC Filtering",
    "Normalization",
    "Variable Features",
    "Scaling",
    "PCA",
    "Embeddings",
    "Variance",
    "Loadings",
    "Neighbors",
    "Clustering",
    "UMAP",
    "Marker Detection"
  ),
  Time_sec = c(
    time_load["elapsed"],
    time_create["elapsed"],
    time_mt["elapsed"],
    time_qc["elapsed"],
    time_norm["elapsed"],
    time_hvg["elapsed"],
    time_scale["elapsed"],
    time_pca["elapsed"],
    time_embed["elapsed"],
    time_var["elapsed"],
    time_loadings["elapsed"],
    time_neighbors["elapsed"],
    time_clusters["elapsed"],
    time_umap["elapsed"],
    time_markers["elapsed"]
  )
)

print(runtime_summary)

write.csv(runtime_summary, "seurat_02/runtime_summary_pcl.csv", row.names = FALSE)
