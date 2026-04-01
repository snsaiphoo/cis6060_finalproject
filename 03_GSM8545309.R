library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# 22128 features across 11742 samples

# Set seed for reproducibility
set.seed(42)


# Load Data
time_load <- system.time({
  data <- Read10X("GSM8545309_GP/")
})

# Create Seurat Object
time_create <- system.time({
  obj <- CreateSeuratObject(
    counts = data,
    project = "growth_plate",
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
# The quality control (QC) violin plots illustrate the distribution of gene counts (nFeature_RNA), transcript counts (nCount_RNA), and mitochondrial gene expression percentage (percent.mt) across all cells. Most cells exhibit a moderate number of detected genes (approximately 1,000–3,000) and transcript counts, indicating generally good sequencing depth and cell quality. However, a subset of cells shows unusually high gene and transcript counts, which may represent doublets or multiplets. Additionally, while the majority of cells have low mitochondrial gene expression (<10%), there is a noticeable tail of cells with elevated mitochondrial percentages, suggesting the presence of stressed or dying cells. Overall, the dataset appears high quality but benefits from filtering to remove low-quality cells and potential doublets prior to downstream analysis.
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# The feature scatter plots show a strong positive correlation (r ≈ 0.91) between nCount_RNA and nFeature_RNA, indicating consistent library complexity across cells. In contrast, percent.mt shows little correlation with sequencing depth (r ≈ −0.01). A subset of cells with low counts and high mitochondrial expression suggests the presence of low-quality or stressed cells, supporting the need for QC filtering.
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# QC Filtering
# In order to remove low-quality tails seen in Vln plots. Remove potential doublets
# found in cells above 7000, the mitochondrial pt, most good cells are <5%, so 10% is safe
time_qc <- system.time({
  obj <- subset(
    obj,
    subset =
      nFeature_RNA > 800 &
      nFeature_RNA < 6000 &
      percent.mt < 5
  )
})

# Normalization
time_norm <- system.time({
  obj <- NormalizeData(obj)
})

# Remove unwanted genes before HVG selection
genes_keep <- rownames(obj)[!grepl("^Hb|^Rpl|^Rps", rownames(obj))]
obj <- obj[genes_keep, ]

# Variable Features
time_hvg <- system.time({
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000
  )
})

length(VariableFeatures(obj))

seurat_hvg <- VariableFeatures(obj)
write.csv(seurat_hvg, "seurat_03/seurat_hvg_03.csv", row.names = FALSE)

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


p <- VizDimLoadings(obj, dims = 1, reduction = "pca")

p +
  ggtitle("PC1 Loadings - Seurat 03") +      
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)  # center the title
  )

p <- DimPlot(obj, reduction = "pca")

# Get variance explained
var_explained <- obj[["pca"]]@stdev^2
var_percent <- round(100 * var_explained / sum(var_explained), 1)

# Add labels
p +  ggtitle("PC1-PC2 - Seurat 03") +   
  xlab(paste0("PC1 (", var_percent[1], "%)")) +
  ylab(paste0("PC2 (", var_percent[2], "%)")) +
  theme(
    plot.title = element_text(hjust = 0.7)  # center the title
  )


# PCA Evaluation Metrics

# Embeddings
time_embed <- system.time({
  pca_embeddings <- Embeddings(obj, "pca")[, 1:2]
})
write.csv(pca_embeddings, "seurat_03/seurat_pca_embeddings_03.csv")

# Variance explained
time_var <- system.time({
  stdev <- obj[["pca"]]@stdev
  var_explained <- (stdev^2) / sum(stdev^2)
})

variance_df <- data.frame(
  PC = paste0("PC", 1:length(var_explained)),
  Variance = var_explained
)

write.csv(variance_df, "seurat_03/seurat_variance_explained_03.csv")

# Loadings
time_loadings <- system.time({
  loadings <- Loadings(obj, "pca")
})

pc1_vals <- sort(abs(loadings[,1]), decreasing = TRUE)[1:10]
pc2_vals <- sort(abs(loadings[,2]), decreasing = TRUE)[1:10]

pc1_df <- data.frame(Gene = names(pc1_vals), Loading = pc1_vals)
pc2_df <- data.frame(Gene = names(pc2_vals), Loading = pc2_vals)

write.csv(pc1_df, "seurat_03/seurat_pc1_genes_03.csv", row.names = FALSE)
write.csv(pc2_df, "seurat_03/seurat_pc2_genes_03.csv", row.names = FALSE)

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
    "UMAP"
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
    time_umap["elapsed"]
  )
)

print(runtime_summary)

write.csv(runtime_summary, "seurat_03/runtime_summary_growth_plate.csv", row.names = FALSE)
