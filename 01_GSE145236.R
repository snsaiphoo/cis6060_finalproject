library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# https://satijalab.org/seurat/articles/pbmc3k_tutorial


# Dataset 1: 4,425 cells × 15,769 genes

# the raw count matrix
# rows --> genes
# columns --> cells
# values --> number of transcripts detected
injured.data <- Read10X("GSM4319249_injured")

# creates a Seurat object
# the raw expression matrix
# metadata about cells
# normalization results
# PCA results
# clustering results 
# min.cells = 3, this removes genes that appear in fewer than 3 cells, these genes are usually noise or sequencing artifacts 
# min.features = 200, this removes cells that express fewer than 200 genes
# to remove extremely rare genes
# to remove low-quality cells
injured <- CreateSeuratObject(counts = injured.data, project="injured", min.cells = 3, min.features = 200)

injured

# why mitochondrial percentage matters
# healthy cells mostly contain nuclear RNA transcripts
# when a cell is damaged or dying, the cytoplasm leaks and mitochondrial RNA becomes overrepresented 
# so cells with very high mitochondrial percentages are usually dying cells, stressed cells, broken cells, debris captured in droplets 
injured[["percent.mt"]] <- PercentageFeatureSet(injured, pattern = "^mt-")


# Visualize QC metrics as a violin plot
injured_qc <- VlnPlot(injured,
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      ncol = 3) +  plot_annotation(
                        title = "Quality Control Metrics for Injured Sample",
                        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
                      )


#ggsave("figures/injured_qc.png", plot = injured_qc, width = 8, height = 6, dpi = 300)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(injured, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(injured, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scatter_injured <- plot1 + plot2 + plot_annotation(
  title = "Scatterplots of the Injured Sample",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

#ggsave("figures/injured_scatter.png", plot = scatter_injured, width = 8, height = 6, dpi = 300)


injured <- subset(
  injured,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    percent.mt < 10 # more common in muscle tissue, due to more mitochondria
)


VlnPlot(injured, features = c("nFeature_RNA","nCount_RNA","percent.mt"))

plot1 <- FeatureScatter(injured, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(injured, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalizing the data 
injured <- NormalizeData(injured)


# Identification of High Variable Features (Feature Selection)
# vst is the variance stabilizing transformation
# this helps account for genes with higher average expression, as these tend to naturally have higher variance 
# it finds genes whose variability is higher than expected given their average expression level 
# those genes are considered highly variable features
set.seed(123)
injured <- FindVariableFeatures(injured, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(injured), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(injured)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
hvg_injured <- plot2 + plot_annotation(
  title = "HVGs of the Injured Sample",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)

ggsave("figures/injured_hvg.png", plot = hvg_injured, width = 10, height = 6, dpi = 300)


# Scaling the Data
all.genes <- rownames(injured)
injured <- ScaleData(injured, features = all.genes)

# Run PCA
injured <- RunPCA(injured, features = VariableFeatures(object = injured))

# Screeplot of both Injured and Uninjured Samples
scree_injured <- ElbowPlot(injured, ndims = 10, reduction = "pca") +
  plot_annotation(
    title = "Screeplot of Injured Sample",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ) + theme_bw() + scale_x_continuous(breaks = 1:10)

#ggsave("figures/screeplot_injured.png", plot = scree_injured, width = 10, height = 6, dpi = 300)


# The results of the Screeplot show that the most variation has been 
# captured in PC1-PC5 for both injured and uninjured 

# Examine and visualize PCA results a few different ways
print(injured[["pca"]], dims = 1:5, nfeatures = 5)

# Gene loadings of injured sample
VizDimLoadings(injured, dims = 1:2, reduction = "pca") + plot_annotation(
  title = "Gene Loadings of PC1 & PC2 - Injured Sample",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
)


DimPlot(injured, reduction = "pca")

DimHeatmap(injured, dims = 1, cells = 500, balanced = TRUE) 

hm_injured <- DimHeatmap(injured, dims = 1:9, cells = 500, balanced = TRUE)

ggsave("figures/heatmap_injured.png", plot = hm_injured, width = 12, height = 12, dpi = 300)

injured <- RunUMAP(injured, dims = 1:10)
DimPlot(injured, reduction = "umap", label = TRUE)
