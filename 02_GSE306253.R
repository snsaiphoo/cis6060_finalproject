library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# https://satijalab.org/seurat/articles/pbmc3k_tutorial

# 22240 features across 6045 samples 

# the raw count matrix
# rows --> genes
# columns --> cells
# values --> number of transcripts detected
pcl.data <- Read10X("GSM9195629_rep1")

# creates a Seurat object
pcl <- CreateSeuratObject(counts = pcl.data, project="GSM9195629_rep1", min.cells = 3, min.features = 200)

pcl
