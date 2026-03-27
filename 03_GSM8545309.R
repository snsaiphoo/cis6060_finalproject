library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
# https://satijalab.org/seurat/articles/pbmc3k_tutorial


# the raw count matrix
# rows --> genes
# columns --> cells
# values --> number of transcripts detected
gp.data <- Read10X("GSM8545309_GP/")

# creates a Seurat object
gp <- CreateSeuratObject(counts = gp.data, project="GSM8545309_GP", min.cells = 3, min.features = 200)

gp
