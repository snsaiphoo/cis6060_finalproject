library(anndataR)
library(Seurat)

adata <- read_h5ad("GSE308009_ASW02/GSE308009_Aged_SW_02_anndata.h5ad")


agedsw <- CreateSeuratObject(counts = adata$X, project="AgedSW", min.cells = 3, min.features = 200)

agedsw

