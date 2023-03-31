if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")

## load matrix
ck.data <- Read10X(data.dir = "E:/苏州研究院/单细胞测序结果/result/data/surCK_matrix/")
ck <- CreateSeuratObject(counts = ck.data, project = "surCK")

## QC and selecting cells for further analysis
ck[["percent.mt"]] <- PercentageFeatureSet(ck, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(ck, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ck <- subset(ck, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 75)
## Normalizing the data
ck <- NormalizeData(ck, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features
ck <- FindVariableFeatures(ck, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ck), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ck)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2
##Scaling the data 
all.genes <- rownames(ck)
ck <- ScaleData(ck, features = all.genes)

##Run PCA
ck <- RunPCA(ck, features = VariableFeatures(object = ck))
print(ck[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ck, dims = 1:2, reduction = "pca")
DimPlot(ck, reduction = "pca")
DimHeatmap(ck, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(ck)

## Cluster the cells
ck <- FindNeighbors(ck, dims = 1:12)
ck <- FindClusters(ck, resolution = 0.2)

## UMAP
ck <- RunUMAP(ck, dims = 1:12)
DimPlot(ck, reduction = "umap")
# ck <- RunTSNE(ck, dims = 1:15)
# DimPlot(ck, reduction = "tsne")

saveRDS(ck, "../RDS data/surCK.rds")

