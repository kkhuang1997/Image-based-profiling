library("Rfast")
library("MASS")
library("LMGene")
library("onlinePCA")
library("RColorBrewer")
library(paletteer)
library(ggbiplot)
library(dplyr)
library(feather)

## load processed features
df_features_meta <- arrow::read_feather("../python code/intermediate_data/StageOne_feature_meta.feather")

## feature modification and extract non_numeric features
df_features = apply((df_features_meta[,0:401]), MARGIN = 2, as.numeric)
df_meta = df_features_meta[, 402:410]
## feature selection
fea_mask = colMads(df_features) == 0 ##  bad feature
df_features = df_features[,which(fea_mask==FALSE)]
df_features = as.data.frame(df_features)

## Seurat pca
library(Seurat)
## Step 1: Normalize and pre-process the data
row.names(df_features) = paste("organoid", rownames(df_features), sep = "_")
fea_seu_obj  = CreateSeuratObject(t(df_features),project = "Process_features")
fea_seu_obj <- NormalizeData(fea_seu_obj, normalization.method = "LogNormalize", scale.factor = 10000) 
all.features <- rownames(fea_seu_obj)
fea_seu_obj <- ScaleData(fea_seu_obj, features = all.features)
fea_PCA <- RunPCA(fea_seu_obj, features = all.features, approx=FALSE)
print(fea_PCA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(fea_PCA, dims = 1:2, reduction = "pca")
DimPlot(fea_PCA, reduction = "pca", cols = "gray12")
DimHeatmap(fea_PCA, dims =1:6, cells = 200, balanced = TRUE)
ElbowPlot(fea_PCA)

## umap
fea_PCA <- FindNeighbors(fea_PCA, dims = 1:10) 
fea_PCA <- FindClusters(fea_PCA, resolution = 0.2)
fea_PCA <- RunUMAP(fea_PCA, dims = 1:10)
DimPlot(fea_PCA, reduction = "umap", label = T,
        cols = c("#e79070","#8cab4e","#d58ec9","#fad95e","#2c858d",
                 "#9dd93c","#8dbeb1","#47c297","#8da1dd","#f0615e","#b76cdc",
                 "#82a2ff","#c66b85","#6992b0","#726FFB"))
# fea_PCA <- RunTSNE(fea_PCA, dims = 1:15)
# DimPlot(fea_PCA, reduction = "tsne")
fea_PCA$v_1 <- fea_PCA@reductions[["umap"]]@cell.embeddings[,1]
fea_PCA$v_2 <- fea_PCA@reductions[["umap"]]@cell.embeddings[,2]
fea_umap <- c(fea_PCA$v_1, fea_PCA$v_2) %>% matrix(.,ncol = 2)
colnames(fea_umap) <- c('UMAP1', 'UMAP2')
df_umap_meta = cbind(fea_umap, df_meta)

## save umap data
write_feather(df_umap_meta, "../python code/intermediate_data/StageFour_umap_meta.feather")
