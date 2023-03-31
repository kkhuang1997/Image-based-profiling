library("Rfast")
library("MASS")
library("LMGene")
library("onlinePCA")
library("RColorBrewer")
library(paletteer)
library(ggbiplot)
library(dplyr)

## delete null element
features_for_jpg = readRDS("../rds data/features_for_1960_jpg.rds")
for(ii in seq_len(length(features_for_jpg))){
  if(dim(features_for_jpg[[ii]])[1] == 0 ){
    features_for_jpg[[ii]][2] = 404
  }
}

## feature modification and extract non_para features
features = features_for_jpg %>%
  do.call(what = rbind, args = .) %>%
  as_tibble(features_for_jpg) %>% 
  mutate(., "Org_size" = .["x.0.s.area"]) 
non_para_features = features[,402:405]
features = apply((features[,1:401]), MARGIN = 2, as.numeric)

## feature selection
fea_mask = colMads(features) == 0 ##  bad feature
features = features[,which(fea_mask==FALSE)]
rownames(features) = seq_len(dim(features)[1])
features = as.data.frame(features)

## delete inf and na value
features <- lapply(features, function(x) {
  replace(x, is.na(x) | is.infinite(x), 0)
}) %>% as.data.frame(.)

## Seurat pca
library(Seurat)
## Step 1: Normalize and pre-process the data
row.names(features) = paste("organoid", rownames(features), sep = "_")
fea_seu_obj  = CreateSeuratObject(t(features),project = "extracted_features")
fea_seu_obj <- NormalizeData(fea_seu_obj, normalization.method = "LogNormalize", scale.factor = 10000) 
all.genes <- rownames(fea_seu_obj)
fea_seu_obj <- ScaleData(fea_seu_obj, features = all.genes)
fea_org <- RunPCA(fea_seu_obj, features = all.genes)
print(fea_org[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(fea_org, dims = 1:2, reduction = "pca")
DimPlot(fea_org, reduction = "pca", cols = "gray12")
DimHeatmap(fea_org, dims =1:6, cells = 200, balanced = TRUE)
ElbowPlot(fea_org)

## umap
fea_org <- FindNeighbors(fea_org, dims = 1:9) 
fea_org <- FindClusters(fea_org, resolution = 0.2)
fea_org <- RunUMAP(fea_org, dims = 1:9)
DimPlot(fea_org, reduction = "umap", label = T,
        cols = c("#e79070","#8cab4e","#d58ec9","#fad95e","#2c858d",
                 "#9dd93c","#8dbeb1","#47c297","#8da1dd","#f0615e","#b76cdc",
                 "#82a2ff","#c66b85","#6992b0","#726FFB"))

fea_org$v_1 <- fea_org@reductions[["umap"]]@cell.embeddings[,1]
fea_org$v_2 <- fea_org@reductions[["umap"]]@cell.embeddings[,2]
rd_umap <- c(fea_org$v_1, fea_org$v_2) %>% matrix(.,ncol = 2)
colnames(rd_umap) <- c('UMAP1', 'UMAP2')
umap_df = cbind(rd_umap, non_para_features)

saveRDS(umap_df, "../rds data/umap_df.rds")
saveRDS(fea_org, "../rds data/fea_seurat.rds")
