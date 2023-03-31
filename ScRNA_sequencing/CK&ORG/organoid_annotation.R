if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(clustree))install.packages("clustree")
if(!require(viridis))install.packages("viridis")

## load data
org = readRDS("../RDS data/organoid.rds")

## 9 cluster
org <- FindNeighbors(org, dims = 1:13, k.param = 50, n.trees = 100)
org <- FindClusters(org, resolution = 0.4)
# DimPlot(org, reduction = "tsne")
# DimPlot(org, reduction = "umap", cols ='Paired')
DimPlot(org, reduction = "umap", cols = c("#003049","#d62828","#A098D6","#fcbf49","#eae2b7","#81b29a","#209A90","#D4B7AF", "#4A3F33") )

## top10 markers
diff.wilcox = FindAllMarkers(org)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# write.csv(top20, "../../output/cell_identity/top20_diff_genes_9cluster.csv", row.names = F)
DoHeatmap(org, features = top10$gene)

## tsne feature plot
# FeaturePlot(org, features = c("OLFM4","AREG","PCNA","FABP1","LYZ","TFF3","TOP2A","SOX4","MALAT1"),reduction = "tsne",cols = c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))
# FeaturePlot(org, features = c("CLDN2","ASCL2","AXIN2","RNF43","CDX2","EPHB2"), reduction = "tsne", cols = c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))

## umap feature plot
DotPlot(org, features = c("PCNA", "TOP2A", "SOX4","AREG","LYZ","CA2","CEACAM5","OLFM4","FABP1","TFF3", "S100A4","ODAM","CDK1","MT-TE"), cols = "RdBu" ) + RotatedAxis()
FeaturePlot(org, features =c("PCNA", "TOP2A", "SOX4","AREG","LYZ","CA2","CEACAM5","OLFM4","FABP1","TFF3", "S100A4","ODAM","CDK1","MT-TE"), reduction = "umap", cols =c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))

##  tsne cell annotation
# new.cluster.ids <- c("TA","Stem/progenitor","Goblet","Enterocyte","Stem","Enterocyte","TA","TA","MT")
# names(new.cluster.ids) <- levels(org)
# org <- RenameIdents(org, new.cluster.ids)
# DimPlot(org, reduction = "tsne", label = F)

## umap annotation
new.cluster.ids <- c("TA","Progenitor","Paneth","Enterocyte","Stem","Goblet","S100A4+","Cell cycle","MT")
names(new.cluster.ids) <- levels(org)
org <- RenameIdents(org, new.cluster.ids)
DimPlot(org, reduction = "umap", label = F, 
        cols = c("#003049","#d62828","#A098D6","#fcbf49","#eae2b7","#81b29a","#209A90","#D4B7AF", "#4A3F33"))

## remove mitochondrion cluster
org = subset(org, 
             idents = unique(Idents(org))[!unique(Idents(org))=="MT"])
DimPlot(org, reduction = "umap", label = F, 
        cols = c("#003049","#d62828","#A098D6","#fcbf49","#eae2b7","#81b29a","#209A90","#D4B7AF"))
DotPlot(org, features = c("PCNA", "TOP2A", "SOX4","AREG","LYZ","CEACAM5","OLFM4","FABP1","TFF3", "S100A4","ODAM","CDK1"), cols = "RdBu" ) + RotatedAxis()
FeaturePlot(org, features =c("PCNA", "TOP2A", "SOX4","AREG","LYZ","CEACAM5","OLFM4","FABP1","TFF3", "S100A4","ODAM","CDK1","CHGA","CHGB"), reduction = "umap", cols =c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))

saveRDS(org,file = "../RDS data/organoid_annotation.rds")


## using slingshot
library(slingshot)
sds <- slingshot(Embeddings(org, "tsne"), clusterLabels = org$seurat_clusters, 
                 start.clus = 4, stretch = 0)
## Assign a color to each cell based on some value
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}
cell_colors <- cell_pal(org$cell_type, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(org$seurat_clusters, hue_pal())

