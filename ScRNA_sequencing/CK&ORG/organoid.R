if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(clustree))install.packages("clustree")
if(!require(viridis))install.packages("viridis")
if(!require(glue))install.packages("glue")
if(!require(knitr))install.packages("knitr")
if(!require(BiocParallel))install.packages("BiocParallel")

##load organoid data
org.data <- Read10X(data.dir = "E:/苏州研究院/单细胞测序结果/result/data/surCK1201_matrix/")
org <- CreateSeuratObject(counts = org.data, project = "org")

##QC and selecting cells for further analysis
org[["percent.mt"]] <- PercentageFeatureSet(org, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(org, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
org <- subset(org, subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & percent.mt < 40)
## Normalizing the data
org <- NormalizeData(org, normalization.method = "LogNormalize", scale.factor = 10000)
## Scaling the data
all.genes <- rownames(org)
org <- ScaleData(org, features = all.genes)
## Run PCA
org <- FindVariableFeatures(org, selection.method = "vst", nfeatures = 2000)
org <- RunPCA(org, features = VariableFeatures(object = org))
DimHeatmap(org, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(org)
## Cluster the cells
org <- FindNeighbors(org, dims = 1:13)
org <- FindClusters(org, resolution = 0.4)
## Run umap
org <- RunUMAP(org, dims = 1:13)
DimPlot(org, reduction = "umap")
# org <- RunTSNE(org, dims = 1:13)
# DimPlot(org, reduction = "tsne")

## Cluster annotation
diff.wilcox = FindAllMarkers(org)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "../table/org_diff_genes.csv", row.names = F)
write.csv(top10, "../table/org_top10_diff_genes.csv", row.names = F)

## Top10 genes
top10_genes <- read.csv("../../output/cell_identify/top10_diff_genes_wilcox_modified.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(org))
plot1 = DoHeatmap(org, features = top10_genes, slot = "data", group.by = "seurat_clusters", group.bar = T, size = 6)
ggsave("../../output/image/top10_markers_modified.pdf", plot=plot1, width=12, height=8.13) 

## SingleR identify cell types
## HumanPrimaryCellAtlasData as ref dataset 
if(!require(SingleR))install.packages("SingleR")
refdata <- celldex::HumanPrimaryCellAtlasData()
testdata <- GetAssayData(org, slot="data")
clusters <- org@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts",  assay.type.ref = "logcounts"
                    )
table(cellpred$labels) ##only Epithelial_cells 

## Human colorectal datasets from the scRNAseq package
if(!require(scRNAseq))install.packages("scRNAseq")
sce <- HeOrganAtlasData()
sce <- sce[,!is.na(sce$label)]
# SingleR() expects reference datasets to be normalized and log-transformed
library(scuttle)
sceM <- logNormCounts(sce)
cellpred2 <- SingleR(test = testdata, ref = sceM, labels = sceM$label, 
                     method = "cluster", clusters = clusters,
                     assay.type.test = "logcounts",  assay.type.ref = "logcounts"
                     )
table(cellpred2$labels) 

## MonacoImmuneData as ref dataset
refdata <- MonacoImmuneData()
cellpred3 <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                     method = "cluster", clusters = clusters,
                     assay.type.test = "logcounts",  assay.type.ref = "logcounts"
                     )
table(cellpred3$labels) 

## DICE as ref dataset
refdata <- DatabaseImmuneCellExpressionData()
cellpred4 <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main,
                     method = "cluster", clusters = clusters,
                     assay.type.test = "logcounts",  assay.type.ref = "logcounts"
)
table(cellpred4$labels) 

saveRDS(org, file = "../RDS data/organoid.rds")


