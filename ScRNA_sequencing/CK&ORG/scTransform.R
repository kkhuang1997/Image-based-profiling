if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")
library(Seurat)
library(patchwork)

## load data
org <- readRDS("../RDS data/organoid.rds")
ck <- readRDS("../RDS data/surCK.rds")
seurat_list <- list(org, ck)
seurat_list <- lapply(seurat_list, SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
## find anchors
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",
                                  anchor.features = features)
# saveRDS(anchors, "../../output/anchors.rds")
anchors <- readRDS("E:/RDS_data/anchors.rds")
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
# saveRDS(combined.sct, "../../output/combined_sct.rds")

## cluster
combined.sct <- readRDS("E:/RDS_data/20220209/combined_sct.rds")
combined.sct <- RunPCA(combined.sct)
ElbowPlot(combined.sct)
combined.sct <- FindNeighbors(combined.sct, dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 0.4)
combined.sct <- RunTSNE(combined.sct, dims = 1:30)
DimPlot(combined.sct, reduction = "tsne", group.by = "orig.ident", cols =c("#D8383A","#14517c"))
DimPlot(combined.sct, reduction = "tsne", repel = T, label = T, label.size = 4)
DimPlot(combined.sct, reduction = "tsne", split.by = "orig.ident")

## top10 markers
diff.wilcox = FindAllMarkers(combined.sct, min.pct = 0.25)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "../table/combined_top10_diff_genes.csv", row.names = F)
write.csv(all.markers, "../table/combined_diff_genes.csv", row.names = F)
DoHeatmap(combined.sct, features = top10$gene)
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(combined.sct, features = top5$gene)

## the number of cells in each cluster
plot.data <- as.data.frame(combined.sct@meta.data) %>%
  select(orig.ident, seurat_clusters) %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(count = n()) %>%
  mutate(clust_total = sum(count)) %>%
  mutate(clust_prop = count / clust_total)  %>%
  group_by(orig.ident)  %>%
  mutate(dataset_total = sum(count)) %>%
  ungroup() %>%
  mutate(dataset_prop = count / dataset_total)

ggplot(plot.data, aes(x = seurat_clusters, y = count, fill = orig.ident)) +
  geom_col()
ggplot(plot.data, aes(x = seurat_clusters, y = clust_prop, fill = orig.ident)) +
  geom_col()
ggplot(plot.data, aes(x = seurat_clusters, y = count, fill = orig.ident)) +
  geom_col(position = position_dodge(0.9))

## Conserved markers
DefaultAssay(combined.sct) <- "RNA"
ident_0 <- FindConservedMarkers(combined.sct, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
head(ident_0)
ident_1 <- FindConservedMarkers(combined.sct, ident.1 = 1, grouping.var = "orig.ident", verbose = T)
ident_2 <- FindConservedMarkers(combined.sct, ident.1 = 2, grouping.var = "orig.ident", verbose = T)
ident_3 <- FindConservedMarkers(combined.sct, ident.1 = 3, grouping.var = "orig.ident", verbose = T)
ident_4 <- FindConservedMarkers(combined.sct, ident.1 = 4, grouping.var = "orig.ident", verbose = T)
ident_5 <- FindConservedMarkers(combined.sct, ident.1 = 5, grouping.var = "orig.ident", verbose = T)
ident_6 <- FindConservedMarkers(combined.sct, ident.1 = 6, grouping.var = "orig.ident", verbose = T)
ident_7 <- FindConservedMarkers(combined.sct, ident.1 = 7, grouping.var = "orig.ident", verbose = T)
ident_8 <- FindConservedMarkers(combined.sct, ident.1 = 8, grouping.var = "orig.ident", verbose = T)
ident_9 <- FindConservedMarkers(combined.sct, ident.1 = 9, grouping.var = "orig.ident", verbose = T)
ident_10 <- FindConservedMarkers(combined.sct, ident.1 = 10, grouping.var = "orig.ident", verbose = T)
ident_11 <- FindConservedMarkers(combined.sct, ident.1 = 11, grouping.var = "orig.ident", verbose = T)
ident_12 <- FindConservedMarkers(combined.sct, ident.1 = 12, grouping.var = "orig.ident", verbose = T)
ident_13 <- FindConservedMarkers(combined.sct, ident.1 = 13, grouping.var = "orig.ident", verbose = T)
ident_14 <- FindConservedMarkers(combined.sct, ident.1 = 14, grouping.var = "orig.ident", verbose = T)

## feature plot
DefaultAssay(combined.sct) <- "RNA"
FeaturePlot(combined.sct, features = c("TOP2A", "PCNA", "MKI67", 
                                       "FABP1","KRT20", "TFF3" ,"TFF1",
                                       "OLFM4", "LGR5",
                                       "LYZ", "CD24",
                                       "CD79A", "MZB1", "CD27",
                                       "PECAM1", "VWF", "CD34", "CLDN5",
                                       "COL1A1","MMP2","DCN","COL1A2"),
            reduction = "tsne",cols = c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))


## annotation
new.cluster.ids <- c("TA","Enterocyte/Goblet","Stem","Paneth","Plasma B","Plasma B",
                     "Stem","Plasma B","Stem","MT","MT","KRT20+","Endothelial","CAF","Endothelial")
names(new.cluster.ids) <- levels(combined.sct)
combined.sct <- RenameIdents(combined.sct, new.cluster.ids)
DimPlot(combined.sct, label = F)

## remove mt cells
combined.sct = subset(combined.sct, 
             idents = unique(Idents(combined.sct))[!unique(Idents(combined.sct))=="MT"])
DimPlot(combined.sct, reduction = "tsne", label = F, 
        cols = c("#003049","#fcbf49","#eae2b7","#A098D6",
                 "#F25F53","#D3CF38","#92CD8B","#93AED5"))


## view conserved cell type markers across conditions
markers.to.plot <- c("TOP2A", "PCNA", "MKI67", 
                     "FABP1","KRT20", "TFF3" ,"TFF1",
                     "OLFM4", "LGR5",
                     "LYZ", "CD24",
                     "CD79A", "MZB1", "CD27",
                     "PECAM1", "VWF", "CD34", "CLDN5",
                     "COL1A1","MMP2","DCN","COL1A2")
DotPlot(combined.sct, features = markers.to.plot, cols = c("#D8383A","#14517c"), dot.scale = 7, split.by = "orig.ident") +
  RotatedAxis()
saveRDS(combined.sct, "../RDS data/combined_sct_annotation.rds")

## cell group comparison
plot.data <- combined.sct@meta.data %>%
  select(seurat_clusters, orig.ident) %>%
  mutate(tissue = fct_collapse(seurat_clusters,
                               TA = c("0"),
                               Enterocyte_Goblet = c("1"),
                               Stem = c("2","6","8"),
                               Paneth = c("3"),
                               Plasma_B = c("4","5","7"),
                               KRT20 = c("11"),
                               Endothelial = c("12","14"),
                               CAF = c("13"),
                               Others = c("9","10"))) %>%
  mutate(tissue = fct_relevel(tissue,
                              "TA", "Enterocyte/Goblet", "Stem",
                              "Paneth", "Plasma B", "KRT20+",
                              "Endothelial","CAF",
                              "Others")) %>%
  group_by(orig.ident, tissue) %>%
  summarise(Count = n()) %>%
  mutate(Prop = Count / sum(Count))
 
f2E <- ggplot(plot.data, aes(x = orig.ident, y = Prop, fill = tissue)) +
  geom_col() +
  labs(y = "Proportion of cells") +
  guides(fill = guide_legend(direction = "vertical", ncol = 2)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid =element_blank(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.text.x = element_text( hjust = 1,angle = 45),
        axis.title = element_text(color = "black", size = 15),
        plot.title = element_text(size=15,hjust=0.5, color = 'black'))+
  scale_fill_manual(values=c("#003049","#eae2b7","#A098D6","#92CD8B",
                             "#93AED5","#fcbf49","#F25F53","#D3CF38"))
f2E

