if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")

## load data
solid.data <- Read10X(data.dir = "E:/苏州研究院/单细胞测序结果/result/data/solid_Matrix/")
solid <- CreateSeuratObject(counts = solid.data, project = "solid")

##QC and selecting cells for further analysis
solid[["percent.mt"]] <- PercentageFeatureSet(solid, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(solid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
solid <- subset(solid, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 50)
## Normalizing the data
solid <- NormalizeData(solid, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features
solid <- FindVariableFeatures(solid, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(solid), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(solid)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

## Scaling the data
all.genes <- rownames(solid)
solid <- ScaleData(solid, features = all.genes,
                   model.use = "negbinom")

## PCA
solid <- RunPCA(solid, features = VariableFeatures(object = solid))
DimHeatmap(solid, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(solid)

## Cluster the cells
solid <- FindNeighbors(solid, dims = 1:20, k.param = 70)
solid <- FindClusters(solid, resolution = 0.8)

## plot tsne
solid <- RunTSNE(solid, dims = 1:20)
DimPlot(solid, reduction = "tsne")
        # cols = c("#bfc2e1","#609278","#8e787b","#b1c75f","#f2d900","#e6a1c6","#68b0d6"))

## cluster identity
diff.wilcox = FindAllMarkers(solid)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(all.markers, "../table/SOLID&CYSTIC/diff_genes_solid.csv", row.names = F)
write.csv(top100, "../table/SOLID&CYSTIC/diff_genes_top100_solid.csv", row.names = F)
# saveRDS(solid,"../../output/solid_cystic/solid.rds")

FeaturePlot(cp, features =c("DEFA6", "DEFA5", ## Paneth
                            "TFF1","TFF3",  ## Goblet
                            
                            ),
            reduction = "tsne", cols =c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))

FeaturePlot(solid, features ="BEST",
            reduction = "tsne", cols =c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"))

## cell annotation
library(ComplexHeatmap)
new.cluster.ids <- c("CRC_related Sig", "Enterocyte/goblet", "MT", "EPCAM+", "Paneth", 
                     "TA1","TA2")
names(new.cluster.ids) <- levels(solid)
solid <- RenameIdents(solid, new.cluster.ids)
DoHeatmap(solid, features = top100$gene,
          group.colors = c("#bfc2e1","#609278","#8e787b","#b1c75f","#f2d900","#e6a1c6","#68b0d6"),
          size = 2) + 
  scale_fill_gradientn(colours =c("navy","white","firebrick3"))

# saveRDS(solid,"../../output/solid_cystic/solid_with_annotation.rds")


## gesa
library(tibble)
library(ggplot2)
library(purrr)
library(fgsea)
gene_diff <- read.delim('../../output/cell_identity/diff_genes_solid.csv', sep = ',', check.names = FALSE)
geneList = gene_diff$avg_log2FC
names(geneList) = as.character(gene_diff$gene)
geneList=sort(geneList,decreasing = T)
load("./csc_signatures.rda")
custom_barcode_plot <- function(df, sig){
  ## named vector of gene-level stats
  stat_vector <- setNames(df$avg_log2FC, df$gene)
  ## genes in signature
  sig_genes <- signatures[[sig]]
  
  ## generate barcode plot
  bc_plot <- plotEnrichment(sig_genes, stat_vector)
  
  ## remove unwanted layers
  bc_plot$layers <- list()
  
  ## add barcode at the bottom
  lowest_pos <- min(bc_plot$data[,2])
  dash_length <- abs(reduce(range(bc_plot$data[,2]), `-`)*0.1)
  middle <- which.min(abs(sort(df$avg_log2FC, decreasing=T)))
  
  bc_plot_custom <- bc_plot + geom_segment(aes(x=x, xend=x), y=lowest_pos,
                                           yend=lowest_pos-dash_length) + 
    geom_line(colour='#4daf4a') + 
    geom_hline(yintercept=lowest_pos, colour='#cccccc') + 
    geom_hline(yintercept=0, colour='#cccccc') + xlab('') +
    theme_classic() +
    geom_tile(data=tibble(rank=1:length(stat_vector), 
                          y=lowest_pos-(1.25*dash_length)), 
              aes(x=rank, y=y, fill=rank),
              width=1,
              height=0.5*dash_length) +
    scale_fill_gradient2(low ='#b2182b', high='#2166ac', 
                         mid='#f7f7f7', midpoint = middle) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none') + 
    ggtitle(paste(sig, 'signature')) +
    ylab('Enrichment score')
  
  return(bc_plot_custom)
}
solid_gsea <- fgsea(pathways = signatures,
                  stats = geneList, nperm=100000)

bc_plots_affy <- map(1:nrow(solid_gsea), function(j){
  bcp <- custom_barcode_plot(gene_diff, solid_gsea$pathway[j]) + 
    annotate('text', x=Inf , y=Inf, hjust=1, vjust=1, 
             label=paste('NES =', round(solid_gsea$NES[j], 2), 
                         '\nFDR =', solid_gsea$padj[j]))
  
  return(bcp)
})
## plot to canvas
reduce(rev(bc_plots_affy), `+`) 
