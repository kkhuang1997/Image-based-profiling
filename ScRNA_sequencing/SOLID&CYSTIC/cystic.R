if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")

## read  genes, barcode, matrix
cp.data <- Read10X(data.dir = "../../data/surCP_matrix/")
cp <- CreateSeuratObject(counts = cp.data, project = "surCP")


##QC and selecting cells for further analysis
cp[["percent.mt"]] <- PercentageFeatureSet(cp, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(cp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cp <- subset(cp, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 70)
# cp <- subset(cp, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 75)
## Normalizing the data
cp <- NormalizeData(cp, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features
cp <- FindVariableFeatures(cp, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cp), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cp)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot1 + plot2

## Scaling the data
all.genes <- rownames(cp)
cp <- ScaleData(cp, features = all.genes, 
                model.use = "negbinom", 
                vars.to.regress = "percent.mt")

cp <- FindVariableFeatures(cp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cp), 10)
## PCA
cp <- RunPCA(cp, features = VariableFeatures(object = cp))
DimHeatmap(cp, dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(cp)

## Cluster the cells
cp <- FindNeighbors(cp, dims = 1:20, k.param = 150)
cp <- FindClusters(cp, resolution = 0.3)

## plot umap,tsne
# cp <- RunUMAP(cp, dims = 1:15)
# DimPlot(cp, reduction = "umap")
cp <- RunTSNE(cp, dims = 1:25)
DimPlot(cp, reduction = "tsne", 
        cols = c("#8e787b","#b1c75f","#f2d900","#e6a1c6","#68b0d6"))

## cluster identity
diff.wilcox = FindAllMarkers(cp)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top100 = all.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(all.markers, "../../output/cell_identity/diff_genes_cys.csv", row.names = F)
write.csv(top100, "../../output/cell_identity/diff_genes_top100_cys.csv", row.names = F)

## cell annotation
library(ComplexHeatmap)
new.cluster.ids <- c("Immune-related Sig", "Stem" ,"Enterocyte/goblet", "Paneth", 
                     "Wnt-related Sig")
names(new.cluster.ids) <- levels(cp)
cp <- RenameIdents(cp, new.cluster.ids)
DoHeatmap(cp, features = top100$gene,
          group.colors = c("#8e787b","#b1c75f","#f2d900","#e6a1c6","#68b0d6"),
          size = 2) + 
  scale_fill_gradientn(colours =c("navy","white","firebrick3"))

saveRDS(cp,"../../output/solid_cystic/cystic.rds")


## gesa
library(tibble)
library(ggplot2)
library(purrr)
library(fgsea)
gene_diff <- read.delim('../table/SOLID&CYSTIC/diff_genes_cys.csv', sep = ',', check.names = FALSE)
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

cys_gsea <- fgsea(pathways = signatures,
                    stats = geneList, nperm=10000)

## generate plots
bc_plots_affy <- map(1:nrow(cys_gsea), function(j){
  bcp <- custom_barcode_plot(gene_diff, cys_gsea$pathway[j]) + 
    annotate('text', x=Inf , y=Inf, hjust=1, vjust=1, 
             label=paste('NES =', round(cys_gsea$NES[j], 2), 
                         '\nFDR =', cys_gsea$padj[j]))
  
  return(bcp)
})
## plot to canvas
reduce(rev(bc_plots_affy), `+`) 
