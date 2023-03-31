## Trajectory Inference using velocyto
## Convert data from Seurat to anndata format
## load data
org = readRDS("../RDS data/organoid_annotation.rds")
org$barcode <- colnames(org)
org$umap_1 <- org@reductions$umap@cell.embeddings[,1]
org$umap_2 <- org@reductions$umap@cell.embeddings[,2]
org$celltype <- org@active.ident
write.csv(org@meta.data, file = "../table/metadata.csv",quote=F, row.names=F)

##write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(org, assay='RNA', slot='counts')
writeMM(counts_matrix, file=('../table/counts.mtx'))

## write dimesnionality reduction matrix (pca & umap)
write.csv(org@reductions$pca@cell.embeddings, file='../table/pca.csv', quote=F, row.names=F)
write.csv(org@reductions$umap@cell.embeddings, file='../table/umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='../table/gene_names.csv',
  quote=F,row.names=F,col.names=F
)
