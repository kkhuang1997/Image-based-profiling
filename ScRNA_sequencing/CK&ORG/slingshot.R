if(!require(slingshot))install.packages("slingshot")
if(!require(BUSpaRse))install.packages("BUSpaRse")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(tidymodels))install.packages("tidymodels")
if(!require(Seurat))install.packages("Seurat")
if(!require(scales))install.packages("scales")
if(!require(viridis))install.packages("viridis")
if(!require(Matrix))install.packages("Matrix")
# library(dyno)

##load seurat obj
org <- readRDS("../RDS data/organoid_annotation.rds")
# org <- FindNeighbors(org, dims = 1:13)
# org <- FindClusters(org, resolution = 0.4)
cols = c("#003049","#d62828","#A098D6","#fcbf49","#eae2b7","#81b29a","#209A90","#D4B7AF")
DimPlot(org, reduction = "umap", cols = cols)

## dyno Slingshot wrapper
dataset <- wrap_expression(
  counts = t(org@assays$RNA@counts),
  expression = t(org@assays$RNA@data)
)

## set stem cell as starting point
start_cells = names(org@active.ident[org@active.ident == "Stem"])
dataset <- add_prior_information(
  dataset,
  start_id = start_cells
)

## add cluster info
dataset <- add_grouping(
  dataset,
  org@active.ident
)

## add dimensionality reduction
dataset <- dynwrap::add_dimred(
  dataset,
  org@reductions$umap@cell.embeddings
)

## trajectory inferring
guidelines <- guidelines_shiny(dataset)
methods_selected <- c("slingshot", "paga_tree", "paga","scorpius")
model_slingshot <- infer_trajectory(dataset, methods_selected[1])
plot_dimred(
  model_slingshot,
  expression_source = dataset$expression, 
  dimred = dataset$dimred,
  grouping = dataset$grouping,
  color_density = "grouping")

plot_heatmap(
  model_slingshot,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)

plot_dimred(model_slingshot, 
            "pseudotime", 
            dimred = dataset$dimred,
            expression_source = dataset$expression,
            grouping = dataset$grouping,
            pseudotime = calculate_pseudotime(model_slingshot))
saveRDS(model_slingshot, "../RDS data/model_slingshot.rds")


## single slingshot trajectory inferring
org_sce <- as.SingleCellExperiment(org)
org$umap_1 <- org@reductions$umap@cell.embeddings[,1]
org$umap_2 <- org@reductions$umap@cell.embeddings[,2]
rd_umap <- c(org$umap_1, org$umap_2)
rd_umap <- matrix(rd_umap,ncol = 2)
colnames(rd_umap) <- c('UMAP1', 'UMAP2')
plot(rd_umap, col = rgb(0,0,0,.5), pch=20, asp = 1)
reducedDims(org_sce) <- SimpleList(UMAP = rd_umap)

## GMM cluster (9 clusters)
library(mclust)
library(RColorBrewer)
cl1 <- Mclust(rd_umap)$classification
colData(org_sce)$GMM <- cl1
plot(rd_umap, col = brewer.pal(12,"Set3")[cl1], pch=20, asp = 1)
## k-means cluster (8 clusters)
cl2 <- kmeans(rd_umap, centers = 8)$cluster
colData(org_sce)$KMEANS <- cl2
plot(rd_umap, col = brewer.pal(12,"Set3")[cl2], pch=20, asp = 1)
## original seurat clusters
colData(org_sce)$SEU_CLUSTER <- as.double(org@meta.data[["seurat_clusters"]])
cl3 = colData(org_sce)$SEU_CLUSTER
plot(rd_umap, col = brewer.pal(12,"Set3")[cl3], pch=20, asp = 1)

##using slingshot
sds <- slingshot(org_sce, clusterLabels = "SEU_CLUSTER", reducedDim = "UMAP")

library(grDevices)
library(cowplot)
## global trajectory
plot(reducedDims(sds)$UMAP, col = cols[sds$SEU_CLUSTER], pch=20, asp = 1.5)
lines(SlingshotDataSet(sds), lwd=2, type = 'lineages', col = 'black')
## lineage trajectory (3 lineages)
umap_df = as.data.frame(reducedDims(sds)$UMAP)
umap_df %>%
  ggplot(aes(UMAP1, UMAP2, color = sds$slingPseudotime_3)) + 
  geom_point(alpha = 0.5, size = 0.35) +
  scale_color_viridis_c(option = "C") +
  theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "Pseudotime")
## save slingshot object
saveRDS(sds, "../RDS data/slingshot_object.rds")
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sds$slingPseudotime_1, breaks=100)]
# plot(reducedDims(sds)$UMAP, col = plotcol, pch=20, asp = 1.5)
# lines(SlingshotDataSet(sds), lwd=2, col = 'black')


sds = readRDS("../RDS data/slingshot_object.rds")
## Find temporally expressed genes using TradeSeq
library(gam)
Y <- log2(counts(sds) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
t <- sds$slingPseudotime_2
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]  

# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.
library(Polychrome)
library(RColorBrewer)
heatdata <- assays(sds)$counts[topgenes, order(t, na.last = NA)]
heatdata <- as.matrix(heatdata)
#heatclus <- (sds@colData@listData[["UMAP_CLUSTER"]])[order(t, na.last = NA)]
heatclus <- sds$GMM[order(t, na.last = NA)]

#my_color <- createPalette(8, c("#14517c","#E7EFFA","#96C37D","#F3D266","#D8383A"), M=1000)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(12,"Set3")[heatclus],col = coul)


sds = readRDS("../RDS data/slingshot_object.rds")
org <- readRDS("../RDS data/organoid_annotation.rds")
## Find temporally expressed genes using a random forest model
# Get top highly variable genes
top_hvg <- HVFInfo(SCTransform(org)) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(300, residual_variance) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(org, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,1], dat_use) # using lineage 1
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
# train vs. validation
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)
# validation using RMSE, RSQ, and MAE
val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)
metrics(data = val_results, truth, estimate)

## plot genes deemed the most important to predicting pseudotime
var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[81:100]
par(mfrow = c(4, 5)) ## set layout
pal <- plasma(100, end = 0.95) ## set color
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(reducedDim(sds), col = colors, 
       pch = 16, cex = 0.5, main = top_genes[i])
  lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')
}

