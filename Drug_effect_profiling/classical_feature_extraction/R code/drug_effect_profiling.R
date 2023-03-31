library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(here)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)
library(ggridges)
library(dplyr)
library(feather)
library(arrow)
library(readxl)
library(pheatmap)
source("./utils.R")
library(SCOPEAnalysis)

## load and prepare svm data
root = "../table/drug_effect_red/"
Lines = c("CK33P02","CK38P01","CK40P01","CK44P01","CK51P06")
drug_effect_metadata = vector(mode = "list", length = length(Lines))
drug_effect_profiles = vector(mode = "list", length = length(Lines))
for (idx in seq_along(Lines)) {
  drug_effect_metadata[[idx]] <- read_csv(file.path(root,paste0("SVM_Accuracies_Ipca_",Lines[idx],".csv")))
  drug_effect_profiles[[idx]] <- read_csv(file.path(root,paste0("SVM_Profiles_mean_Ipca_",Lines[idx],".csv")))
  
  colnames(drug_effect_metadata[[idx]])[1] = 'Drug'
  colnames(drug_effect_profiles[[idx]])[1] = 'Drug'
  
  drug_effect_metadata[[idx]] = mutate(drug_effect_metadata[[idx]], "Line" = Lines[idx])
  drug_effect_profiles[[idx]] = mutate(drug_effect_profiles[[idx]], "Line" = Lines[idx])
  
  rownames(drug_effect_metadata[[idx]]) = paste(drug_effect_metadata[[idx]]$Line,
                                              drug_effect_metadata[[idx]]$Drug, 
                                              sep = ".")
  rownames(drug_effect_profiles[[idx]]) = paste(drug_effect_profiles[[idx]]$Line,
                                                drug_effect_profiles[[idx]]$Drug, 
                                                sep = ".")
  
}

drug_effect_metadata = do.call(what = rbind, args = drug_effect_metadata)
drug_effect_profiles = do.call(what = rbind, args = drug_effect_profiles)

## save processed data
write.csv(drug_effect_metadata, "../table/drug_effect_red/drug_effect_metadata_Ipca.csv")
write.csv(drug_effect_profiles, "../table/drug_effect_red/drug_effect_profiles_Ipca.csv")


## reload data
drug_effect_metadata <- read.csv("../table/drug_effect_blue/drug_effect_metadata_Ipca.csv",row.names = 1)
drug_effect_profiles <- read.csv("../table/drug_effect_blue/drug_effect_profiles_Ipca.csv",row.names = 1)

####
#### Define Active Drugs
####

ggplot(data = drug_effect_metadata, mapping = aes(x = AUC_Mean)) + 
  geom_histogram(bins = 50) + 
  # theme_vignette() +
  ggtitle("AUROC Histogram") + 
  xlab("Mean AUROC") + ylab("")

## distance between drugs and DMSO correlates with the AUROC
ggplot(data = drug_effect_metadata, 
       mapping = aes(x = AUC_Mean, y = Distance)) + 
  geom_point(size = 0.5) + 
  # theme_vignette() + 
  labs(
    caption = paste("Spearman rho =", round(
      cor(drug_effect_metadata$AUC_Mean, 
          drug_effect_metadata$Distance, method = "spearman"), 2))) + 
  xlab("AUROC") + labs(color = "") + 
  geom_smooth()

ggplot(data = drug_effect_metadata, 
       mapping = aes(x = AUC_Mean, y = AUC_Std)) + 
  geom_point(size = 0.5) + 
  labs(x = "AUROC Mean", y = "AUROC Standard Deviation", color = "") + 
  geom_smooth()


## active vs. inactive drugs
auc_thresh = 0.85
profiles_active = drug_effect_profiles[drug_effect_metadata$AUC_Mean >= auc_thresh, ]
metadata_active = drug_effect_metadata[drug_effect_metadata$AUC_Mean >= auc_thresh, ]
## No. of active drugs varies per line
ggplot_df = as.data.frame(table(metadata_active$Line))
ggplot_df$TextLoc = ggplot_df$Freq - 0.8
ggplot(data = ggplot_df, mapping = aes(x = Var1, y = Freq)) + 
  geom_col() + 
  geom_text(mapping = aes(x = Var1, y = TextLoc, label = Freq), 
            color = "white", fontface = "bold", size = 3) + 
  # theme_vignette() + 
  theme(axis.line = element_blank(), axis.text.x = element_blank(), 
        axis.ticks = element_blank()) + 
  xlab("") + ylab("") + ggtitle("Number of Active Drugs") + 
  coord_flip()


#####
##### heatmap of all drug effect profiles (include inactive) across all lines
#####
## cosine similarity
drug_effect_profiles = data.matrix(drug_effect_profiles)
angles = get_angles(drug_effect_profiles)
d = acos(angles) * 180 / pi
hc = hclust(as.dist(d), method = "median")

## annotation
layout <- read_csv("../table/layout.csv")
pathway_target = as.data.frame(t(get_target_pathway(drug_effect_profiles)))
pathway_target=as.data.frame(lapply(pathway_target,as.character))
row.names(pathway_target) = rownames(drug_effect_profiles)

annotation <- read.csv("../table/heatmap_annotation.csv", row.names = 1)


## set color 
hm_colorscale = colorRampPalette(c("navy", "white", "firebrick3"))(200)

anno_colorScale = list(Target = sort(unique(annotation$Target)),
                       Pathway = sort(unique(annotation$Pathway)),
                       Morphology = sort(unique(annotation$Morphology)))

anno_colorScale$Target = setNames(
  object = c("#C02C38", "#DE7897", "#7E1671", "#61649F", "#619AC3", 
             "#2BAE85", "#BACF65", "#FBCD31", "#9A8878", "#F05A46", 
             "#008080", "#e6beff", "#aa6e28", "#fffac8", "#800000", 
             "#aaffc3"), 
  nm = c(anno_colorScale$Target))

anno_colorScale$Pathway = setNames(
  object = c("#FAD2D2", "#EC6183", "#B466CB", "#53A8A8", "#99D698"),
             # "#E0F77A", "#54AAD7"),
  nm = c(anno_colorScale$Pathway)
)

anno_colorScale$Morphology = setNames(
  object = c("#E77882", "#EBC94E", "#417D9F"),
  nm = c(anno_colorScale$Morphology)
)

pheatmap(d, annotation_row = annotation,
         color = hm_colorscale,
         annotation_colors = anno_colorScale,
         show_rownames = F, show_colnames = F,
         cluster_rows = hc, cluster_cols = hc,
         border_color = NA)

