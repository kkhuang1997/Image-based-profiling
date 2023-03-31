library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(tibble)
library(readr)
library(here)
library(ggrastr)
library(cowplot)
library(princurve)
library(scico)
library(ggridges)
library(dplyr)
library(viridis)

# live_dead_df = read.csv("../live_dead_profiling.csv")
# organoids_morphology_drug =  read.csv("../../image_based profiling/table/organoids_morphology_drug.csv") 
# umap_df_for_1960_jpg = read.csv("../../image_based profiling/table/umap_df_for_1960_jpg.csv")
  
# 
# live_dead_df = live_dead_df %>% mutate(success_prob = success_score / success_num)   %>%
#   mutate(failure_prob = failure_score / failure_num) %>%
#   mutate(live_prob = success_score - failure_score) 
# live_dead_df[is.na(live_dead_df)] <- 0
# ## data normalization
# live_dead_df$live_prob =  (live_dead_df$live_prob - min(live_dead_df$live_prob)) / (max(live_dead_df$live_prob) - min(live_dead_df$live_prob))
#   
# 
# df <- umap_df_for_1960_jpg %>%
#   left_join(organoids_morphology_drug,by = 'Field') %>%
#   left_join(live_dead_df, by = "Field")
# 
# df <- mutate(df, Size_log10 = log10(Org_size))
# df <- mutate(df, Size_ln = log(Org_size))
# 
# ## UMAP live_dead probability distribution
# df %>%
#   ggplot(aes(UMAP1, UMAP2, color = live_prob)) + 
#   geom_point(alpha = 0.8, size = 1) +
#   scale_color_viridis_c() +
#   theme_cowplot() + 
#   labs(x = "UMAP 1",
#        y = "UMAP 2",
#        color = "live_prob")
# 
# 
# ## UMAP live_dead probability distribution with drug treatment
# df_drug <- df  %>% filter(Line_no_passage =="CK07")
# df_drug %>%
#   ggplot(aes(UMAP1, UMAP2, color = live_prob)) + 
#   geom_point(alpha = 0.8, size = 1) +
#   scale_color_viridis_c() +
#   theme_cowplot() + 
#   labs(x = "UMAP 1",
#        y = "UMAP 2",
#        color = "live_prob")
# 
# ## organoid viability vs organoid sizes
# library(ggstatsplot)
# org_via_size <- df %>% group_by(df$Field) %>%
#   summarise(count = n(),
#             Median_size = median(Size_log10),
#             Mean_live_prob = mean(live_prob)) %>% na.omit()
# 
# ggscatterstats(
#   data = org_via_size, ##%>% filter(Drug == "FOLFOX"),
#   x = Median_size,
#   y = Mean_live_prob,
#   bf.message = FALSE,
#   ##title = "FOLFOX"
# )
# cor(org_via_size$Median_size, org_via_size$Mean_live_prob)


library(readxl)
apoptosis_morphology <- read_excel("E:/classical feature extraction/table/meta_apoptosis_normed_intensity_morphology_norm_intensity_per_line_withoutNA.xls")

org_apop_morph = apoptosis_morphology %>% group_by(apoptosis_morphology$Field) %>%
  summarise(count = n(),
            Org_morphology = Morphology,
            Apoptosis_intensity = Norm_mean,
            Apoptosis_intensity_perline = Norm_mean_per_line) %>%
  unique() %>% na.omit()
quantile(org_apop_morph$Apoptosis_intensity)
org_apop_morph <- mutate(org_apop_morph, Quantile = case_when(0.112107<=Apoptosis_intensity & Apoptosis_intensity< 2.174184 ~"Q1",
                                                              2.174184<=Apoptosis_intensity & Apoptosis_intensity<2.697279 ~ "Q2",
                                                              2.697279<=Apoptosis_intensity & Apoptosis_intensity<3.221439 ~ "Q3",
                                                              3.221439<=Apoptosis_intensity & Apoptosis_intensity<=7.820739 ~ "Q4"))

write.csv(org_apop_morph,"../table/org_apop_morph_withoutNA.csv")
org_apop_morph_statistics <- read_csv("../table/org_apop_morph_withoutNA_statistics.csv") 
org_apop_morph_statistics<-org_apop_morph_statistics[,-1]
rownames(org_apop_morph_statistics)<-c("Q1","Q2","Q3","Q4")
library(pheatmap)
col = colorRampPalette(c("navy", "white", "firebrick3"))(200)

anno_row = data.frame(
  Apoptosis_intensity = rep(c("Quartile_1", "Quartile_2","Quartile_3","Quartile_4"), each = 1),
  row.names = rownames(org_apop_morph_statistics))
anno_colors = list(
  Apoptosis_intensity = c(Quartile_1 = "#97D16A", 
                Quartile_2 = "#C4E2EE",
                Quartile_3 = "#EFA709",
                Quartile_4 = "#D40000"))

pheatmap(scale(org_apop_morph_statistics,center=FALSE),
         treeheight_row = 0, treeheight_col = 0,
         color = col,
         angle_col = 45,
         annotation_row = anno_row,
         annotation_colors = anno_colors,
         border_color = "white",
         cellwidth = 20, cellheight = 20,
         cluster_rows = FALSE, cluster_cols = FALSE)
