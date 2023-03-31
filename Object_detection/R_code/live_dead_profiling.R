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

## load data
live_dead_df = read.csv("../table/live_dead_profiling.csv")
organoids_morphology_drug =  read.csv("../table/organoids_morphology_drug.csv") 
umap_df_for_1960_jpg = read.csv("../table/umap_df_for_1960_jpg.csv")
  
## data normalization
live_dead_df = live_dead_df %>% mutate(success_prob = success_score / success_num)   %>%
  mutate(failure_prob = failure_score / failure_num) %>%
  mutate(live_prob = success_score - failure_score) 
live_dead_df[is.na(live_dead_df)] <- 0
live_dead_df$live_prob =  (live_dead_df$live_prob - min(live_dead_df$live_prob)) / (max(live_dead_df$live_prob) - min(live_dead_df$live_prob))
  

df <- umap_df_for_1960_jpg %>%
  left_join(organoids_morphology_drug,by = 'Field') %>%
  left_join(live_dead_df, by = "Field")

df <- mutate(df, Size_log10 = log10(Org_size))
df <- mutate(df, Size_ln = log(Org_size))

## UMAP live_dead probability distribution
df %>%
  ggplot(aes(UMAP1, UMAP2, color = live_prob)) + 
  geom_point(alpha = 0.8, size = 1) +
  scale_color_viridis_c() +
  theme_cowplot() + 
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "live_prob")


## UMAP live_dead probability distribution with 3 drugs in CK07
df_drug <- df  %>% filter(Line_no_passage =="CK07") %>% 
  filter(Drug == "FOLFOX") %>%
  filter(Morphology != "intermediate")
gg_drug = df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 2, color = "#f1f1f1") + 
  geom_point_rast(data = df_drug,
                  aes(color = live_prob),size =2, alpha = 0.8, shape =16) + 
  scale_color_viridis_c() +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "live_prob of FOLFOX") +
  theme_cowplot(font_size = 8) 
gg_drug

## organoid viability vs organoid sizes
library(ggstatsplot)
org_via_size <- df %>% group_by(df$Field) %>%
  summarise(count = n(),
            Mean_size = mean(Size_log10),
            Median_size = median(Size_log10),
            Mean_live_prob = mean(live_prob)) %>% na.omit()

ggscatterstats(
  data = org_via_size, ##%>% filter(Drug == "FOLFOX"),
  x = Median_size,
  y = Mean_live_prob,
  bf.message = FALSE,
  ##title = "FOLFOX"
)
cor(org_via_size$Median_size, org_via_size$Mean_live_prob)


## organoid viability vs organoid morphology
org_via_morph <- df %>% group_by(df$Field) %>%
  summarise(count = n(),
            Org_morphology = Morphology,
            Mean_live_prob = mean(live_prob),
            .groups = "drop" ) %>%
  unique() %>% na.omit() 
quantile(org_via_morph$Mean_live_prob)
org_via_morph <- mutate(org_via_morph, Quantile = case_when(0.0<=Mean_live_prob & Mean_live_prob< 0.2467192 ~"Q1",
                                                         0.2467192<=Mean_live_prob & Mean_live_prob<0.2673416 ~ "Q2",
                                                         0.2673416<=Mean_live_prob & Mean_live_prob<0.3072741 ~ "Q3",
                                                         0.3072741<=Mean_live_prob & Mean_live_prob<=1.0000000 ~ "Q4"
                                                         ))

# write.csv(org_via_morph,"../table/org_via_morph.csv")
# library(readr)
org_via_morph_statistics <- read_csv("../table/org_via_morph_statistics.csv") 
org_via_morph_statistics<-org_via_morph_statistics[,-1]
rownames(org_via_morph_statistics)<-c("Q1","Q2","Q3","Q4")
library(pheatmap)
col = colorRampPalette(c("navy", "white", "firebrick3"))(200)

anno_row = data.frame(
  Live_prob = rep(c("Quartile_1", "Quartile_2","Quartile_3","Quartile_4"), each = 1),
  row.names = rownames(org_via_morph_statistics))
anno_colors = list(
  Live_prob = c(Quartile_1 = "#97D16A", 
                Quartile_2 = "#C4E2EE",
                Quartile_3 = "#EFA709",
                Quartile_4 = "#D40000"))

pheatmap(scale(org_via_morph_statistics,center=FALSE),
         treeheight_row = 0, treeheight_col = 0,
         color = col,
         angle_col = 45,
         annotation_row = anno_row,
         annotation_colors = anno_colors,
         border_color = "white",
         cellwidth = 20, cellheight = 20,
         cluster_rows = FALSE, cluster_cols = FALSE)
