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
library(feather)

##  load umap data
umap_df = arrow::read_feather("../python code/intermediate_data/StageFour_umap_meta.feather")

## umap glance
umap_df %>% 
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 0.5, size = 0.35) + 
  scale_color_brewer(type = "qual", palette = "Set2") +
  theme_cowplot() +
  labs(x = "UMAP 1",
       y = "UMAP 2") + 
  theme(legend.position = "bottom") + 
  coord_fixed()

# Organoid Size Distributions
gg_size_dist <- umap_df %>% 
  ggplot(aes(Size)) + 
  geom_histogram(bins = 30) + 
  theme_cowplot() 
gg_size_dist_log <- gg_size_dist + 
  scale_x_log10()
print(gg_size_dist_log)


## Size_dist according lines
umap_df <- mutate(umap_df, Size_log10 = log10(Size))
umap_df <- mutate(umap_df, Size_ln = log(Size))
organoid_size_factor_09 <- umap_df %>% group_by(Line) %>% 
  summarise(x = quantile(Size_ln, 0.9)) %>% 
  #summarise(x = mean(size_log)) %>% 
  arrange(x) %>% .$Line

size_dist_morph_ridge <- umap_df %>% 
  mutate(Line = factor(Line, levels = organoid_size_factor_09)) %>% 
  ggplot() +
  geom_density_ridges_gradient(aes(y = Line, x = Size_log10, fill = stat(x)), scale = 1) +
  #geom_density(aes(x = size_log, group = replicate, color = morphological_class)) + 
  #facet_wrap(~ line) + 
  scale_fill_viridis_c() +
  labs(x = "log(size)",
       fill = "size") + 
  theme(legend.position = "bottom") +
  theme_cowplot() + 
  coord_fixed(ratio = 1)
size_dist_morph_ridge


## UMAP size distribution
umap_df %>%
  ggplot(aes(UMAP1, UMAP2, color = Size_log10)) + 
  geom_point(alpha = 0.5, size = 0.35) +
  scale_color_viridis_c() +
  theme_cowplot() + 
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "log10(size)")

## 5 Organoid line distribution (no evidence)
loi <- c("CK40P01", "CK51P06", "CK38P01", "CK33P02", "CK44P01")
gg_line <- umap_df %>% dplyr::select(-Line) %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  geom_point_rast(data = umap_df %>%
                    filter(Line %in% loi) %>% 
                    mutate(Line = factor(Line, levels = loi)) %>% 
                  sample_frac(0.5),
                  aes(color = Line),alpha = .4, size = 0.35, shape=16) + 
  facet_wrap( ~ Line, ncol = 2) +
  # scale_color_brewer(type = "qual", palette = "Set2") +
  scale_color_manual(values = c("#d43232","#aaa4d6","#81b29a","#657f8d","#726FFB")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  theme_cowplot(font_size = 8) + 
  theme(legend.position = "nothing")
gg_line

## organoid_morphology
library(readxl)
organoid_morphology <- read_excel("../table/meta_apoptosis_normed_intensity_morphology_norm_intensity_per_line_withoutNA.xls")
df <- umap_df %>%
  left_join(organoid_morphology,by = 'Field') %>% na.omit()

gg_morph <- df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>% 
                  filter(Morphology != "intermediate") %>%
                  sample_frac(1),
                  aes(color = Morphology), size = 0.35, shape=16) + 
  # scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = c("#fa7f6f", "#beb8dc", "#82b0d2")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  theme_cowplot(font_size = 8)
gg_morph


## organoid_morphology in 5 lines 
## random sampling 2000 points
df_CK51P06 <- df %>% filter(Line.x == "CK51P06") %>% .[sample(nrow(.), size = 2000), ]
df_CK33P02 <- df %>% filter(Line.x == "CK33P02") %>% .[sample(nrow(.), size = 2000), ]
df_CK38P01 <- df %>% filter(Line.x == "CK38P01") %>% .[sample(nrow(.), size = 2000), ]
df_CK40P01 <- df %>% filter(Line.x == "CK40P01") %>% .[sample(nrow(.), size = 2000), ]
df_CK44P01 <- df %>% filter(Line.x == "CK44P01") %>% .[sample(nrow(.), size = 2000), ]

df_sampled_lines <- rbind(df_CK51P06, df_CK33P02, df_CK38P01,
                          df_CK40P01, df_CK44P01)

gg_morph <- df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 2, color = "#f1f1f1") + 
  geom_point_rast(data = df_sampled_lines  %>% 
                  filter(Morphology != "intermediate"),
                  # sample_frac(0.5),
                  aes(color = Morphology) ,alpha = 0.8 ,size =2, shape=16) + 
  # scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = c("#fa7f6f", "#beb8dc", "#82b0d2")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  theme_cowplot(font_size = 8) 
gg_morph


## organoid_morphology with 27 drugs
gg_drug <- df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 1, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>% 
                  filter(Drug == "Trametinib") %>%
                  filter(Morphology != "intermediate") %>%
                  sample_frac(2, replace = TRUE),
                  aes(color = Morphology),size =1.5, alpha = 0.8, shape =16) + 
  # scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = c("#F28C61","#6CC3AD","#93A3CB")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "Trametinib")+
  theme_cowplot(font_size = 8) 
gg_drug

