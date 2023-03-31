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

## load data
umap_df = read_delim("../table/umap_df_for_1960_jpg.csv",escape_double = FALSE, trim_ws = TRUE)
umap_df <- mutate(umap_df, Size_log10 = log10(Org_size))
umap_df <- mutate(umap_df, Size_ln = log(Org_size))

# Organoid Size Distributions
gg_size_dist <- umap_df %>% 
  ggplot(aes(Org_size)) + 
  geom_histogram(bins = 30) + 
  theme_cowplot() 
gg_size_dist_log <- gg_size_dist + 
  scale_x_log10()
print(gg_size_dist_log)

## size_dist according to lines
organoid_size_factor_09 <- umap_df %>% group_by(Line_no_passage) %>% 
  filter(Line_no_passage == "CK02"| Line_no_passage == "CK04" | 
           Line_no_passage == "CK06" | Line_no_passage == "CK07" |
           Line_no_passage == "CK11" | Line_no_passage == "CK13" | 
           Line_no_passage == "CK08" | Line_no_passage == "CK10") %>%
  summarise(x = quantile(Size_log10, 0.9)) %>% 
  #summarise(x = mean(size_log)) %>% 
  arrange(x) %>% .$Line_no_passage 

org_size_by_line <- umap_df %>% 
  mutate(Line = factor(Line_no_passage, levels = organoid_size_factor_09)) %>% 
  na.omit() %>%
  ggplot() +
  geom_density_ridges_gradient(aes(y = Line, x = Size_log10, fill = stat(x)), scale = 1) +
  scale_fill_viridis_c() +
  labs(x = "log10(size)",
       fill = "size") + 
  theme(legend.position = "bottom") +
  theme_cowplot() + 
  coord_fixed(ratio = 1)
org_size_by_line

## size_dist according date
organoid_size_by_date <- umap_df %>% group_by(Date) %>% 
  summarise(x = quantile(Size_ln, 0.9)) %>%
  .$Date
organoid_size_by_date_1 <- organoid_size_by_date[1:7] 
# organoid_size_by_date_2 <- organoid_size_by_date[36:40]
# organoid_size_by_date_3 <- organoid_size_by_date[43:47]

size_dist_morph_by_date_ridge <- umap_df %>% 
  mutate(Date = factor(Date, levels = organoid_size_by_date_1)) %>% 
  na.omit(.) %>%
  ggplot() +
  geom_density_ridges_gradient(aes(y = Date, x = Size_log10, fill = stat(x)), scale = 1, alpha=0.5) +
  scale_fill_viridis_c(option = "C") +
  labs(y = "Day",
       x = "log10(size)",
       fill = "size") + 
  theme(legend.position = "bottom") +
  theme_cowplot() + 
  coord_fixed(ratio = 1)
size_dist_morph_by_date_ridge

## UMAP size distribution
umap_df %>%
  ggplot(aes(UMAP1, UMAP2, color = Size_log10)) + 
  geom_point(alpha = 0.5, size = 0.35) +
  scale_color_viridis_c() +
  theme_cowplot() + 
  labs(x = "UMAP 1",
       y = "UMAP 2",
       color = "log10(size)")

## four organoid line differences
loi <- c("CK02","CK08","CK10", "CK12")
gg_line <- umap_df %>% dplyr::select(-Line_no_passage) %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") + 
  geom_point_rast(data = umap_df %>%
                  filter(Line_no_passage %in% loi) %>% 
                  mutate(Line_no_passage = factor(Line_no_passage, levels = loi)) %>% 
                  sample_frac(1),
                  aes(color = Line_no_passage), alpha = .65,size = 1, shape=16) + 
  facet_wrap( ~ Line_no_passage, ncol =2) +
  scale_color_manual(values = c("#d43232","#aaa4d6","#81b29a","#657f8d")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  theme_cowplot(font_size = 8) + 
  theme(legend.position = "nothing")
gg_line


## organoid morphology across all lines
organoid_morphology <- read_delim(("../table/organoids_morphology.csv"),escape_double = FALSE, trim_ws = TRUE) 
df <- umap_df %>%
  left_join(organoid_morphology,by = 'Field')

gg_cystic <- umap_df %>% 
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 0.35, color = "#f1f1f1") +
  geom_point_rast(data = df %>%
                  filter(Morphology != "intermediate") %>%
                  sample_frac(1),
                  aes(color = Morphology), size = 0.35, shape=16) +
  scale_color_manual(values = c("#fa7f6f", "#beb8dc")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  #caption = "control treated organoids") + 
  theme_cowplot(font_size = 8) + 
  #theme(legend.position = "nothing")  + 
  coord_fixed()
gg_cystic

## organoid morphology across 4 lines
# set.seed(123)
# organoid_morphology <- read_delim(("../table/organoids_morphology.csv"),escape_double = FALSE, trim_ws = TRUE)
# df <- umap_df %>%
#   left_join(organoid_morphology,by = 'Field')

## random sampling 1000 point
df_CK01 <- df %>% filter(Line_no_passage == "CK01") %>% .[sample(nrow(.), size = 1000), ]
df_CK02 <- df %>% filter(Line_no_passage == "CK02")
df_CK03 <- df %>% filter(Line_no_passage == "CK03") %>% .[sample(nrow(.), size = 1000), ]
df_CK04 <- df %>% filter(Line_no_passage == "CK04")
df_CK05 <- df %>% filter(Line_no_passage == "CK05") %>% .[sample(nrow(.), size = 1000), ]
df_CK06 <- df %>% filter(Line_no_passage == "CK06") %>% .[sample(nrow(.), size = 1000), ]
df_CK07 <- df %>% filter(Line_no_passage == "CK07") %>% .[sample(nrow(.), size = 1000), ]
df_CK08 <- df %>% filter(Line_no_passage == "CK08") 
df_CK09 <- df %>% filter(Line_no_passage == "CK09") %>% .[sample(nrow(.), size = 1000), ]
df_CK10 <- df %>% filter(Line_no_passage == "CK10")
df_CK11 <- df %>% filter(Line_no_passage == "CK11") 
df_CK12 <- df %>% filter(Line_no_passage == "CK12") 
df_CK13 <- df %>% filter(Line_no_passage == "CK13") %>% .[sample(nrow(.), size = 1000), ]
df_CK24 <- df %>% filter(Line_no_passage == "CK24") %>% .[sample(nrow(.), size = 1000), ]

## select 4 lines
df <- rbind(df_CK02, df_CK03, df_CK07,df_CK08)

gg_morph <- umap_df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 2, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>% 
                  filter(Morphology != "intermediate"),
                  # sample_frac(0.5),
                  aes(color = Morphology) ,alpha = 0.8 ,size =2, shape=16) + 
  scale_color_manual(values = c("#fa7f6f", "#beb8dc", "#82b0d2")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2")+
  theme_cowplot(font_size = 8) 
gg_morph
  

## organoid morphological statistics across 8 lines
organoid_morphology_sta <- read_csv(("../table/organoid_morphology_statistics.csv"))

gg_morph_statistics_by_lines <- organoid_morphology_sta %>%
  filter(Line=="CK02"| Line=="CK04" | Line=="CK06" | Line=="CK07" | Line=="CK08" |
           Line=="CK10" | Line=="CK11" | Line=="CK13") %>%
  ggplot(aes(x = Line, y = Num, fill=Morphology)) + 
  geom_bar(stat = "identity", position = position_fill(), width = 0.5) +
  labs(x = "Line", y ="Num_ratio") +
  theme_cowplot() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#fa7f6f", "#beb8dc", "#82b0d2")) 
gg_morph_statistics_by_lines


## organoid morphology with 3 drugs in CK07
organoid_morphology_drug <- read_delim(("../table/organoids_morphology_drug.csv"),escape_double = FALSE, trim_ws = TRUE) 
df <- umap_df %>%
  left_join(organoid_morphology_drug,by = 'Field')%>%
  filter(Line_no_passage=="CK07")

gg_drug <- umap_df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 2, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>% 
                  filter(Drug == "DMSO") %>%
                  filter(Morphology != "intermediate"),
                  aes(color = Morphology),size =2, alpha = 0.8, shape =16) + 
  # scale_color_brewer(type = "qual", palette = "Set1") +
  scale_color_manual(values = c("#F28C61","#6CC3AD","#93A3CB")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "DMSO")+
  theme_cowplot(font_size = 8) 
gg_drug

## organoid sizes with 3 drugs 7 conc in CK07
gg_conc <- umap_df %>%
  ggplot(aes(UMAP1, UMAP2)) + 
  geom_point_rast(alpha = 1, size = 2, color = "#f1f1f1") + 
  geom_point_rast(data = df  %>% 
                  filter(Drug == "DMSO")%>%
                  sample_frac(0.5),
                  aes(color = sort(as.character(Conc))),size =2, alpha =1, shape =16) + 
  scale_color_manual(values = c("#001959","#B28C32","#F28C61","#6CC3AD",
                                "#F9CCF9","#93A3CB","#577646")) +
  theme_classic() +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "DMSO") +
  theme_cowplot(font_size = 8)+
  theme(legend.position = "bottom")
gg_conc


## line chart of organoid_sizes with 3 drugs 7 conc in CK07
## 5fu
df_5fu_200 <- df %>%
  filter(Drug == "5-FU" & Conc == "200")
df_5fu_200_size_log10 = median(df_5fu_200$Size_log10)
df_5fu_100 <- df %>%
  filter(Drug == "5-FU" & Conc == "100")
df_5fu_100_size_log10 = median(df_5fu_100$Size_log10)
df_5fu_50 <- df %>%
  filter(Drug == "5-FU" & Conc == "50")
df_5fu_50_size_log10 = median(df_5fu_50$Size_log10)
df_5fu_10 <- df %>%
  filter(Drug == "5-FU" & Conc == "10")
df_5fu_10_size_log10 = median(df_5fu_10$Size_log10)
df_5fu_5 <- df %>%
  filter(Drug == "5-FU" & Conc == "5")
df_5fu_5_size_log10 = median(df_5fu_5$Size_log10)
df_5fu_1 <- df %>%
  filter(Drug == "5-FU" & Conc == "1")
df_5fu_1_size_log10 = median(df_5fu_1$Size_log10)
df_5fu_0_5 <- df %>%
  filter(Drug == "5-FU" & Conc == "0.5")
df_5fu_0_5_size_log10 = median(df_5fu_0_5$Size_log10)

##Oxaliplatin
df_oxa_200 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "200")
df_oxa_200_size_log10 = median(df_oxa_200$Size_log10)
df_oxa_100 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "100")
df_oxa_100_size_log10 = median(df_oxa_100$Size_log10)
df_oxa_50 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "50")
df_oxa_50_size_log10 = median(df_oxa_50$Size_log10)
df_oxa_10 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "10")
df_oxa_10_size_log10 = median(df_oxa_10$Size_log10)
df_oxa_5 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "5")
df_oxa_5_size_log10 = median(df_oxa_5$Size_log10)
df_oxa_1 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "1")
df_oxa_1_size_log10 = median(df_oxa_1$Size_log10)
df_oxa_0_5 <- df %>%
  filter(Drug == "Oxaliplatin" & Conc == "0.5")
df_oxa_0_5_size_log10 = median(df_oxa_0_5$Size_log10)

##FOLFOX
df_fol_200 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "200")
df_fol_200_size_log10 = median(df_fol_200$Size_log10)
df_fol_100 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "100")
df_fol_100_size_log10 = median(df_fol_100$Size_log10)
df_fol_50 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "50")
df_fol_50_size_log10 = median(df_fol_50$Size_log10)
df_fol_10 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "10")
df_fol_10_size_log10 = median(df_fol_10$Size_log10)
df_fol_5 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "5")
df_fol_5_size_log10 = median(df_fol_5$Size_log10)
df_fol_1 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "1")
df_fol_1_size_log10 = median(df_fol_1$Size_log10)
df_fol_0_5 <- df %>%
  filter(Drug == "FOLFOX" & Conc == "0.5")
df_fol_0_5_size_log10 = median(df_fol_0_5$Size_log10)


## correlation of cell_viability vs. organoid size
# organoid_morphology_drug <- read_delim(("../table/organoids_morphology_drug.csv"),escape_double = FALSE, trim_ws = TRUE) 
df <- umap_df %>%
  left_join(organoid_morphology_drug,by = 'Field') %>%
  filter(Line_no_passage =="CK07") %>%
  filter(Date == "20220121" | Date == "20220122")
 
org_via_size <- df %>% group_by(df$Field) %>%
  summarise(count = n(),
            Median_size = median(Size_log10),
            Mean_viability = mean(Cell_viability)) %>%
  mutate(Normal_viability = (Mean_viability - min(Mean_viability)) / (max(Mean_viability) - min(Mean_viability)))

library(ggstatsplot)
ggscatterstats(
  data = org_via_size, ##%>% filter(Drug == "FOLFOX"),
  x = Median_size,
  y = Normal_viability,
  bf.message = FALSE
  ##title = "FOLFOX"
)
cor(org_via_size$Median_size, org_via_size$Normal_viability)


## cell_viability vs. organoid morphology
library(plotrix)
org_via_morph <- df %>% group_by(df$Morphology) %>%
  summarise(count = n(),
            mean_viability = mean(Cell_viability),
            median_viability = median(Cell_viability),
            sem = std.error(Cell_viability),
            sd = sd(Cell_viability))

org_via_morph_5fu <- df %>% filter(Drug == "5-FU") %>% 
  group_by(Morphology) %>%
  summarise(count = n(),
            mean_viability = mean(Cell_viability),
            median_viability = median(Cell_viability),
            sem = std.error(Cell_viability),
            sd = sd(Cell_viability))
    
org_via_morph_oxa <- df %>% filter(Drug == "Oxaliplatin") %>% 
  group_by(Morphology) %>%
  summarise(count = n(),
            mean_viability = mean(Cell_viability),
            median_viability = median(Cell_viability),
            sem = std.error(Cell_viability),
            sd = sd(Cell_viability))

org_via_morph_fol <- df %>% filter(Drug == "FOLFOX") %>% 
  group_by(Morphology) %>%
  summarise(count = n(),
            mean_viability = mean(Cell_viability),
            median_viability = median(Cell_viability),
            sem = std.error(Cell_viability),
            sd = sd(Cell_viability))
