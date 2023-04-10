library(reshape2)
library(ggplot2)
library(ggridges)
library(forcats)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)

## load data
org_conf <- read_excel("../table/org_confidence.xlsx", 
                             col_types = c("skip", "text", "numeric"))

org_conf %>%
  mutate(text = fct_reorder(Date, Live_prob)) %>%
  ggplot(aes(y=Date, x=Live_prob,  fill=Date)) +
  geom_density_ridges_gradient(alpha=0.5, scale = 2, rel_min_height = 0.011) +
  theme_ridges() +
  scale_fill_viridis_d() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
  ) +
  xlab("Live_prob") +
  ylab("Date")
