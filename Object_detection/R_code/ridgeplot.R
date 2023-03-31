library(reshape2)
library(ggplot2)
library(ggridges)
library(forcats)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(readxl)

## load data
org_conf = read.csv("../table/org_confidence.csv", row.names = 1)
org_conf[,1] = as.character(org_conf[,1])
org_conf %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot(aes(y=text, x=value,  fill=text)) +
  geom_density_ridges_gradient(alpha=0.5, scale = 2, rel_min_height = 0.011) +
  theme_ridges() +
  scale_fill_viridis_d() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8),
  ) +
  xlab("Success probability(%)") +
  ylab("Date")
