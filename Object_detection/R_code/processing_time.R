library(ggplot2)
library(cowplot)

## load data
data = read.csv("../table/processing_time.csv", row.names = 1)

# Barplot
ggplot(data, aes(x = Model, y = Time)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() +
  labs(y = "Time(hours)") + 
  coord_flip()
