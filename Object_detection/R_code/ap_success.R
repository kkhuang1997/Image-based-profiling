library(ggplot2)
library(cowplot)

## load data
data = read.csv("../table/AP_success.csv", row.names = 1)

# Barplot
ggplot(data, aes(x=models, y=AP_success)) + 
  geom_bar(stat = "identity") +
  theme_cowplot() +
  labs(x = "Models",
       y = "AP_success(%)") + 
  coord_flip()
