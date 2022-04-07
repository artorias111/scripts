#manhattan plots

library("tidyverse")
library("plotly")
library("DT")
library("ggbeeswarm")
library("knitr")
library("ggrepel")
library("genetics")
library("ggnewscale")
library("cowplot")

load('../processed_data/ce_len_isotype.Rdata')


eleg_pop_plot <- ce_len_isotype %>% ggplot(aes(x="",y=length)) + 
  geom_boxplot() +
  ylim(0,50) + 
  theme(axis.title = element_blank(), legend.position = "none",axis.text.y = element_blank(),axis.ticks.y = element_blank())


eleg_pop_plot

save(brig_pop_plot,file="../processed_data/eleg_pop_boxplot.Rdata")