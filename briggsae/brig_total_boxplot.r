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

load('../processed_data/cb_len_isotype.Rdata')


brig_pop_plot <- cb_len_isotype %>% ggplot(aes(x="",y=length)) + 
  geom_boxplot() +
  ylim(0,15) + 
  theme(axis.title = element_blank(), legend.position = "none",axis.text.y = element_blank(),axis.ticks.y = element_blank())


brig_pop_plot

save(brig_pop_plot,file="../processed_data/brig_pop_boxplot.Rdata")

