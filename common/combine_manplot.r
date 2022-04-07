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

load("/Users/shrirambhat/Desktop/telomere/processed_data/eleg_manplot.Rdata")
load('/Users/shrirambhat/Desktop/telomere/processed_data/brig_manplot.Rdata')
load('/Users/shrirambhat/Desktop/telomere/processed_data/trop_manplot.Rdata')

p <- plot_grid(eleg_manplot,brig_manplot,trop_manplot,ncol=1,align = "v",labels = c("A","B","C"))
p
ggsave("../plots/manplot.png",width=7.5,height=8,units="in",dpi=300)