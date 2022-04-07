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

#load brig geno
load('../processed_data/brig_geno.Rdata')
load('../processed_data/brig_fgeno.Rdata')
load('../processed_data/cb_len_isotype.Rdata')
#load('../processed_data/cbraex51_model.Rdata')


cbraex5 <- brig_geno %>% filter(POS>=14669282 & POS<=14679748) %>% filter(CHROM=="I")
#cbraex5 <- brig_fgeno %>% filter(POS==14670406 | POS==14672327 | POS ==14673503 | POS ==14673745| POS ==14674803| POS ==14676169) %>% filter(CHROM==1)
brig_temp1 <- cbraex5 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)
brig_temp <- brig_temp1 %>% inner_join(cb_len_isotype) %>% filter(altref==-1) %>% dplyr::select(length)


cbraex5_plot_2 <- brig_temp %>% ggplot(aes(x="",y=length)) + 
  geom_boxplot() + 
  theme(axis.title = element_blank(), legend.position = "none",axis.text.y = element_blank(),axis.ticks.y = element_blank())


cbraex5_plot_2

save(cbraex5_plot_2,file="../processed_data/cbraex5_boxplot_2.Rdata")
