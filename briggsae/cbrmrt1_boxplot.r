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
load('../processed_data/cb_len_isotype.Rdata')
load('../processed_data/cbrmrt1_model.Rdata')


cbrmrt1 <- brig_geno %>% filter(POS>=15329043 & POS<=15332022) %>% filter(CHROM=="I")
brig_temp1 <- cbrmrt1 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)
# -1 is ref, 1 is alt
brig_temp <- brig_temp1 %>% inner_join(cb_len_isotype) %>% filter(altref==1) %>% dplyr::select(POS,length,isotype)

cbrmrt1_plot <- brig_temp %>% ggplot(aes(x=POS,y=length,group=POS)) + 
  geom_boxplot() + 
  #geom_jitter(color="black", alpha=0.1) + 
  xlim(15329043,15332022) + 
  theme(axis.title = element_blank(), legend.position = "none",)+
  geom_vline(xintercept = c(15329043,15332022,15329877,15330271,15331113,15331537),linetype="dashed",alpha=0.5) 

cbrmrt1_plot

p <- plot_grid(cbrmrt1_model,cbrmrt1_plot,ncol=1, align='v',rel_widths=c(1,1),rel_heights = c(1,12))
p
#ggsave("../plots/cbrmrt1_boxplot.png",dpi=300,height=4,width = 8, units = "in")
