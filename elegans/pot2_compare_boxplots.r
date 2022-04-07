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

#load eleg geno
load('../processed_data/eleg_fgeno.Rdata')
load('../processed_data/ce_len_isotype.Rdata')
load('../processed_data/pot2_model.Rdata')

#celpot2 <- eleg_fgeno %>% filter(POS>=14524173 & POS<=14525108) %>% filter(CHROM==2)
celpot2 <- eleg_fgeno %>% filter(POS==14524804 | POS==14524396) %>% filter(CHROM==2)
eleg_temp1 <- celpot2 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)


eleg_temp <- eleg_temp1 %>% inner_join(ce_len_isotype) %>% filter(altref==1) %>% dplyr::select(POS,length,isotype)

   celpot2_plot <- eleg_temp %>% ggplot(aes(x=POS,y=length,group=POS)) + 
   geom_boxplot() + 
   #geom_jitter(color="black", alpha=0.1) + 
   xlim(14524173,14525108) + 
   theme(axis.title = element_blank(), legend.position = "none",)+
   geom_vline(xintercept = c(14524173,14525108,14524396,14524804),linetype="dashed",alpha=0.5) 


celpot2_plot
p <- plot_grid(pot2_model,celpot2_plot,ncol=1, align='v',rel_widths=c(1,1),rel_heights = c(1,12))
p
#ggsave("../plots/pot2_boxplot_high_impact.png",dpi=300,height=4,width = 8, units = "in")
