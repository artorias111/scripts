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
load('../processed_data/brig_pop_boxplot.Rdata')
load('../processed_data/cbrqx9184_model.Rdata')

#cbraex5 <- brig_geno %>% dplyr::filter(POS>14761948 & POS<14772711) %>% dplyr::filter(CHROM=="I")
cbraex5 <- brig_fgeno %>% filter(POS==14804441 | POS==14798602 | POS==14800045) %>% filter(CHROM==1)
brig_temp1 <- cbraex5 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)
brig_temp <- brig_temp1 %>% inner_join(cb_len_isotype) %>% filter(altref==1) %>% dplyr::select(POS,length,isotype)

cbraex5_plot <- brig_temp %>% ggplot(aes(x=POS,y=length,group=POS)) + 
  geom_boxplot() + 
  #scale_x_continuous(breaks=c(14670406,14672327,14673372,14673503,14673745,14674803,14676169),lim=c(14669282,14679748))+
  scale_x_continuous(lim=c(14797708,14806589))+
  #geom_jitter(color="black", alpha=0.1) + 
  #theme(axis.title = element_blank(), legend.position = "none",axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))+
  theme(axis.title = element_blank(), legend.position = "none")+
  ylim(0,15)
  #geom_hline(yintercept=3.68,color='blue')


  #geom_vline(xintercept = c(14669282,14679748,14669282,14671667,14676069,14676749,14677416,14677725,14679067,14679211,14670011,14671846,14676265,14676905,14677619,14677823,14679164,14679748),linetype="dashed",alpha=0.5)
  #geom_vline(xintercept = c(14670406,14672327,14673372,14673503,14673745,14674803,14676169),linetype="dashed",alpha=0.5)

cbraex5_plot

p <- plot_grid(cbrqx9184_model,cbraex5_plot,ncol=1, align='v',rel_widths=c(1,1),rel_heights = c(1,12))
x <- plot_grid(NULL,brig_pop_plot,ncol=1,align='v',rel_widths=c(1,1),rel_heights = c(1,12))




#ggsave("../plots/cbraex5_boxplot.png",dpi=300,height=4,width = 8, units = "in")
#ggsave("../plots/cbraex5_boxplot_high_impact.png",dpi=300,height=4,width = 8, units = "in")

q <- plot_grid(p,x,rel_widths = c(3,1))
q
ggsave("../plots/cbrqx9184_boxplot_highimpact.png",dpi=300,height=4,width = 12, units = "in")
