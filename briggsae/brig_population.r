#population plots

library("tidyverse")
library("plotly")
library("DT")
library("ggbeeswarm")
library("knitr")
library("ggrepel")
library("genetics")
library("ggnewscale")
library("cowplot")

cb_len <- as_tibble(read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/lengths/brig_strain.tsv')) %>% arrange(desc(length))
cb_len_isotype <- read_tsv('//Users/shrirambhat/Desktop/telomere/data/briggsae/lengths/brig_isotype.tsv')


pb <- cb_len_isotype %>% mutate(isotype=fct_reorder(isotype,length)) %>%
      ggplot(aes(x=isotype,y=length)) + 
        geom_bar(stat='identity') +
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none") +
        xlab("briggsae isotypes")
pb
ggsave("../plots/brig_population.png",width=7.5,height=3,units="in",dpi=300)
