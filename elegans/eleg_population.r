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

ce_len <- as_tibble(read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/lengths/eleg_strain.tsv')) %>% arrange(desc(length))
ce_len_isotype <- read_tsv('//Users/shrirambhat/Desktop/telomere/data/elegans/lengths/eleg_isotype.tsv')


pe <- ce_len_isotype %>% mutate(isotype=fct_reorder(isotype,length)) %>%
  ggplot(aes(x=isotype,y=length)) + 
  geom_bar(stat='identity') +
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none") +
  xlab("elegans isotypes")

pe
ggsave("../plots/eleg_population.png",width=7.5,height=3,units="in",dpi=300)
