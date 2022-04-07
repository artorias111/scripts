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
load('../processed_data/eleg_geno.Rdata')
load('../processed_data/ce_len_isotype.Rdata')
load('../processed_data/eleg_fgeno.Rdata')

celpot2 <- eleg_fgeno %>% filter(POS==14524396) %>% filter(CHROM==2)
eleg_temp1 <- celpot2 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)

iso_list <- eleg_temp1 %>% filter(altref==-1) %>% dplyr::select(isotype) %>% distinct() %>% add_column(status=1)
cel_isotypes <- ce_len_isotype %>% mutate(isotype=fct_reorder(isotype,length))

x <- as_tibble(setdiff(cel_isotypes$isotype,iso_list$isotype)) 
y <- x %>% add_column(status=-1)
y <- y %>% dplyr::rename(isotype=value)

z <- full_join(iso_list,y)
k <- inner_join(ce_len_isotype,z)
k <- k %>% mutate(status=as.character(status)) %>% mutate(status=str_replace(status,"-1","alt")) %>% mutate(status=str_replace(status,"1","ref"))

celpot2_bar <- k %>% mutate(isotype=fct_reorder(isotype,length)) %>%
      ggplot(aes(x=isotype,y=length,fill=as.factor(status))) + 
        labs(fill="isotype") + 
        geom_bar(stat='identity') +
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
        xlab(expression(paste(italic('C. elegans')," isotypes")))+
        ggtitle("I:14524396")

celpot2_bar
ggsave("../plots/celpot2_barplot_finemapping2.png",dpi=300,height=4,width = 8, units = "in")