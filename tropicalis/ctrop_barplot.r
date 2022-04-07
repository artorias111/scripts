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
load('../processed_data/trop_geno.Rdata')
load('../processed_data/ct_len_isotype.Rdata')


ctrop <- trop_geno %>% filter(POS==15331537) %>% filter(CHROM=="I")
brig_temp1 <- cbrmrt1 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)

iso_list <- brig_temp1 %>% filter(altref==-1) %>% dplyr::select(isotype) %>% distinct() %>% add_column(status=1)
cbrmrt_strains <- cb_len_isotype %>% mutate(isotype=fct_reorder(isotype,length))

x <- as_tibble(setdiff(cbrmrt_strains$isotype,iso_list$isotype)) 
y <- x %>% add_column(status=-1)
y <- y %>% dplyr::rename(isotype=value)

z <- full_join(iso_list,y)
k <- inner_join(cb_len_isotype,z)
k <- k %>% mutate(status=as.character(status)) %>% mutate(status=str_replace(status,"-1","alt")) %>% mutate(status=str_replace(status,"1","ref"))

cbrmrt_bar <- k %>% mutate(isotype=fct_reorder(isotype,length)) %>%
      ggplot(aes(x=isotype,y=length,fill=as.factor(status))) + 
        labs(fill="isotype") + 
        geom_bar(stat='identity') +
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
        xlab(expression(paste(italic('C. briggsae')," isotypes")))+
  ggtitle("I:15331537")

cbrmrt_bar
ggsave("../plots/cbrmrt1_barplot4.png",dpi=300,height=4,width = 8, units = "in")