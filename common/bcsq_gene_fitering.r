#bcsq gene filtering

library("tidyverse")
library("plotly")
library("DT")
library("ggbeeswarm")
library("knitr")
library("ggrepel")
library("genetics")
library("ggnewscale")
library("cowplot")

eleg_finemapping <- as_tibble(read.delim("~/Desktop/telomere/data/elegans/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length_II_13272357-15257187_bcsq_genes.tsv"))
brig_finemapping <- as_tibble(read.delim("~/Desktop/telomere/data/briggsae/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length_I_14312709-15147409_bcsq_genes.tsv"))

  
pot2 <- eleg_finemapping %>% filter(CHROM=="II") %>% filter(POS>=14524173 & POS<=14525108) %>% arrange(VARIANT_IMPACT)
aex5 <- brig_finemapping  %>% filter(CHROM=="I") %>% filter(POS>=14669282 & POS<=14679748) %>% arrange(VARIANT_IMPACT) 

p <- pot2 %>% dplyr::select(POS,REF,ALT,TRANSCRIPT_BIOTYPE,CONSEQUENCE,VARIANT_IMPACT) %>% dplyr::distinct()
a <- aex5 %>% dplyr::filter(VARIANT_IMPACT=="HIGH") %>% dplyr::distinct()
a
write_tsv(a,'aex_bcsq.tsv')
