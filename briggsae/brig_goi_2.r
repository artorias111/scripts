#shortlisting goi round 2
library('tidyverse')
library('ape')

brig_finemapping <- as_tibble(read.delim("~/Desktop/telomere/data/briggsae/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length_I_14312709-15147409_bcsq_genes.tsv"))
brig_finemapping <- brig_finemapping %>% arrange(desc(VARIANT_LOG10p)) 
brig <- brig_finemapping %>% filter(VARIANT_IMPACT=='HIGH') %>% select(POS,REF,ALT,WBGeneID,TRANSCRIPT_BIOTYPE,CONSEQUENCE,VARIANT_IMPACT,VARIANT_LOG10p) %>% filter(VARIANT_LOG10p>17) %>% distinct()


brig_gff <- ape::read.gff('../data/briggsae/gffs/QX1410.SB.final_merger.v2.1.renamed.csq.gff')

brig_goi <- brig_gff %>% filter(type=='gene') %>% filter(seqid=='I') %>% arrange(start)

positions <- brig [,1]
#14764654 14764729 14804441 14798602 14800045

a <- brig_goi %>% filter(14764654>start & 14764654<end)

