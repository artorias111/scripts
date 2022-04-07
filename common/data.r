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

#genotype matrix
eleg_geno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/nemascan_runs/Analysis_Results-20220308/Genotype_Matrix/Genotype_Matrix.tsv')
save(eleg_geno,file="../processed_data/eleg_geno.Rdata")
brig_geno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/nemascan_runs/Analysis_Results-20220308/Genotype_Matrix/Genotype_Matrix.tsv')
save(brig_geno,file="../processed_data/brig_geno.Rdata")
trop_geno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/tropicalis/nemascan_runs/Analysis_Results-20220308/Genotype_Matrix/Genotype_Matrix.tsv')
save(trop_geno,file="../processed_data/trop_geno.Rdata")

#fine-mapping genotype matrix
eleg_fgeno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length.II:13272357-15257187.ROI_Genotype_Matrix.tsv')
save(eleg_fgeno,file="../processed_data/eleg_fgeno.Rdata")
brig_fgeno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length.I:14312709-15147409.ROI_Genotype_Matrix.tsv')
save(brig_fgeno,file="../processed_data/brig_fgeno.Rdata")
#trop_fgeno <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/tropicalis/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length.I:775059-2137061.ROI_Genotype_Matrix.tsv')


#telomere lengths
ce_len <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/lengths/eleg_strain.tsv') %>% arrange(desc(length))
ce_len_isotype <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/lengths/eleg_isotype.tsv')
save(ce_len_isotype,file="../processed_data/ce_len_isotype.Rdata")

cb_len <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/lengths/brig_strain.tsv') %>% arrange(desc(length))
cb_len_isotype <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/lengths/brig_isotype.tsv')
save(cb_len_isotype,file="../processed_data/cb_len_isotype.Rdata")

ct_len <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/tropicalis/lengths/trop_strain.tsv') %>% arrange(desc(length))
ct_len_isotype <- read_tsv('/Users/shrirambhat/Desktop/telomere/data/tropicalis/lengths/trop_isotype.tsv')
save(ct_len_isotype,file="../processed_data/ct_len_isotype.Rdata")