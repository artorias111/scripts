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

load("../functions/gene_model.Rdata")
load("../processed_data/goi_list.Rdata")
goi_list <- goi_list %>% add_row(
  gene='cbrmrt1',
  chr='I',
  strand='+',
  txstart=15327767,
  txend=15332867,
  codingstart=15327767,
  codingend=15332867,
  numexons=13,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='15327767,15328194,15328397,15329043,15329367,15329585,15329710,15330731,15331104,15331280,15331843,15332165,15332648',
  exonends='15328147,15328349,15328540,15329324,15329540,15329665,15329931,15331051,15331235,15331789,15332022,15332603,15332867',
  type='transcript'
)

cbrmrt1_model <- gene_model(goi_list,'cbrmrt1')
cbrmrt1_model
save(cbrmrt1_model,file="../processed_data/cbrmrt1_model.Rdata")
