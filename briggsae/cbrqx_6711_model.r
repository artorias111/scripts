#orthofinder - cbrmrt1.2

library("tidyverse")
library("plotly")
library("DT")
library("ggbeeswarm")
library("knitr")
library("ggrepel")
library("genetics")
library("ggnewscale")
library("cowplot")

#

load("../functions/gene_model.Rdata")
load("../processed_data/goi_list.Rdata")

goi_list <- goi_list %>% add_row(
  gene='qx6711',
  chr='I',
  strand='+',
  txstart=2996606,
  txend=2999968,
  codingstart=txstart,
  codingend=txend,
  numexons=8,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='2996606,2996837,2997851,2998136,2999071,2999352,2999483,2999710',
  exonends='2996784,2997343,2997982,2998456,2999292,2999432,2999656,2999968',
  type='transcript'
)

cbrqx6711_model <- gene_model(goi_list,'qx6711')
cbrqx6711_model
save(cbrqx6711_model,file="../processed_data/cbrqx6711_model.Rdata")
