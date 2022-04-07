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

#

load("../functions/gene_model.Rdata")
load("../processed_data/goi_list.Rdata")

goi_list <- goi_list %>% add_row(
  gene='qx9180',
  chr='I',
  strand='+',
  txstart=14761948,
  txend=14772711,
  codingstart=txstart,
  codingend=txend,
  numexons=8,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='14761948,14764636,14769040,14769719,14770386,14770696,14772039,14772183',
  exonends='14762109,14764815,14769236,14769875,14770589,14770794,14772136,14772711',
  type='transcript'
)

cbrqx9180_model <- gene_model(goi_list,'qx9180')
cbrqx9180_model
save(cbrqx9180_model,file="../processed_data/cbrqx9180_model.Rdata")
