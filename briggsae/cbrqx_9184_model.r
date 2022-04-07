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
  gene='qx9184',
  chr='I',
  strand='+',
  txstart=14797708,
  txend=14806589,
  codingstart=txstart,
  codingend=txend,
  numexons=12,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='14797708,14797896,14798398,14800092,14800739,14802352,14803481,14804101,14804449,14805953,14806242,14806439',
  exonends='14797849,14798009,14799392,14800233,14801491,14802598,14803767,14804273,14804712,14806195,14806393,14806589',
  type='transcript'
)

cbrqx9184_model <- gene_model(goi_list,'qx9184')
cbrqx9184_model
save(cbrqx9184_model,file="../processed_data/cbrqx9184_model.Rdata")
