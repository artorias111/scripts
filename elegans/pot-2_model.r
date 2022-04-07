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

load("../processed_data/goi_list.Rdata")

#14524173 .. 14525108

goi_list <- goi_list %>% add_row(
  gene='pot2',
  chr='II',
  strand='+',
  txstart=14524173,
  txend=14525108,
  codingstart=14524173,
  codingend=14525108,
  numexons=3,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='14524173,14524585,14524884',
  exonends='14524539,14524832,14525108',
  type='transcript'
)

pot2_model <- gene_model(goi_list,'pot2')
pot2_model
save(pot2_model,file="../processed_data/pot2_model.Rdata")

