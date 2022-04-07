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
  gene='qx9296',
  chr='I',
  strand='+',
  txstart=15440188,
  txend=15446390,
  codingstart=txstart,
  codingend=txend,
  numexons=10,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='15440188,15440420,15440974,15441158,15442278,15442544,15442669,15442885,15446061,15446190',
  exonends='15440366,15440929,15441105,15441478,15442499,15442624,15442842,15443056,15446144,15446390',
  type='transcript'
)

cbr_model <- gene_model(goi_list,'qx9296')
cbr_model
save(cbr_model,file="../processed_data/cbrqx9296_model.Rdata")
