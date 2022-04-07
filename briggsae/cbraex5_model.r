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

#14669282 .. 14679748

load("../functions/gene_model.Rdata")
load("../processed_data/goi_list.Rdata")

goi_list <- goi_list %>% add_row(
  gene='cbraex5',
  chr='I',
  strand='+',
  txstart=14669282,
  txend=14679748,
  codingstart=14669282,
  codingend=14679748,
  numexons=8,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='14669282,14671667,14676069,14676749,14677416,14677725,14679067,14679211',
  exonends='14670011,14671846,14676265,14676905,14677619,14677823,14679164,14679748',
  type='transcript'
)

cbraex5_model <- gene_model(goi_list,'cbraex5')
save(cbraex5_model,file="../processed_data/cbraex51_model.Rdata")
