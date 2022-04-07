args <- commandArgs(trailingOnly = TRUE)


library("Biostrings")
library(stringr)
library(dplyr)
library(tidyr)
library(ape)


#set output file suffix
suffix <- "CE-CB.txt"
#set output directory
outdir <- "/projects/b1059/projects/Nicolas/collabs/forShriram/blast_test/"

#read C. briggsae (C.b.) protein FASTA into AA stringset object
pred <- readAAStringSet("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/prot/QX1410.SB.final_anno.dedup.renamed.prot.fa")

#read C. elegans (C.e.) protein FASTA into stringset object
elegansprot <- readAAStringSet("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.prot.fa")

#read blast results (forward and reciprocal)
blast = read.table("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/blast_out/forward/QX1410.SB.final_anno.dedup.renamed.prot.pb.out", header = FALSE, sep = "", dec = ".")
recipro = read.table("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/Plots/blast_out/reciprocals/QX1410.SB.final_anno.dedup.renamed.prot.recipro.pb.out", header = FALSE, sep = "", dec = ".")
colnames(blast) <- c("query","subject","ident","length","mismatch","gapop","qstart","qend","sstart","send","eval","bitscore")
colnames(recipro) <- c("query","subject","ident","length","mismatch","gapop","qstart","qend","sstart","send","eval","bitscore")

#read C.b. GFF 
predGFF <- read.gff("/projects/b1059/projects/Nicolas/c.briggsae/gene_predictions/final_annotations/check_dup_prot/QX1410.SB.final_anno.dedup.renamed.gff")
predGFF <- predGFF %>% dplyr::filter(type=="mRNA") %>%
                       tidyr::separate(attributes,into=c("Transcript_ID","Transcript_parent","Name","Other"),sep=";",extra="merge")
predGFF$Transcript_ID <- gsub(".*=","",predGFF$Transcript_ID)
predGFF$Transcript_parent <- gsub(".*=","",predGFF$Transcript_parent)

#read C.e. GFF 
celeGFF <- read.gff("/projects/b1059/projects/Nicolas/c.elegans/N2/wormbase/WS279/c_elegans.PRJNA13758.WS279.protein_coding.gff")
celeGFF <- celeGFF %>% dplyr::filter(type=="mRNA") %>% 
                       tidyr::separate(attributes,into=c("Subject_ID","Subject_parent","Name","Other"),sep=";",extra="merge")
celeGFF$Subject_ID <- gsub(".*=","",celeGFF$Subject_ID)
celeGFF$Subject_parent <- gsub(".*=","",celeGFF$Subject_parent)


#transform AA stringset object to 2-column dataframe
prednames <- names(pred)
predseq <- paste(pred)
preddf <- data.frame(prednames,predseq)
eprotnames <- names(elegansprot)
eprotseq <- paste(elegansprot)
eprotdf <- data.frame(eprotnames,eprotseq)

#calculate protein length for each FASTA entry
preddf$qlen <- str_count(preddf$predseq)
eprotdf$slen <- str_count(eprotdf$eprotseq)


#select best hit(s) per transcript for C.b.
orderedhits <- blast[order(blast$query,-blast$bitscore, blast$eval),]
zeroeval <- orderedhits[orderedhits$eval == 0,]
nonzeroeval <- orderedhits[orderedhits$eval > 0,]
filterednzhits <- nonzeroeval[!duplicated(nonzeroeval$query),]
filteredhits <- bind_rows(zeroeval,filterednzhits)
filteredhits <- filteredhits[order(filteredhits$query),]
filteredhits_qse <- filteredhits[,c(1:2,11,12)]
colnames(filteredhits_qse) <- c('Transcript_ID','Subject_ID','eval','bitscore')


#select best hit(s) per transcript for C.e.
rorderedhits <- recipro[order(recipro$query, -recipro$bitscore, recipro$eval),]
rzeroeval <- rorderedhits[rorderedhits$eval == 0,]
rnonzeroeval <- rorderedhits[rorderedhits$eval > 0,]
rfilterednzhits <- rnonzeroeval[!duplicated(rnonzeroeval$query),]
rfilteredhits <- bind_rows(rzeroeval,rfilterednzhits)
rfilteredhits <- rfilteredhits[order(rfilteredhits$query),]
rfilteredhits_qse <- rfilteredhits[,c(2,1,11,12)]
colnames(rfilteredhits_qse) <- c('Transcript_ID','Subject_ID','eval','bitscore')

#append gene-transcript relationships for every BLAST hit
forwardHits<- dplyr::left_join(filteredhits_qse,predGFF %>% dplyr::select("Transcript_ID","Transcript_parent"),by="Transcript_ID",keep=F)
forwardHits<- dplyr::left_join(forwardHits,celeGFF %>% dplyr::select("Subject_ID","Subject_parent"),by="Subject_ID",keep=F)
forwardHits <- forwardHits[,c(6,2,5,1,3,4)]
reciproHits<- dplyr::left_join(rfilteredhits_qse,predGFF %>% dplyr::select("Transcript_ID","Transcript_parent"),by="Transcript_ID",keep=F)
reciproHits<- dplyr::left_join(reciproHits,celeGFF %>% dplyr::select("Subject_ID","Subject_parent"),by="Subject_ID",keep=F)
reciproHits <- reciproHits[,c(6,2,5,1,3,4)]

#merge forward and reciprocal hits
aggregates <- dplyr::bind_rows(forwardHits,reciproHits, .id = "set")
aggregates <- aggregates[order(aggregates$Transcript_parent),]

#identify which transcripts have matching forward and reciprocal hits
common <- aggregates %>% 
  group_by(Transcript_ID,Subject_ID) %>% 
  mutate(dupe = n()>1)

#identify best hit per locus
common <- common[order(common$Transcript_parent,-common$bitscore,common$eval),]
common_dupes <- common[common$dupe == TRUE,]
common_ndg <- common_dupes[!duplicated(common_dupes$Transcript_parent),]
commonNoDup <- common_ndg[,2:7]

#append protein lengths to best hit per locus and estimate protein length ratio (C.b. protein length / C.e. protein length)
predhits <- dplyr::left_join(commonNoDup,preddf, by=c("Transcript_ID"="prednames"))
predhits <- dplyr::left_join(predhits,eprotdf, by=c("Subject_ID"="eprotnames")) %>% dplyr::select(-eprotseq,-predseq)
predhits$plen_acc <- (predhits$qlen / predhits$slen)
predhits <- predhits[order(predhits$Transcript_parent,-predhits$bitscore,predhits$eval),]
colnames(predhits) <- c("CE_gene","CE_transcript","Query_Gene","Query_transcript","eval","bitscore","Query_ProtLen","CE_ProtLen","ProtLen_ratio")

#set output file dir and name
filenameRatio <- paste0(outdir,"/BestHitPerLocus-",suffix)

#write output tables
write.table(predhits, sep="\t",quote = FALSE,  col.names=FALSE, row.names = FALSE,file = filenameRatio)
