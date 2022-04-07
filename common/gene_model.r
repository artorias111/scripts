#gene models

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

gene_model <- function(df, genename = NA, WBID = NA, gene_color = "blue", intron_color = "black",
                       utr3_color = "gray60", utr5_color = "gray60", gene_alpha = 0.5) {
    if(is.na(genename)){
        if(is.na(WBID)) {
            stop("No gene provided. Give either gene name or WB ID")
        } else {
            selection <- df %>%
                dplyr::filter(wbgene == WBID)
        }
    } else {
        selection <- df %>%
            dplyr::filter(gene == genename)
    }
    
    if(nrow(selection) > 1){
        stop("multiple entries with given gene name.")
    }
    
    if(selection$strand == "+"){
        
        selection$txend <- as.numeric(selection$txend)
        selection$codingend <- as.numeric(selection$codingend)
        exonstarts <- matrix(as.numeric(unlist(strsplit(selection$exonstarts, ",")), ncol = selection$numexons))
        
        exonends <- matrix(as.numeric(unlist(strsplit(selection$exonends, ",")), ncol = selection$numexons))
        
        allexons <- data.frame("starts" = exonstarts, "ends" = exonends)
        
        intronstarts <- exonends[1:(nrow(exonends)-1)]
        
        intronends <- exonstarts[2:(nrow(exonends))]
        
        allintrons <- data.frame("starts" = intronstarts, "ends" = intronends) %>%
            dplyr::mutate(midpoint = (starts+ends)/2)
        
        endUTR <- data.frame("x" = c(selection$codingend, selection$codingend, selection$txend), "y" = c(1, -1, 0))
        
        plot <- ggplot(allexons)+
            geom_segment(aes(x = txstart, xend = txend, y = 0, yend = 0), data = selection, color = "black")+
            geom_rect( aes(xmin =  starts, xmax = ends, ymin = -1 , ymax = 1), fill = gene_color, color = "black", alpha = gene_alpha)+
            geom_segment(aes(x = starts, y = 1, xend = midpoint, yend = 2), data = allintrons, color = intron_color)+
            geom_segment(aes(x = midpoint, y = 2, xend = ends, yend = 1), data = allintrons, color = intron_color)+
            geom_rect(aes(xmin = txstart, xmax = codingstart, ymin = -1, ymax = 1), data = selection, fill = utr3_color)+
            geom_rect(aes(xmin = codingend, xmax = txend, ymin = -1, ymax = 1), data = selection, fill= "white", color = "white", lwd = 1.2)+
            geom_polygon(aes(x = x, y = y), data = endUTR, fill = utr5_color, color = "black")+
            theme_void()
        
        return(plot)
        
    } else {
        
        selection$txend <- as.numeric(selection$txend)
        selection$codingend <- as.numeric(selection$codingend)
        exonstarts <- matrix(as.numeric(unlist(strsplit(selection$exonstarts, ",")), ncol = selection$numexons))
        
        exonends <- matrix(as.numeric(unlist(strsplit(selection$exonends, ",")), ncol = selection$numexons))
        
        allexons <- data.frame("starts" = exonstarts, "ends" = exonends)
        
        intronstarts <- exonends[1:(nrow(exonends)-1)]
        
        intronends <- exonstarts[2:(nrow(exonends))]
        
        allintrons <- data.frame("starts" = intronstarts, "ends" = intronends) %>%
            dplyr::mutate(midpoint = (starts+ends)/2)
        
        endUTR <- data.frame("x" = c(selection$codingstart, selection$codingstart, selection$txstart), "y" = c(1, -1, 0))
        
        plot <- ggplot(allexons)+
            geom_segment(aes(x = txend, xend = txstart, y = 0, yend = 0), data = selection, color = "black")+
            geom_rect( aes(xmin =  starts, xmax = ends, ymin = -1 , ymax = 1), fill = gene_color, color = "black", alpha = gene_alpha)+
            geom_segment(aes(x = starts, y = 1, xend = midpoint, yend = 2), data = allintrons, color = intron_color)+
            geom_segment(aes(x = midpoint, y = 2, xend = ends, yend = 1), data = allintrons, color = intron_color)+
            geom_rect(aes(xmin = txend, xmax = codingend, ymin = -1, ymax = 1), data = selection, fill = utr3_color)+
            geom_rect(aes(xmin = codingstart, xmax = txstart, ymin = -1, ymax = 1), data = selection, fill= "white", color = "white", lwd = 1.2)+
            geom_polygon(aes(x = x, y = y), data = endUTR, fill = utr5_color, color = "black")+
            theme_void()
        
        return(plot)
    }
}

save(gene_model,file="../functions/gene_model.Rdata")

goi_list <- tibble(
  gene='',
  chr='',
  strand='+',
  txstart=0,
  txend=0,
  codingstart=txstart,
  codingend=txend,
  numexons=0,
  wbgene='Wb000',
  gene_name=gene,
  biotype='protein_coding',
  exonstarts='',
  exonends='',
  type='Transcript'
)
save(goi_list,file="../processed_data/goi_list.Rdata")
#what needs to be added: gene, chr, txtstart,txtend,codingstart,codingend,numexons,exonstarts,exonends,