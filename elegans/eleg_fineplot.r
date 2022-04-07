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


#C. elegans fine mapping gene plot


genes_in_region <- data.table::fread(glue::glue("/Users/shrirambhat/Desktop/telomere/data/elegans/nemascan_runs/Analysis_Results-20220308/Fine_Mappings/Data/length_II_13272357-15257187_bcsq_genes.tsv"))
#genes_in_region <- data.table::fread(glue::glue("{trait_name}_{QTL_chrom}_{QTL_start}-{QTL_end}_bcsq_genes.tsv"))
    
gene_df <- genes_in_region %>%
    dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
    dplyr::distinct(WBGeneID, GENE_NAME, PEAK_MARKER, CHROM, STRAND, TRANSCRIPTION_START_POS, 
                    TRANSCRIPTION_END_POS, START_POS, END_POS, VARIANT_LOG10p) %>%
    dplyr::group_by(WBGeneID) %>%
    dplyr::mutate(VARIANT_LOG10p = max(VARIANT_LOG10p, na.rm = T)) %>%
    dplyr::distinct()

peak_variant <- as.numeric(strsplit(unique(gene_df$PEAK_MARKER), split = ":")[[1]][2])

variant_df <- genes_in_region %>%
    dplyr::filter(START_POS == unique(genes_in_region$START_POS)) %>%
    dplyr::distinct(CHROM, POS, GENE_NAME, WBGeneID, VARIANT_LOG10p, VARIANT_IMPACT)

variant_df$VARIANT_IMPACT[is.na(variant_df$VARIANT_IMPACT)] <- "Intergenic"

xs <- unique(gene_df$START_POS)
xe <- unique(gene_df$END_POS)
cq <- unique(gene_df$CHROM)

max_logp <- unique(max(variant_df$VARIANT_LOG10p, na.rm = T))/150

gene_plot <- ggplot(gene_df) +
    aes(text = paste0(GENE_NAME, "\n", WBGeneID)) +
    geom_vline(aes(xintercept = peak_variant/1e6),
               linetype=1, color = "blue")+
    geom_segment(aes(x = ifelse(STRAND == "+", TRANSCRIPTION_START_POS/1e6, TRANSCRIPTION_END_POS/1e6),
                     xend = ifelse(STRAND == "+", TRANSCRIPTION_END_POS/1e6, TRANSCRIPTION_START_POS/1e6),
                     y = VARIANT_LOG10p,
                     yend = VARIANT_LOG10p),
                 arrow = arrow(length = unit(5, "points")), size = 1) +
    geom_segment(aes(x = POS/1e6,
                     xend = POS/1e6,
                     y = VARIANT_LOG10p+max_logp,
                     yend = VARIANT_LOG10p-max_logp,
                     color = VARIANT_IMPACT), data = variant_df) +
    # if snpeff, need to add back moderate and modifier
    scale_color_manual(values = c("LOW" = "gray30",
                                  "HIGH" = "red",
                                  "Linker" = "gray80", 
                                  "Intergenic" = "gray80"),
                       breaks = c("HIGH", "LOW", "Intergenic"),
                       name = "EFFECT")+
    labs(x = "Genomic Position (Mb)",
         y = expression(-log[10](italic(p)))) +
         # )+
    theme_bw(18)+
    xlim(c(xs/1e6, xe/1e6)) +
    theme(legend.position = "top",
          panel.grid = element_blank()) 
gene_plot
ggsave("../plots/eleg_fineplot2.png",dpi=300,height=2,width = 7.5, units = "in")