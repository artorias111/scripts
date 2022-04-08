#population plots

library("tidyverse")

ce_meanlen <- as_tibble(read_tsv('/Users/shrirambhat/Desktop/telomere/data/elegans/lengths/eleg_isotype.tsv'))

pe_mean <- ce_meanlen %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth = 1) +
  theme(text = element_text(size=7)) + 
  geom_vline(xintercept = 16,color="red",alpha=0.3)

pe_mean


ggsave("../plots2/eleg_population_count.png",width=2,height=2,units="in",dpi=300)
