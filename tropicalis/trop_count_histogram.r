#population plots

library("tidyverse")

pt_meanlen <- as_tibble(read_tsv('/Users/shrirambhat/Desktop/telomere/data/tropicalis/lengths/trop_isotype.tsv'))

pt_mean <- pt_meanlen %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth = 0.5) +
  theme(text = element_text(size=7)) + 
  geom_vline(xintercept = 15,color="red",alpha=0.3)

pt_mean


ggsave("../plots2/trop_population_count.png",width=2,height=2,units="in",dpi=300)
