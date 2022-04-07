#population plots

library("tidyverse")

cb_meanlen <- as_tibble(read_tsv('/Users/shrirambhat/Desktop/telomere/data/briggsae/lengths/brig_isotype.tsv'))

pb_mean <- cb_meanlen %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth = 0.5) +
  theme(text = element_text(size=5)) + 
  geom_vline(xintercept = 4,color="red",alpha=0.3)

pb_mean


ggsave("../plots/brig_population_count.png",width=2,height=2,units="in",dpi=600)
