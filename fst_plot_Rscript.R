
library(ggplot2)
library(tidyverse)
library(dplyr)

###FST Wild Norway - Farmed NORWAY

#fst_Euro <- read.table("Fst_stats_Euro.weir.fst", header=TRUE) %>% as_tibble %>% rename_all(~c("CHROM","POS","FST"))

fst_Euro_window <- read.table("~/ATLANTIDES/Fst/Fst_stats_Euro_window5000.windowed.weir.fst", header=TRUE) %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))

#data preparation pour plot

data_cum <- fst_Euro_window %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(BIN_START)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)

fst_Euro_window <- fst_Euro_window %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(bp_cum = BIN_START + bp_add)

axis_set <- fst_Euro_window %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))

quantile(fst_Euro_window$WEIGHTED_FST, c(.01, .99))
quantile_fst_euro=0.286489820

fst_plot <- ggplot(fst_Euro_window, aes(x = bp_cum, y = WEIGHTED_FST,color = as.factor(CHROM))) +
  geom_point(alpha = 0.75,size=0.80) + ylim(c(0,1)) + 
  scale_x_continuous(label =axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#419c62", "#70d495"), unique(length(axis_set$CHROM)))) +
  labs(y="FST",x="Chromosome") + 
  theme(legend.position ="none", axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = quantile_fst_euro, color = "red", linetype = "dashed") +
  ggtitle("FST - European samples")

fst_plot

###FST-American

fst_american_window <- read.table("~/ATLANTIDES/Fst/Fst_stats_American.windowed.weir.fst", header=TRUE) %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))

#data preparation pour plot

data_cum2 <- fst_american_window %>% 
  group_by(CHROM) %>% 
  summarise(max_bp2 = max(BIN_START)) %>% 
  mutate(bp_add2 = lag(cumsum(as.numeric(max_bp2)), default = 0)) %>% 
  select(CHROM, bp_add2)

fst_american_window <- fst_american_window %>% 
  inner_join(data_cum2, by = "CHROM") %>% 
  mutate(bp_cum2 = BIN_START + bp_add2)

axis_set2 <- fst_american_window %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum2))

quantile(fst_american_window$WEIGHTED_FST, c(.01, .99))
quantile_fst_american=0.4794965000

fst_plot_american <- ggplot(fst_american_window, aes(x = bp_cum2, y = WEIGHTED_FST,color = as.factor(CHROM))) +
  geom_point(alpha = 0.75,size=0.80) + ylim(c(0,1)) + 
  scale_x_continuous(label =axis_set2$CHROM, breaks = axis_set2$center) +
  scale_color_manual(values = rep(c("#419c62", "#70d495"), unique(length(axis_set2$CHROM)))) +
  labs(y="FST",x="Chromosome") + 
  theme(legend.position ="none", axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = quantile_fst_american, color = "red", linetype = "dashed") +
  ggtitle("FST - North American samples")

fst_plot_american

###FST - plot the two datatsets
library(gridExtra)

grid.arrange(fst_plot, fst_plot_american)
