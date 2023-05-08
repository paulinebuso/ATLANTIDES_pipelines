

library(ggplot2)
library(tidyverse)
library(dplyr)

###FST Wild Norway - Farmed NORWAY

fst_euro <- read.table("~/ATLANTIDES/Fst/Fst_stats_Euro.windowed.weir.fst", header=TRUE) %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))

#data preparation pour plot
#remove ssa and ssa0

fst_euro <- fst_euro %>%
  dplyr::mutate(CHROM=str_remove_all(CHROM,"ssa")) %>% 
  dplyr::mutate_at('CHROM',as.numeric)

quantval<-quantile(fst_euro$WEIGHTED_FST, c(0.01, 0.99)) # calculate quantiles to have significance lines
quartile1 <- quantval[1]
quartile2 <- quantval[2]

#add space between chromosome + draw a beautiful axis with cumulative position
#chromosome lenght 
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+20000000)[1:28])
chromname<-unique(fst_euro$CHROM)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
fst_euro$newpos<-fst_euro$BIN_START+listval[match(fst_euro$CHROM,listval$chromname),3]
##in this plot:
## color chrom alternatively
## label x axis with chrom name
## remove legend using theme(legend.position="none")
euro_plot <- ggplot(fst_euro, aes(x = newpos, y = WEIGHTED_FST,color = as.factor(CHROM))) + 
  scale_color_manual(values = rep(c("#1c1616", "#8a8484"), length(unique(fst_euro$CHROM))))+
  geom_point(alpha = 0.90,size=1)+ 
  theme_classic() +
  geom_hline(yintercept = quartile2, color = "red", linetype = 1, linewidth = 0.6)+
  scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+
  theme(legend.position="none")+
  labs(y="Fst",x="Chromosome")+
  ylim(c(0,1))
 
euro_plot 

ggsave("euro_fst_plot.png",width=12,height=5,dpi=300)
dev.off()

###American

fst_american <- read.table("~/ATLANTIDES/Fst/Fst_stats_American.windowed.weir.fst", header=TRUE) %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","BIN_END","N_VARIANTS","WEIGHTED_FST","MEAN_FST"))
fst_american <- fst_american %>%
  dplyr::mutate(CHROM=str_remove_all(CHROM,"ssa")) %>% 
  dplyr::mutate_at('CHROM',as.numeric)

quantval<-quantile(fst_american$WEIGHTED_FST, c(0.01, 0.99)) # calculate quantiles to have significance lines
quartile3 <- quantval[1]
quartile4 <- quantval[2]

#add space between chromosome + draw a beautiful axis with cumulative position
#chromosome lenght 
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+20000000)[1:28])
chromname<-unique(fst_american$CHROM)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
fst_american$newpos<-fst_american$BIN_START+listval[match(fst_american$CHROM,listval$chromname),3]
##in this plot:
## color chrom alternatively
## label x axis with chrom name
## remove legend using theme(legend.position="none")
american_plot <- ggplot(fst_american, aes(x = newpos, y = WEIGHTED_FST,color = as.factor(CHROM))) + 
  scale_color_manual(values = rep(c("#1c1616", "#8a8484"), length(unique(fst_american$CHROM))))+
  geom_point(alpha = 0.90,size=1)+ 
  theme_classic() +
  geom_hline(yintercept = quartile4, color = "red", linetype = 1, linewidth = 0.6)+
  scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+
  theme(legend.position="none")+
  labs(y="Fst",x="Chromosome")+
  ylim(c(0,1))


american_plot 

ggsave("american_fst_plot.png",width=12,height=5,dpi=300)
dev.off()

###other method : 

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
  scale_color_manual(values = rep(c("#1c1616", "#8a8484"), unique(length(axis_set2$CHROM)))) +
  labs(y="FST",x="Chromosome") + 
  theme(legend.position ="none", axis.text.x = element_text(angle = 60, hjust = 1)) +
  geom_hline(yintercept = quantile_fst_american, color = "red", linetype = "dashed") +
  ggtitle("FST - North American samples")

fst_plot_american

###FST - plot the two datatsets
library(gridExtra)

grid.arrange(fst_plot, fst_plot_american)
