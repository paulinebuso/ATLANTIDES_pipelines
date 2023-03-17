#Call some libraries we need

library(ggplot2)
library(tidyverse)
library(dplyr)

#Download the output from bash pipeline. It contains chromosome name, start position of the window, the number of SNPs per window, tajima's values.
tajima_euro <- read.table("~/ATLANTIDES/Tajima/Europe/tajima_farmed_euro.Tajima.D") %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","N_SNPS","TajimaD"))

#Data preparation
#install.packages("stringr")

library(stringr) #library for data editing, with tidyverse package.

tajima_euro <- tajima_euro %>%
  dplyr::mutate(CHROM=str_remove_all(CHROM,"ssa")) %>% 
  dplyr::mutate_at('CHROM',as.numeric) %>%
  dplyr::mutate_at('TajimaD', ~tidyr::replace_na(.,0))
#here we used dplyrs package attached to stringr and tidyverse packages. 
#line 1 : it allows to remove the ssa in chromosome name colum. 
#line 2 : allow to transform the character in chromosome nam colum as numeric values. Both lines are useful for the plot legend.
#line 3 : transform the NA to 0 because we'll calculate later the quantiles in order to fix a significant treshold. Quantiles calculation does not support NA and 0 values.

tajima_euro<-tajima_euro[tajima_euro$TajimaD!=0,] #remove "0"

which(is.na(tajima_euro), arr.ind=TRUE) #verify we have no NA

#Q-Q plot 
library(qqman)
qq(tajima_euro$TajimaD, main = "Q-Q plot Tajima D values", cex = 0.8)

###Manhattan plot with ggplot2 package 

#Create a column with the accumulated position for the x-axis. "= bp_cum"
data_cum <- tajima_euro %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(BIN_START)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)

tajima_euro <- tajima_euro %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(bp_cum = BIN_START + bp_add)

#Get the center of chromosomes
axis_set <- tajima_euro %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))

ylim <- tajima_euro %>% 
  filter(TajimaD == min(TajimaD)) %>% 
  mutate(ylim = abs(floor(TajimaD)) + 2) %>% 
  pull(ylim)

#Calcul of the treshold
quantile(tajima_euro$TajimaD, c(.01, .99)) # calculate quantiles to have significance lines 
#Here we define anything larger than 1% of the distribution as “extraordinarily large D” (i.e., potential regions under non-neutral evolution). 
#Other thresholds, for example, 0.1%, are also acceptable (nobody knows the “correct answer”)... 
#but we need to do multiple analyses to find a region that is “possibly truly” adaptive. 
#We define the top 1% as “candidate adaptive regions” for the time being, and will find overlap with candidate regions from other analyses later.

#create an object with quantiles 0.01 and 0.99 values for the significant lines in the plot
quartile1_tajima <- -1.489670
quartile2_tajima <- 4.899038

#Plot
euro_plot <- ggplot(tajima_euro, aes(x = bp_cum, y =TajimaD,color = as.factor(CHROM))) +
  geom_point(alpha = 0.75,size=0.80) + ylim(c(-4,6)) + 
  geom_hline(yintercept = quartile1_tajima, color = "red", linetype = "dashed") +
  geom_hline(yintercept = quartile2_tajima, color = "red", linetype = "dashed") +
  scale_x_continuous(label =axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#3699ad", "#ab0c36"), unique(length(axis_set$CHROM)))) +
  scale_size_continuous(range = c(0.4,2)) +
  labs(y="Tajima's D p-value",x="Chromosome") + 
  theme(legend.position="none") +
  ggtitle("Tajima's D - Farmed European samples")

euro_plot

###Marie's plot : allow us to see Tajima's distribution on each ch independantly
#Full plot
ggplot(tajime_euro, aes(x= BIN_START, y=TajimaD)) +
  geom_point(alpha = 0.75,size=0.80, color = "purple") +
  facet_wrap(~CHROM,scales="free_x") + 
  theme_bw() +
  labs(y="Tajima's D value",x="Position") +
  geom_hline(yintercept = quartile1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = quartile2, color = "red", linetype = "dashed")

#on 3 ch only
Marie_3ch<-tajima_euro[tajima_euro$CHROM%in%c("19", "20", "21"),]

ggplot(Marie_3ch, aes(x= BIN_START, y=TajimaD)) +
  geom_point(alpha = 0.75,size=0.80, color = "purple") +
  facet_grid(rows = vars(CHROM)) + 
  theme_bw() +
  labs(y="Tajima's D value",x="Position") 

ggplot(Marie_3ch, aes(x= BIN_START, y=TajimaD)) +
  geom_point(alpha = 0.75,size=0.80, color = "purple") +
  facet_wrap(~CHROM,scales="free_x") + 
  theme_bw() +
  labs(y="Tajima's D value",x="Position") +
  geom_hline(yintercept = quartile1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = quartile2, color = "red", linetype = "dashed")

###TajimaD Histogram distribution
ggplot(tajima_euro, aes(TajimaD)) + geom_histogram(color = "blue", fill="white") + 
  ylab("Frequency") + 
  xlab("TajimaD") +
  theme_bw() 

###ChromoMAP : to see the chromosomal regions with the higest/lowest values of Tajima. Facultative.
#Need two files : a chromosome file (name, start size, end size) and file with chromosomes annotations (chromosome annotation (ssa+bin start pasted), chromosome number (ssa..),2* Bin start, values of Tajima) :
#load the chromosome file (ssa01...+ start size and end size).

setwd("~/ATLANTIDES/Tajima_test")
CH_file <- read.table("chromosome_file.csv") %>% as_tibble %>% rename_all(~c("CHROM","SIZE_START","SIZE_END"))

#Remove "0" values of TajimaD file by creating a new file (useful to not crash the initial one :D). If not, ChromoMap does not work
TajimaD2<-TajimaD[TajimaD$N_SNPS!=0,]

#Make an unique name of each sample by pasting chromosome number and bin start of TajimaD2 file
TajimaD3<-paste(TajimaD2$CHROM,TajimaD2$BIN_START, sep="_", collaspse=NULL)

#Write a new table as the "chromosome annotations file : 1.Chromosome name (TajimaD3), 2.Chromosome number (ssa01), 3.binstart 4. binstart 5.binsend 6.N-SNP 7.TajimaD values
Annotations_tab<-cbind.data.frame(TajimaD3,TajimaD2$CHROM,TajimaD2$BIN_START,TajimaD2$BIN_START,TajimaD2$TajimaD)
write.table(Annotations_tab, file = "~/ATLANTIDES/Tajima_test/Annotations_tab.csv", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

#Draw the plot 
install.packages("chromoMap")
library(chromoMap)
chromoMap("chromosome_file.csv","Annotations_tab.csv", data_based_color_map = T, heat_map = TRUE, data_type = "numeric",  plots = "bar")

#Pareil mais avec que deux chromosomes pour pas que ca rame 
CH_file <- read.table("chromosome_file.csv") %>% as_tibble %>% rename_all(~c("CHROM","SIZE_START","SIZE_END"))

Two_CH<-CH_file[CH_file$CHROM%in%c("ssa01", "ssa02"),]
write.table(Two_CH, file = "~/ATLANTIDES/Tajima_test/Two_CH.csv", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

TajimaD_chromomap<-TajimaD[TajimaD$N_SNPS!=0,] 

Tajima_2CH<-TajimaD_chromomap[TajimaD_chromomap$CHROM%in%c("ssa01", "ssa02"),]
write.table(Tajima_2CH, file = "~/ATLANTIDES/Tajima_test/Tajima_2CH.csv", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

Chr_names<-paste(Tajima_2CH$CHROM,Tajima_2CH$BIN_START, sep="_", collaspse=NULL)

Annotations_tab_2CH<-cbind.data.frame(Chr_names,Tajima_2CH$CHROM,Tajima_2CH$BIN_START,Tajima_2CH$BIN_START,Tajima_2CH$TajimaD)
write.table(Annotations_tab_2CH, file = "~/ATLANTIDES/Tajima_test/Annotations_tab_2CH.csv", append = FALSE, quote = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

install.packages("chromoMap")
library(chromoMap)
#chromoMap("Two_CH.csv","Annotations_tab_2CH.csv", data_based_color_map = T, heat_map = TRUE, data_type = "numeric",  plots = "bar")
#not enough contrast 

chromoMap("Two_CH.csv","Annotations_tab_2CH.csv", data_based_color_map = T, data_colors = list(c("red","yellow")), data_type = "numeric", plots = "scatter", plot_filter = list(c("lt",0,"red")))

####-----------------AMERICAN-----------------####Same as previously with american samples.

library(ggplot2)
library(tidyverse)
library(dplyr)

tajima_american <- read.table("~/ATLANTIDES/Tajima/American/tajima_farmed_american.Tajima.D") %>% as_tibble %>% rename_all(~c("CHROM","BIN_START","N_SNPS","TajimaD"))

###Data preparation

#install.packages("stringr")
library(stringr)
tajima_american <- tajima_american %>%
  dplyr::mutate(CHROM=str_remove_all(CHROM,"ssa")) %>% 
  dplyr::mutate_at('CHROM',as.numeric) %>%
  dplyr::mutate_at('TajimaD', ~tidyr::replace_na(.,0))

tajima_american<-tajima_american[tajima_american$TajimaD!=0,] #remove "0"

which(is.na(tajima_american), arr.ind=TRUE) #verify we have no NA

###Q-Q plot 
library(qqman)
qq(tajima_american$TajimaD, main = "Q-Q plot Tajima D values", cex = 0.8)

###Manhattan plot with ggplot2 package 

data_cum <- tajima_american %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(BIN_START)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHROM, bp_add)

tajima_american <- tajima_american %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(bp_cum = BIN_START + bp_add)

axis_set <- tajima_american %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))

ylim <- tajima_american %>% 
  filter(TajimaD == min(TajimaD)) %>% 
  mutate(ylim = abs(floor(TajimaD)) + 2) %>% 
  pull(ylim)

quantile(tajima_american$TajimaD, c(.01, .99)) # calculate quantiles to have significance lines 
#Here we define anything larger than 1% of the distribution as “extraordinarily large D” (i.e., potential regions under non-neutral evolution). 
#Other thresholds, for example, 0.1%, are also acceptable (nobody knows the “correct answer”)... 
#but we need to do multiple analyses to find a region that is “possibly truly” adaptive. 
#We define the top 1% as “candidate adaptive regions” for the time being, and will find overlap with candidate regions from other analyses later.

quartile1 <- -2.120700
quartile2 <- 3.900521


#Plot
us_plot <- ggplot(tajima_american, aes(x = bp_cum, y =TajimaD,color = as.factor(CHROM))) +
  geom_point(alpha = 0.75,size=0.80) + ylim(c(-4,6)) + 
  geom_hline(yintercept = quartile1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = quartile2, color = "red", linetype = "dashed") +
  scale_x_continuous(label =axis_set$CHROM, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("#3699ad", "#ab0c36"), unique(length(axis_set$CHROM)))) +
  scale_size_continuous(range = c(0.4,2)) +
  labs(y="Tajima's D p-value",x="Chromosome") + 
  theme(legend.position="none") +
  ggtitle("Tajima's D - Farmed American samples")

us_plot

####--------------Europe - American - comparison--------------####
library(gridExtra)
grid.arrange(euro_plot, us_plot)
#here we plot the two previous plot showing Tajima's D distribution in each farmed population and compare the distribution
