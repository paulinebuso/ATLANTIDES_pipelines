setwd("~/ATLANTIDES/SNP_quality/celian")

install.packages("readr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)

##Genotype quality
#read data
var_qual <- read_delim("wild_european.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)
summary(var_qual)
Q1_qual<-summary(var_qual$qual)[2]

#PLOT

ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light() +
  geom_vline(xintercept=Q1_qual, size=1.0, color="orange")+ 
  xlim(0,100)+
  xlab("quality")+
  ylab("density") + 
  annotate(x=22,y=0.05,label=paste0("Q1 = ", Q1_qual),vjust=2,geom="label") 

## Mean read depth
#read data

var_depth <- read_delim("wild_european.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
summary(var_depth$mean_depth)

Q1_depth<-summary(var_depth$mean_depth)[2]

#NOTE: the x-axis (xlim) is limited because the extreme values are too high, so the graph is not readable if I do not limit it 

#depth ch1
depth1<-var_depth[var_depth$chr == "ssa01",]
summary(depth1$mean_depth)

Q1_depth_ch1<-summary(depth1$mean_depth)[2]
ggplot(depth1, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+ 
  geom_vline(xintercept=Q1_depth_ch1, size=1.0, color="orange")+
  xlim(0,50) +
  xlab("mean depth ch1")+
  ylab("density")+
  annotate(x=4,y=0.6,label=paste0("Q1 = ", Q1_depth_ch1),vjust=2,geom="label")


##depth ch2
depth2<-var_depth[var_depth$chr == "ssa02",]
summary(depth2$mean_depth)

Q1_depth_ch2<-summary(depth2$mean_depth)[2]

ggplot(depth2, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+ 
  geom_vline(xintercept=Q1_depth_ch2, size=1.0, color="orange")+
  xlim(0,50) +
  xlab("mean depth ch2")+
  ylab("density")+
  annotate(x=4,y=0.6,label=paste0("Q1 = ", Q1_depth_ch2),vjust=2,geom="label")

#ch3

depth3<-var_depth[var_depth$chr == "ssa03",]
summary(depth3$mean_depth)
Q1_depth_ch3<-summary(depth3$mean_depth)[2]
ggplot(depth3, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+ 
  geom_vline(xintercept=Q1_depth_ch3, size=1.0, color="orange")+
  xlim(0,50) +
  xlab("mean depth ch3")+
  ylab("density")+
  annotate(x=4,y=0.6,label=paste0("Q1 = ", Q1_depth_ch3),vjust=2,geom="label")

#ch4
depth4<-var_depth[var_depth$chr == "ssa04",]
summary(depth4$mean_depth)

#ch5
depth5<-var_depth[var_depth$chr == "ssa05",]
summary(depth5$mean_depth)

#ch6
depth6<-var_depth[var_depth$chr == "ssa06",]
summary(depth6$mean_depth)

#ch7
depth7<-var_depth[var_depth$chr == "ssa07",]
summary(depth7$mean_depth)

#ch8
depth8<-var_depth[var_depth$chr == "ssa08",]
summary(depth8$mean_depth)

#ch9
depth9<-var_depth[var_depth$chr == "ssa09",]
summary(depth9$mean_depth)

#ch10
depth10<-var_depth[var_depth$chr == "ssa10",]
summary(depth10$mean_depth)

#Stats depth --> juste pour savoir le seuil de filtration du premier quartile
depth11<-var_depth[var_depth$chr == "ssa11",]
summary(depth11$mean_depth) #Q1=6.16

depth12<-var_depth[var_depth$chr == "ssa12",]
summary(depth12$mean_depth) #Q1=6.97

depth13<-var_depth[var_depth$chr == "ssa13",]
summary(depth13$mean_depth) #Q1=7.45

depth14<-var_depth[var_depth$chr == "ssa14",]
summary(depth14$mean_depth) #Q1=6.91

depth15<-var_depth[var_depth$chr == "ssa15",]
summary(depth15$mean_depth) #Q1=6.41

depth16<-var_depth[var_depth$chr == "ssa16",]
summary(depth16$mean_depth) #Q1=6.62

depth17<-var_depth[var_depth$chr == "ssa17",]
summary(depth17$mean_depth) #Q1=4.08
Q1_depth_ch17<-summary(depth17$mean_depth)[2]
ggplot(depth17, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+ 
  geom_vline(xintercept=Q1_depth_ch17, size=1.0, color="orange")+
  xlim(0,100) +
  xlab("mean depth")+
  ylab("density")+
  annotate(x=4,y=0.6,label=paste0("Q1 = ", Q1_depth_ch17),vjust=2,geom="label")

depth18<-var_depth[var_depth$chr == "ssa18",]
summary(depth18$mean_depth) #Q1=6.63

depth19<-var_depth[var_depth$chr == "ssa19",]
summary(depth19$mean_depth) #Q1=7.0089

depth20<-var_depth[var_depth$chr == "ssa20",]
summary(depth20$mean_depth) #Q1=7.25

depth21<-var_depth[var_depth$chr == "ssa21",]
summary(depth21$mean_depth) #Q1=7.54

depth22<-var_depth[var_depth$chr == "ssa22",]
summary(depth22$mean_depth) #Q1=8.054

depth23<-var_depth[var_depth$chr == "ssa23",]
summary(depth23$mean_depth) #Q1=7.955

depth24<-var_depth[var_depth$chr == "ssa24",]
summary(depth24$mean_depth) #Q1=7.768

depth25<-var_depth[var_depth$chr == "ssa25",]
summary(depth25$mean_depth) #Q1=7.06

depth26<-var_depth[var_depth$chr == "ssa26",]
summary(depth26$mean_depth) #Q1=6.39

depth27<-var_depth[var_depth$chr == "ssa27",]
summary(depth27$mean_depth) #Q1=6.94

depth28<-var_depth[var_depth$chr == "ssa28",]
summary(depth28$mean_depth) #Q1=7.23

depth29<-var_depth[var_depth$chr == "ssa29",]
summary(depth29$mean_depth) #Q1=6.67

## Missingness Variant

#read data
var_miss <- read_delim("wild_european.lmiss", delim = "\t", col_names = c("chr", "pos", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

summary(var_miss$fmiss)

Q3_miss<-summary(var_miss$fmiss)[5]
    
ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+
  geom_vline(xintercept=Q3_miss, size=1.0, color="orange")+
  xlim(0,0.4) +
  xlab("percent missing")+
  ylab("density")+
  annotate(x=0.1,y=25,label=paste0("Q3 = ", Q3_miss),vjust=2,geom="label")

#ch1
miss1<-var_miss[var_miss$chr == "ssa01",]

summary(miss1$fmiss)

Q3_miss_ch1<-summary(miss1$fmiss)[5]

ggplot(miss1, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+
  geom_vline(xintercept=Q3_miss_ch1, size=1.0, color="orange")+
  xlim(0,0.4) +
  xlab("percent missing")+
  ylab("density")+
  annotate(x=0.1,y=25,label=paste0("Q3 = ", Q3_miss_ch1),vjust=2,geom="label")

#ch2
miss2<-var_miss[var_miss$chr == "ssa02",]

summary(miss2$fmiss)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.00000 0.08935 0.03571 0.99107 

Q3_miss_ch2<-summary(miss2$fmiss)[5]

ggplot(miss2, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light()+
  geom_vline(xintercept=Q3_miss_ch2, size=1.0, color="orange")+
  xlim(0,0.4) +
  xlab("percent missing")+
  ylab("density")+
  annotate(x=0.1,y=25,label=paste0("Q3 = ", Q3_miss_ch2),vjust=2,geom="label")

## Minor allele frequency

var_freq <- read_delim("wild_european.frq", delim = "\t", col_names = c("chr", "pos", "snpcount", "variant.kb", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

# find minor allele frequency

summary(var_freq$maf)

ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light() +
  xlab("minor allele frequency")+
  ylab("density") +
  xlim(0,1)

ggplot(var_freq, aes(maf)) + geom_histogram(bins=30, color="black", fill="white") +
  theme_light() +
  xlim(0,1) 

#optional : add xlim(0.01,0.5)

#SNP density maf > 0.01
    
var_freq_0_01<-var_freq[var_freq$maf>0.01, ] #subset of only > 0.01 values
ggplot(var_freq_0_01, aes(pos)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+
  theme_light() +
  xlab("position")+
  ylab("density")
  
#SNP density scatterplot (100kb windows)
#read data

var_density <- read_delim("wild_european.snpden", delim = "\t", col_names = c("chr", "pos_start", "SNP_count", "vars_per_kb"), skip = 1)
    
summary(var_density$vars_per_kb)

ggplot(var_density, aes(pos_start, vars_per_kb))+geom_point(size=1)+
  theme_light() +
  xlab("position")+
  ylab("SNPs/kb") +
  facet_wrap(~chr,scales="free_x")

#SNP_density maf 0.01 scatterplot (100kb windows)

var_density_maf_0.01 <- read_delim("farmed_american_quality_maf0.01.snpden", delim = "\t", col_names = c("chr", "pos_start", "SNP_count", "vars_per_kb"), skip = 1)

summary(var_density_maf_0.01$vars_per_kb)
    
ggplot(var_density_maf_0.01, aes(pos_start, vars_per_kb))+geom_point(size=1, color="blue")+ 
  theme_light() +
  xlab("position")+
  ylab("SNPs/kb") +
  facet_wrap(~chr,scales="free_x")

## Singletons


var_singleton <- read_delim("wild_european.singletons", delim = "\t", col_names = c("chr", "pos", "singleton/doubleton", "allele", "indv"), skip = 1)

ggplot(var_singleton, aes(pos)) + geom_density(fill = "dodgerblue1", colour = "blue", alpha = 0.3) +
  theme_light() +
  xlab("position")+
  ylab("density") +
  facet_wrap(~chr,scales="free_x")

#ch1

s1<-var_singleton[var_singleton$chr == "ssa01",]

ggplot(s1, aes(pos)) + geom_density(fill = "dodgerblue1", colour = "blue", alpha = 0.3) +
  theme_light() +
  xlab("position")+
  ylab("density")
#facet_wrap(~chr,scales="free_x")

summary(s1)

#ch2
s2<-var_singleton[var_singleton$chr == "ssa02",]

ggplot(s2, aes(pos)) + geom_density(fill = "dodgerblue1", colour = "blue", alpha = 0.3) +
  theme_light() +
  xlab("position")+
  ylab("density")
#facet_wrap(~chr,scales="free_x")

summary(s1)

