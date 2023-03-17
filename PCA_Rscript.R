setwd("~/ATLANTIDES/PCA/test")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)


####-----------------------PCA-FINALE----------------------------####

pca <-data.table::fread("~/ATLANTIDES/PCA/test/PCA_4pops.eigenvec") %>% as_tibble %>% rename_all(~c("samples","samples2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) #for entire PCA plot
pca <- pca[,-2]
pca <- pca[,-1]

pops_file <- read.table("~/ATLANTIDES/PCA/test/salmon_pops.tsv", header=FALSE, sep="\t")  %>% as_tibble %>% rename_all(~c("V1","Populations"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
pca <- cbind(pca,pops_file)

#Plot

ggplot(pca, aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), alpha = 1, size = 1.5, show.legend = TRUE) +
  labs(x = "PC1 (56% variation)", y = "PC2 (19% variation)") + 
  theme_bw() +
  scale_color_manual(values = c("#fc2634", "#092fd9", "#fc79a7","#8aceeb")) 
  #stat_ellipse(geom = "polygon", aes(fill = Populations), alpha = 0.20) 

ggsave("PCA_4pops.pdf",width=7,height=5,dpi=600)
dev.off()

#ed4779","#ed0e0e","#75bbc7","#430aad"
#scale_color_manual(values = c("#f2ce05","#f26805","#810996","#dc51f5")) 
#scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#PCA variation
pca_variation <-data.table::fread("~/ATLANTIDES/PCA/PCA_4pops.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 

#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(pca_variation,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1)) +
  labs(x = "PC", y = "Variance proportion") +
  theme_bw()
  
#PC1=56%
#PC2=19%

###PCA - 4 pop - by ch###

MOWI_pca_ssa02 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa02.eigenvec")  #for entire PCA plot

MOWI_pca_ssa02_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa02.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa02_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=66%
#PC2=15%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa02 <- merge(MOWI_pca_ssa02, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa02,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (66% variation)", y = "PC2 (15% variation") + 
  theme_bw() + 
  ggtitle("Ch2 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-----------------------------------

MOWI_pca_ssa03 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa03.eigenvec")  #for entire PCA plot

MOWI_pca_ssa03_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa03.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa03_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=71%
#PC2=12%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa03 <- merge(MOWI_pca_ssa03, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa03,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (71% variation)", y = "PC2 (12% variation") + 
  theme_bw() + 
  ggtitle("Ch3 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#------------------------------------

MOWI_pca_ssa04 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa04.eigenvec")  #for entire PCA plot

MOWI_pca_ssa04_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa04.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa04_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=9%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa04 <- merge(MOWI_pca_ssa04, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa04,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (9% variation") + 
  theme_bw() + 
  ggtitle("Ch4 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------------------------
MOWI_pca_ssa05 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa05.eigenvec")  #for entire PCA plot

MOWI_pca_ssa05_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa05.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa05_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=8%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa05 <- merge(MOWI_pca_ssa05, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa05,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch5 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#--------------------------------
MOWI_pca_ssa12 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa12.eigenvec")  #for entire PCA plot

MOWI_pca_ssa12_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa12.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa12_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=71%
#PC2=9%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa12 <- merge(MOWI_pca_ssa12, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa12,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (71% variation)", y = "PC2 (9% variation") + 
  theme_bw() + 
  ggtitle("Ch12 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#----------------------

MOWI_pca_ssa17 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa17.eigenvec")  #for entire PCA plot

MOWI_pca_ssa17_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa17.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa17_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=68%
#PC2=12%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa17 <- merge(MOWI_pca_ssa17, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa17,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (68% variation)", y = "PC2 (12% variation") + 
  theme_bw() + 
  ggtitle("Ch17 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------------

MOWI_pca_ssa18 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa18.eigenvec")  #for entire PCA plot

MOWI_pca_ssa18_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa18.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa18_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=65%
#PC2=10%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa18 <- merge(MOWI_pca_ssa18, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa18,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (65% variation)", y = "PC2 (10% variation") + 
  theme_bw() + 
  ggtitle("Ch18 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 


#----------------------------------

MOWI_pca_ssa19 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa19.eigenvec")  #for entire PCA plot

MOWI_pca_ssa19_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa19.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa19_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=70%
#PC2=7%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa19 <- merge(MOWI_pca_ssa19, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa19,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (70% variation)", y = "PC2 (7% variation") + 
  theme_bw() + 
  ggtitle("Ch19 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#----------------------------

MOWI_pca_ssa24 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa24.eigenvec")  #for entire PCA plot

MOWI_pca_ssa24_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa24.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa02_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=75%
#PC2=7%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa24 <- merge(MOWI_pca_ssa24, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa24,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (75% variation)", y = "PC2 (7% variation") + 
  theme_bw() + 
  ggtitle("Ch24 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#---------------------------


MOWI_pca_ssa25 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa25.eigenvec")  #for entire PCA plot

MOWI_pca_ssa25_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa25.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa25_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=8%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa25 <- merge(MOWI_pca_ssa25, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa25,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch25 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------

MOWI_pca_ssa26 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa26.eigenvec")  #for entire PCA plot

MOWI_pca_ssa26_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa26.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa26_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=71%
#PC2=10%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa26 <- merge(MOWI_pca_ssa26, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa26,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (71% variation)", y = "PC2 (10% variation") + 
  theme_bw() + 
  ggtitle("Ch26 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#--------------------

MOWI_pca_ssa27 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa27.eigenvec")  #for entire PCA plot

MOWI_pca_ssa27_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa27.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa27_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=8%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa27 <- merge(MOWI_pca_ssa27, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa27,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch27 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#---------------------------------

MOWI_pca_ssa28 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa28.eigenvec")  #for entire PCA plot

MOWI_pca_ssa28_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa28.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa28_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=78%
#PC2=6%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa28 <- merge(MOWI_pca_ssa28, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa28,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (78% variation)", y = "PC2 (6% variation") + 
  theme_bw() + 
  ggtitle("Ch28 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#----------------------------


MOWI_pca_ssa29 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa29.eigenvec")  #for entire PCA plot

MOWI_pca_ssa29_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa29.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa29_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=75%
#PC2=7%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa29 <- merge(MOWI_pca_ssa29, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa29,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (75% variation)", y = "PC2 (7% variation") + 
  theme_bw() + 
  ggtitle("Ch29 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#---------------------


MOWI_pca_ssa08 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa08.eigenvec")  #for entire PCA plot

MOWI_pca_ssa08_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa08.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa08_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=59%
#PC2=15%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa08 <- merge(MOWI_pca_ssa08, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa08,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (59% variation)", y = "PC2 (15% variation") + 
  theme_bw() + 
  ggtitle("Ch08 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------


MOWI_pca_ssa07 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa07.eigenvec")  #for entire PCA plot

MOWI_pca_ssa07_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa07.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa07_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=8%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa07 <- merge(MOWI_pca_ssa07, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa07,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch07 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------

MOWI_pca_ssa13 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa13.eigenvec")  #for entire PCA plot

MOWI_pca_ssa13_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa13.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa13_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=74%
#PC2=7%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa13 <- merge(MOWI_pca_ssa13, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa13,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (74% variation)", y = "PC2 (4% variation") + 
  theme_bw() + 
  ggtitle("Ch13 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------

MOWI_pca_ssa16 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa16.eigenvec")  #for entire PCA plot

MOWI_pca_ssa16_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa16.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa16_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=71%
#PC2=8%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa16 <- merge(MOWI_pca_ssa16, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa16,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (71% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch16 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

#-------------------------

MOWI_pca_ssa23 <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa23.eigenvec")  #for entire PCA plot

MOWI_pca_ssa23_var <-data.table::fread("~/ATLANTIDES/PCA/pca_ssa23.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 
#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(MOWI_pca_ssa23_var,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1))
#PC1=73%
#PC2=8%
#PC3=5%
#PC4=3%

#input the populations file 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("#IID","Location","Populations","Origin"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged_ssa23 <- merge(MOWI_pca_ssa23, pops_file, by ="#IID")

#Plot
ggplot(PCA_file_merged_ssa23,aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC1 (73% variation)", y = "PC2 (8% variation") + 
  theme_bw() + 
  ggtitle("Ch23 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 

ggplot(PCA_file_merged_ssa23,aes(x = PC4, y = PC3)) + geom_point((aes(colour = Populations)), size = 1, show.legend = TRUE) +
  labs(x = "PC4 (3% variation)", y = "PC3 (5% variation") + 
  theme_bw() + 
  ggtitle("Ch23 - PCA - 4pops") +
  scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 
