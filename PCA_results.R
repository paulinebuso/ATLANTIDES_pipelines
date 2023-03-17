setwd("~/ATLANTIDES/PCA/test")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)


####-----------------------PCA------------------------####
#PCA variation
#read the first output containing the dregree of genetic variaiton per PCA axis. 
pca_variation <-data.table::fread("~/ATLANTIDES/PCA/PCA_4pops.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 

#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(pca_variation,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1)) +
  labs(x = "PC", y = "Variance proportion") +
  theme_bw()
  
#We'll use the two first PCA axis for the plot 
#PC1=56%
#PC2=19%

#Final plot
#read the second output containing the genetic variation (values per PCA axis + per sample)
pca <-data.table::fread("~/ATLANTIDES/PCA/test/PCA_4pops.eigenvec") %>% as_tibble %>% rename_all(~c("samples","samples2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) #for entire PCA plot

#remove the two first colunms because the samples are not written the way I want 
pca <- pca[,-2]
pca <- pca[,-1]

#read the file containing the right name of samples + population (farmed/wild euro/american)
pops_file <- read.table("~/ATLANTIDES/PCA/test/salmon_pops.tsv", header=FALSE, sep="\t") %>% as_tibble %>% rename_all(~c("Samples","Populations"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
pca <- cbind(pca,pops_file)

#Plot using th two first PCA variation axis
ggplot(pca, aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations)), alpha = 1, size = 1.5, show.legend = TRUE) +
  labs(x = "PC1 (56% variation)", y = "PC2 (19% variation)") + 
  theme_bw() +
  scale_color_manual(values = c("#fc2634", "#092fd9", "#fc79a7","#8aceeb")) 
  #stat_ellipse(geom = "polygon", aes(fill = Populations), alpha = 0.20) 

#save the plot
ggsave("PCA_4pops.pdf",width=7,height=5,dpi=600)
dev.off()

#colour ideas :
#ed4779","#ed0e0e","#75bbc7","#430aad"
#scale_color_manual(values = c("#f2ce05","#f26805","#810996","#dc51f5")) 
#scale_color_manual(values = c("#32a8a4","#32a850","#eb349b","#db72e0")) 


