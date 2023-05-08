library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)


####-----------------------PCA----------------------------####

pca <-data.table::fread("~/ATLANTIDES/PCA/PCA_4pops.eigenvec") %>% as_tibble %>% rename_all(~c("samples","samples2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")) #for entire PCA plot
pca <- pca[,-2]
pca <- pca[,-1]

pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops_pca2.tsv", header=FALSE, sep="\t")  %>% as_tibble %>% rename_all(~c("samples","Populations2","Location","Populations","Type"))

#merge the two tables to have location info on the pca + legend by pop or type of samples
pca <- cbind(pca,pops_file)

#Plot
ggplot(pca, aes(x = PC1, y = PC2)) + geom_point((aes(colour = Populations, shape = Type)), alpha = 1, size = 2, show.legend = TRUE) +
  labs(x = "PC1 (53.85% variation)", y = "PC2 (10.74% variation)") + 
  theme_bw() +
  scale_color_manual(values = c("#ffb703","#092fd9","#6fc441","#fc2634","#fc79a7","#8aceeb")) 

ggsave("PCA_4pops.png",width=7,height=5,dpi=600)
dev.off()

#PCA variation
pca_variation <-data.table::fread("~/ATLANTIDES/PCA/PCA_4pops.eigenval") %>% rename (eigenval = V1) %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 

#eigenval file tell us the 10 principal components + second ligne "rename.. mutate(var.prop) allow to add a colomn and tell the number on the component

#Histogram
ggplot(pca_variation,aes(PC,var.prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1,10,1)) +
  labs(x = "PC", y = "Variance proportion") +
  theme_bw()


