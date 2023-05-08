setwd("~/ATLANTIDES/admixture")

library(stringr)
library(ggplot2)
library(dplyr)

cv <- read.table("cv_error.txt", header=TRUE) 

ggplot(cv,aes(x=K,y=CV_error)) + geom_line()+scale_x_continuous(breaks=seq(1,16,1))+
  labs(title="Cross-validation error - ADMIXTURE", xlab="K value", ylab="CV error") +
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))+
  theme_classic()

ggsave("Admixture_cross-validation.png",width=7,height=5,dpi=600)
dev.off()

#pophelper
install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
remotes::install_github('royfrancis/pophelper')
library(pophelper)
setwd("~/ATLANTIDES/admixture")
slist= readQ(files=c("input_admixture.2.Q"))
slist2= readQ(files=c("input_admixture.3.Q"))
slist3= readQ(files=c("input_admixture.4.Q"))
slist4= readQ(files=c("input_admixture.5.Q"))
slist5= readQ(files=c("input_admixture.6.Q"))

twolabset <- data.frame(Population= c(rep("MOWI", 112), 
                                      rep("Wild Norwegian", 98), 
                                      rep("Gaspe New Brunswick", 1), 
                                      rep("St John River", 13),
                                      rep("Gaspe New Brunswick", 6),
                                      rep("Penobscot River", 7),
                                      rep("St John River", 2),
                                      rep("Penobscot River", 2),
                                      rep("St John River", 1),
                                      rep("Gaspe New Brunswick", 2),
                                      rep("St John River", 29),
                                      rep("Penobscot River", 2),
                                      rep("Gaspe New Brunswick", 2),
                                      rep("St John River", 1),
                                      rep("Gaspe New Brunswick", 5),
                                      rep("St John River", 7),
                                      rep("Wild Canadian", 80)),
                                Loc=c(rep("Farmed Norwegian", 112),
                                      rep("Wild Norwegian", 98),
                                      rep("Farmed Canadian", 80),
                                      rep("Wild Canadian", 80)),
                                      stringsAsFactors = FALSE)
                            

p <- plotQ(
  qlist = c(slist3,slist4,slist5), 
  imgoutput = "join",
  clustercol = c(c("#fc2634","#fc79a7","#f7dc6f","#8aceeb","#092fd9","#fae8e9")),
  sortind = NA,
  grplab = twolabset,
  selgrp = "Loc",
  ordergrp = TRUE,
  panelspacer = 0.1,
  showsp = FALSE,
  splab = c("K=4","K=5","Ancestry K=6"),
  splabsize = 3,
  grplabspacer = -0.6,
  grplabheight = 0,
  grplabpos = 0.8,
  grplabsize = 0.7,
  grplabangle = 30,
  showindlab = FALSE,
  useindlab = TRUE,
  indlabheight = 0.7,
  indlabsize = 1.5,
  indlabangle = 90,
  indlabvjust = 1,
  indlabhjust = 1,
  indlabcol = "grey30",
  pointsize = 1.5,
  linepos = 1,
  linesize = 0.2,
  divgrp = c("Population", "Loc"),
  divcol = "Black",
  divtype = 2,
  divsize = 0.2,
  barsize = 1,
  barbordersize = 0,
  barbordercolour = "Black",
  showyaxis = TRUE,
  outputfilename = "admixture_4pops",
  imgtype = "pdf",
  height = 1.2,
  exportplot = TRUE,
  returnplot = TRUE,
  exportpath = getwd()
)
