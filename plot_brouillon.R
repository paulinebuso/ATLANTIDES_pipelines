setwd("~/ATLANTIDES/admixture")

library(stringr)
library(ggplot2)
library(dplyr)

cv <- read.table("cross_validation.txt", header=TRUE) 

ggplot(cv,aes(x=K,y=cv_error)) + geom_line()+scale_x_continuous(breaks=seq(1,16,1))+
  labs(title="Cross-validation error - ADMIXTURE", xlab="K value", ylab="CV error") +
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))+
  theme_classic()

ggsave("Admixture_cross-validation.pdf",width=7,height=5,dpi=600)
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

twolabset <- data.frame(Population= c(rep("Farmed European", 112), 
                                      rep("Wild European", 98), 
                                      rep("Farmed North American", 80), 
                                      rep("Wild American", 79)), 
                                Loc=c(rep("Norway", 112),
                                      rep("Norway", 98),
                                      rep("Canada", 80),
                                      rep("Canada", 79)),
                                      stringsAsFactors = FALSE)

p <- plotQ(
  qlist = c(slist,slist2,slist3,slist4,slist5), 
  imgoutput = "join",
  clustercol = c(c("#fc2634", "#092fd9", "#fc79a7","#8aceeb","#3ec964","#fac132")),
  sortind = NA,
  grplab = twolabset,
  selgrp = "Loc",
  ordergrp = TRUE,
  panelspacer = 0.1,
  showsp = FALSE,
  splab = c("K=2","K=3","K=4","K=5","K=6"),
  splabsize = 3,
  grplabspacer = -0.6,
  grplabheight = 0,
  grplabpos = 0.8,
  grplabsize = 0.8,
  grplabangle = 30,
  showindlab = FALSE,
  useindlab = TRUE,
  indlabheight = 0.7,
  indlabsize = 1.2,
  indlabangle = 90,
  indlabvjust = 1,
  indlabhjust = 1,
  indlabcol = "grey30",
  pointsize = 1.5,
  linepos = 1,
  linesize = 0.2,
  divgrp = c("Population","Loc"),
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

#Loc=c(rep("Gaspe New Brunswick", 1), 
# rep("St John River", 13),
# rep("Gaspe New Brunswick", 6),
# rep("Penobscot River", 7),
#rep("St John River", 2),
#rep("Penobscot River", 2),
# rep("St John River", 1),
# rep("Gaspe New Brunswick", 2),
# rep("St John River", 29),
# rep("Penobscot River", 2),
# rep("Gaspe New Brunswick", 2),
# rep("St John River", 1),
# rep("Gaspe New Brunswick", 5),
# rep("St John River", 7),
# rep("Bonaventure", 10),
# rep("De La Chaloupe", 10),
# rep("Jupiter", 10),
# rep("Laval", 10),
# rep("Malbaie", 10),
# rep("Cascapedia", 10),
# rep("Saint Paul", 10),
# rep("Du vieux Fort", 10),
# rep("MOWI - Norwegian West cost", 112), 
# rep("Aaroey", 10), 
# rep("Daleelva", 2),
# rep("Eidfjordvassdraget", 2),
# rep("Etneelva", 2),
# rep("Flaamselva", 10),
# rep("Flekeelva", 2),
# rep("Gloppenelva", 2),
# rep("Joelstra", 10),
# rep("Laerdalselva", 2),
# rep("Loneelva Osteroey", 2),
#rep("Nausta", 10),
#rep("Oselva", 2),
#rep("Oselvvassdraget", 10),
# rep("Ryggelva", 10),
# rep("Suldalslaagen", 10),
# rep("Vikedalselva", 10),
# rep("Vorma", 2)) ,
#stringsAsFactors = FALSE)

#autre methode

library(reshape2)
library(plyr)
library(tidyverse)
library(RColorBrewer)

admixture <- read.table("input_admixture.7.Q") %>% as_tibble %>% rename_all(~c("K1","K2","K3","K4","K5","K6","K7"))
id <- read.table("samples_list.txt") 
pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("V1","Location","Populations","Origin"))

#merge admixture ouput + populations and samples in the same order as the admixture output 
admixture <- cbind(id,admixture)
admixture <- merge(admixture, pops_file, by="V1")

#remove column I don't need : loc and origin
admixture <- admixture[,-9]
admixture <- admixture[,-10]

#rename the column
colnames(admixture) <- c("IND","K1","K2","K3","K4","K5","K6", "K7", "POP")

#Gather the colum "K1" and "K2" into the column "ANCESTRY"
admixture_long <- melt(admixture,id.vars=c("IND","POP"),variable.name="ANCESTRY",value.name="PERC")
names(admixture_long)
class(admixture_long$ANCESTRY)
levels(admixture_long$ANCESTRY)

# subset only 50%
#admixture_long_50 <- subset(admixture_long, subset=admixture_long$PERC>=0.50)

#graph admixture results

col <- c("#3250a8","#eb1a90", "#55e066", "#f09d18", "#9620ba", "#f5df1b", "#58dbab")

ggplot(admixture_long,aes(x=POP,y=PERC,fill= ANCESTRY)) + 
  geom_bar(stat="identity", position="stack") +
  xlab("Populations") + ylab("Ancestry") +
  scale_fill_manual(values=col, name= "K")

#autre methode

all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:7){
  data <- read_delim(paste0("input_admixture.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- id$V1
  data$k <- k
  
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}
all_data


pops_file <- read.table("~/ATLANTIDES/PCA/salmon_pops.csv", header = TRUE, sep=",")  %>% as_tibble %>% rename_all(~c("sample","Location","Populations","Origin"))
all_data_plot <- merge(all_data, pops_file, by="sample")

#plot
all_data_plot %>%
  filter(k == 4) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set3",name="K",
                    labels=c("1","2","3","4","5","6","7"))
