setwd("~/ATLANTIDES/admixture")

library(stringr)
library(ggplot2)
library(dplyr)

#Read th table containing the Cross-validation errors to see which value of K has the lowest CV-error. 
cv <- read.table("cross_validation.txt", header=TRUE) 

#plot
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

#K=6 has the lowest value of K

#PLot of the Admixture result with pophelper package

install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
remotes::install_github('royfrancis/pophelper')
library(pophelper)

#help for the package on this link : http://www.royfrancis.com/pophelper/articles/#plotting-1
#Create a slist containing each admixture file from K=2 to K=6 
slist= readQ(files=c("input_admixture.2.Q"))
slist2= readQ(files=c("input_admixture.3.Q"))
slist3= readQ(files=c("input_admixture.4.Q"))
slist4= readQ(files=c("input_admixture.5.Q"))
slist5= readQ(files=c("input_admixture.6.Q"))

#Create an object containing the populations : 2 choices (define the population by type "wild/farmed" + geographic origin OR segregate them by river of origin)
#Here I chose to whrite only the type and the general location because there are too much different rivers for wild populations = unreadable on the plot
#the labset with the river is at the end of the script.

twolabset <- data.frame(Population= c(rep("Farmed European", 112), 
                                      rep("Wild European", 98), 
                                      rep("Farmed North American", 80), 
                                      rep("Wild American", 79)), 
                                Loc=c(rep("Norway", 112),
                                      rep("Norway", 98),
                                      rep("Canada", 80),
                                      rep("Canada", 79)),
                                      stringsAsFactors = FALSE)

#plot 
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

#Labset by river of origin
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

#Other method  : 

all_data <- tibble(sample=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

for (k in 1:6){
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
  filter(k == 6) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_brewer(palette="Set3",name="K",
                    labels=c("2","3","4","5","6"))
