### Celian Diblasi
### 04/02/2025

### read Atlantic salmon duplicates gene data 
genetab<-read.table(file="ss4r_dups_and_singletons_ENSrapid_convPipeline.tsv",header=TRUE)
### available here: https://gitlab.com/sandve-lab/salmonid_synteny

### get numnber of singletons and number of pairs
nbsingletons<-length(genetab[genetab$type=="singleton","type"])
nbpairs<-length(genetab[genetab$type=="ss4r","type"])


nbgenesfound<-153 ### from Supplementary table (norway only) - nb of genes selected in norway populations



singletons_elem<-rep("singletons",nbsingletons)
duplicated_elem <- unlist(lapply(1:nbpairs, function(i) rep(paste0("duplicate_", i), 2)))

n_iterations <- 10000

### make a list of all gene possible to sample
elements <- c(singletons_elem, duplicated_elem)

### expectations within a population (select ohnologs in one population)

table_expectations<-c()
for( i in 1:n_iterations){
  ## pick random gene
  picked_elem<-sample(elements,nbgenesfound,replace=FALSE)
  picked_elem_nosing<-picked_elem[picked_elem!="singletons"]
  table_elem<-table(picked_elem_nosing)
  ### check if there is any duplicate
  nb_duplicate_pick<-sum(table_elem==2)
  newrow<-c(i,nb_duplicate_pick)
  table_expectations<-rbind.data.frame(table_expectations,newrow,stringsAsFactors = FALSE)
}

colnames(table_expectations)<-c("run","nbofdup")

library(ggplot2)
ggplot(table_expectations,aes(x=nbofdup))+geom_histogram(color="black",fill="dodgerblue3")+theme_bw(15)+
  ylab("Number of iterations")+xlab("Number of duplicates pairs found by chance")

### check how many rows have more than 2 duplicates (as we found in our results)
length(table_expectations[table_expectations$nbofdup>=2,"nbofdup"])



### expectation between population (one duplicate in one population and the other duplicate in another population)

nbgenesfoundbis <- 120 ### from Supplementary table (american population)


table_expectations<-c()
for( i in 1:n_iterations){
  picked_elem_sp1<-sample(elements,nbgenesfound,replace=FALSE)
  picked_elem_sp2<-sample(elements,nbgenesfoundbis,replace=FALSE)
  picked_elem_sp1_nosing<-picked_elem_sp1[picked_elem_sp1!="singletons"]
  picked_elem_sp2_nosing<-picked_elem_sp2[picked_elem_sp2!="singletons"]
  picked_elem_all<-c(picked_elem_sp1_nosing,picked_elem_sp2_nosing)
  table_elem<-table(picked_elem_all)
  nb_duplicate_pick<-sum(table_elem==2)
  newrow<-c(i,nb_duplicate_pick)
  table_expectations<-rbind.data.frame(table_expectations,newrow,stringsAsFactors = FALSE)
}

colnames(table_expectations)<-c("run","nbofdup")

library(ggplot2)
ggplot(table_expectations,aes(x=nbofdup))+geom_histogram(color="black",fill="dodgerblue3")+theme_bw(15)+
  ylab("Number of iterations")+xlab("Number of duplicates found by chance")

### check how many iterations have more than 4 (found in our results) duplicates (one gene in one population and the duplicate in the other population)
length(table_expectations[table_expectations$nbofdup>=4,"nbofdup"])




