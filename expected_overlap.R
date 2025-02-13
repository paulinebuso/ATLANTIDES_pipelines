### Celian Diblasi
### 13/02/2025

### importot expected overlap data from expectations_overlap.slurm
expectedoverlap<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/pauline_data/Overlap_expected.txt",header=TRUE)
expectedoverlap$type<-"expected"

### then create a table with real observations
EU_data<-data.frame(Run=0,overlap_FST_TAJ=8,overlap_XPEHH_TAJ=0,overlap_FST_XPEHH=104,type="observed")
AM_data<-data.frame(Run=0,overlap_FST_TAJ=2,overlap_XPEHH_TAJ=0,overlap_FST_XPEHH=339,type="observed")


library(ggplot2)
expectedoverlap<-expectedoverlap[1002:2001,]
expectedoverlap$overlap_FST_TAJ<-as.numeric(as.character(expectedoverlap$overlap_FST_TAJ))
expectedoverlap$overlap_FST_XPEHH<-as.numeric(as.character(expectedoverlap$overlap_FST_XPEHH))
expectedoverlap$overlap_XPEHH_TAJ<-as.numeric(as.character(expectedoverlap$overlap_XPEHH_TAJ))



ggplot(expectedoverlap, aes(x = overlap_FST_TAJ))+
  geom_histogram(binwidth = 1, alpha = 0.7,position="stack",fill="lightblue4") +
  theme_bw(15) +geom_vline(data = AM_data, mapping = aes(xintercept = overlap_FST_TAJ),color="#C7190B",size=1)+
  geom_vline(data = EU_data, mapping = aes(xintercept = overlap_FST_TAJ),color="#050FA3",size=1)


ggplot(expectedoverlap, aes(x = overlap_XPEHH_TAJ))+
  geom_histogram(binwidth = 1, alpha = 0.7,position="stack",fill="lightblue4") +
  theme_bw(15) +geom_vline(data = AM_data, mapping = aes(xintercept = overlap_XPEHH_TAJ),color="#C7190B",size=1)+
  geom_vline(data = EU_data, mapping = aes(xintercept = overlap_XPEHH_TAJ),color="#050FA3",size=1)

ggplot(expectedoverlap, aes(x = overlap_FST_XPEHH))+
  geom_histogram(binwidth = 1, alpha = 0.7,position="stack",fill="lightblue4") +
  theme_bw(15) +geom_vline(data = AM_data, mapping = aes(xintercept = overlap_FST_XPEHH),color="#C7190B",size=1)+
  geom_vline(data = EU_data, mapping = aes(xintercept = overlap_FST_XPEHH),color="#050FA3",size=1)
