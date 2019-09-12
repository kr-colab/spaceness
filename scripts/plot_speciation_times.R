setwd("~/Dropbox/speciation_cyclical_migration/simple_model_output/")
library(tidyverse);library(data.table)
as0 <- fread("allopatry_s0.txt")

ggplot(data=as0,aes(x=V1))+
  theme_classic()+
  geom_histogram(fill=NA,col="black",bins=100)+
  geom_vline(xintercept=mean(as0$V1),col="red")+
  annotate(geom="text",x=mean(as0$V1)+1.2*sd(as0$V1),y=50,label=paste("mean =",round(mean(as0$V1))))

hist(as0$V1,breaks=50)
mean(as0$V1)
