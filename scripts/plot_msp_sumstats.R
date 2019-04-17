#msp grid sims summary stat viz 
setwd("~/spaceness/")
library(ggplot2);library(plyr);library(magrittr);library(data.table)

ss <- fread("sumstats/ss_msp_grid.txt",data.table = F)
ss <- ss[,!grepl("sfs",colnames(ss))]
ss$sigma2 <- ss$sigma*4
mss <- melt(ss,id.vars=c("sigma","sigma2"))
ggplot(data=mss,aes(x=sigma2*10,y=value))+
  theme_classic()+
  xlab("sigma")+
  facet_wrap(~variable,scales="free")+
  geom_point(shape=21,stroke=0.3,size=1)+geom_smooth(span=.2,fill=NA,col="orange",lwd=0.5)


