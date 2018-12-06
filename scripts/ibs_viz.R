setwd("~/spaceness/")
library(plyr);library(ggplot2);library(magrittr);library(data.table);library(wesanderson)

theme_set(theme_classic()+theme(axis.text=element_text(size=7),
                                axis.title=element_text(size=7),
                                strip.text = element_text(size=7),
                                title = element_text(size=8,hjust=0.5),
                                legend.title = element_text(size=7),
                                legend.text = element_text(size=7),
                                strip.background = element_blank()))

ss <- fread("sumstats/ss_spatial_W35.txt",data.table=F)
rm <- fread("sumstats/ss_randmates_W35.txt",data.table=F)
ss$model <- "spatial"
rm$model <- "random mating"
ss <- rbind(ss,rm)
sfs <- ss[,15:(ncol(ss)-1)]
#ibsbins <- ss[,113:(ncol(ss)-1)]
ss <- ss[,c(1:14,ncol(ss))]
colnames(ss) <- c("sigma","segsites","pi","thetaW","tajD","het_o","fis",
                  "gen_sp_corr","ibs_mean","ibs_var","ibs_skew","ibs_95p","ibs_num_blocks_per_pair","ibs_mean_spat_corr","model")
#ss <- subset(ss,sigma != 0)
popsize <- fread("sumstats/popsize_spatial_W35.txt")
rmpopsize <- fread("sumstats/popsize_randmates_W35.txt")
popsize <- rbind(popsize,rmpopsize)
names(popsize) <- c("sigma","census_n")
ss <- merge(ss,popsize,by="sigma")

ibsbins$sigma <- ss$sigma
tmp <- melt(ibsbins,id.vars=c("sigma"))
tmp$variable <- unlist(lapply(tmp$variable,function(e) as.numeric(gsub("V","",e))-111))
tmp <- subset(tmp,sigma != 0)

######### IBS block size plots #############
ggplot(data=subset(tmp,sigma==min(tmp$sigma)| sigma == max(tmp$sigma)),
       aes(x=variable,y=log(value)/1e8,col=factor(sigma)))+
  xlim(0,1000)+
  geom_line()

ggplot(data=tmp,aes(x=variable,y=log(value),col=sigma,z=sigma))+
  scale_fill_viridis_c()+
  xlim(0,10000)+
  xlab("IBS Block Size")+ylab("log(count)")+
  stat_summary_hex(fun="mean",bins = 100)

############ mixed summary stat blocks ###########
ss$census_n_2 <- ss$census_n
ss$theta_4Nmu <- 4*ss$census_n*4e-8
mss <- melt(ss,id.vars=c("sigma","model","census_n_2"))
mss <- subset(mss,sigma > 0.5)

ggplot(data=mss,aes(x=sigma,y=value,col=model,fill=model))+
  facet_wrap(~variable,scales="free")+
  xlab("sigma")+ylab("")+
  scale_color_manual(values = c("cornflowerblue","chartreuse4"))+
  scale_fill_manual(values = c("cornflowerblue","chartreuse4"))+
  #geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5)+
  geom_smooth(lwd=0.75,fill=NA,span=0.5)

ggplot(data=mss,aes(x=((sigma^2*pi)/35^2)*census_n_2,y=value,col=model,fill=model))+
  facet_wrap(~variable,scales="free")+
  xlab("Mean Neighborhood Size")+ylab("")+
  scale_color_manual(values = c("cornflowerblue","chartreuse4"))+
  scale_fill_manual(values = c("cornflowerblue","chartreuse4"))+
  #geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5)+
  geom_smooth(lwd=0.75,fill=NA,span=0.5)

ggplot(data=mss,aes(x=((((sigma*3)^2*pi)/35^2)*census_n_2)/census_n_2,y=value,col=model,fill=model))+
  facet_wrap(~variable,scales="free")+
  xlab("Potential Mates / Total Population")+ylab("")+
  scale_color_manual(values = c("cornflowerblue","chartreuse4"))+
  scale_fill_manual(values = c("cornflowerblue","chartreuse4"))+
  geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5)+
  geom_smooth(lwd=0.75,fill=NA,span=0.5)

ggplot(data=mss,aes(x=sigma/35*100,y=value,col=model,fill=model))+
  facet_wrap(~variable,scales="free")+
  xlab("Interaction Radius (% Total Range)")+ylab("")+
  scale_color_manual(values = c("cornflowerblue","chartreuse4"))+
  scale_fill_manual(values = c("cornflowerblue","chartreuse4"))+
  geom_point(shape=21,stroke=0.5)+
  geom_smooth(lwd=0.75)
