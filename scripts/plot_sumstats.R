setwd("~/spaceness/")
library(data.table);library(ggplot2);library(magrittr);library(cowplot)
theme_set(theme_classic()+theme(axis.text=element_text(size=8),
                                axis.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                title = element_text(size=8,hjust=0.5)))

k5 <- fread("sims/sumstats.txt")
k5labs <- fread("sims/mutated/labels.txt")
k5$labs <- k5labs$V1
k5$k <- 5
k10 <- fread("sims/sumstats_k10.txt")
k10labs <- fread("sims/mutated/k10/labels.txt")
k10$labs <- k10labs$V1
k10$k <- 10
k3 <- fread("sims/sumstats_k3.txt")
k3labs <- fread("sims/mutated/k3/labels.txt")
k3$labs <- k3labs$V1
k3$k <- 3

pd <- rbind(k5,k10,k3)
pd$neighbors <- 2*pi*pd$labs^2*pd$k


names(pd) <- c("segsites","mpd","thetaW","tajD","het_o","Fis",paste0("SFS",0:100),"dispersal","k","neighbors")
sfs <- pd[,7:(ncol(pd)-2)]
pd <- melt(pd[,c("segsites","mpd","thetaW","tajD","het_o","Fis","dispersal","k","neighbors")],id.vars=c("dispersal","k","neighbors"))

#summary stat plot
ggplot(data=subset(pd,dispersal<=2),aes(x=neighbors,y=value,col=factor(k)))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        legend.box.background = element_blank())+
  scale_color_brewer(palette = "Accent")+
  xlab("Neighbors")+
  facet_wrap(~variable,scales="free")+
  geom_point(size=0.5)+
  geom_smooth(lwd=0.5,fill=NA)

ggplot(data=subset(pd,dispersal<=2),aes(x=dispersal/35*100,y=value,col=factor(k)))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        legend.box.background = element_blank())+
  scale_color_brewer(palette = "Accent")+
  xlab("Dispersal Kernal (% total range)")+
  facet_wrap(~variable,scales="free")+
  geom_point(size=0.5)+
  geom_smooth(lwd=0.5,fill=NA)

ggplot(data=subset(pd,k==5),aes(x=neighbors,y=value))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        legend.box.background = element_blank())+
  scale_color_brewer(palette = "Accent")+
  xlab("Neighbors")+
  facet_wrap(~variable,scales="free")+
  geom_point(size=0.5)+
  geom_smooth(lwd=0.5,fill=NA)



#sim vs pred plot
pred <- fread("predictions.txt")
pred$sim <- pd$dispersal[151:200]

ggplot(data=pred,aes(x=sim,y=V1,fill=V1-sim))+
  theme_classic()+
  scale_fill_distiller(palette = "RdYlBu",name="Residual")+
  geom_point(shape=21,size=2)+
  annotate(geom="segment",x=0,xend=2,y=0,yend=2)


# 

#SFS plot
sfs$dispersal_class <- factor(.bincode(sfs$dispersal,seq(0.1,2,0.3)))
sfs$sim <- 1:nrow(sfs)
sfs <- melt(sfs,id.vars=c("dispersal","dispersal_class","sim"))
sfs <- na.omit(sfs)

ddply(sfs,.(dispersal_class),summarize,nsingletons=mean())
ggplot(data=sfs,aes(x=variable,y=value,fill=dispersal_class))+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  facet_wrap(~dispersal_class)+
  geom_bar(stat="identity")


#plot of proportion uncoalesced trees for 10k outs
coal <- fread("uncoal_prop.txt",data.table=F)
trees <- fread("uncoal_tsnames.txt",header=F,data.table=F)

coal$tree <- trees$V1
coal$dispersal <- unlist(lapply(coal$tree,function(e){
  strsplit(e,"_") %>% unlist() %>% .[2] %>% as.numeric()
  }))
coal$gen <- unlist(lapply(coal$tree,function(e){
  strsplit(e,".trees") %>% unlist() %>% .[2] %>% as.numeric()
}))

ss <- fread("sumstats/10kouts_sumstats.txt",data.table = F)[,1:7]
colnames(ss) <- c("dispersal","segsites","pi","thetaW","tajD","het_o","Fis")
ss$gen <- coal$gen

coal_by_gen <- ggplot(data=coal,aes(x=gen,y=V1,col=dispersal))+
  ggtitle("Coalescence by Generation\nNWF Spatial Model, n~=1000")+
  ylab("Proportion Uncoalesced Gene Trees")+xlab("Generation")+
  theme(legend.position = "none")+
  scale_color_distiller(palette = "RdYlBu")+
  geom_point(position=position_jitter(width=2000),size=0.75)

tajd_by_gen <- ggplot(data=ss,aes(x=gen,y=tajD,col=dispersal))+
  ggtitle("Tajima's D by Generation\nNWF Spatial Model, n~=1000")+
  scale_color_distiller(palette = "RdYlBu")+
  geom_point(position=position_jitter(width=2000),size=0.75)

ss_final <- subset(ss,gen==200000)
mss <- melt(ss_final,id.vars=c("dispersal","gen"))
ggplot(data=mss,aes(x=dispersal,y=value))+
  theme(strip.background = element_blank())+
  facet_wrap(~variable,scales="free")+
  geom_point(shape=1)+
  geom_smooth(fill=NA,wd=0.75,col="forestgreen")

pdf("~/spaceness/figures/10kouts_coalescence_by_generation.pdf",width=7.5,height=3.5)
ggdraw()+
  draw_plot(coal_by_gen,0,0,0.5,1)+
  draw_plot(tajd_by_gen,0.5,0,.5,1)
dev.off()




