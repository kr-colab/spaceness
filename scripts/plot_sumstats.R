setwd("~/spaceness/")
library(data.table);library(ggplot2);library(magrittr);library(cowplot);library(plyr);library(wesanderson)
theme_set(theme_classic()+theme(axis.text=element_text(size=7),
                                axis.title=element_text(size=7),
                                strip.text = element_text(size=7),
                                title = element_text(size=8,hjust=0.5),
                                legend.title = element_text(size=7),
                                legend.text = element_text(size=7),
                                strip.background = element_blank()))
discrete_palette <- c("steelblue1","forestgreen","chartreuse3","slateblue3")


spat <- fread("sumstats/ss_spatial.txt")
spat$type <- "spatial\nmate choice\n"
spat$popsize <- fread("sumstats/popsize_spatial.txt",header = F)$V1
rm <- fread("sumstats/ss_random_mating.txt")
rm$type <- "random\nmating\n"
rm$popsize <- fread("sumstats/popsize_random_mating.txt",header = F)$V1

pd <- rbind(spat,rm)
names(pd) <- c("dispersal","segsites","pi","thetaW","tajD","het_o","Fis",paste0(0:100),"type","Census N")
sfs <- pd[,8:(ncol(pd)-1)]
sfs$dispersal <- pd$dispersal
sfs$segsites <- pd$segsites
sfs$type <- pd$type
pd$mean_neighbors <- ((pi*pd$dispersal^2)/16^2)*1000
pd$max_neighbors <- (3*((pi*pd$dispersal^2)/16^2))*1000
pd <- melt(pd[,c("segsites","pi","tajD","het_o","Fis","dispersal","type","Census N","mean_neighbors","max_neighbors")],
           id.vars=c("dispersal","type","mean_neighbors","max_neighbors"))

#summary stat plot
pd$variable <- factor(pd$variable,levels=c("Census N","segsites","pi","tajD","het_o","Fis"))
ss_by_dispersal <- ggplot(data=pd,aes(x=dispersal/16*100,y=value,fill=type,col=type))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        strip.text = element_text(size=8))+
  #scale_color_brewer(palette = "Dark2",name="Model")+
  scale_fill_manual(values=discrete_palette,name="Model")+
  scale_color_manual(values=discrete_palette,guide=F)+
  xlab("Interaction Radius (% range width)")+
  facet_wrap(~variable,scales="free")+
  geom_point(size=1.4,shape=21,stroke=0.02,col="black",alpha=0.7)+
  geom_smooth(lwd=0.5,fill=NA)

pdf("~/spaceness/figures/sumstat_by_dispersal_spat_v_rm.pdf",width=7.5,height=3.5,useDingbats = F)
ss_by_dispersal <- ss_by_dispersal + guides(fill = guide_legend(override.aes = list(size=4)))
print(ss_by_dispersal)
dev.off()

ss_by_neighbors <- ggplot(data=pd,aes(x=mean_neighbors,y=value,fill=type,col=type))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        strip.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=9),
        strip.text = element_text(size=8))+
  #scale_color_brewer(palette = "Dark2",name="Model")+
  scale_fill_manual(values=discrete_palette,name="Model")+
  scale_color_manual(values=discrete_palette,guide=F)+
  xlab("Mean Neighborhood Size")+
  facet_wrap(~variable,scales="free")+
  geom_point(size=1.4,shape=21,stroke=0.02,col="black",alpha=0.7)+
  geom_smooth(lwd=0.5,fill=NA)

pdf("~/spaceness/figures/sumstat_by_neighbors_spat_v_rm.pdf",width=7.5,height=3.5,useDingbats = F)
ss_by_neighbors <- ss_by_neighbors + guides(fill = guide_legend(override.aes = list(size=4)))
print(ss_by_neighbors)
dev.off()


#SFS plots for binned runs
spat_sfs <- fread("sumstats/ss_spatial_bins.txt",data.table=F)
rm_sfs <- fread("sumstats/ss_randmates_bins.txt",data.table=F)
spat_sfs$model <- "Spatial Mate Choice"
rm_sfs$model <- "Random Mating"
sfs <- rbind(spat_sfs,rm_sfs)

sfs <- sfs[,c(1,2,8:(ncol(sfs)))]
colnames(sfs) <- c("sigma","segsites",paste(0:100),"model")
sfs$dispersal_class <- round(sfs$sigma,2)
sfs$sim <- 1:nrow(sfs)
sfs <- subset(sfs,sigma<=1)
sfs <- melt(sfs,id.vars=c("sigma","dispersal_class","sim","segsites","model"))
sfs$sfs_class <- as.numeric(as.character(sfs$variable))
sfs <- na.omit(sfs)
sfs <- subset(sfs,sfs_class>0 & sfs_class<100)

sfsbins <- ddply(sfs,.(variable,dispersal_class,model),summarize,meanCount=mean(value),meanSegsites=mean(segsites))
sfsbins$sfs_class <- as.numeric(as.character(sfsbins$variable))
sfsbins$log_segscaled_n <- log(sfsbins$meanCount/sfsbins$meanSegsites)

sfsbins <- subset(sfsbins,dispersal_class != 0.5) #can remove once rm 0.5 sims are done

sfs_diffs <- ddply(sfsbins,.(sfs_class,dispersal_class),function(e){
  diff <- e$log_segscaled_n[e$model=="Spatial Mate Choice"]-e$log_segscaled_n[e$model=="Random Mating"]
  return(diff)
})

sfsplot <- ggplot(data=sfsbins,aes(x=sfs_class,y=log_segscaled_n,col=factor(dispersal_class)))+
  facet_grid(.~model)+
  theme(legend.position = "none",
        )+
  scale_color_manual(values=wes_palette("Zissou1",nlevels(factor(sfsbins$dispersal_class)),"continuous"),
                     name="Interaction\nRadius")+
  ylab("log(count per segregating site)")+xlab("SFS Class")+
  geom_line(lwd=0.25)

sfsplot_diffs <- ggplot(data=sfs_diffs,aes(x=sfs_class,col=factor(dispersal_class),fill=factor(dispersal_class),y=V1))+
  scale_color_manual(values=wes_palette("Zissou1",nlevels(factor(sfs_diffs$dispersal_class)),"continuous"),
                     name="Interaction\nRadius",guide=F)+
  scale_fill_manual(values=wes_palette("Zissou1",nlevels(factor(sfs_diffs$dispersal_class)),"continuous"),
                     name="Interaction\nRadius")+
  xlab("SFS Class")+ylab(expression(italic(SFS[spatial])-italic(SFS[random~mating])))+
  geom_point(shape=21,col="black",stroke=0.05,size=1.05,alpha=0.7)+
  geom_smooth(fill=NA,lwd=0.5)

sfsplot_diffs <- sfsplot_diffs+guides(fill=guide_legend(override.aes = list(size=4),
                                                        keyheight = 3,
                                                        default.unit = "mm"))

pdf("~/spaceness/figures/sfs_spatial_v_rm.pdf",width=6.5,height=2)
ggdraw()+
  draw_plot(sfsplot,0,0,0.5,1)+
  draw_plot(sfsplot_diffs,0.5,0,0.5,1)
dev.off()




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

coal_by_gen <- ggplot(data=coal,aes(x=gen,y=V1,fill=dispersal))+
  #ggtitle("Coalescence by Time Step\nNWF Spatial Model, n~=1000")+
  ylab("Proportion Uncoalesced Gene Trees")+xlab("")+
  theme(legend.position = c(0.75,0.65))+
  scale_fill_distiller(palette = "RdYlBu",name="Interaction\nDistance")+
  geom_point(position=position_jitter(width=2000),shape=21,stroke=0.01,col="grey40")
coal_by_gen <- coal_by_gen + guides(fill=guide_colorbar(barwidth=unit(4,"mm"),barheight=unit(22,"mm"),ticks.colour = "black"))

tajd_by_gen <- ggplot(data=ss,aes(x=gen,y=tajD,fill=dispersal))+
  theme(legend.position = "none")+
  #ggtitle("Tajima's D by Time Step\nNWF Spatial Model, n~=1000")+
  ylab("Tajima's D")+xlab("")+
  scale_fill_distiller(palette = "RdYlBu")+
  geom_point(position=position_jitter(width=2000),shape=21,stroke=0.01,col="grey40")

pdf("~/spaceness/figures/10kouts_coalescence_by_generation.pdf",width=6.5,height=2.5)
ggdraw()+
  draw_plot(tajd_by_gen,0,0,0.5,1)+
  draw_plot(coal_by_gen,0.5,0,.5,1)+
  draw_text("Time Steps",0.52,0.035,9)
dev.off()

ss_final <- subset(ss,gen==200000)
mss <- melt(ss_final,id.vars=c("dispersal","gen"))
ss_by_dispersal <- ggplot(data=mss,aes(x=dispersal,y=value))+
  theme(strip.background = element_blank())+ylab("")+
  facet_wrap(~variable,scales="free")+
  xlab("Time Steps")+
  geom_point(shape=1)+
  geom_smooth(fill=NA,lwd=0.75,col="forestgreen")

pdf("~/spaceness/figures/10kouts_final_gen_ss_by_dispersal.pdf",width=7.5,height=3.5)
print(ss_by_dispersal)
dev.off()


