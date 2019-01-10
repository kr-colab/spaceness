setwd("~/spaceness/")
library(plyr);library(ggplot2);library(magrittr);library(data.table);library(wesanderson)

theme_set(theme_classic()+theme(axis.text=element_text(size=7),
                                axis.title=element_text(size=7),
                                strip.text = element_text(size=7),
                                title = element_text(size=8,hjust=0.5),
                                legend.title = element_text(size=7),
                                legend.text = element_text(size=7),
                                strip.background = element_blank()))

ss_m <- fread("sumstats/ss_spatial_midpoint_W35.txt",data.table=F)
rm_m <- fread("sumstats/ss_randmates_midpoint_W35.txt",data.table=F)
ss_r <- fread("sumstats/ss_spatial_W35.txt",data.table=F)
rm_r <- fread("sumstats/ss_randmates_W35.txt",data.table=F)
ss_p <- fread("sumstats/ss_spatial_point_W35.txt",data.table=F)
rm_p <- fread("sumstats/ss_randmates_point_W35.txt",data.table=F)
ss_m$model <- "spatial"
ss_m$sampling <- "midpoint"
#rm_m$model <- "random mating"
#rm_m$sampling <- "midpoint"
ss_r$model <- "spatial"
ss_r$sampling <- "random"
rm_r$model <- "random mating"
rm_r$sampling <- ""
ss_p$model <- "spatial"
ss_p$sampling <- "point"
#rm_p$model <- "random mating"
#rm_p$sampling <- "point"
ss <- rbind(ss_m,ss_r,ss_p,rm_r)
sfs <- ss[,grepl("sfs",colnames(ss))]
ss <- ss[,!grepl("sfs",colnames(ss))]

#ss <- subset(ss,sigma != 0)
popsize <- fread("sumstats/ss_spatial_W35.txt.popsizes")
rmpopsize <- fread("sumstats/ss_randmates_W35.txt.popsizes")
popsize <- rbind(popsize,rmpopsize)
names(popsize) <- c("sigma","census_n")
ss <- merge(ss,popsize,by="sigma",all.x=T,all.y=F)

############ mixed summary stat blocks ###########
ss$census_n_2 <- ss$census_n
ss$theta_4Nmu <- 4*ss$census_n*1e-8
ss$mean_neighbors <- ((ss$sigma^2*pi)/35^2)*ss$census_n
ss$max_neighbors <- (((ss$sigma*3)^2*pi)/35^2)*ss$census_n
ss$col <- paste(ss$model,ss$sampling)
mss <- melt(ss,id.vars=c("sigma","model","sampling","census_n_2","mean_neighbors","max_neighbors","col"))
mss <- subset(mss,sigma > 0.25 & sigma < 3)

p <- ggplot(data=mss,aes(x=sigma,y=value,col=col,fill=col))+
  facet_wrap(~variable,scales="free")+
  theme(legend.position = c(0.8,0.08),legend.background = element_blank())+
  xlab(expression(sigma))+ylab("")+
  scale_color_manual(values=c("steelblue","darkgoldenrod1","darkorange","firebrick4"),name="Model")+
  scale_fill_manual(values=c("steelblue","darkgoldenrod1","darkorange","firebrick4"))+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5,size=.75)+
  geom_smooth(lwd=0.5,fill=NA,span=0.35)

ggplot(data=mss,aes(x=mean_neighbors,y=value,col=col,fill=col))+
  facet_wrap(~variable,scales="free")+
  xlab("Mean Neighborhood Size")+ylab("")+
  scale_color_manual(values=c("steelblue1","steelblue","midnightblue","darkgoldenrod1","darkorange","darkorange3"))+
  scale_fill_manual(values=c("steelblue1","steelblue","midnightblue","darkgoldenrod1","darkorange","darkorange3"))+
  geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5,size=1)+
  geom_smooth(lwd=0.75,fill=NA,span=0.35)

ggplot(data=mss,aes(x=max_neighbors/census_n_2,y=value,col=col,fill=col))+
  facet_wrap(~variable,scales="free")+
  xlab("Potential Mates / Total Population")+ylab("")+
  scale_color_manual(values=c("steelblue1","steelblue","midnightblue","darkgoldenrod1","darkorange","darkorange3"))+
  scale_fill_manual(values=c("steelblue1","steelblue","midnightblue","darkgoldenrod1","darkorange","darkorange3"))+
  geom_point(shape=21,stroke=0.25,fill=NA,alpha=0.5,size=1)+
  geom_smooth(lwd=0.75,fill=NA,span=0.35)

pdf("~/spaceness/figures/sumstats_by_sigma.pdf",width=6.5,height=4)
p+guides(color=guide_legend(ncol=1,byrow=F,keyheight = 0.15,default.unit = "in"))
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


