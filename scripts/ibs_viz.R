setwd("~/spaceness/")
library(plyr);library(ggplot2);library(magrittr);library(data.table);library(car);library(wesanderson)
library(cowplot)
source("~/R/colordict.R")
pal <- grDevices::colorRampPalette(color=c("steelblue4","skyblue","lightgoldenrod1","orangered","red4"),
                                   bias=1,space="rgb",interpolate="linear")
theme_set(theme_classic()+theme(axis.text=element_text(size=8),
                                axis.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                title = element_text(size=9,hjust=0.5),
                                legend.title = element_text(size=8),
                                legend.text = element_text(size=8),
                                strip.background = element_blank()))

ss_m <- fread("sumstats/ss_spatial_midpoint_W50.txt",data.table=F)
ss_r <- fread("sumstats/ss_spatial_random_W50.txt",data.table=F)
ss_p <- fread("sumstats/ss_spatial_point_W50.txt",data.table=F)
rm_r <- fread("sumstats/ss_randmates_random_W50.txt",data.table=F)
rm_m <- fread("sumstats/ss_randmates_midpoint_W50.txt",data.table=F)
rm_p <- fread("sumstats/ss_randmates_point_W50.txt",data.table=F)
ss_m$model <- "spatial"
ss_m$sampling <- "midpoint"
ss_r$model <- "spatial"
ss_r$sampling <- "random"
rm_r$model <- "random mating"
rm_r$sampling <- "random"
rm_m$model <- "random mating"
rm_m$sampling <- "midpoint"
rm_p$model <- "random mating"
rm_p$sampling <- "point"
ss_p$model <- "spatial"
ss_p$sampling <- "point"
ss <- rbind(ss_m,ss_r,ss_p,rm_r,rm_m,rm_p)
popsize <- fread("sumstats/ss_spatial_random_W50.txt.popsizes")
rmpopsize <- fread("sumstats/ss_randmates_random_W50.txt.popsizes")
popsize <- rbind(popsize,rmpopsize)
names(popsize) <- c("sigma","census_n")
ss <- merge(ss,popsize,by="sigma",all.x=T,all.y=F)
sfs <- ss[,grepl("sfs",colnames(ss))] #split off sfs for separate analysis
ss <- ss[,!grepl("sfs",colnames(ss))]
ss$neighborhood_size <- 4*pi*ss$sigma^2*5
ss$col <- paste(ss$model,ss$sampling)

#test for differences by sampling in spatial and random mating summary stat distributions
mss <- melt(ss,id.vars=c("sigma","model","sampling","neighborhood_size","col"))
mss <- subset(mss,variable!="census_n")
sampling_diffs <- ddply(mss,.(variable,model),function(e){
  meantest <- summary(aov(value~sampling,data=e))[[1]][,5][1]
  vartest <- leveneTest(value~sampling,data=e)[[3]][1] 
  return(c(round(meantest,5),round(vartest,5)))
}) #40 warnings for forcing sampling to factor, which is fine
sampling_diffs <- arrange(sampling_diffs,model,variable)
colnames(sampling_diffs) <- c("variable","model","p(different\nmeans)","p(different\nvariance)")
xtable(sampling_diffs,digits=6)

#drop alternate sampling schemes for random mating simulations
ss <- rbind(ss_m,ss_r,ss_p,rm_r)
popsize <- fread("sumstats/ss_spatial_random_W50.txt.popsizes")
rmpopsize <- fread("sumstats/ss_randmates_random_W50.txt.popsizes")
popsize <- rbind(popsize,rmpopsize)
names(popsize) <- c("sigma","census_n")
ss <- merge(ss,popsize,by="sigma",all.x=T,all.y=F)
sfs <- ss[,grepl("sfs",colnames(ss))] 
ss <- ss[,!grepl("sfs",colnames(ss))]
ss$neighborhood_size <- 4*pi*ss$sigma^2*5
ss$col <- paste0(ss$model,"\n",ss$sampling," sampling")

############ plot all summary stats ###########
mss <- melt(ss,id.vars=c("sigma","model","sampling","neighborhood_size","col"))
mss <- subset(mss,variable!="gen_dist_mean")
mss$variable <- factor(mss$variable,levels=c("segsites","pi","thetaW","tajD",
                                             "het_o","fis","gen_dist_var","gen_dist_skew",
                                             "ibs_mean","ibs_var","ibs_skew",
                                             "ibs_blocks_per_pair","ibs_blocks_over_1e6_per_pair",
                                             "gen_sp_corr","ibs_mean_spat_corr","ibs_blocks_spat_corr",
                                             "ibs_1e6blocks_spat_corr","ibs_skew_spat_corr",
                                             "census_n"))
mss$variable <- mapvalues(mss$variable,from=c("segsites","pi","thetaW","tajD",
                                              "het_o","fis","gen_dist_var","gen_dist_skew",
                                              "ibs_mean","ibs_var","ibs_skew",
                                              "ibs_blocks_per_pair","ibs_blocks_over_1e6_per_pair",
                                              "gen_sp_corr","ibs_mean_spat_corr","ibs_blocks_spat_corr",
                                              "ibs_1e6blocks_spat_corr","ibs_skew_spat_corr",
                                              "census_n"),
                          to=c("segregating\nsites","pi","Watterson's\ntheta","Tajima's\nD",
                             "observed\nheterozygosity","Fis","var(dxy)","skew(dxy)",
                             "mean(IBS)","var(IBS)","skew(IBS)","nIBS","mean(IBS>1e6)",
                             "corr(dxy,dist)","corr(mean(IBS),dist)","corr(nIBS,dist)",
                             "corr(IBS>1e6,dist)","corr(skew(IBS),dist)",
                             "census N"))

mss <- subset(mss,sigma<=4)
p <- ggplot(data=mss,aes(x=neighborhood_size,y=value,col=col,fill=col))+
  facet_wrap(~variable,scales="free",ncol=4)+
  theme(legend.position = "bottom",legend.direction = "horizontal")+
  xlab("Neighborhood Size")+ylab("")+
  scale_color_manual(values=pal(4),name="Model")+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  scale_x_log10()+
  geom_point(shape=21,stroke=0.4,fill=NA,size=.9)+
  geom_smooth(lwd=0.3,fill=NA,span=0.2)
p <- p+guides(color=guide_legend(keyheight = 0.15,nrow = 1,title.position = "top",
                                 default.unit = "in",override.aes = list(size=5,shape=16)))


pdf("~/spaceness/figures/sumstats_by_neighbors_allstats.pdf",width=6,height=6)
print(p)
dev.off()

#subset of sumstats for main figure 
mss <- subset(mss,variable %in% c("pi","Watterson's\ntheta","observed\nheterozygosity","Fis",
                                  "Tajima's\nD","var(dxy)","nIBS","mean(IBS>1e6)",
                                  "corr(dxy,dist)","corr(mean(IBS),dist)",
                                  "corr(skew(IBS),dist)","corr(IBS>1e6,dist)"))
mss$variable <- factor(mss$variable,levels=c("pi","Watterson's\ntheta","observed\nheterozygosity","Fis",
                                             "Tajima's\nD","var(dxy)","nIBS","mean(IBS>1e6)",
                                             "corr(dxy,dist)","corr(mean(IBS),dist)",
                                             "corr(skew(IBS),dist)","corr(IBS>1e6,dist)"))
p <- ggplot(data=mss,aes(x=neighborhood_size,y=value,col=col,fill=col))+
  facet_wrap(~variable,scales="free",ncol=4)+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,-30),)+
  xlab("Neighborhood Size")+ylab("")+
  scale_color_manual(values=getcolordict(list('dark_medici_blue',
                                              'burnt_sienna',
                                              'isabella_color',
                                              'red_orange')),name="Model")+
  #scale_color_manual(values=c("black","steelblue4","steelblue2","forestgreen"),name="Model")+
  #scale_color_brewer(palette = "Paired",name="Model")+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  scale_x_log10()+
  geom_point(shape=21,stroke=0.4,fill=NA,size=.9)+
  geom_smooth(lwd=0.3,fill=NA,span=0.2)
p <- p+guides(color=guide_legend(keyheight = 0.15,nrow = 1,title.position = "top",
                                 default.unit = "in",override.aes = list(size=5,shape=16)))

#SFS plots
sfs <- cbind(sfs,ss[,c("sigma","sampling","model","segsites")])
sfs$neighbors <- 4*pi*sfs$sigma^2*5
colnames(sfs) <- sapply(colnames(sfs),function(e) gsub("sfs_","",e))
msfs <- melt(sfs,id.vars=c("sigma","neighbors","sampling","model","segsites"))
msfs$value <- msfs$value/msfs$segsites
msfs$variable <- as.numeric(as.character(msfs$variable))
msfs$variable <- msfs$variable/120
msfs <- subset(msfs,sigma<4)
msfs <- msfs[msfs$neighbors %in% unique(msfs$neighbors)[sample(1:length(unique(msfs$neighbors)),100)],]
p2 <- ggplot(data=msfs,aes(x=variable,y=value,z=neighbors))+
  theme(axis.text.x=element_text(size=7,angle=45,hjust=1),plot.background = element_blank())+
  facet_wrap(paste0(model,"\n",sampling," sampling")~.,nrow=1)+
  #facet_grid(sampling~model)+
  scale_y_log10()+
  #scale_color_distiller(palette = "RdYlBu",name="Neighborhood\nSize")+
  scale_color_gradientn(colors=pal(100),name="Neighborhood\nSize")+
  xlab("Allele Frequency")+ylab("Proportion of\nSegregating Sites")+
  geom_smooth(fill=NA,lwd=0.25,alpha=0.5,span=0.3,aes(group=neighbors,col=neighbors))
  #geom_line(aes(col=neighbors,group=neighbors))
  #stat_summary_hex(fun=median,bins=100)
p2 <- p2+guides(color=guide_colorbar(barwidth=unit(4,"mm"),barheight=unit(22,"mm")))


pdf("~/spaceness/figures/sfs_w_sumstats.pdf",width=6.5,height=6,useDingbats = F)
ggdraw()+
  draw_plot(p,0,0,1,0.7)+
  draw_plot(p2,0,0.7,1,0.3)+
  draw_label("A",0.07,0.96)+
  draw_label("B",0.07,0.72)
dev.off()


#example sampling maps
rs <- fread("~/Desktop/spaceviz/random_sampling_locs.txt")[sample(1:nrow(rs),60),]
ps <- fread("~/Desktop/spaceviz/point_sampling_locs.txt")[sample(1:nrow(ps),60),]
ms <- fread("~/Desktop/spaceviz/midpoint_sampling_locs.txt")[sample(1:nrow(ms),60),]
rs$sampling <- "random"
ps$sampling <- "point"
ms$sampling <- "midpoint"
sampling <- rbind(rs,ps,ms)

samplemaps <- ggplot(data=sampling,aes(x=V1,y=V2))+coord_fixed()+
  theme(axis.title=element_blank())+
  facet_wrap(~sampling)+
  geom_point(size=0.5)
pdf("~/spaceness/figures/sampling_maps.pdf",width=5.5,height=1.5,useDingbats = F)
print(samplemaps)
dev.off()

