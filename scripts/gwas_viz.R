#visualizing GWAS on SLiM simulations
library(ggplot2);library(magrittr);library(plyr);library(data.table)
library(qvalue);library(cowplot);library(wesanderson)
setwd("~/projects/spaceness/gwas/")
setwd("~/spaceness/gwas")
pal <- grDevices::colorRampPalette(color=c("steelblue4","skyblue","lightgoldenrod1","orangered","red4"),
                                   bias=1,space="rgb",interpolate="linear")
theme_set(theme_classic()+theme(
          axis.text=element_text(size=6),
          axis.title=element_text(size=6),
          strip.text=element_text(size=6),
          legend.text=element_text(size=6),
          legend.title = element_text(size=6),
          strip.background = element_blank(),
          title=element_text(size=6)))

#####################################################
###################  qq plots  ######################
#####################################################
#helper function for qqplot confidence interval given binomial sampling for alpha=0.05
panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) { 
  require(grid)
  conf.points = min(conf.points, n-1);
  mpts<-matrix(nrow=conf.points*2, ncol=2)
  for(i in seq(from=1, to=conf.points)) {
    mpts[i,1]<- -log10((i-.5)/n)
    mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
    mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
    mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
  }
  #grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  return(mpts)
}
# setwd("random_sampling")
# files <- list.files("out_corner",full.names = T) %>% grep(".assoc.linear",.,value = T)
# qq_corner <- data.table(P=NA,sigma=NA,observed=NA,expected=NA,neighbors=NA)[-1]
# for(f in files){
#   sigma <- strsplit(f,"_") %>% unlist() %>% .[3] %>% as.numeric()
#   if(sigma < 4){
#     gwas <- fread(f)
#     gwas <- subset(gwas,TEST=="ADD")
#     expected <- -log10(1:length(gwas$P)/length(gwas$P))
#     observed <- -log10(sort(gwas$P))
#     qq_corner <- rbind(qq_corner,data.table(P=gwas$P,sigma=sigma,observed=observed,expected=expected,neighbors=4*pi*sigma^2*5))
#     print(f)
#   }
# }
# qq_corner$phenotype <- "corner"
# 
# files <- list.files("out_clinal",full.names = T) %>% grep(".assoc.linear",.,value = T)
# qq_clinal <- data.table(P=NA,sigma=NA,observed=NA,expected=NA,neighbors=NA)[-1]
# for(f in files){
#   sigma <- strsplit(f,"_") %>% unlist() %>% .[3] %>% as.numeric()
#   if(sigma < 4){
#     gwas <- fread(f)
#     gwas <- subset(gwas,TEST=="ADD")
#     expected <- -log10(1:length(gwas$P)/length(gwas$P))
#     observed <- -log10(sort(gwas$P))
#     qq_clinal <- rbind(qq_clinal,data.table(P=gwas$P,sigma=sigma,observed=observed,expected=expected,neighbors=4*pi*sigma^2*5))
#     print(f)
#   }
# }
# qq_clinal$phenotype <- "clinal"
# 
# files <- list.files("out_patchy",full.names = T) %>% grep(".assoc.linear",.,value = T)
# qq_patchy <- data.table(P=NA,sigma=NA,observed=NA,expected=NA,neighbors=NA)[-1]
# for(f in files){
#   sigma <- strsplit(f,"_") %>% unlist() %>% .[3] %>% as.numeric()
#   if(sigma < 4){
#     gwas <- fread(f)
#     gwas <- subset(gwas,TEST=="ADD")
#     expected <- -log10(1:length(gwas$P)/length(gwas$P))
#     observed <- -log10(sort(gwas$P))
#     qq_patchy <- rbind(qq_patchy,data.table(P=gwas$P,sigma=sigma,observed=observed,expected=expected,neighbors=4*pi*sigma^2*5))
#     print(f)
#   }
# }
# qq_patchy$phenotype <- "patchy"
# 
# files <- list.files("out_normal",full.names = T) %>% grep(".assoc.linear",.,value = T)
# qq_nonspatial <- data.table(P=NA,sigma=NA,observed=NA,expected=NA,neighbors=NA)[-1]
# for(f in files){
#   sigma <- strsplit(f,"_") %>% unlist() %>% .[3] %>% as.numeric()
#   if(sigma < 4){
#     gwas <- fread(f)
#     gwas <- subset(gwas,TEST=="ADD")
#     expected <- -log10(1:length(gwas$P)/length(gwas$P))
#     observed <- -log10(sort(gwas$P))
#     qq_nonspatial <- rbind(qq_nonspatial,data.table(P=gwas$P,sigma=sigma,observed=observed,expected=expected,neighbors=4*pi*sigma^2*5))
#     print(f)
#   }
# }
# qq_nonspatial$phenotype <- "nonspatial"
# 
# qq <- rbind(qq_nonspatial,qq_clinal,qq_corner,qq_patchy)
# setwd("~/projects/spaceness/gwas/")
# 
# save(qq,file="qqplot_data.Rdata")
load("qqplot_data.Rdata")

conf <- data.frame(panel.qqconf(2e5,10000)) 

qqp <- ggplot(data=qq,aes(x=expected,y=observed,z=neighbors))+
  theme(strip.background = element_blank(),
        axis.text=element_text(size=6),
        axis.title.x=element_text(vjust = 2),
        axis.ticks=element_line(size=0.4),
        axis.title=element_text(size=7),
        strip.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text=element_text(size=6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,-30),
        axis.line.x=element_line(size=0.4),
        axis.line.y=element_line(size=0.4))+
  xlim(0,6)+ylim(0,10)+
  facet_wrap(~phenotype,nrow=2,ncol=2)+
  xlab(Expected~-log10(italic(p)))+ylab(Observed~-log10(italic(p)))+
  scale_fill_gradientn(colors=pal(50),name="Neighborhood\nSize")+
  #scale_fill_distiller(palette = "YlGnBu",name="Neighborhood\nSize")+
  stat_summary_hex(fun=median,bins=100)+
  #geom_line(data=conf[1:10000,],aes(x=X1,y=X2,z=NULL),linetype=2,lwd=0.35)+
  #geom_line(data=conf[10001:20000,],aes(x=X1,y=X2,z=NULL),linetype=2,lwd=0.35)+
  annotate(geom="segment",x=0,xend=6,y=0,yend = 6,lwd=0.35,col="black")
qqp <- qqp+guides(fill=guide_colorbar(barwidth = unit(27,"mm"),barheight = unit(4,"mm")))

###################################################################################################
################################# n significant SNPs by sigma #####################################
###################################################################################################
#summary of n significant snps by sigma (~20min)
# phenotypes <- c("out_normal","out_clinal","out_corner","out_patchy")
# pd <- data.frame(sigma=numeric(),
#                  n_sig=integer(),
#                  n_snps=integer(),
#                  prop_sig=numeric(),
#                  phenotype=character(),
#                  pc_corr=logical())
# for(x in phenotypes){
#   files <- list.files(x,full.names = T) %>% grep("assoc.linear",.,value = T)
#   for(i in files){
#     sigma <- strsplit(i,"_") %>% unlist() %>% .[3] %>% as.numeric()
#     if(sigma<4){
#       print(i)
#       gw <- fread(i)
#       gw <- subset(gw,TEST=="ADD")
#       gw$p_adj <- p.adjust(gw$P,method="fdr")
#       n_sig <- nrow(subset(gw,p_adj<=0.05))
#       n_snps <- nrow(gw)
#       prop_sig <- n_sig/n_snps
#       sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
#       pccorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig)
#       pccorr$phenotype <- x
#       pccorr$pc_corr <- "PC covariates"
#       pd <- rbind(pd,pccorr)
#     }
#   }
# 
#   files <- list.files(x,full.names = T) %>% grep("qassoc",.,value = T)
#   for(i in files){
#     sigma <- strsplit(i,"_") %>% unlist() %>% .[3] %>% as.numeric()
#     if(sigma<4){
#       print(i)
#       gw <- fread(i)
#       gw$p_adj <- p.adjust(gw$P,method="fdr")
#       n_sig <- nrow(subset(gw,p_adj<=0.05))
#       n_snps <- nrow(gw)
#       prop_sig <- n_sig/n_snps
#       sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
#       nocorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig,pc_corr=F)
#       nocorr$phenotype <- x
#       nocorr$pc_corr <- "Uncorrected"
#       pd <- rbind(pd,nocorr)
#     }
#   }
# }
# write.csv(pd,"gwas_summary_scaledmuts.csv",row.names = F)
pd <- read.csv("gwas_summary_scaledmuts.csv")

pd$neighbors <- 4*pi*pd$sigma^2*5
pd$max_neighbors <- 9*pi*pd$sigma^2*5
pd$log10_n_sig <- log(pd$n_sig,10)
pd$log10_n_sig[is.infinite(pd$log10_n_sig)] <- 0
pd$log10_prop_sig <- log(pd$prop_sig,10)
pd$log10_prop_sig[is.infinite(pd$log10_prop_sig)] <- 0
pd$phenotype <- mapvalues(pd$phenotype,
                          c("out_normal","out_clinal","out_corner","out_patchy"),
                          c("nonspatial","clinal","corner","patchy"))
pd$phenotype <- factor(pd$phenotype,levels=c("nonspatial","clinal","corner","patchy"))
nsig_by_sigma <- ggplot(data=pd,aes(x=neighbors,y=prop_sig,col=phenotype))+
  theme(strip.background = element_blank(),
        axis.text=element_text(size=6),
        axis.ticks=element_line(size=0.4),
        axis.title=element_text(size=7),
        axis.title.x=element_text(vjust = 2),
        strip.text = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text=element_text(size=6),
        legend.position = "bottom",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,-30),
        legend.key.height = unit(1,"mm"),
        axis.line.x=element_blank(),
        axis.line.y=element_line(size=0.4))+
  coord_cartesian(clip="off")+
  ylab("Proportion Significant SNPs (FDR < 0.05)")+xlab("Neighborhood Size")+
  scale_y_continuous(trans="log10")+
  scale_x_log10()+
  facet_wrap(~pc_corr,ncol=1,scales="free_y")+
  #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values=c("deepskyblue","green4","goldenrod3","red4"),name="Environment")+
  #geom_smooth(fill=NA,size=0.5)+
  geom_point(shape=21,stroke=0.35,size=0.4)
nsig_by_sigma <- nsig_by_sigma+guides(col=guide_legend(nrow=2,byrow=TRUE,override.aes = list(size=3,shape=16)))

##########################################################################################
################################ PC variation by sigma ###################################
##########################################################################################
pcfiles <- list.files("random_sampling/out_normal",full.names = T) %>% grep("pca_var_explained",.,value=T)
propvar <- data.frame(sigma=NA,pvx=NA)
for(f in pcfiles){
  sigma <- strsplit(f,"_") %>% unlist() %>% .[4] %>% as.numeric()
  a <- fread(f)
  propvar <- rbind(propvar,c(sigma,sum(a$V1)))
}
propvar$neighbors <- 4*5*pi*propvar$sigma^2
pvxplot <- ggplot(data=subset(propvar,sigma<4),aes(x=neighbors,y=pvx))+
                  scale_x_log10()+
                  xlab("Neighborhood Size")+
                  ylab("Proportion Variance\nExplained, PC 1-10")+
                  geom_line(lwd=0.5)
           
#########################################################################################
############################ sampling/phenotype maps ####################################
#########################################################################################
load("mapsdata.Rdata")
levels(l$type) <- c("nonspatial","clinal","corner","patchy")
maps <- ggplot(data=l,aes(x=x,y=y,col=phenotype,z=phenotype))+
  coord_fixed()+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())+
  facet_wrap(~type,nrow=1)+
  #scale_fill_gradientn(colors=pal(1000))+
  scale_color_distiller(palette = "RdYlBu",direction = 1)+
  #scale_color_distiller(palette = "RdYlBu",direction = -1,name="Phenotype",guide=F)+
  #scale_fill_distiller(palette = "RdYlBu",direction = -1,name="Phenotype",guide=F)+
  geom_point(size=0.25)
  #stat_summary_2d(fun="mean",bins = 40)
maps <- maps+guides(col=guide_colorbar(barwidth=unit(3,"mm"),barheight=unit(16,"mm")))

n <- ggplot(data=subset(l,type=="nonspatial"),aes(x=x,y=y,col=phenotype,alpha=phenotype^6))+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())+
  ggtitle("nonspatial")+coord_fixed()+
  scale_color_gradient2(low="white",high="deepskyblue3",guide=F)+
  scale_alpha_continuous(guide=F,range=c(0,1))+
  geom_point(size=0.35)
cl <- ggplot(data=subset(l,type=="clinal"),aes(x=x,y=y,col=phenotype,alpha=phenotype^8))+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())+
  ggtitle("clinal")+coord_fixed()+
  scale_color_gradient2(low="white",high="green4",guide=F)+
  scale_alpha_continuous(guide=F,range=c(0,1))+
  geom_point(size=0.35)
cr <- ggplot(data=subset(l,type=="corner"),aes(x=x,y=y,col=phenotype,alpha=phenotype^6))+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())+
  ggtitle("corner")+coord_fixed()+
  scale_color_gradient2(low="white",high="darkgoldenrod3",guide=F)+
  scale_alpha_continuous(guide=F,range=c(0,1))+
  geom_point(size=0.35)
pt <- ggplot(data=subset(l,type=="patchy"),aes(x=x,y=y,col=phenotype,alpha=phenotype^6))+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank())+
  ggtitle("patchy")+coord_fixed()+
  scale_color_gradient2(low="white",high="red4",guide=F)+
  scale_alpha_continuous(guide=F,range=c(0,1))+
  geom_point(size=0.35)


###################################################################################
###################### phenotype:distance correlations ############################
###################################################################################
# library(pbapply)
# pairs=combn(1:1000,2)
# clinal <- subset(l,type=="clinal")
# phendist <- c();spdist <- c();xdist <- c()
# dists <- pbapply(pairs[,1:10000],2,function(e){
#   xdist <- append(xdist,abs(clinal$x[e[1]]-clinal$x[e[2]]))
#   gdist <- append(spdist,spDists(x=as.matrix(clinal[e[1],c("x","y")]),y=as.matrix(clinal[e[2],c("x","y")])))
#   phendist <- append(phendist,abs(clinal$phenotype[e[1]]-clinal$phenotype[e[2]]))
#   c(gdist,phendist,xdist)
# })
# dists <- data.frame(t(dists))
# # spd <- clinal %>% as.matrix() %>% spDists() %>% as.dist() %>% c()
# # dists <- data.frame(phenotype=dists,spatial=spd)
# ggplot(dists,aes(x=X3,y=X2))+
#   scale_fill_distiller(palette = "YlOrBr")+
#   stat_bin_2d()+
#   xlab("Spatial Distance")+ylab("Phenotype Distance")+
#   geom_smooth(method="lm",fill=NA,col="black")

###################################################################################
############################## output summary pdf's ###############################
###################################################################################

#summary plots for n significant snps and qqplots
pdf("gwas_summary_loglog.pdf",width=6,height=5,useDingbats = F)
ggdraw()+
  #draw_plot(maps,0.06,0.72,0.6,.3)+
  draw_plot(n,0.06,0.72,0.15,0.3)+
  draw_plot(cl,0.2,0.72,0.15,0.3)+
  draw_plot(cr,0.35,0.72,0.15,0.3)+
  draw_plot(pt,0.5,0.72,0.15,0.3)+
  draw_plot(pvxplot,0.67,0.7,0.3,0.27)+
  draw_plot(nsig_by_sigma,0,0,0.4,0.73)+
  draw_plot(qqp,0.4,0,0.6,0.7)+
  draw_label("A",0.05,0.97)+
  draw_label("B",0.66,0.97)+
  draw_label("C",0.05,0.7)+
  draw_label("D",0.45,0.7)
dev.off()


#individual manhattan plots
sim <- "sigma_3.5976011606040497_.trees_4924299.assoc.linear"
palette <- wes_palette("Darjeeling1",4)
palette <- RColorBrewer::brewer.pal(4,"RdYlBu")

normal <- fread(paste0("random_sampling/out_normal/",sim))
normal <- subset(normal,TEST=="ADD")
normal$phenotype <- "nonspatial"
normal$q <- p.adjust(normal$P,"fdr")#
normal$cutoff <- (max(subset(normal,q<=0.05)$P)+min(subset(normal,q>0.05)$P))/2
if(is.infinite(normal$cutoff[1])) normal$cutoff <- 0.05/nrow(normal) #... must be a better way to do this

clinal <- fread(paste0("random_sampling/out_clinal/",sim))
clinal <- subset(clinal,TEST=="ADD")
clinal$phenotype <- "clinal"
clinal$q <- p.adjust(clinal$P,"fdr")
clinal$cutoff <- (max(subset(clinal,q<=0.05)$P)+min(subset(clinal,q>0.05)$P))/2
if(is.infinite(clinal$cutoff[1])) clinal$cutoff <- 0.05/nrow(normal) 

corner <- fread(paste0("random_sampling/out_corner/",sim))
corner <- subset(corner,TEST=="ADD")
corner$phenotype <- "corner"
corner$q <- p.adjust(corner$P,"fdr")
corner$cutoff <- (max(subset(corner,q<=0.05)$P)+min(subset(corner,q>0.05)$P))/2
if(is.infinite(corner$cutoff[1])) corner$cutoff <- 0.05/nrow(normal) 

patchy <- fread(paste0("random_sampling/out_patchy/",sim))
patchy <- subset(patchy,TEST=="ADD")
patchy$phenotype <- "patchy"
patchy$q <- p.adjust(patchy$P,"fdr")
patchy$cutoff <- (max(subset(patchy,q<=0.05)$P)+min(subset(patchy,q>0.05)$P))/2
if(is.infinite(patchy$cutoff[1])) patchy$cutoff <- 0.05/nrow(normal) 
 
# rsnps <- fread(paste0("random_snps/",sim))
# rsnps$phenotype <- "random_snps"
# rsnps$q <- p.adjust(rsnps$P,"fdr")
# rsnps$cutoff <- (max(subset(rsnps,q<=0.05)$P)+min(subset(rsnps,q>0.05)$P))/2
# if(is.infinite(rsnps$cutoff[1])) rsnps$cutoff <- 0.05/nrow(normal) 
# pheno_snps <- fread(paste0("random_snps/",gsub(".trees.assoc.linear",".treesphenotype_snp_indices.txt",sim)))

gwas <- rbind(normal,clinal,corner,patchy)
gwas$p_minus_log_10 <- -log(gwas$P,10)
gwas$phenotype <- factor(gwas$phenotype,levels=c("nonspatial","clinal","corner","patchy"))

manhattan_plots <- ggplot(data=gwas,aes(x=BP,y=p_minus_log_10))+
  theme(#axis.text.x=element_blank(),
        plot.title = element_text(size=6),
        #axis.ticks.x=element_blank(),
        #axis.title.x=element_blank(),
        strip.background = element_blank(),
        axis.ticks.y=element_line(size=0.4),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=7),
        #strip.text = element_blank(),#element_text(size=6),
        legend.title = element_text(size=6),
        legend.text=element_text(size=6),
        axis.line=element_line(size=0.4))+
  ggtitle("Neighborhood Size = 800, PCA covariates")+
  facet_grid(phenotype~.)+
  #scale_color_manual(values=pal(4),guide=F)+
  geom_point(size=0.2,stroke=0.2,shape=21,col="grey")+
  geom_point(data=subset(gwas,q<0.05),size=0.5,stroke=0.2,shape=21,col="red")+
  geom_hline(yintercept=-log(0.05/nrow(normal),10),linetype=3,size=0.35)+
  geom_hline(aes(yintercept=-log10(cutoff)),linetype=2,size=0.35)+
  #geom_line(aes(y=min(-log(cutoff,10)),size=0.5))+
  ylab(expression(italic(-log[10](p))))

png("~/projects/spaceness/figures/manhattan_plots_NS800.png",width=3,height=4,units="in",res=300)
print(manhattan_plots)
dev.off()

#maps
locs_name <- tools::file_path_sans_ext(sim) %>% tools::file_path_sans_ext() %>% paste0("_locs.txt")
phen_name <- tools::file_path_sans_ext(sim) %>% tools::file_path_sans_ext() %>% paste0(".phenotypes")
pc_name <- tools::file_path_sans_ext(sim) %>% tools::file_path_sans_ext() %>% paste0(".pca")
l <- data.frame(x=numeric(),y=numeric(),ind=character(),ind2=character(),
                phenotype=numeric(),PC1=numeric(),PC2=numeric(),type=character())
for(i in c("out_normal","out_clinal","out_corner","out_patchy")){
  locs <- fread(paste0(i,"/",locs_name))
  phen <- fread(paste0(i,"/",phen_name))
  pc <- fread(paste0(i,"/",pc_name))[,3:4]
  locsdf <- cbind(locs,phen,pc)
  locsdf$type <- i
  colnames(locsdf) <- c("x","y","ind","ind2","phenotype","PC1","PC2","type")
  l <- rbind(l,locsdf)
}
l$type <- gsub("out_","",l$type)
l$type <- factor(l$type,levels=c("normal","clinal","corner","patchy"))
maps <- ggplot(data=l,aes(x=x,y=y,z=phenotype,fill=phenotype))+
  facet_wrap(~type,nrow=1)+
  ggtitle("Spatial Position")+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)+
  geom_point(shape=21,stroke=0.05,size=0.5)
pcplots <- ggplot(data=l,aes(x=PC1,y=PC2,fill=phenotype))+
  facet_wrap(~type,ncol=1,scales="free")+
  theme_minimal()+
  theme(axis.text = element_blank(),
        plot.title = element_text(size=6),
        axis.ticks = element_blank(),
        axis.title.y=element_text(size=5,margin=margin(r=-1)),
        axis.title.x=element_text(size=5,margin=margin(t=-1)),
        legend.title = element_text(size=6),
        legend.text=element_text(size=6),
        legend.position = "right",
        legend.box.margin = margin(l=-11),
        strip.text = element_blank(),#element_text(size=6,color="white"),
        #panel.grid.major=element_line(size=0.15,color="black"),
        #panel.grid.minor=element_line(size=0.15,color="black"),
        strip.background = element_blank())+
  ggtitle("PC position")+
  scale_fill_distiller(palette = "RdYlBu",direction = -1)+
  geom_point(shape=21,stroke=0.05,size=0.5)
pcplots <- pcplots+guides(fill=guide_colorbar(barwidth = unit(3,"mm"),barheight = unit(20,"mm")))

#full plot w example sim
png("gwas_summary_sigma023.png",res=400,width=6.5,height=2.25,unit="in")
ggdraw()+
  draw_plot(nsig_by_sigma,0,0,0.3,1)+
  draw_plot(manhattan_plots,0.28,0,0.42,1)+
  draw_plot(maps,0.69,-0.01,0.125,1)+
  draw_plot(pcplots,0.81,-0.01,0.2,1)+
  draw_label("A",0.025,0.96,size = 9)+
  draw_label("B",0.325,0.96,size = 9)+
  draw_label("C",0.69,0.96,size = 9)
dev.off()


  
  
