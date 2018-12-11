#visualizing GWAS on SLiM simulations
library(ggplot2);library(magrittr);library(plyr);library(data.table)
library(qvalue);library(cowplot);library(wesanderson)
setwd("~/spaceness/gwas/")
source("~/R/ggthemes.R")
#individual manhattan plots
sim <- "out_normal/sigma_0.613949565110263_.trees1500000.trees.assoc.linear" %>% gsub("out_normal/","",.)
palette <- wes_palette("Darjeeling1",4)
palette <- RColorBrewer::brewer.pal(4,"RdYlBu")

normal <- fread(paste0("out_normal/",sim))
normal$phenotype <- "nonspatial"
normal$q <- p.adjust(normal$P,"fdr")#qvalue(normal$P)$qvalues
normal$cutoff <- (max(subset(normal,q<=0.05)$P)+min(subset(normal,q>0.05)$P))/2
if(is.infinite(normal$cutoff[1])) normal$cutoff <- 0.05/nrow(normal) #... must be a better way to do this

clinal <- fread(paste0("out_clinal/",sim))
clinal$phenotype <- "clinal"
clinal$q <- p.adjust(clinal$P,"fdr")#qvalue(clinal$P)$qvalues
clinal$cutoff <- (max(subset(clinal,q<=0.05)$P)+min(subset(clinal,q>0.05)$P))/2
if(is.infinite(clinal$cutoff[1])) clinal$cutoff <- 0.05/nrow(normal) 

corner <- fread(paste0("out_corner/",sim))
corner$phenotype <- "corner"
corner$q <- p.adjust(corner$P,"fdr")#qvalue(corner$P)$qvalues
corner$cutoff <- (max(subset(corner,q<=0.05)$P)+min(subset(corner,q>0.05)$P))/2
if(is.infinite(corner$cutoff[1])) corner$cutoff <- 0.05/nrow(normal) 

patchy <- fread(paste0("out_patchy/",sim))
patchy$phenotype <- "patchy"
patchy$q <- p.adjust(patchy$P,"fdr")#qvalue(patchy$P)$qvalues
patchy$cutoff <- (max(subset(patchy,q<=0.05)$P)+min(subset(patchy,q>0.05)$P))/2
if(is.infinite(patchy$cutoff[1])) patchy$cutoff <- 0.05/nrow(normal) 

rsnps <- fread(paste0("random_snps/",sim))
rsnps$phenotype <- "random_snps"
rsnps$q <- p.adjust(rsnps$P,"fdr")
rsnps$cutoff <- (max(subset(rsnps,q<=0.05)$P)+min(subset(rsnps,q>0.05)$P))/2
if(is.infinite(rsnps$cutoff[1])) rsnps$cutoff <- 0.05/nrow(normal) 
pheno_snps <- fread(paste0("random_snps/",gsub(".trees.assoc.linear",".treesphenotype_snp_indices.txt",sim)))

gwas <- rbind(normal,clinal,corner,patchy)
gwas$p_minus_log_10 <- -log(gwas$P,10)
gwas$phenotype <- factor(gwas$phenotype,levels=c("nonspatial","clinal","corner","patchy"))

manhattan_plots <- ggplot(data=gwas,aes(x=BP,y=p_minus_log_10,col=phenotype))+
  theme(axis.text.x=element_blank(),
        plot.title = element_text(size=6),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        strip.background = element_blank(),
        axis.ticks.y=element_line(size=0.4),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=6),
        strip.text = element_blank(),#element_text(size=6),
        legend.title = element_text(size=6),
        legend.text=element_text(size=6),
        axis.line=element_line(size=0.4))+
  ggtitle("Neighbors = 5, PCA covariates")+
  facet_grid(phenotype~.)+
  scale_color_manual(values=palette,guide=F)+
  geom_point(size=0.2,stroke=0.2,shape=21)+
  geom_point(data=subset(gwas,q<0.05),size=0.2,stroke=0.2,shape=21,color="black")+
  geom_hline(yintercept=-log(0.05/nrow(normal),10),linetype=2,size=0.35)+
  #geom_line(aes(y=min(-log(cutoff,10)),size=0.5))+
  ylab(expression(italic(-log[10](p))))
  
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
  facet_wrap(~type,ncol=1)+
  theme_minimal()+
  theme(axis.text = element_blank(),
        plot.title = element_text(size=6),
        axis.ticks = element_blank(),
        axis.title.y=element_text(size=5,margin=margin(r=-1)),
        axis.title.x=element_text(size=5,margin=margin(t=-1)),
        legend.title = element_text(size=6),
        legend.text=element_text(size=6),
        legend.position = "none",
        strip.text = element_blank(),#element_text(size=6,color="white"),
        #panel.grid.major=element_line(size=0.15,color="black"),
        #panel.grid.minor=element_line(size=0.15,color="black"),
        strip.background = element_blank())+
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
  
  
#summary of n significant snps by sigma (~20min)
#load("gwas_summary.Rdata")
phenotypes <- c("out_normal","out_clinal","out_corner","out_patchy")
pd <- data.frame(sigma=numeric(),
                 n_sig=integer(),
                 n_snps=integer(),
                 prop_sig=numeric(),
                 phenotype=character(),
                 pc_corr=logical())
for(x in phenotypes){
  files <- list.files(x,full.names = T) %>% grep("assoc.linear",.,value = T)
  for(i in files){
    gw <- fread(i)
    gw$p_adj <- p.adjust(gw$P,method="fdr")
    n_sig <- nrow(subset(gw,p_adj<=0.05))
    n_snps <- nrow(gw)
    prop_sig <- n_sig/n_snps
    sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    pccorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig)
    pccorr$phenotype <- x
    pccorr$pc_corr <- "PC covariates"
    pd <- rbind(pd,pccorr)
  }
  
  files <- list.files(x,full.names = T) %>% grep("qassoc",.,value = T)
  for(i in files){
    gw <- fread(i)
    gw$p_adj <- p.adjust(gw$P,method="fdr")
    n_sig <- nrow(subset(gw,p_adj<=0.05))
    n_snps <- nrow(gw)
    prop_sig <- n_sig/n_snps
    sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    nocorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig,pc_corr=F)
    nocorr$phenotype <- x
    nocorr$pc_corr <- "Uncorrected"
    pd <- rbind(pd,nocorr)
  }
}

pd$neighbors <- ((pd$sigma^2*pi)/50^2)*1.1e4
pd$max_neighbors <- (((pd$sigma*3)^2*pi)/50^2)*1.1e4
pd$log10_n_sig <- log(pd$n_sig,10)
pd$log10_n_sig[is.infinite(pd$log10_n_sig)] <- 0
pd$phenotype <- mapvalues(pd$phenotype,
                          c("out_normal","out_clinal","out_corner","out_patchy"),
                          c("nonspatial","clinal","corner","patchy"))
pd$phenotype <- factor(pd$phenotype,levels=c("nonspatial","clinal","corner","patchy"))
nsig_by_sigma <- ggplot(data=subset(pd,sigma>0.3),aes(x=neighbors,y=n_sig,col=phenotype))+
  theme(strip.background = element_blank(),
                        axis.text=element_text(size=6),
                        axis.ticks=element_line(size=0.4),
                        axis.title=element_text(size=6),
                        axis.title.x=element_text(vjust = 3),
                        strip.text = element_text(size=6),
                        legend.title = element_text(size=6),
                        legend.text=element_text(size=6),
                        legend.position = "bottom",
                        legend.margin=margin(0,0,0,0),
                        legend.box.margin=margin(-20,0,0,-30),
                        legend.key.height = unit(1,"mm"),
                        axis.line.x=element_blank(),
                        axis.line.y=element_line(size=0.4))+
  coord_cartesian(clip="off")+
  ylab("significant SNPs (FDR < 0.05)")+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~pc_corr,ncol=1,scales="free_y")+
  scale_color_manual(values=palette)+
  geom_smooth(fill=NA,size=0.5)+
  geom_point(shape=21,stroke=0.35,size=0.4)
nsig_by_sigma <- nsig_by_sigma+guides(col=guide_legend(nrow=2,byrow=TRUE))

#booya
png("gwas_summary.png",res=1000,width=6.5,height=2.25,unit="in")
ggdraw()+
  draw_plot(nsig_by_sigma,0,0,0.3,1)+
  draw_plot(manhattan_plots,0.28,0,0.42,1)+
  draw_plot(maps,0.69,-0.01,0.125,1)+
  draw_plot(pcplots,0.81,-0.01,0.2,1)+
  draw_label("A",0.025,0.96,size = 9)+
  draw_label("B",0.325,0.96,size = 9)+
  draw_label("C",0.69,0.96,size = 9)
dev.off()
 

  
  
