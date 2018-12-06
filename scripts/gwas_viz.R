#visualizing GWAS on SLiM simulations
library(ggplot2);library(magrittr);library(plyr);library(data.table);library(fdrtool)
setwd("~/spaceness/gwas/")

#individual manhattan plot
gwas <- fread("~/spaceness/gwas/out_phen_corner/sigma_0.9050194305800039_.trees1500000.trees.assoc.linear")
gwas <- subset(gwas,TEST=="ADD")
#gwas$p_adj <- p.adjust(gwas$P,method="fdr")
gwas$p_adj <- qvalue(gwas$P)$qvalues
gwas$p_minus_log_10 <- -log(gwas$P,10)
cutoff <- (max(subset(gwas,p_adj<=0.05)$P)+min(subset(gwas,p_adj>0.05)$P))/2
if(is.infinite(cutoff)) cutoff <- 0.05/nrow(gwas) #... can you calculate an FDR qval if it's outside the range of the data? 

png("gwas_sim_spatial_phen_corner_PCA_corr.png",width=6.5,height=1,units="in",res=1200)
ggplot(data=gwas,aes(x=BP,y=p_minus_log_10))+
  theme_classic()+theme(axis.text=element_text(size=7),
                        axis.title = element_text(size=7),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.x=element_blank())+
  geom_point(size=0.2,alpha=0.8)+
  geom_hline(yintercept = -log(cutoff,10),linetype=2,col="grey",size=0.35)+
  #geom_point(data=subset(gwas,p_adj<0.05),size=0.2,shape=21,stroke=0.4,col="red")+
  ylab(expression(italic(-log[10](p))))
  #annotate(geom="text",label=expression(italic(FDR==0.05)),
  #         x=4e6,y=-log(.05/nrow(gwas),10)-0.07*max(gwas$p_minus_log_10),size=1.5,col="white")
dev.off()

#how does the number of genome-wide significant SNPs change with dispersal? 
phenotypes <- c("out_gaussian","out","out_phen_corner")
pd <- data.frame(sigma=numeric(),n_sig=integer(),n_snps=integer(),prop_sig=numeric(),phenotype=character(),pc_corr=logical())
for(x in phenotypes){
  files <- list.files(x,full.names = T) %>% grep("assoc.linear",.,value = T)
  for(i in files){
    gwas <- fread(i)
    gwas$p_adj <- p.adjust(gwas$P,method="fdr")
    n_sig <- nrow(subset(gwas,p_adj<=0.05))
    n_snps <- nrow(gwas)
    prop_sig <- n_sig/n_snps
    sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    #sigma <- strsplit(i,"_") %>% unlist() %>% .[2] %>% as.numeric()
    pccorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig)
  }
  pccorr$phenotype <- x
  pccorr$pc_corr <- T
  
  files <- list.files(x,full.names = T) %>% grep("qassoc",.,value = T)
  for(i in files){
    gwas <- fread(i)
    gwas$p_adj <- p.adjust(gwas$P,method="fdr")
    n_sig <- nrow(subset(gwas,p_adj<=0.05))
    n_snps <- nrow(gwas)
    prop_sig <- n_sig/n_snps
    sigma <- strsplit(i,"sigma_") %>% unlist() %>% .[2] %>% strsplit("_") %>% unlist() %>% .[1] %>% as.numeric()
    nocorr <- data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig,pc_corr=F)
  }
  nocorr$phenotype <- x
  nocorr$pc_corr <- F
  pd <- rbind(pd,pccorr)
  pd <- rbind(pd,nocorr)
}

pd$neighbors <- ((pd$sigma^2*pi)/50^2)*1.1e4
pd$max_neighbors <- (((pd$sigma*3)^2*pi)/50^2)*1.1e4
ggplot(data=pd,aes(x=max_neighbors,y=log(n_sig,10),col=phenotype))+
  theme_classic()+theme(strip.background = element_blank())+
  ylab(expression(italic(log[10]("n significant SNPs"))))+
  facet_wrap(~pc_corr,ncol=1,scales="free_y")+
  geom_smooth(fill=NA,method="glm",formula=y~log(x))+
  geom_point()



  
  
