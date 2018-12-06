#visualizing GWAS on SLiM simulations
library(ggplot2);library(magrittr);library(plyr);library(data.table)
setwd("~/spaceness/gwas/")

#individual manhattan plot
gwas <- fread("out/sigma_0.3707859307192668_.trees1500000.trees.assoc.linear")

gwas <- subset(gwas,TEST=="ADD")
nrow(gwas)

gwas$p_adj <- p.adjust(gwas$P,method="fdr")
gwas$p_minus_log_10 <- -log(gwas$P,10)
  
  
png("gwas_sim_spatial_phen_PCA_corr_sigma0905.png",width=6.5,height=1,units="in",res=1200)
ggplot(data=gwas,aes(x=BP,y=p_minus_log_10))+
  theme_classic()+theme(axis.text=element_text(size=7),
                        axis.title = element_text(size=7),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title.x=element_blank())+
  geom_point(size=0.2,shape=21,stroke=0.4)+
  ylab(expression(italic(-log[10](p))))+
  geom_hline(yintercept=-log(.05/nrow(gwas),10),col="grey",linetype=2)+
  annotate(geom="text",label=expression(italic(p[corr]==0.05)),
           x=4e6,y=-log(.05/nrow(gwas),10)-0.07*max(gwas$p_minus_log_10),size=1.5,col="white")
dev.off()


#how does the number of genome-wide significant SNPs change with dispersal? 
files <- list.files("out",full.names = T) %>% grep("assoc.linear",.,value = T)
pd <- data.frame(sigma=numeric(),n_sig=integer(),n_snps=integer(),prop_sig=numeric())
for(i in files){
  gwas <- fread(i)
  gwas$p_adj <- p.adjust(gwas$P,method="fdr")
  n_sig <- nrow(subset(gwas,p_adj<=0.05))
  n_snps <- nrow(gwas)
  prop_sig <- n_sig/n_snps
  sigma <- strsplit(i,"_") %>% unlist() %>% .[2] %>% as.numeric()
  pd <- rbind(pd,data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig))
}
pd$pc_corr <- T

files <- list.files("out",full.names = T) %>% grep("qassoc",.,value = T)
for(i in files){
  gwas <- fread(i)
  gwas$p_adj <- p.adjust(gwas$P,method="fdr")
  n_sig <- nrow(subset(gwas,p_adj<=0.05))
  n_snps <- nrow(gwas)
  prop_sig <- n_sig/n_snps
  sigma <- strsplit(i,"_") %>% unlist() %>% .[2] %>% as.numeric()
  pd <- rbind(pd,data.frame(sigma=sigma,n_sig=n_sig,n_snps=n_snps,prop_sig=prop_sig,pc_corr=F))
}

ggplot(data=pd,aes(x=sigma,y=n_sig,col=pc_corr))+
  theme_classic()+
  #facet_wrap(~pc_corr,ncol=1,scales="free_y")+
  geom_point()


