#exploratory analysis for gwas on phenotypes created from 100 additive SNPs

#indices of snps affecting phenotype
pheno_snps <- fread("~/Desktop/sigma_0.6544936332655427_.trees1500000.treesphenotype_snp_indices.txt")
pheno_snps <- pheno_snps$V1+1 #bc python 0-indexing

rsnps <- fread("~/Desktop/sigma_0.6544936332655427_.trees1500000.trees.qassoc",data.table=F)
rsnps$phenotype <- "random_snps"
rsnps$q <- p.adjust(rsnps$P,"fdr")
rsnps$fdr_cutoff <- (max(subset(rsnps,q<=0.05)$P)+min(subset(rsnps,q>0.05)$P))/2
rsnps$bonferroni_cutoff <- 0.05/nrow(rsnps)
rsnps$PCcovar <- F
rsnps$pheno_snp <- F
rsnps$pheno_snp[pheno_snps] <- T

rsnps2 <- fread("~/Desktop/sigma_0.6544936332655427_.trees1500000.trees.assoc.linear",data.table=F)
rsnps2$phenotype <- "random_snps"
rsnps2$q <- p.adjust(rsnps2$P,"fdr")
rsnps2$fdr_cutoff <- (max(subset(rsnps2,q<=0.05)$P)+min(subset(rsnps2,q>0.05)$P))/2
rsnps2$bonferroni_cutoff <- 0.05/nrow(rsnps2)
rsnps2$PCcovar <- T
rsnps2$pheno_snp <- F
rsnps2$pheno_snp[pheno_snps] <- T

rsnps3 <- rbind(rsnps[,c("BP","P","fdr_cutoff","bonferroni_cutoff","PCcovar","q","pheno_snp")],
               rsnps2[,c("BP","P","fdr_cutoff","bonferroni_cutoff","PCcovar","q","pheno_snp")])
rsnps3$minus_log_10_P <- -log(rsnps$P,10)

ggplot(data=rsnps3,aes(x=BP,y=minus_log_10_P))+
  facet_wrap(~PCcovar,ncol=1)+
  geom_point(shape=1,alpha=0.6)+
  geom_point(data=subset(rsnps3,pheno_snp==T),col="red")+
  geom_hline(yintercept=-log(.05/nrow(rsnps2),10),linetype=2,col="orange")

#how similar is the p distribution for phenotype vs neutral snps?
wilcox.test(rsnps$P,rsnps[pheno_snps,]$P)

#n significant hits by correction and true snp status
ddply(rsnps3,.(PCcovar,pheno_snp),function(e){
  nrow(subset(e,P<bonferroni_cutoff))/nrow(e)
})



