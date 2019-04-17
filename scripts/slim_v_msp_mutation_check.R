#plots to compare msprime and slim mutation effects on summary stats
setwd("~/spaceness/sumstats")

slim_spat <- fread("ss_spatial_random_W16muts.txt")
slim_spat$model <- "spatial"
slim_spat$muts <- "slim"
msp_spat <- fread("ss_spatial_random_W16_msp_muts.txt")
msp_spat$model <- "spatial"
msp_spat$muts <- "msp"
slim_rand <- fread("ss_random_random_W16muts.txt")
slim_rand$model <- "random"
slim_rand$muts <- "slim"
msp_rand <- fread("ss_random_random_W16_msp_muts.txt")
msp_rand$model <- "random"
msp_rand$muts <- "msp"
ss <- rbind(slim_spat,msp_spat,slim_rand,msp_rand)
ss <- data.frame(ss)[,!grepl("sfs",colnames(ss))]

popsize_spat <- fread("ss_spatial_random_W16muts.txt.popsizes")
popsize_rand <- fread("ss_random_random_W16muts.txt.popsizes")
popsize <- rbind(popsize_spat,popsize_rand)
colnames(popsize) <- c("sigma","census_n")
ss <- merge(ss,popsize,by="sigma")
ss$theta <- 4*(ss$census_n/2)*1e-8
ss$census_n2 <- ss$census_n

mss <- melt(ss,id.vars=c("model","muts","sigma","census_n2"))
pd <- subset(mss,variable %in% c("pi","theta"))
p <- ggplot(data=pd,aes(x=sigma,y=value,col=paste(variable,model)))+
  ggtitle("Ne estimates for constant mu applied in msprime v SLiM")+
  facet_wrap(~muts)+ylab("")+
  scale_color_brewer(palette = "Paired")+
  geom_point(shape=21,size=0.8,stroke=0.2)+
  geom_smooth(fill=NA,lwd=0.5,span=0.6)

pdf("~/spaceness/figures/msp_v_slim_mutations.pdf",width=6,height=2.5,useDingbats = F)
print(p)
dev.off()



mss$variable <- factor(mss$variable,levels=c("segsites","pi","thetaW","theta","tajD",
                                             "het_o","fis","gen_dist_var","gen_dist_skew",
                                             "ibs_mean","ibs_var","ibs_skew","ibs_95p",
                                             "gen_sp_corr","ibs_mean_spat_corr","ibs_1e6blocks_spat_corr","ibs_skew_spat_corr",
                                             "census_n"))
mss$variable <- mapvalues(mss$variable,from=c("segsites","pi","thetaW","theta","tajD",
                                              "het_o","fis","gen_dist_var","gen_dist_skew",
                                              "ibs_mean","ibs_var","ibs_skew","ibs_95p",
                                              "gen_sp_corr","ibs_mean_spat_corr","ibs_1e6blocks_spat_corr","ibs_skew_spat_corr",
                                              "census_n"),
                          to=c("segregating\nsites","pi","Watterson's\ntheta","`4Nmu`","Tajima's\nD",
                               "observed\nheterozygosity","Fis","var(dxy)","skew(dxy)",
                               "mean(IBS)","var(IBS)","skew(IBS)","ibs_95p",
                               "corr(dxy,dist)","corr(mean(IBS),dist)","corr(IBS>1e6,dist)","corr(skew(IBS),dist)",
                               "census N"))

ggplot(data=subset(mss,muts=="slim"),aes(x=sigma,y=value,col=model))+
  facet_wrap(~variable,scales="free")+
  geom_point(shape=21,size=0.8,stroke=0.2)+
  geom_smooth(fill=NA,lwd=0.5)

