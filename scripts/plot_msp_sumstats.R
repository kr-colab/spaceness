#msp grid sims summary stat viz 
setwd("~/spaceness/")
library(ggplot2);library(plyr);library(magrittr);library(data.table)

theme_set(theme_classic()+theme(axis.text=element_text(size=8),
                                axis.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                title = element_text(size=9,hjust=0.5),
                                legend.title = element_text(size=8),
                                legend.text = element_text(size=8),
                                strip.background = element_blank()))

ss <- fread("sumstats/ss_msp_grid.txt",data.table = F)

ss <- ss[,!grepl("sfs",colnames(ss))]
mss <- melt(ss,id.vars=c("sigma"))

ggplot(data=mss,aes(x=sigma,y=value))+
  xlab("sigma")+
  scale_x_log10(labels = function(x) format(x, scientific = FALSE))+
  facet_wrap(~variable,scales="free")+
  geom_point(shape=21,stroke=0.3,size=1)+geom_smooth(span=.2,fill=NA,col="orange",lwd=0.5)

spaceness <- fread("sumstats/ss_spatial_random_W50.txt")
spaceness <- spaceness[,c("sigma","segsites","pi","thetaW","tajD","het_o","gen_sp_corr")]
spaceness$width <- Inf

ss5 <- fread("sumstats/ss_msp_grid_Ne6100_w5.txt",data.table=F)
ss5 <- ss5[,!grepl("sfs",colnames(ss5))]
ss5$width <- 5
ss10 <- fread("sumstats/ss_msp_grid_Ne6100_w10.txt",data.table=F)
ss10 <- ss10[,!grepl("sfs",colnames(ss10))]
ss10$width <- 10
ss20 <- fread("sumstats/ss_msp_grid_Ne6100_w20.txt",data.table=F)
ss20 <- ss20[,!grepl("sfs",colnames(ss20))]
ss20$width <- 20
#ss50 <- fread("sumstats/ss_msp_grid_Ne6100_w50.txt",data.table=F)
#ss50 <- ss50[,!grepl("sfs",colnames(ss50))]
#ss50$width <- 50
ss <- rbind(ss5,ss10,ss20)
ss <- ss[,names(ss) != "fis"]
ss <- rbind(ss,spaceness)
ss <- ss[ss$sigma<=4,]
mss <- melt(ss,id.vars=c("sigma","width"))

sscurves <- ggplot(data=mss,aes(x=4*pi*sigma^2*5,y=value,color=factor(width)))+
  #theme(legend.position="bottom")+
  scale_color_manual(values=c("steelblue3","khaki3","red3","black"),name="Demes Per\nSide")+
  xlab("sigma")+
  scale_x_log10(labels = function(x) format(x, scientific = FALSE))+
  facet_wrap(~variable,scales="free")+
  #geom_point(shape=21,stroke=0.15,size=1)+
  geom_smooth(span=.2,fill=NA,lwd=0.5)

pdf("figures/grid_ss_curves.pdf",useDingbats = F,width=6,height=3)
print(sscurves)
dev.off()

get_grid_coords <- function(width){
  x <- c();y <- c()
  for(i in 1:width){
    for(j in 1:width){
      x <- append(x,i)
      y <- append(y,j)
    }
  }
  out <- data.frame(x,y)
  return(out)
}
a <- get_grid_coords(5)
a$width <- 5
b <- get_grid_coords(10)
b$width <- 10
c <- get_grid_coords(20)
c$width <- 20
d <- rbind(a,b,c)

maps <- ggplot(data=d,aes(x=x,y=y,fill=factor(width)))+
  theme(strip.text = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank(),axis.title = element_blank(),
        axis.line = element_blank())+
  facet_wrap(~width,scales="free",ncol=1)+
  scale_fill_manual(values=c("steelblue3","khaki3","red3","black"),guide=F)+
  #coord_equal()+
  geom_point(shape=21,size=1.25,stroke=0.01,color="grey")

pdf("figures/msp_grid_maps.pdf",width=1,height=3,useDingbats = F)
print(maps)
dev.off()





