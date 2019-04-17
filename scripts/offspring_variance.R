#estimating and plotting genealogical parameters (var(nOff),census_n,gen_time)
library(data.table);library(magrittr);library(plyr);library(progress)
library(foreach);library(doMC)
setwd("~/spaceness/")
source("~/R/colordict.R")
registerDoMC(cores=10)

colornames <- list('red_orange','dark_medici_blue')
theme_set(theme_classic()+theme(axis.text=element_text(size=8),
                                axis.title=element_text(size=8),
                                strip.text = element_text(size=8),
                                title = element_text(size=9,hjust=0.5),
                                legend.title = element_text(size=8),
                                legend.text = element_text(size=8),
                                strip.background = element_blank(),
                                axis.title.y=element_blank()))

#get variance in offspring by model and sigma
files <- list.files("~/Desktop/spaceviz/10kpeds/",full.names = T)
getvaroff <- function(f) {
  ped <- fread(f)
  sigma <- strsplit(basename(f),"_") %>% unlist() %>% .[3] %>% as.numeric()
  model <- strsplit(basename(f),"_") %>% unlist() %>% .[1]
  ped <- ddply(ped,.(ind_ID),function(e) e[e$timestep==max(e$timestep),]) #trim to last generation entry for each individual
  ped <- subset(ped,timestep<150)
  data.frame(sigma=sigma,model=model,varoff=var(ped$n_offspring))
}
varoff <- foreach(f=files,.combine = rbind) %dopar% getvaroff(f)
varoff$model <- mapvalues(varoff$model,from=c("rm","sp"),to=c("random mating","spatial"))
varoff$stat <- "var(n offspring)"
varoff <- varoff[,c(1,3,2,4)]

#estimate number of individuals born per timestep
ped <- fread(files[1])
nbirths <- c()
for(i in 130:170){
  t <- subset(ped,timestep==i)
  nbirths <- append(nbirths,sum(!t$ind_ID %in% ped$ind_ID[ped$timestep<i]))
}
print(mean(nbirths))


#census population size
popsizes <- fread("~/spaceness/sumstats/ss_spatial_random_W50.txt.popsizes")
popsizes$model <- "spatial"
rmpops <- fread("~/spaceness/sumstats/ss_randmates_random_W50.txt.popsizes")
rmpops$model <- "random mating"
popsizes <- rbind(popsizes,rmpops)
popsizes$stat <- "Census N"
names(popsizes) <- names(varoff)

gen_times <- fread("~/spaceness/W50sp_gentimes.txt")
gen_times$model <- "spatial"
gen_times2 <- fread("~/spaceness/W50rm_gentimes.txt")
gen_times2$model <- "random mating"
gen_times <- rbind(gen_times,gen_times2)
gen_times <- gen_times[,c(2,1,3)]
gen_times$stat <- "Generation Time"
names(gen_times) <- names(popsizes)

pd <- rbind(popsizes,gen_times,varoff)
pd$neighbors <- 4*pi*pd$sigma^2*5
names(pd) <- c("sigma","value","model","stat","neighbors")
pd <- subset(pd,sigma<4)
p <- ggplot(data=pd,aes(x=neighbors,y=value,col=model))+
  xlab("Neighborhood Size")+
  # theme(legend.position="bottom",legend.direction="horizontal",
  #       legend.box.margin = margin(-10,0,0,-30),
  #       axis.text.x=element_text(size=7,hjust=0.7))+
  facet_wrap(~stat,scales="free")+
  scale_x_log10()+
  scale_color_manual(values=getcolordict(colornames))+
  geom_point(shape=21,stroke=0.3,size=0.8)+
  geom_smooth(lwd=0.5,fill=NA,span=.4)
p <- p+guides(col=guide_legend(override.aes = list(size=4)))
pdf("figures/pop_params.pdf",width=5.5,height=1.5)
print(p)
dev.off()

save(pd,file="~/Desktop/spaceviz/gen_params_sumary.Rdata")
