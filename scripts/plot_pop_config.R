#msprime population configuration plotting script
library(ggplot2);library(argparse)
parser <- ArgumentParser(description='plot population configuration',python_cmd="python")
parser$add_argument('--migration_matrix',metavar='migmat',type="character",
                    help='path to model params file (comma-delimited text with columns x,y,Ne,mright,mleft,mup,mdown)')
parser$add_argument('--outpath',metavar="outpath",type="character",
                    help='output file path')
args <- parser$parse_args()

#debug params
args <- list(migmat="~/Desktop/mig_matrix.txt",outpath="~/Desktop/test_pop_config.pdf")

##### convert migration matrix to param table for plotting #####
mig_matrix_to_plot_params <- function(mig_matrix_path){
  mig_matrix <- read.table(mig_matrix_path)
  npops <- ncol(mig_matrix)
  width <- sqrt(npops)
  params <- data.frame(matrix(ncol=6,nrow=npops))
  for(i in 1:npops){
    x <- i-width*(ceiling(i/width)-1)
    y <- ceiling(i/width)
    mr <- mig_matrix[i,i+1]
    ml <- mig_matrix[i,i-1]
    mu <- mig_matrix[i,i+width]
    md <- mig_matrix[i,i-width]
    row <- sapply(list(x,y,mr,ml,mu,md),function(e) if(class(e)!="numeric") 0 else e)
    params[i,] <- row
  }
  colnames(params) <- c("x","y","mr","ml","mu","md")
  return(params)
}

plot_params_to_mig_matrix <- function(plot_params,to_file=T,outpath){
  npops <- nrow(plot_params)
  width <- sqrt(npops)
  migmat <- matrix(data=0,nrow=npops,ncol=npops)
  for(i in 1:nrow(plot_params)){
    row <- plot_params[i,]
    if(i+1<=npops) migmat[i,i+1] <- row$mr
    if(i-1>0) migmat[i,i-1] <- row$ml
    if(i+width<=npops) migmat[i,i+width] <- row$mu
    if(i-width>0) migmat[i,i-width] <- row$md
  }
  if(to_file) write.table(migmat,outpath,col.names = F,row.names = F)
  return(migmat)
}

new_grid_params <- function(width,rate){
  return(data.frame(x=rep(1:width,width),
                    y=c(sapply(1:width,function(e) rep(e,width))),
                    mr=rate,ml=rate,mu=rate,md=rate)
  )
}

params <- new_grid_params(10,1/(1e5/100))
#params$mu[params$y==4] <- 0
#params$md[params$y==5] <- 0
params$mr[params$x==4] <- 0
params$ml[params$x==5] <- 0
mig_matrix_out <- plot_params_to_mig_matrix(params,T,"~/Desktop/mig_matrix_var.txt")




#plot
ggplot(data=params,aes(x=x,y=y))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())+
  geom_point(shape=1,size=6,aes())+
  xlim(1,sqrt(nrow(params)))+ylim(1,sqrt(nrow(params)))+
  #scale_color_viridis_c()+
  scale_color_distiller(palette="RdYlBu",name="Migration Rate")+
  geom_segment(aes(x=x+.25,xend=x+.75,y=y,yend=y,col=mr),arrow=arrow(length=unit(2,"mm")),size=.75)+
  geom_segment(aes(x=x-.25,xend=x-.75,y=y,yend=y,col=ml),arrow=arrow(length=unit(2,"mm")),size=.75)+
  geom_segment(aes(x=x,xend=x,y=y+.25,yend=y+.75,col=mu),arrow=arrow(length=unit(2,"mm")),size=.75)+
  geom_segment(aes(x=x,xend=x,y=y-.25,yend=y-.75,col=md),arrow=arrow(length=unit(2,"mm")),size=.75)

pdf(args$outpath,useDingbats = F,width=6,height=5)
print(p)
dev.off()






