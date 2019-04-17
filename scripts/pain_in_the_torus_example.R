#toy example of Felsenstein's pain in the torus problem
library(ggplot2)
theme_set(theme_classic()+theme(axis.line=element_blank(),
                                axis.title = element_blank(),
                                axis.ticks=element_blank(),
                                axis.text.y=element_blank()))

#infinite landscape
positions <- runif(100,1,100)

ggplot(data=data.frame(positions),aes(x=positions))+
  coord_polar()+
  xlim(1,100)+
  geom_point(aes(y=0))

 #ngen <- 100

#for(gen in 1:ngen){
  offlocs <- c()
  for(i in 1:length(positions)){
    nOff <- rpois(1,1)
    if(nOff){
      for(j in 1:nOff){
        newpos <- positions[i]+rnorm(1,0,1)
        if(newpos > 100) newpos-100
        if(newpos < 1) newpos+100
        offlocs <- append(offlocs,newpos)
      }
    }
  }
  positions <- offlocs
  
  ggplot(data=data.frame(positions),aes(x=positions))+
    #ggtitle(gen)+
    xlim(1,100)+
    coord_polar()+
    geom_point(aes(y=0))
#  
#  png(paste0("~/Desktop/spaceviz/torus_pain_ex/",gen,".png"),width=4,height=4,unit="in",res=400)
#  print(p)
#  dev.off()
#}
system("cd ~/Desktop/spaceviz/;
        convert -delay 20 torus_pain_ex/{1..99}.png -delay 200 torus_pain_ex/100.png torus_pain_ex.gif")


#finite landscape with reflected offspring positions
positions <- runif(100,1,100)

ggplot(data=data.frame(positions),aes(x=positions))+
  xlim(1,100)+
  geom_point(aes(y=0))+
  geom_density()

ngen <- 100
for(gen in 1:ngen){
  offlocs <- c()
  for(i in 1:length(positions)){
    nOff <- rpois(1,1)
    if(nOff){
      for(j in 1:nOff){
        newpos <- positions[i]+rnorm(1,0,1)
        if(newpos > 100) newpos <- 100-(newpos-100)
        if(newpos < 1) newpos <- abs(newpos)
        offlocs <- append(offlocs,newpos)
      }
    }
  }
  positions <- offlocs
  
  p <- ggplot(data=data.frame(positions),aes(x=positions))+
    ggtitle(gen)+
    xlim(1,100)+
    geom_point(aes(y=0))+
    geom_density()
  
  png(paste0("~/Desktop/spaceviz/torus_pain_ex/",gen,".png"),width=4,height=4,unit="in",res=400)
  print(p)
  dev.off()
}
system("cd ~/Desktop/spaceviz/;
        convert -delay 20 torus_pain_ex/{1..99}.png -delay 200 torus_pain_ex/100.png torus_pain_ex_finite.gif")
