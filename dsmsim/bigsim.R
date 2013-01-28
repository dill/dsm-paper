#set.seed(100)

library(mgcv)
library(soap) #for inSide
library(mrds)
library(wisp)
library(dsm)

source("line2seg.R")
source("dist2seg.R")
source("countsegs.R")
source("sim.data.R")
source("runsim.R")

set.seed(1011)

# generate the region
myreg <- generate.region(x.length =100, y.width = 50)

# some settings go here
samp.sizes<-500#c(100,250,500,1000)
n.sims<-100
density.plots<-list()
n.transects<- 14
n.boot <- 200

# 1 - two humps
density.plots[[1]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 1, southeast = 1, northwest = 1)
# two splodges one at either side
density.plots[[1]] <- add.hotspot(density.plots[[1]], 10,25,40,15)
density.plots[[1]] <- add.hotspot(density.plots[[1]], 90,25,40,15)

# 2 - gradient with transect length
density.plots[[2]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 100, southeast = 1, northwest = 100)

# 3 - gradient perpendicular transect length
density.plots[[3]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 100, southeast = 100, northwest = 1)

## 4 - lots of splodges
#density.plots[[4]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 1, southeast = 1, northwest = 1)
#for(i in 1:50){
#   x<-runif(1,0,100)
#   y<-runif(1,0,50)
#   density.plots[[4]] <- add.hotspot(density.plots[[4]], x,y,10,3)
#}


#par(mfrow=c(2,2))
#for(i in 1:length(density.plots)){
#   image(density.plots[[i]]$matrix,asp=1)
#   contour(density.plots[[i]]$matrix,add=T)
#}


for(samp.size in samp.sizes){
  i <- 1
  for(density.plot in density.plots){

    big.res<-c()
    # prediction grid
    pred.data<-expand.grid(x=seq(0,density.plot$reg$length,
                                 len=density.plot$n.interval.x),
                           y=seq(0,density.plot$reg$width,
                                 len=density.plot$n.interval.y))

    this.res<-run.sim(density.plot,samp.size,n.sims,
                      pred.data,preddata.varprop,n.transects,n.boot)

    big.res<-rbind(big.res,
                   cbind(this.res,
                         rep(samp.size,nrow(this.res)),
                         rep(i,nrow(this.res))))

    write.csv(big.res,file=paste("dsmres-",i,"-",samp.size,".csv",sep=""))
    i<-i+1
  }
}

