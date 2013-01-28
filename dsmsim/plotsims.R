# plot simulation scenarios


library(mgcv)
library(mrds)
library(wisp)
library(dsm)

set.seed(1011)

# generate the region
myreg <- generate.region(x.length =100, y.width = 50)

density.plots <- list()

# 1 - two humps
density.plots[[1]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 1, southeast = 1, northwest = 1)
# two splodges one at either side
density.plots[[1]] <- add.hotspot(density.plots[[1]], 10,25,40,15)
density.plots[[1]] <- add.hotspot(density.plots[[1]], 90,25,40,15)

# 2 - gradient with transect length
density.plots[[2]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 100, southeast = 1, northwest = 100)

# 3 - gradient perpendicular transect length
density.plots[[3]] <- generate.density(myreg,nint.x = 100, nint.y = 50, southwest = 100, southeast = 100, northwest = 1)



pdf("sim-plots.pdf",width=12,height=2.8)

par(mfrow=c(1,3))
for(i in 1:length(density.plots)){
   image(density.plots[[i]]$matrix,asp=1,las=1,
         x=seq(0,100,1),y=seq(0,50,1),xlab="x",ylab="y")
   contour(z=density.plots[[i]]$matrix,x=seq(1,100,1),y=seq(1,50,1),add=T)
}

dev.off()
