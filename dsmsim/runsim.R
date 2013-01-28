run.sim<-function(mydens,n.groups,n.sims,pred.data,preddata.varprop,n.transects,n.boot){
   bigres<-c()

   for(i in 1:n.sims){
      simdat<-sim.data(mydens,n.groups,n.transects)
      lt.spatdat<-simdat$lt.spatdat
      transects<-simdat$transects
      lt.samp<-simdat$lt.samp

      # make the segment data
      segments<-line2seg(lt.spatdat,transects)

      #### PLOTTING
      #par(mfrow=c(1,2))
      #image(z=mydens$matrix,x=seq(0,mydens$reg$length,len=mydens$n.interval.x),
      #                      y=seq(0,mydens$reg$width,len=mydens$n.interval.y),
      #      xlab="x",ylab="y",
      #      zlim=c(0,100),col=heat.colors(1000),main="True surface",asp=1)
      #contour(z=mydens$matrix,
      #        x=seq(0,mydens$reg$length,len=mydens$n.interval.x),
      #        y=seq(0,mydens$reg$width,len=mydens$n.interval.y),
      #        zlim=c(0,100),col="green",add=TRUE)
      #### /PLOTTING

      # fit detection function 
      dist.dat<-lt.spatdat
      names(dist.dat)[1]<-"object"

      dist<-ddf(dsmodel=~mcds(key="hn",formula="~1"),data=dist.dat,
                   method="ds",meta.data=list(width=transects$width[1]))

      # setup segment data
      segdata <- data.frame(x = segments$x,
                            y = segments$y,
                            Effort = segments$Effort ,
                            Transect.Label = segments$transect.id,
                            Sample.Label = segments$Sample.Label)

      # setup obsdata
      obsdata <- dist2seg(lt.spatdat,segments)

      # fit the density surface
      dens <- dsm(N~s(x,y),dist,segdata,obsdata)

      #### PLOTTING
      #vis.gam(dens$result,plot.type="contour",asp=1,type="response",
      #        color="heat",nCol=1000,main="Estimated surface")
      #points(lt.spatdat$x,lt.spatdat$y,pch=19,cex=0.5)
      #### /PLOTTING

      off.set <- 1.020408*1.010101 #*predict(dist)$fitted[1]

      #cat("True abundance =",sum(lt.samp$detected==1,na.rm=T)+
      #                       sum(is.na(lt.samp$detected)),"\n")
      #cat("mrds N in covered region =",dist$Nhat,"\n")
      #cat("Volume under density surface =",sum(dsm.predict(dens,pred.data,off=off.set)))
      #cat("\n")

      # MVB varprop
      var.prop <- dsm.var.prop(dens,pred.data,off.set)
      res.varprop <- summary(var.prop)
      #print(summary(var.prop))

      # GAM-style
      var.gam <- dsm.var.gam(dens,pred.data,off.set)
      res.vargam <- summary(var.gam)
      #print(summary(var.gam))

      # moving block
      var.movblk <- dsm.var.movblk(dens, pred.data, n.boot=n.boot,
                              block.size=5,
                              off.set=off.set,bar=FALSE,
                              ds.uncertainty=FALSE)
      res.movblk <- summary(var.movblk)

      var.movblk.du <- dsm.var.movblk(dens, pred.data, n.boot=n.boot,
                              block.size=5,
                              off.set=off.set,bar=FALSE,
                              ds.uncertainty=TRUE)
      res.movblk.du <- summary(var.movblk.du)


      ## order: true, mrds, gam
      #res<-c(sum(lt.samp$detected==1,na.rm=T)+sum(is.na(lt.samp$detected)),
      #       dist$Nhat,
      #       sum(dsm.predict(dens,pred.data,off=off.set)))

      cat("Predicted abundance :",
                    sum(predict(dens,pred.data,off=off.set)),"\n")
      cat("True abundance :",sum(lt.samp$detected==1,na.rm=T)+
                             sum(is.na(lt.samp$detected)),"\n")

      bigres <- rbind(bigres,
              c("varprop",res.varprop$cv,res.varprop$se),
              c("vargam",res.vargam$cv,res.vargam$se),
              c("movblk",res.movblk$cv,res.movblk$se),
              c("movblk.du",res.movblk.du$cv,res.movblk.du$se),
              c("pointest",sum(predict(dens,pred.data,off=off.set)),NA))

      write.table(rbind(
                    c("varprop",res.varprop$cv,res.varprop$se),
                    c("vargam",res.vargam$cv,res.vargam$se),
                    c("movblk",res.movblk$cv,res.movblk$se),
                    c("movblk.du",res.movblk.du$cv,res.movblk.du$se),
                    c("pointest",sum(predict(dens,pred.data,off=off.set)),NA)),
                 file="working.csv",append=TRUE,sep=",")

   }
   return(bigres)
}
