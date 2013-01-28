sim.data<-function(mydens,n.groups,n.transects){

   ## take a sample using line transects
   # at the moment no groups
   lt.poppars<-setpars.population(density.pop=mydens, number.groups=n.groups,
                                  size.method="poisson",size.min=1,size.max=10,
                                  size.mean=3, exposure.method="beta",
                                  exposure.min=0,exposure.max=1,
                                  exposure.mean=0.8, exposure.shape=0.5)
   lt.pop<-generate.population(lt.poppars)
   lt.despars<-setpars.design.lt(myreg,n.transects=n.transects,
                                 n.units=n.transects, visual.range=1,
                                 percent.on.effort=100)
   lt.des<-generate.design.lt(lt.despars)#, seed=3)
   lt.survpars<-setpars.survey.lt(lt.pop, lt.des, disthalf.min=0.15,
                                  disthalf.max=0.5)
   lt.samp<-generate.sample.lt(lt.survpars)
   #summary(lt.samp)
   #plot.density.sample.3d(mydens,lt.samp, scale.fact=1.5)

   ## create a data.frame with the information that we'd like
   ## to make the spatial model

   # 1 == detected
   # 0 == undetected, but in covered region
   # NA == out of covered region
   ind<-as.logical(lt.samp$detected)&!is.na(lt.samp$detected)

  cat(sum(ind), "detections\n")

   lt.spatdat<-data.frame(id=c(1:length(lt.samp$distance))[ind],
                          distance=lt.samp$distance[ind],
                          x=lt.samp$population$posx[ind],
                          y=lt.samp$population$posy[ind],
                          transect=lt.samp$transect[ind]
                         )
   transects<-data.frame(id=1:lt.samp$design$n.units,
                         start.x=lt.samp$design$pos.x,
                         end.x=lt.samp$design$pos.x,   # same here
                         start.y=lt.samp$design$start.y,
                         end.y=lt.samp$design$end.y,
                         width=rep(lt.samp$design$visual.range,
                                   lt.samp$design$n.units)
                        )

   return(list(lt.spatdat=lt.spatdat,transects=transects,lt.samp=lt.samp))
}

