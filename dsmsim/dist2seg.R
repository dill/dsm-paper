# which observations are in which segments?
dist2seg<-function(lt.spatdat,segments){

   # observations
   xydat<-data.frame(id=lt.spatdat$id,
                     x=lt.spatdat$x,
                     y=lt.spatdat$y,
                     distance=lt.spatdat$distance)

   # storage
   segobj<-data.frame(object=NA,distance=NA,Sample.Label=NA,size=NA)

   for(i in 1:nrow(segments)){
      # were there any observations in this segment?
      inout<-soap:::inSide(list(x=c(segments[i,1],segments[i,3],segments[i,5], 
                                     segments[i,7],segments[i,1]),
                                 y=c(segments[i,2],segments[i,4],segments[i,6],
                                     segments[i,8],segments[i,2])),
                            x=xydat$x,y=xydat$y)

      # if there were, attribute that segment to them
      if(sum(inout)>0){
         ids<-xydat$id[inout]
         segobj<-rbind(segobj,
                       cbind(object=ids,
                             distance=xydat$distance[inout],
                             Sample.Label=rep(segments$Sample.Label[i],sum(inout)),
                             size=1))
      }

      # leave only those observations not yet matched in xydat
      xydat<-xydat[!inout,]
   }

   # take out the first line of NAs
   segobj<-segobj[-1,]

   return(segobj)
}
