# function to turn lines into segments
line2seg<-function(lt.spatdat,transects){

   width<-transects$width[1]
   # need to calculate the number of splits
   # from the lengths, keep it simple here!
   seg.len<-2*transects$width[1] # move inside loop eventually?

   # do something smart here!
   n.segs<-ceiling(abs(transects$end.y[1]-transects$start.y[1])/seg.len)
   segments<-c()
   seg.id<-1

   # split the transects
   for(i in transects$id){
      # start of the transect vector 
      cl1<-c(transects$start.x[i],transects$start.y[i])
      # unit vector in the direction we want to go
      cl2<-(c(transects$end.x[i],transects$end.y[i])-
            c(transects$start.x[i],transects$start.y[i]))/
           sqrt((transects$end.x[i]-transects$start.x[i])^2+
                (transects$end.y[i]-transects$start.y[i])^2)

      # if we finish at the stop, we start again at the the top on the
      # next transect. So, need to take that into account otherwise
      # the segments screw up...
      if(transects$end.y[i]>transects$start.y[i]) for.i<-1:n.segs
      if(transects$end.y[i]<transects$start.y[i]) for.i<-0:(n.segs-1)

      # actually make the segments
      for(j in for.i){
         # top of the centre line for that segment
         cl<-cl1+j*seg.len*cl2 
         # tl tr br bl, centrepoint x, centrepoint y, 
         # transect id, segment within transect id,
         # overall segment label (id)
         segments<-rbind(segments,
                         c(cl+c(-width,0),
                           cl+c(width,0),
                           cl+c(width,-seg.len),
                           cl+c(-width,-seg.len),
                           cl+c(0,-seg.len/2),
                           i,j,seg.id,seg.len)
                        )
         seg.id<-seg.id+1
      }
   }

   # make a data frame of this
   segments<-as.data.frame(segments)
   names(segments)<-c("tl.x","tl.y","tr.x","tr.y",
                      "br.x","br.y","bl.x","bl.y","x","y",
                      "transect.id","tsegment.id","Sample.Label","Effort")
   return(segments)
}
