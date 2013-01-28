## count in each segment
countsegs<-function(lt.spatdat,segments){
   counts<-double(nrow(segments))
   for(i in 1:nrow(segments)){
      counts[i]<-sum(inSide(list(x=c(segments[i,1], segments[i,3], segments[i,5], 
                                     segments[i,7],segments[i,1]),
                                 y=c(segments[i,2], segments[i,4], segments[i,6],
                                     segments[i,8],segments[i,2])),
                            x=lt.spatdat$x,y=lt.spatdat$y))
   }

   return(counts)
}
