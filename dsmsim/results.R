# let's look at those results...

library(reshape2)
library(ggplot2)
library(RColorBrewer)


type<-1#2

#samp.size <- 5000
#ind<-1


for(samp.size in c(500,5000)){
  for (ind in 1:3){

    # load the CVs and variances
    cvs <- read.csv(paste("dsmres-",ind,"-",samp.size,".csv",sep=""),
                    header=TRUE)
    cvs <- cvs[,-c(1,5,6)]
    names(cvs) <- c("varmethod","CV","var")

    # point estimates
    pes <- as.numeric(cvs[,2][cvs$varmethod=="pointest"])
    cvs <- cvs[cvs$varmethod!="pointest",]

    # using asymptotic CI
    asymp.ci <- exp(1.96*sqrt(log(1+as.numeric(cvs$CV))))
    preds<-pes
    cis <- data.frame(
                      lower=c(preds/asymp.ci[seq(1,400,by=4)],
                              preds/asymp.ci[seq(2,400,by=4)],
                              preds/asymp.ci[seq(3,400,by=4)],
                              preds/asymp.ci[seq(4,400,by=4)]),
                      upper=c(preds*asymp.ci[seq(1,400,by=4)],
                              preds*asymp.ci[seq(2,400,by=4)],
                              preds*asymp.ci[seq(3,400,by=4)],
                              preds*asymp.ci[seq(4,400,by=4)]),
                      method=c(rep("varprop",100),
                               rep("GAM",100),
                               rep("MBB",100),
                               rep("MBB (DU)",100)))

    # plot type
    if(type==1){
      # plot in order
      cis <- cbind(cis, y=c(t(matrix(1:400,4,100))))
    }else{
      # plot out of order
      cis <- cbind(cis, y=1:400)
    }

    truth <- samp.size

    # is truth bigger than lower? ...
    bigger <- rep(truth,4) > cis$lower
    # and smaller than upper?
    smaller <- rep(truth,4) < cis$upper
    # > cvs$varmethod[bigger&smaller]
    # [1] movblk.du movblk.du movblk.du vargam    movblk    movblk.du
    # Levels: movblk movblk.du vargam varprop


    p <- ggplot(cis)
    # cis
    p <- p + geom_segment(aes(x=lower,xend=upper,y=y,yend=y,col=method))
    # truth
    p <- p + geom_line(aes(x=rep(samp.size,2),
                           y=c(-10,410)))
    # point estimates
    if(type==1){
      pred.df <- data.frame(y=seq(1,400,4),t=preds)
      p <- p + geom_linerange(aes(x=t,ymin=y,ymax=y+4),
                            col="blue",data=pred.df)
    }else{
      pred.df <- data.frame(y=seq(1,400,1),t=c(t(matrix(preds,ncol=4))))
      p <- p + geom_point(aes(x=t,y=y),size=0.75,
                            col="blue",data=pred.df)
    }
    # make it look nice
    p <- p + coord_cartesian(ylim=c(1,400))
    p <- p + coord_cartesian()
    p <- p + theme(panel.grid.major=element_blank(),
                   panel.background=element_blank(),
                   axis.ticks.y=element_blank(),
                   legend.key=element_blank(),
                   axis.text.y=element_blank(),
                   panel.grid.minor=element_blank())
    p <- p + labs(x="Abundance",y="Simulation",col="Method")
    p <- p + scale_colour_manual(values=brewer.pal(name="Dark2",n=4))
    p <- p + ggtitle(paste("n=",samp.size,",","sim=",ind,sep=""))

    print(p)

    ggsave(paste("res-",samp.size,"-",ind,"-",type,".pdf",sep=""))

  }
}



