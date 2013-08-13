# let's look at those results...

library(reshape2)
library(ggplot2)
library(RColorBrewer)


type<-2

plot.it <- TRUE
#plot.it <- FALSE

for(samp.size in c(5000,500)){
cat(samp.size,"\n")
  for (ind in 1:3){
cat(ind,"\n")

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
                      method=c(rep("VARPROP",100),
                               rep("GAMU",100),
                               rep("MBB",100),
                               rep("MBB+SDU",100)))

    # plot type
    if(type==1){
      # plot in order
      cis <- cbind(cis, y=c(t(matrix(1:400,4,100))))
    }else{
      # plot out of order
      cis <- cbind(cis, y=1:400)
    }

    truth <- samp.size

    # how many times was truth not covered by the CI
    bigger <- rep(truth,4) > cis$upper
    smaller <- rep(truth,4) < cis$lower

    cat("varprop ",sum(bigger[seq(1,400,by=4)]|smaller[seq(1,400,by=4)]),
        "\n",sep="")
    cat("gamu ",sum(bigger[seq(2,400,by=4)]|smaller[seq(2,400,by=4)]),
        "\n",sep="")
    cat("mbb ",sum(bigger[seq(3,400,by=4)]|smaller[seq(3,400,by=4)]),
        "\n",sep="")
    cat("mbb+du ",sum(bigger[seq(4,400,by=4)]|smaller[seq(4,400,by=4)]),
        "\n",sep="")














    if(plot.it){

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

    }# end plotting





  }
cat("\n\n")
}



