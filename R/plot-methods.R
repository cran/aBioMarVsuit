
##-------
# Pushpike Thikarathne
# I-BioStat, University Of Hasselt, Belgium.

setMethod("plot", signature(x="FurtherValidation", y="missing"),
  function(x,  y, ptype=1, ...) {
 if (class(x)!="FurtherValidation") stop("Invalid class object")


if (ptype==1) {
    dotsCall <- substitute(list(...))
    ll <- eval(dotsCall)
    if(!hasArg("xlab")) ll$xlab <- "Outer Iteration"
    if(!hasArg("ylab")) ll$ylab <- "HR estimate"
    ll$main <- "Distribution of HR on Test Data \n for low risk group"
    if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
    if(!hasArg("cex.main")) ll$cex.main <- 1
    if(!hasArg("col")) ll$col <- 3
    ll$x<-sapply(1:x@n.cv,function(j) x@HRInner[,1,j])
    do.call(boxplot,args=ll)
  points(1: x@n.cv,x@HRTE[,1],pch=19,col="red",cex = 1.5)
  for (i in 1: x@n.cv) lines(c(i+0.05,i+0.05),
  c(x@HRTE[i,2],x@HRTE[i,3]),col="red")

}

if (ptype==2) {

 dotsCall <- substitute(list(...))
    ll <- eval(dotsCall)
    if(!hasArg("xlab")) ll$xlab <- "Estimated HR"
    if(!hasArg("ylab")) ll$ylab <- ""
    ll$main <- "Estimated HR Density for low Risk Group "
    if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
    if(!hasArg("cex.main")) ll$cex.main <- 1
    if(!hasArg("col")) ll$col <- 2
    
        ll$x <- density(x@HRTE[,1]) 
        
         do.call(plot,args=ll)
        qq<-quantile(sort(x@HRTE[,1]),prob=c(0.05,0.5,0.95))
        abline(v=qq[1],col=3)
        abline(v=qq[2],col=3)
        abline(v=median(x@HRTE[,1]),col=3,lwd=5)     
}

})


#-----------------------------------------------------------------

setMethod("plot", signature(x="cvelasticnetcox", y="missing"),
  function(x,  y, ptype=1, ...) {
 if (class(x)!="cvelasticnetcox") stop("Invalid class object")


if (ptype==1) {

    DistHR<-data.frame(HRTest=x@HRTE[,1],HRTrain=x@HRT[,1])
    
    colnames(DistHR)<-c("Test","Train")
    dotsCall <- substitute(list(...))
    ll <- eval(dotsCall)
    if(!hasArg("xlab")) ll$xlab <- ""
    if(!hasArg("ylab")) ll$ylab <- "HR estimate"
    ll$main <- "Distribution of HR on Test and Train Data \n for low risk group"
    if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
    if(!hasArg("cex.main")) ll$cex.main <- 1
    if(!hasArg("col")) ll$col <- 2:3
    ll$x<-DistHR
   do.call(boxplot,args=ll)

}


if (ptype==2) {

    HRTest=x@HRTE[,1]
    dotsCall <- substitute(list(...))
    ll <- eval(dotsCall)
    if(!hasArg("xlab")) ll$xlab <- "Estimated HR on Test Data"
    if(!hasArg("ylab")) ll$ylab <- "Number of non zero coef."
    ll$main <- "HR vs number of genes"
    if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
    if(!hasArg("cex.main")) ll$cex.main <- 1
    if(!hasArg("col")) ll$col <- 2
    ll$x<-HRTest
    ll$y<-x@n.g
   do.call(plot,args=ll)

}


if (ptype==3) {

    Freq=colSums(x@gene.mat)
    names(Freq)<-rownames(x@ReduGdata)
    sFreq<-sort(Freq,decreasing = TRUE)
   sFreq<-sFreq[sFreq>0]
   
   maxG<-length(sFreq)
   if (maxG>30) maxG<-30
  
  
    dotsCall <- substitute(list(...))
    lll <- eval(dotsCall)

         lll$height<-sFreq[1:maxG]
      if(!hasArg("xlab" ))   lll$xlab<-""
      if(!hasArg("ylab" ))   lll$ylab<-"Frequency"
      if(!hasArg("main"))    lll$main<-"Mostly Selected Genes"
          lll$col<-rainbow(maxG)
          lll$names<-names(sFreq)[1:maxG]
          if(!hasArg("cex.lab")) lll$cex.lab <- 1
       if(!hasArg("las"))   lll$las<-2
       lll$cex.names<-0.65
          do.call(barplot,args=lll) 

}



}
)
#-----------------------------------------------------------------

setMethod("plot", signature(x="Permutation", y="missing"),
  function(x,  y, ...) {
 if (class(x)!="Permutation") stop("Invalid class object")

HR<-x@HRlowPerm[,1]
HR<-na.exclude(HR)
n<-x@n.perm
vv=x@HRlowObs[1]

 pvalue<-sum(vv>HR)/n


dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- paste("Estimated HR \n Emperical p-value: ",round(pvalue,4),sep="")
           if(!hasArg("ylab")) ll$ylab <- ""
           ll$main <- "Null Distribution of HR on Permuted Data \n for low risk group"
           if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
           if(!hasArg("cex.main")) ll$cex.main <- 1
          if(!hasArg("col")) ll$col <- 1
            if(!hasArg("ylim")) ll$ylim <- c(0,4)
            
        ll$x <- density(HR,from=0,to=(max(HR)+0.25)) 
        do.call(plot,args=ll)
        abline(v=vv,col=2)
        #CI for permuated cases
        qq<-quantile(sort(HR),prob=c(0.05,0.95))
        abline(v=qq[1],col=3)
        abline(v=qq[2],col=3)
        abline(v=median(HR),col=3,lwd=3)

})








#-----------------------------------------------------------------

setMethod("plot", signature(x="CVMV", y="missing"),
  function(x,  y, ...) {
 if (class(x)!="CVMV") stop("Invalid class object")
HRp.test<-x@HRp.test
HRp.train<-x@HRp.train
 nCV<-x@nCV
 dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "HR" 
           if(!hasArg("ylab")) ll$ylab <- "Cross Validation index"
           ll$main <- "Estimated HR on Test Data \n for low risk group"
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
           if(!hasArg("cex.main")) ll$cex.main <- 1.3
          if(!hasArg("col")) ll$col <- 2
           
           ll$x<-HRp.test[,1]
  if(!hasArg("ylim")) ll$ylim <- c(0,2)
  
  
par(mfrow=c(1,2))
t1 <- which(HRp.test[,1]<1)
do.call(plot,args=ll)
#plot(HRp.test[,1],ylim=c(0,2),ylab="HR",main="")
for(i in 1:nCV){
lines(c(i,i),HRp.test[i,2:3])
}
for(i in t1){
lines(c(i,i),HRp.test[i,2:3],col=2)
}
abline(h=1)


Results<-data.frame(HRpTrain=HRp.train[,1],HRpTest=as.numeric(HRp.test[,1]))
ll$x<-Results
ll$names<-c("Train ","Test ")
ll$main <- "Estimated HR on Train and Test Data \n for low risk group"
 if(!hasArg("col")) ll$col <- 2:3
do.call(boxplot,args=ll)

#boxplot(Results,names=c("Train ","Test "),main="",ylab="HR",col=c("green","red"))


})








#-----------------------------------------------------


setMethod("plot", signature(x="CVSeqInc", y="missing"),
  function(x,  y, ...) {
  
par(mfrow=c(1,2))
nn<-dim(x@HRPC)[3]
PC.HRp<-x@HRPC[,1,1:nn]
colnames(PC.HRp)<-x@top

PL.HRp<-x@HRPL[,1,1:nn]
colnames(PL.HRp)<-x@top
 
 dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "Top K Genes" 
           if(!hasArg("ylab")) ll$ylab <- "Cross Validated HR"
           ll$main <- "Estimated HR on Test Data \n for Top K Genes (PCA)"
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
           if(!hasArg("cex.main")) ll$cex.main <- 1.3
          if(!hasArg("col")) ll$col <- 1:nn
            if(!hasArg("ylim")) ll$ylim <- c(0,5)
           ll$x<-PC.HRp
           do.call(boxplot,args=ll)

#boxplot(PC.HRp,main="HR (PC1)",xlab="Top K Genes",ylab="HR",col=1:nn)

 ll$x<-PL.HRp
 ll$main <- "Estimated HR Test Data \n for Top K Genes (PLS)"
 do.call(boxplot,args=ll)
#boxplot(PL.HRp,main="HR (PLS1)",xlab="Top K Genes",ylab="HR",col=1:nn)

return(invisible())
}

)

#-----------------------------------------------------


setMethod("plot", signature(x="CVGbyG", y="missing"),
  function(x,  y, Which=1, ...) {
  HRTEST  <- x@HRTE[Which,,][1,]
  HRTrain <- x@HRT[Which,,][1,]
  Results<-data.frame(HRTrain,HRTEST)
  dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "" 
           if(!hasArg("ylab")) ll$ylab <- "Cross Validated HR"
           if(!hasArg("main")) ll$main <- paste("Estimated HR of low risk for Gene ", Which, "\n Number of CVs : ",x@n.cv,sep="")
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
           if(!hasArg("cex.main")) ll$cex.main <- 1.3
           if(!hasArg("ylim")) ll$ylim <- c(0,15)
           if(!hasArg("col")) ll$col <- 2:3
           if(!hasArg("names"))  ll$names=c("Train","Test")
           ll$x<-Results
           do.call(boxplot,args=ll)
return(invisible())
}

)


#-----------------------------------------------------


setMethod("plot", signature(x="CVfordimMethods", y="missing"),
  function(x,  y, ...) {
  dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "Risk groups" 
           if(!hasArg("ylab")) ll$ylab <- "Cross Validated HR"
           if(!hasArg("main")) ll$main <- paste("HR for Low risk group\n by ", x@PCAorPLS,sep="")
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
           if(!hasArg("ylim")) ll$ylim <- c(0,15)
           if(!hasArg("col")) ll$col <- c(2,3)
           ll$x<-x@Results
           do.call(boxplot,args=ll)
  #boxplot(x@Results,ylim=c(0,max(x@Results,na.rm=T)),names=c("Train ","Test"),main=mtitle,ylab="HR",col=c("green","red"))
    return(invisible())
  }
)


#-----------------------------------------------------


setMethod("plot", signature(x="GeneSpecific", y="missing"),
  function(x,  y, ...) {
  dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "" 
           if(!hasArg("ylab")) ll$ylab <- "Estimated HR"
           if(!hasArg("main")) ll$main <- "Distribution of \n Gene-Specific HR for low risk"
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
           if(!hasArg("ylim")) ll$ylim <- c(0,1.5)
           if(!hasArg("col")) ll$col <- c(2,3)
           ll$x<-x@HRp[,1]
           
           par(mfrow=c(1,2))
           do.call(boxplot,args=ll)
  #boxplot(x@HRp[,1],ylim=c(0,max(x@HRp[,1],na.rm=T)),names=c("HR"),main=mtitle,ylab="HR",col=c("green"))
  #CountHighRisk<-sapply(1:ncol(x@gr), function(k) sum(x@gr[,k]==1))
          
          lll<-list()
          lll$height<-sapply(1:ncol(x@gr), function(k) sum(x@gr[,k]==1)/nrow(x@gr))
          lll$xlab<-"Patient Index"
          lll$ylab<-"Proportion"
          lll$main<-"Proportion of being \n Classified as Low Risk"
          lll$col<-rainbow(ncol(x@gr))
          lll$names<-1:ncol(x@gr)
          do.call(barplot,args=lll) 
    return(invisible())
  }
)

#-----------------------------------------------------
