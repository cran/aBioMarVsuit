####**********************************************************************
####**********************************************************************
####
####  cross validated gene by gene analysis
####
####  Copyright 2012, CenStat, Uhasselt
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Pushpike J Thilakarathne, PhD.
####    Interuniversity Institute for Biostatistics and statistical Bioinformatics (I-BioStat) 
####    University of Hasselt and Catholic university of Leuven
####    3000 Leuven
####    Belgium
####
####    email:  pushpike@gmail.com
####    URL:    http://www.uhasselt.be/fiche?voornaam=Pushpike&naam=Thilakarathne
####
####    
####**********************************************************************
####**********************************************************************

########################################################################


GridAnalysis<-function(SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150, ProgFact=NULL, Plots = FALSE,DimMethod=c("PLS","PCA")){

  DimMethod <- match.arg(DimMethod)
  
if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
if (missing(Gdata)) stop("Argument 'Gdata' is missing...")
if (missing(Censor)) stop("Argument 'Censor' is missing...")


if (ReduceDim) {
    DataForReduction<-list(x=Gdata,y=SurvTime, censoring.status=Censor, featurenames=rownames(Gdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:NuFeToBeSel]
    TentativeList
    
    ReduGdata<-Gdata[TentativeList,]
    } else {
    ReduGdata<-Gdata
}


n.genes<-nrow(ReduGdata)
n.patients<-ncol(ReduGdata)

Object<-SurFitPcaClasif(SurvTime,ReduGdata,Censor,ReduceDim=FALSE, NuFeToBeSel=150, ProgFact=ProgFact, Plots = FALSE, MedianCut = NULL )
cut.off<-quantile(Object$p1, probs = seq(0.10, 0.9, 0.05))
grid.HRGSpos121<-matrix(NA,ncol=5,nrow=length(cut.off))

for (i in 1:length(cut.off)) {
            
       if (DimMethod=="PLS") {
       Temp<- SurFitPlsClasif(SurvTime,ReduGdata,Censor,ReduceDim=FALSE,NuFeToBeSel=150,ProgFact=ProgFact, Plots = FALSE,  MedianCut =  cut.off[i])
       } else {
       Temp<-  SurFitPcaClasif(SurvTime,ReduGdata,Censor,ReduceDim=FALSE, NuFeToBeSel=150, ProgFact=ProgFact, Plots = FALSE, MedianCut = cut.off[i] )
       } 
        grid.HRGSpos121[i,]<-c(cut.off[i],summary(Temp$SurFit)[[8]][1,] )
}


if (Plots) { plotCI(x=grid.HRGSpos121[,2],xaxt="n", ui=grid.HRGSpos121[,5],li=grid.HRGSpos121[,4], 
            col="black", barcol="blue", lwd=1,ylab="HR",xlab="CutOFF (Percentiles)",main=paste("Dimension reduction :",DimMethod,sep=""))
            axis(side=1,at=1:length(cut.off),seq(0.10, 0.9, 0.05)*100,cex=0.7)
            }
data1<-grid.HRGSpos121[,-3]
colnames(data1)<-c("CutOff","EstimatedHR","LowerCI","UpperCI")
return(data1)
}
