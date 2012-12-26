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
SeqIncreaseGenes<-function(TopK=15,SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150, ProgFact=NULL, Plots = FALSE,DimMethod=c("PLS","PCA"),...)
{
Decrease=FALSE
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

object<-GeneSpecificCoxPh( SurvTime, ReduGdata,  Censor,  ReduceDim=FALSE, NuFeToBeSel=150,    ProgFact=ProgFact, MedianCut = NULL)

Names.KGenes<-object@gNames
index.Top.KGenes <- order(object@HRp[,1],decreasing =Decrease)
index.Top.KGenes<-index.Top.KGenes[1:TopK]

TopSet<-Names.KGenes[index.Top.KGenes]

Glist.resutls.plus<-matrix(NA,TopK,4)

for (i in 1:length(TopSet[1:TopK]) ){
    
    glist<-1:i 
    
      if (DimMethod=="PLS") {
       Temp<- SurFitPlsClasif(SurvTime,ReduGdata[intersect(rownames(ReduGdata),TopSet[glist]),],Censor,
       ReduceDim=FALSE,NuFeToBeSel=150,ProgFact=ProgFact, Plots = FALSE,  MedianCut =  NULL)
       } else {
       Temp<-  SurFitPcaClasif(SurvTime,ReduGdata[intersect(rownames(ReduGdata),TopSet[glist]),],Censor,
       ReduceDim=FALSE, NuFeToBeSel=150, ProgFact=ProgFact, Plots = FALSE, MedianCut = NULL )
       } 
    
    if (is.null(ProgFact)) Glist.resutls.plus[i,]<-c(i,(summary(Temp$SurFit)[[8]][1,])[-2] )
    if (!is.null(ProgFact)) Glist.resutls.plus[i,]<-c(i,(summary(Temp$SurFit)[[8]][1,])[-2] )
    }

if (Plots) {
  dotsCall <- substitute(list(...))
  ll <- eval(dotsCall)
           if(!hasArg("xlab")) ll$xlab <- "Top K Genes" 
           if(!hasArg("ylab")) ll$ylab <- " HR for low risk"
           if(!hasArg("main")) ll$main <- paste("Estimated HR \n Based on ",DimMethod,sep="")
           if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
           if(!hasArg("cex.main")) ll$cex.main <- 1.3
           #if(!hasArg("ylim")) ll$ylim <- c(0,15)
           if(!hasArg("col")) ll$col <-2
           ll$x<-Glist.resutls.plus[,2]
           ll$xaxt<-"n"
           ll$ui=Glist.resutls.plus[,4]
           ll$li=Glist.resutls.plus[,3]
           
           do.call(plotCI,args=ll)
           axis(1);axis(2);box()
}

colnames(Glist.resutls.plus)<-c("Topk","EstimatedHR","LowerCI","UpperCI")
return(Glist.resutls.plus)
}
