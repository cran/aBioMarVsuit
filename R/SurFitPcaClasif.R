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


#---------------------------------------------------------------------------------------------------
#----------------------------------- function for PCA and survival and classification ---------------



SurFitPcaClasif<-function(                                                        
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        Plots = FALSE,
                        MedianCut = NULL 
                              
){

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
#n.train<-(n.patients-floor(n.patients/fold))
#n.test<-floor(n.patients/fold)
#
#pIndex <- c(1:n.patients)
#TrainIndex <-sort(sample(pIndex,n.train,replace=F) )
#TestIndex  <-c(1:n.patients)[-c(intersect(TrainIndex,c(1:n.patients)))] 
#
#dataTrain<-list(x=Gdata[,TrainIndex],y=SurvTime[TrainIndex], censoring.status=Censor[TrainIndex], featurenames=rownames(Gdata))
#dataTest<-list(x=Gdata[,TestIndex],y=SurvTime[TestIndex], censoring.status=Censor[TestIndex], featurenames= rownames(Gdata))
#
#a<- superpc.train(dataTrain, type="survival")$feature.scores
#
#cv.obj<-superpc.cv(a, dataTrain, max.features=nrow(Gdata), n.components = 1)
#
#TempThre<-cv.obj$threshold[which.max(cv.obj$scor)]
#
##superpc.plotcv(cv.obj)
#
#
#fit.red<- superpc.predict.red(a,dataTrain, dataTest,threshold=2.11,n.components = 1)
#superpc.listfeatures(dataTrain, a,  fit.red, num.features=NuFeToBeSel)

    if (is.matrix(ReduGdata)) {
        pc1 <- f.pca(as.matrix(ReduGdata))[[6]][,1]
        } else {
        pc1<-ReduGdata
        }         
        
    if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime,Censor,pc1)
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ pc1,data=cdata)
    }
    if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime,Censor,pc1,ProgFact)
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
         } else {
         
         stop(" Argument 'ProgFact' is NOT a data frame ")
         }

    }
              
       #risk Score       
                TrtandPC1<-summary(m0)[[7]][c("pc1"),1]
                p1 <- TrtandPC1*pc1
              
   TempRes<-EstimHR(p1,Sdata=cdata, ProgFact=ProgFact,Plots = Plots, MedianCut = MedianCut )
        
tempp<-list(SurFit=TempRes$SurFit,p1=p1,gs=TempRes$gs,pc1=pc1)
class(tempp)<-"SurFitPcaOut"
return(tempp)
}


#---------------------------------------------------------------------------------------------------
#--------------------------------- END OF function for PC and survival and classification ----------
