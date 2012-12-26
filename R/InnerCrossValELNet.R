####**********************************************************************
####**********************************************************************
####
####  Inner Cross Valdiation of LAsso and Elastic net for cox with gene expression
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
#----------------------------------------------------------- Begin InnerCrossValELNet
InnerCrossValELNet<-function(fold=3,n.cv=50,n.innerCV=100,MixParAlpha=0.1,
                             GexprMat,TopGenes,WeightsFixed=FALSE,
                             SurvTime, Censor,ProgFact=NULL)
{

if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
if (missing(GexprMat)) stop("Argument 'GexprMat' is missing...")
if (missing(Censor)) stop("Argument 'Censor' is missing...")

TopK<-TopGenes
if (length(TopK)!=sum(is.element(rownames(GexprMat),TopK))) stop("TopK should be a subset of GexprMat")

#TopK<-rownames(ReduGexprMat)

options(warn=-1)
t1 <- proc.time()
            n.patients<-length(SurvTime)
            n.train<-(n.patients-floor(n.patients/fold))
            n.test<-floor(n.patients/fold) 
            ind.train <-matrix(0,n.cv,n.train)
            ind.test  <-matrix(0,n.cv,n.test)
            
            HRInner<-array(NA,c(n.innerCV,1,n.cv))
            WeightMat<-matrix(NA,ncol=length(TopK),nrow=n.cv)
            
            #HR--------

            HRTE<-matrix(NA,nrow=n.cv,ncol=3)
            
            #----------record majority vote results
            MajorV<-matrix(NA,nrow=n.test,ncol=n.cv)
            

            
            #------------------------  Begin FOR LOOP :1 ---------------------------
   message("Cross Valdiations are being runing ...")
   pIndex <- c(1:n.patients)
   SurvTime[SurvTime<=0]<-quantile(SurvTime,probs =0.02)
   
            for (j in 1:n.cv){
                ind.train[j,] <-sort(sample(pIndex,n.train,replace=F) )
                ind.test[j,] <-c(1:n.patients)[-c(intersect(ind.train[j,] ,c(1:n.patients)))] 
            }



for (j in 1:n.cv){
message(paste(" ", j, " Outer Iteration ",sep=""))
       perProg<-NULL
       if (!is.null(ProgFact))     perProg<-ProgFact[ind.train[j,],]
            #inner cross valdiation     
           ResT<-CVLassoElasticNetCoxPh(n.cv=n.innerCV,fold=3,SurvTime=SurvTime[ind.train[j,]], 
                 Censor=Censor[ind.train[j,]], GexprMat[is.element(rownames(GexprMat),TopK),ind.train[j,]], 
                 NuFeToBeSel=150, ProgFact=perProg, ReduceDim=FALSE, GeneList=NULL,  StZ=TRUE, alpha=MixParAlpha)     
           
           HRInner[,1,j]<-ResT@HRT[,1]

            MeanWeightGene<-sapply(1:length(TopK),function(i) mean(ResT@coef.mat[,i]))
            names(MeanWeightGene)<-TopK
            WeightMat[j,]<-MeanWeightGene
            #end of inner cross validation
            
            #weights are NOT fixed
        
            if (!WeightsFixed) {  
                                    MajorityVotes<-matrix(NA,ncol=length(ind.test[j,]),nrow=n.innerCV)
                                    for (ii in 1:n.innerCV){
                                            WeightsInner<-ResT@coef.mat[ii,]
                                            names(WeightsInner)<-TopK
                                           # Risk score for Test set
                                            ScoreInner<- WeightsInner%*%GexprMat[is.element(rownames(GexprMat),TopK),ind.test[j,]]
                        
                                            MedianScoreInner<-median(ScoreInner)
                                            MajorityVotes[ii,]<-ifelse(ScoreInner<=MedianScoreInner,1,0)  # 1 => Responsive, 0 => non responsive
                                        }
                                    
                                    ResNres<-   sapply(1:length(ind.test[j,]),function(kk) sum(MajorityVotes[,kk]))<(n.innerCV/2) # patients should get more than 50% votes to be classified as Responsive
                                    MajorV[,j]<-ResNres

                                    #---------------------- Create indicator variable ----------------          
                                    ##HR low
                                            gid<- NULL
                                            n.patients<-n.test 
                                    gid[c(1:n.patients)[ResNres==FALSE]]<- "low"
                                    gid[c(1:n.patients)[ResNres==TRUE]]<- "high"
                                    
        
                                    
                    #extract out of bag sample-----------------------------------------------------------------------------
                           
                            if (is.null(ProgFact)) {
        
                            cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid))
                            m1 <- coxph(Surv(SurvTime, Censor==1) ~ gid,data=cdata)
                            }
                            
                            if (!is.null(ProgFact)) {
                             if (is.data.frame(ProgFact)) { 
                                 nPrgFac<-ncol(ProgFact)
                                 cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid),ProgFact[ind.test[j,],])
                                 NameProg<-colnames(ProgFact)
                                eval(parse(text=paste( "m1 <-coxph(Surv(SurvTime, Censor==1) ~ gid",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
                                 } else {
                                 
                                 stop(" Argument 'ProgFact' is NOT a data frame ")
                                 }
                    
                            }
     
                            HRTE[j,]<-(summary(m1)[[8]][1,])[-2]
           }  else {
           
           
           #-------------------------  now weights are fixed ----------------------------------------
           # Risk score for Test set
            
            RiskScore<-MeanWeightGene%*%GexprMat[is.element(rownames(GexprMat),TopK),ind.test[j,]]
            gid<- NULL
            n.patients<-n.patients<-n.test 
            gid[c(1:n.patients)[RiskScore < median(RiskScore)]]<- "low"
            gid[c(1:n.patients)[RiskScore >= median(RiskScore)]]<- "high"
            
            if (is.null(ProgFact)) {
        
                            cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid))
                            m1 <- coxph(Surv(SurvTime, Censor==1) ~ gid,data=cdata)
                            }
                            
                            if (!is.null(ProgFact)) {
                             if (is.data.frame(ProgFact)) { 
                                 nPrgFac<-ncol(ProgFact)
                                 cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid),ProgFact[ind.test[j,],])
                                 NameProg<-colnames(ProgFact)
                                eval(parse(text=paste( "m1 <-coxph(Surv(SurvTime, Censor==1) ~ gid",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
                                 } else {
                                 
                                 stop(" Argument 'ProgFact' is NOT a data frame ")
                                 }
                    
                            }
            HRTE[j,]<-(summary(m1)[[8]][1,])[-2]
   
            }
}

t2 <- proc.time()
RunTime<-as.vector(t2-t1)[3]
#----------- results -------------
return(new("FurtherValidation",RunTime=RunTime, fold=fold,  n.cv=n.cv,  n.innerCV=n.innerCV,TopK=TopK,
            HRInner=HRInner,HRTE=HRTE,WeightMat=WeightMat))
}
#----------------------------------------------------------- END OF EvaluateELNETbyInnerCV
