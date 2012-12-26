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

#-----------------------------------------  -------------------------------------------

InterMediateCox<-function(Genei,ProgFact,SurvTime,Censor,index){
      
        
    if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],Genei=Genei[index])
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ Genei,data=cdata)
    }
    if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],Genei=Genei[index],ProgFact[index,])
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ Genei",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
         } else {
         
         stop(" Argument 'ProgFact' is NOT a data frame ")
         }

    }
    
return(list(m0=m0,Genei=Genei,cdata=cdata))
}

#-----------------------------------------  -------------------------------------------


CVGeneSpecificCoxPh<-function(fold=3,
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL,
                        n.cv=3)
{
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
gNames<-rownames(ReduGdata)  

            n.train<-(n.patients-floor(n.patients/fold))
            n.test<-floor(n.patients/fold) 
            ind.train <-matrix(0,n.cv,n.train)
            ind.test  <-matrix(0,n.cv,n.test)
            
            #HR--------
            HRT<-array(NA,dim=c(n.genes,4,n.cv))
            HRTE<-array(NA,dim=c(n.genes,4,n.cv))
            
   set.seed(123)
   message("Cross Valdiations are being runing ...")
   pIndex <- c(1:n.patients)
                            
            for (k in 1:n.cv){
                ind.train[k,] <-sort(sample(pIndex,n.train,replace=F) )
                ind.test[k,] <-c(1:n.patients)[-c(intersect(ind.train[k,] ,c(1:n.patients)))] 
            }
    
options(warn=-1)
 for (j in 1:n.cv){ #--------------------------loop over CV ---------------------------------------------------
    
    
    for (i in 1:n.genes){  #---------------------------  STRAT FOR LOOP ------------------------ 
                
                #training set ---------------------------------------------------------
                TrainTemp<-InterMediateCox(ReduGdata[i,ind.train[j,]],ProgFact,SurvTime,Censor,ind.train[j,])
                
                m1 <- TrainTemp$m0
                TrtandGene<-summary(m1)[[7]][c("Genei"),1]
                rm(m1)
                p1     <- TrtandGene[1]*ReduGdata[i,ind.train[j,]]
                p1.test<- TrtandGene[1]*ReduGdata[i,ind.test[j,]]              
               
               TempGenei <-EstimHR(p1,Sdata=data.frame(SurvTime=SurvTime[ind.train[j,]],Censor=Censor[ind.train[j,]]), ProgFact=ProgFact[ind.train[j,],],Plots = FALSE, MedianCut = MedianCut )  
               HRT[i,,j]<-summary(TempGenei$SurFit)[[8]][1,]
                  
                 #testing set ---------------------------------------------------------
                
               TempGeneiTE <-EstimHR(p1.test,Sdata=data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]]), ProgFact=ProgFact[ind.test[j,],],Plots = FALSE, MedianCut = median(p1) )  
               HRTE[i,,j]<-summary(TempGeneiTE$SurFit)[[8]][1,]
                           
              
        }#---------------------------  END OF  FOR LOOP over genes ------------------------
    #print(j)
    
    }#--------------------------END OF loop over CV ---------------------------------------------------

return(new("CVGbyG",HRT=HRT,HRTE=HRTE,ind.train=ind.train,ind.test=ind.test,n.genes=n.genes,n.cv=n.cv, ReduGdata=ReduGdata))
} #-------------------------------  END OF
