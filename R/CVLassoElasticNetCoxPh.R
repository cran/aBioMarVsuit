####**********************************************************************
####**********************************************************************
####
####  cross validated LAsso and Elastic net
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


CVLassoElasticNetCoxPh<-function(n.cv=20,fold=3,SurvTime, Censor, Gdata,ReduceDim=TRUE,  NuFeToBeSel=150,  ProgFact=NULL, 
                         GeneList=NULL,  StZ=TRUE, alpha=1){
            options(warn=-1)
            
            
if (ReduceDim) {
    DataForReduction<-list(x=Gdata,y=SurvTime, censoring.status=Censor, featurenames=rownames(Gdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:NuFeToBeSel]
    TentativeList
    
    ReduGdata<-Gdata[TentativeList,]
    } else {
    ReduGdata<-Gdata
}

if (!is.null(GeneList)) ReduGdata<-Gdata[GeneList,]



n.genes<-nrow(ReduGdata)
n.patients<-ncol(ReduGdata)
Run.Time<-rep(NA, n.cv)
n.train<-(n.patients-floor(n.patients/fold))
n.test<-floor(n.patients/fold)
cv.train <-matrix(0,n.cv,n.train)
cv.test  <-matrix(0,n.cv,n.test)
           
            
#no.of genes retain
n.g<-rep(0, n.cv)

#optimum lambda
lambda<-rep(NA, n.cv)
pld<-rep(NA, n.cv)
#HR--------
HRT<-matrix(NA,nrow=n.cv,ncol=3)
HRTE<-matrix(NA,nrow=n.cv,ncol=3)
          
    #set.seed(123)
   message("Cross Valdiations are being runing ...")
   pIndex <- c(1:n.patients)
   SurvTime[SurvTime<=0]<-quantile(SurvTime,probs =0.02)

       if (is.null(ProgFact))  {   AllX<-cbind(t(ReduGdata)) 
                                   penafac<-rep(1, ncol(X))
       
       } else {
                                   AllX<-cbind(t(ReduGdata),ProgFact)
                                   penafac<-c(rep(0,ncol(ProgFact)), rep(1, nrow(ReduGdata)))
       }
     
    perProgFact<-NULL 
    
    coef.mat<-gene.mat<-matrix(0,nrow=n.cv,ncol=nrow(ReduGdata))
    
       
     
    #------------------------  Begin FOR LOOP :1  ---------------------------
    for (i in 1:n.cv){ 
    message(i,appendLF = FALSE)
    #Gindex<-rep(FALSE,nrow(ReduGdata))
    cv.train[i,] <-sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] <-c(1:n.patients)[-c(intersect(cv.train[i,] ,c(1:n.patients)))] 

     Stime=SurvTime[cv.train[i,]]
     sen= Censor[cv.train[i,]]
   
    X<-AllX[cv.train[i,],]
   
    Run.Time[i] <-system.time(
                    COXNET.cv<-cv.glmnet(x=as.matrix(X), y=Surv(Stime,sen==1),family="cox",alpha =alpha,  nlambda = 100,  penalty.factor =penafac )
                    )[3]
                lab<-COXNET.cv$lambda.min
                seq.lambda<-COXNET.cv$lambda 
                ind.lambda<-which(seq.lambda==COXNET.cv$lambda.min)
                pld[i]<-COXNET.cv$cvm[COXNET.cv$lambda==COXNET.cv$lambda.min]  # partial likelihood deviance  
                #select the nozero coef at optimum lambda
                betahat<-coef(COXNET.cv,s="lambda.min")
                betaFP<-betahat[betahat[,1]!=0,]
                BetaFP<-betaFP
                if (class(BetaFP)[1]=="dgCMatrix") {
                betaFP<-attributes(BetaFP)$x
                names(betaFP)<-attributes(BetaFP)$Dimnames[[1]]
                }
                    #if (!is.null(dim(betaFP))) stop("Too Many variables are selected !!! try to increase alpha ")
                    
                    if (is.null(ProgFact)) { gm<-names(betaFP)
                    } else {
                                             gm<-setdiff(names(betaFP),colnames(ProgFact))
                    }
                   
                     n.g[i]<-length(gm)
                     #if (n.g[i]==0) n.gt.zero[i]<-1
                        
                        #--------------if none of the interactions are selected take the next lambda and see the solution
                       
                        bb=0                       
                                while (n.g[i]<1) {
                                bb<-bb+1
                                betahat<-coef(COXNET.cv,s=seq.lambda[ind.lambda+bb])
                                betaFP<-betahat[betahat[,1]!=0,]
                                    BetaFP<-betaFP
                                    if (class(BetaFP)[1]=="dgCMatrix") {
                                    betaFP<-attributes(BetaFP)$x
                                    names(betaFP)<-attributes(BetaFP)$Dimnames[[1]]
                                    }
                                                   # if (!is.null(dim(betaFP))) stop("Too Many variables are selected !!! try to increase alpha ")
                            
                                                        if (is.null(ProgFact)) { gm<-names(betaFP)
                                                        } else {
                                                                                 gm<-setdiff(names(betaFP),colnames(ProgFact))
                                                        }
                                
                                n.g[i]<-length(gm)
                                lab<-seq.lambda[ind.lambda+bb]
                                }
                        
                       
                        #--------------END  OF if none of the interactions are selected take the next lambda and see the solution
                        #print(betaFP)
                      lambda[i]<-lab
                      gene.mat[i,is.element(rownames(ReduGdata),gm)]<-1
                      coef.mat[i,is.element(rownames(ReduGdata),gm)]<-betaFP[gm]

      
      
                                            # Risk score for training set
                                            score2.lassoT<- betaFP[gm]%*%ReduGdata[gm,cv.train[i,]]
                                            # Risk score for Test set
                                            score2.lassoTE<-betaFP[gm]%*%ReduGdata[gm,cv.test[i,]]
                                         
                                            
                                            #######################
                                            ## train set ###########
                                            Sdata<-data.frame(SurvTime=SurvTime[cv.train[i,]],Censor=Censor[cv.train[i,]])
                                            
                                            if (!is.null(ProgFact)) perProgFact<-ProgFact[cv.train[i,],]
                                            
                                            Results1<-EstimHR(score2.lassoT, Sdata, ProgFact = perProgFact, Plots = FALSE, MedianCut = NULL)
                                           
                                           
                                            if (!is.null(ProgFact)) HRT[i,]<-(summary(Results1$SurFit)[[8]])[1,][-2]
                                            if ( is.null(ProgFact)) HRT[i,]<- summary(Results1$SurFit)[[8]][-2]

                                            
                                                        
                                            #######################
                                            ## test set ###########
                                            
                                            Sdata<-data.frame(SurvTime=SurvTime[cv.test[i,]],Censor=Censor[cv.test[i,]])
                                            if (!is.null(ProgFact)) perProgFact<-ProgFact[cv.test[i,],]
                                            
                                            Results2<-EstimHR(score2.lassoTE, Sdata, ProgFact = perProgFact, Plots = FALSE, MedianCut = median(score2.lassoT))
                                           
                                            if (!is.null(ProgFact)) HRTE[i,]<-(summary(Results2$SurFit)[[8]])[1,][-2]
                                            if ( is.null(ProgFact)) HRTE[i,]<- summary(Results2$SurFit)[[8]][-2]

                                            
           
            }#------------------------  END of FOR LOOP :1 ---------------------------

return(new("cvelasticnetcox",coef.mat=coef.mat,Run.Time=Run.Time,lambda=lambda,n.g=n.g,gene.mat=gene.mat,HRT=HRT,HRTE=HRTE,pld=pld,ReduGdata=ReduGdata))

} 

#----------------------------------------------------------------------------------------------------------
