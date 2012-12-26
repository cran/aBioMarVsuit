####**********************************************************************
####**********************************************************************
####
####  LAsso and Elastic net for cox with gene expression
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
LassoElasticNetCoxPh<-function(SurvTime,
                         Censor,
                        ReduGdata,
                        ProgFact=NULL,
                        Plots = FALSE,
                        MedianCut = NULL ,
                        GeneList=NULL,
                        StZ=TRUE,
                        alpha=1){

            if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
            if (missing(ReduGdata)) stop("Argument 'ReduGdata' is missing...")
            if (missing(Censor)) stop("Argument 'Censor' is missing...")
            
            
                       
            
            n.genes<-nrow(ReduGdata)
            n.patients<-ncol(ReduGdata)

             
        
        if (is.null(GeneList)) { ReduGdata<-ReduGdata
                                } else {
                                
                                ReduGdata<-ReduGdata[is.element(rownames(ReduGdata),GeneList),]
                                
                                }
                                
        
        genes.names<-rownames(ReduGdata)
        #rename the col names of Gene*Trt matrix

       if (is.null(ProgFact))  {   X<-t(ReduGdata)
                                   penafac<-rep(1, ncol(X))
       
       } else {
                                   X<-cbind(t(ReduGdata),ProgFact)
                                   penafac<-c(rep(0,ncol(ProgFact)), rep(1, nrow(ReduGdata)))
       }
        
        SurvTime[SurvTime<=0]<-quantile(SurvTime,probs =0.02)

        COXNET.cv<-cv.glmnet(x=as.matrix(X), y=Surv(as.vector(SurvTime),as.vector(Censor)==1),family="cox",alpha =alpha,  nlambda = 100,  penalty.factor =penafac,standardize = StZ )
        COXNET<-glmnet(x=as.matrix(X), y=Surv(SurvTime,Censor==1),family="cox",alpha =alpha,  nlambda = 100,  penalty.factor =penafac,standardize = StZ  )
        
        lambda<-COXNET.cv$lambda.min
        #windows("Shrinkage of coefficients")
          
        seq.lambda<-COXNET.cv$lambda 
        ind.lambda<-which(seq.lambda==COXNET.cv$lambda.min)
        #select the nozero coef at optimum lambda
        betahat<-coef(COXNET.cv,s="lambda.min")
        betaFP<-betahat[betahat[,1]!=0,]
        
        if (!is.null(dim(betaFP))) stop("Too Many variables are selected !!! try to increase alpha ")
        
        if (is.null(ProgFact)) { gm<-names(betaFP)
        } else {
                                 gm<-setdiff(names(betaFP),colnames(ProgFact))
        }
        
        ng<-length(gm)
        
                lab<-COXNET.cv$lambda.min                     
                                #--------------if none of the genes are selected take the next lambda and see the solution
                                bb=0                       
                                while (ng<1) {
                                bb<-bb+1
                                betahat<-coef(COXNET.cv,s=seq.lambda[ind.lambda+bb])
                                
                                betaFP<-betahat[betahat[,1]!=0,]
                                
                                        if (is.null(ProgFact)) { gm<-names(betaFP)
                                            } else {
                                                                 gm<-setdiff(names(betaFP),colnames(ProgFact))
                                            }
                                lab<-seq.lambda[ind.lambda+bb]
                                ng<-length(gm)
                                }
        
        
        
        #lab<-COXNET.cv$lambda.min                     
    # Risk score 
        RiskScore<-betaFP[gm]%*%t(X[,gm])
        
   if (Plots) {
         par(mfrow=c(2,2),mar=c(5,4,4,1))
        #plot(COXNET.cv,xvar="lambda")
        plot(log(COXNET.cv$lambda),COXNET.cv$cvm,main=paste("Alpha = ",alpha,sep=""),ylim=c(min(COXNET.cv$cvlo),max(COXNET.cv$cvup)),
        xlab="log(lambda)",ylab="Partial Likelihood Deviance",pch=19,col="red")
        for (i in 1:length(COXNET.cv$cvm)) lines(log(c(COXNET.cv$lambda[i],COXNET.cv$lambda[i])),c(COXNET.cv$cvlo[i],COXNET.cv$cvup[i]))
       
        plot(COXNET,xvar="lambda",label=TRUE)
        abline(v=log(COXNET.cv$lambda.min),lwd=2,lty=2,col="black")
        abline(v=log(lab),lwd=2,lty=2,col="red")
       
       
        gid<-as.factor(ifelse(RiskScore>median(RiskScore),"high","low"))
        c2data=data.frame(SurvTime, Censor,gid )
        mfit <- survfit(Surv(SurvTime, Censor == 1)~gid, data = c2data)
        plot(mfit,ylab="Survival",xlab="Time",col=c("red","green"),main="Kaplan-Meier curves ")
        legend(1,0.4,col=c("green","red"),c("Low Risk","High Risk"),lty=rep(1,2),bty="n")
                            
        boxplot(SurvTime~gid,ylim=range(SurvTime),ylab="Survival Time",main="Low vs High Risk groups",col=c("red","green"))
   }  
   

        # windows("Survival Fit")
        Sdata<-data.frame(SurvTime,Censor)
        Results<-EstimHR(RiskScore, Sdata, ProgFact = ProgFact, Plots = FALSE, MedianCut = MedianCut)
        betaS<-betaFP
return(list(RiskScore=RiskScore,Results=Results,betaS=betaS,gm=gm))

}
