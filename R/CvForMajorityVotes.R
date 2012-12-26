
# ------------------------------------------------------  Cross validation for MV ---------------------------------------------------------
CvForMajorityVotes<-function(SurvTime,Censor, ProgFact=NULL, Gdata, ReduceDim=TRUE, NuFeToBeSel=150, fold=3, nCV=100){

options( warn = -1)


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
n.train<-(n.patients-floor(n.patients/fold))
n.test<-floor(n.patients/fold)
ind.train <-matrix(0,nCV,n.train)
ind.test  <-matrix(0,nCV,n.test)
res <-  vector("list", n.genes)
    

HRp.train <- matrix(0,nCV,3)  # Training 
HRp.test <-  matrix(0,nCV,3)     # TTsting 
    
set.seed(123)
message("Cross Valdiations are being runing ...")
   pIndex <- c(1:n.patients)
   res <-res1<-res2<-res3<- vector("list", nCV)      
     
     for (j in 1:nCV){ # loop over cross validations:  ------------------------------------------ j=4
     p1<-NA
     p2<-NA    
     gr.train <- matrix(0, n.genes,ncol(ind.train))
     gr.test  <- matrix(0,n.genes,ncol(ind.test)) 
       
     ind.train[j,] <-sort(sample(pIndex,n.train,replace=F) )
     ind.test[j,] <-c(1:n.patients)[-c(intersect(ind.train[j,] ,c(1:n.patients)))] 
     
                 
           for (i in 1:n.genes){  #---------------------------  STRAT FOR LOOP over genes ------------------------ j=1

                #---------------------------------------  Training Set -------------------------------------------
                genei <- ReduGdata[i,ind.train[j,]]
                
                if (is.null(ProgFact)) {
                        
                cdata <- data.frame(SurvTime=SurvTime[ind.train[j,]],Censor=Censor[ind.train[j,]],genei)
                m0 <- coxph(Surv(SurvTime, Censor==1) ~ genei,data=cdata)
                }
                if (!is.null(ProgFact)) {
                    if (is.data.frame(ProgFact)) {
                    nPrgFac<-ncol(ProgFact)
                    cdata <- data.frame(SurvTime=SurvTime[ind.train[j,]],Censor=Censor[ind.train[j,]],genei,ProgFact[ind.train[j,],])
                    NameProg<-colnames(ProgFact)
                    eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ genei",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
                    } else {
                             
                    stop(" Argument 'ProgFact' is NOT a data frame ")
                    }
                
                }
                #risk Score       
                TrtandPC1<-summary(m0)[[7]][c("genei"),1]
                p1 <- TrtandPC1*genei
              
               TempGenei <-EstimHR(p1,Sdata=cdata, ProgFact=ProgFact[ind.train[j,],],Plots = FALSE, MedianCut = NULL )  
               gr.train[i,]<-TempGenei$gs
                
             
                #---------------------------------------  Testing  Set -------------------------------------------
                mp1 <- median(p1)
                geneit <- ReduGdata[i,ind.test[j,]]
                p2 <- TrtandPC1*geneit

                pos <- c(1:length(p2))[p2<=mp1] 
                neg <- c(1:length(p2))[p2>mp1] 
  
                gs <- NULL
                gs[pos] <- 1
                gs[neg] <- -1 
                gr.test[i,] <- gs
  } # END OF LOOP over Genes---------------------------------------------------------------------
                
 # ------------ count majority votes for j th Cross validation and estimate HR --------------
         
        ggr.train <- per.R<-per.NR <- NULL
        for (k in 1:ncol(ind.train)){
        per.R[k]<-sum(gr.train[,k]==1)
        ggr.train[k]<-ifelse((n.genes-per.R[k])>per.R[k],"high","low")
        }
         

        
        #-------------------- HR estimation for Training  ----------------------
          GS<-as.factor(ggr.train)
          if (is.null(ProgFact)) {
                        
                cdata <- data.frame(SurvTime=SurvTime[ind.train[j,]],Censor=Censor[ind.train[j,]],GS)
                mTrain <- coxph(Surv(SurvTime, Censor==1) ~ GS,data=cdata)
                }
          if (!is.null(ProgFact)) {
                    if (is.data.frame(ProgFact)) {
                    nPrgFac<-ncol(ProgFact)
                    cdata <- data.frame(SurvTime=SurvTime[ind.train[j,]],Censor=Censor[ind.train[j,]],GS,ProgFact[ind.train[j,],])
                    NameProg<-colnames(ProgFact)
                    eval(parse(text=paste( "mTrain <-coxph(Surv(SurvTime, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
                    } else {
                      stop(" Argument 'ProgFact' is NOT a data frame ")
                    }
                
                }      
        HRp.train[j,]<-(summary(mTrain)[[8]][1,])[-2]            
                         
        #-------------------- HR estimation for Testing  ----------------------
        ggr.test <- per.R<-per.NR <- NULL
        for (k in 1:ncol(ind.test)){
        per.R[k]<-sum(gr.test[,k]==1)
        ggr.test[k]<-ifelse((n.genes-per.R[k])>per.R[k],"high","low")
        }
        
        GS<-as.factor(ggr.test)
          if (is.null(ProgFact)) {
                        
                cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],GS)
                mTest <- coxph(Surv(SurvTime, Censor==1) ~ GS,data=cdata)
                }
          if (!is.null(ProgFact)) {
                    if (is.data.frame(ProgFact)) {
                    nPrgFac<-ncol(ProgFact)
                    cdata <- data.frame(SurvTime=SurvTime[ind.test[j,]],Censor=Censor[ind.test[j,]],GS,ProgFact[ind.test[j,],])
                    NameProg<-colnames(ProgFact)
                    eval(parse(text=paste( "mTest <-coxph(Surv(SurvTime, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
                    } else {
                      stop(" Argument 'ProgFact' is NOT a data frame ")
                    }
                
                }      
        HRp.test[j,]<-(summary(mTest)[[8]][1,])[-2]    
             
    }#---------------------------  END OF  FOR LOOP over Cross Validations ------------------------

 pFactors<-NA
 if (!is.null(ProgFact)) pFactors <-colnames(ProgFact)
 
return(new("CVMV",HRp.train=HRp.train,HRp.test=HRp.test,nCV=nCV,ReduGdata=ReduGdata, pFactors=pFactors))

}
