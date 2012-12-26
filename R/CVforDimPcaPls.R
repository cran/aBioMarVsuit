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






#---------------------------------------------------------------------------------------
# function for Cross validation for PLS1 and PCA1 
#---------------------------------------------------------------------------------------

InterMediatePca<-function(ReduGdata,ProgFact,SurvTime,Censor,index){
    if (is.matrix(ReduGdata)) {
        pc1 <- f.pca(as.matrix(ReduGdata[,index]))[[6]][,1]
        } else {
        pc1<-ReduGdata[,index]
        }         
        
    if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],pc1)
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ pc1,data=cdata)
    }
    if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],pc1,ProgFact[index,])
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
         } else {
         
         stop(" Argument 'ProgFact' is NOT a data frame ")
         }

    }
    
return(list(m0=m0,pc1=pc1,cdata=cdata))
}


InterMediatePls<-function(ReduGdata,ProgFact,SurvTime,Censor,index){
 if (is.matrix(ReduGdata)) {
        
                        PLSforGSK<-data.frame(1:length(index))
                        PLSforGSK$g<-as.matrix(t(ReduGdata[,index]))
                        colnames(PLSforGSK)[1]<-c("SurvTime")
                        PLSforGSK[,1]<-SurvTime[index]
                        #plsr.1 <- plsr(SurvTime ~ stdize(g), ncomp = 5, scale =TRUE,data = PLSforGSK, validation =  "CV")
                        plsr.1 <- plsr(SurvTime ~ g, method="simpls",ncomp = 10, scale =TRUE,data = PLSforGSK, validation =  "CV")
                        pc1<-scores(plsr.1)[,1] # extract the first com
        } else {
        pc1<-ReduGdata[,index]
        }         
        
    if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],pc1)
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ pc1,data=cdata)
    }
    if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime=SurvTime[index],Censor=Censor[index],pc1,ProgFact[index,])
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
         } else {
         
         stop(" Argument 'ProgFact' is NOT a data frame ")
         }

    }

return(list(m0=m0,pc1=pc1,cdata=cdata))
}



CVforDimPcaPls<-function(fold=3,
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        Plots = FALSE,
                         n=5 ,
                         DR ="PCA", 
                         mtitle=""){
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
    sen<-Censor               # Censoring indicator
    dfi <- SurvTime 
    if (is.data.frame(ProgFact)) {   data1<-data.frame(sen,dfi,ProgFact)
    } else {
                                     data1<-data.frame(sen,dfi)
    }


    HRp.train <- HRn.train <- matrix(0,n,4)  # Training 
    HRp.test <- HRn.test<- matrix(0,n,4)     # Testing 
    
    n.train<-(n.patients-floor(n.patients/fold))
    n.test<-floor(n.patients/fold)
    cv.train <-matrix(0,n,n.train)
    cv.test  <-matrix(0,n,n.test)
 
 
    set.seed(123)
   message("Cross Valdiations are being runing ...")
   pIndex <- c(1:n.patients)
   res <-res1<-res2<-res3<- vector("list", n)
    #------------------------  Begin FOR LOOP :1  --------------------------- i=1
    for (i in 1:n){ 
    p1<-NA
    p2<-NA  
      
    cv.train[i,] <-sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] <-c(1:n.patients)[-c(intersect(cv.train[i,] ,c(1:n.patients)))] 
                        
                    if (DR=="PCA") { #-------------------------------------------------------------------------------------------
                    
                            Temp1<-InterMediatePca(ReduGdata,ProgFact,SurvTime,Censor,as.vector(cv.train[i,]))
                                    pc1  <- Temp1$pc1
                                    cdata <-Temp1$cdata
                                    m2 <- Temp1$m0
                            Temp2<-InterMediatePca(ReduGdata,ProgFact,SurvTime,Censor,cv.test[i,])
                                    pc1test <- Temp2$pc1
                                    ctestdata <- Temp2$cdata
                            
                          
                            # classification method A1
                            TrtandPC1<-summary(m2)[[7]][c("pc1"),1]
                            p1 <-  TrtandPC1[1]*pc1 
                            p2 <-  TrtandPC1[1]*pc1test 
                    }#-------------------------------------------------------------------------------------------
                    
                    
                    if (DR=="PLS") { #-------------------------------------------------------------------------------------------
                            
                            Temp3<-InterMediatePls(ReduGdata,ProgFact,SurvTime,Censor,cv.train[i,])
                            pls.comp1  <- Temp3$pc1
                            cdata <-Temp3$cdata
                            m3 <- Temp3$m0
                            
                            Temp4<-InterMediatePls(ReduGdata,ProgFact,SurvTime,Censor,cv.test[i,])
                            pls.comp1.test <- Temp4$pc1
                            ctestdata <- Temp4$cdata
                            
                            # classification method A1
                            TrtandPLSc1<-summary(m3)[[7]][c("pc1"),1]
                            p1 <- TrtandPLSc1[1]*pls.comp1
                            p2 <- TrtandPLSc1[1]*pls.comp1.test 
                    
                    }#-------------------------------------------------------------------------------------------
                    
                    
                    #######################
                    ## training set ###########
                    
                    TrainSet<-EstimHR(p1,Sdata=cdata, ProgFact=ProgFact[cv.train[i,],],Plots = FALSE, MedianCut = NULL )
                    
                    HRp.train[i,] <- summary( TrainSet$SurFit)[[8]][1,]
                    
                                      
                                        
                    #######################
                    ## test set ###########
                    
                    mp1 <- median(p1)
                    TestSet<-EstimHR(p2,Sdata=ctestdata, ProgFact=ProgFact[cv.test[i,],],Plots = FALSE, MedianCut = mp1 )
                    
                    HRp.test[i,] <- summary( TestSet$SurFit)[[8]][1,]
                    
                   
        }
        #------------------------  END of FOR LOOP :1


Results<-data.frame(HRpTrain=HRp.train[,1],HRpTest=as.numeric(HRp.test[,1]))
#boxplot(Results,ylim=c(0,5),names=c("Train ","Test"),main=mtitle,ylab="HR",col=c("green","red"))

return(new("CVfordimMethods",Results=Results,n=n,PCAorPLS=DR,cv.train=cv.train,NuFeToBeSel=NuFeToBeSel))
}


#--------------------------------------
