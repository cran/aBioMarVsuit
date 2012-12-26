####**********************************************************************
####**********************************************************************
####
#### 
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
                        plsr.1 <- plsr(SurvTime ~ g, method="simpls",ncomp = 2, scale =TRUE,data = PLSforGSK, validation =  "CV")
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



#------------------- main function --------------------
CVSeqIncreaseGenes<-function(ObjectCv,top=seq(5,100,by=5),SurvTime,Censor, ProgFact=NULL){

  Decrease=FALSE
if (class(ObjectCv)!="CVGbyG") stop("Invalid class object.")
if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
if (missing(Censor)) stop("Argument 'Censor' is missing...")


gene.mat<-matrix(0,ObjectCv@n.cv,ObjectCv@n.genes)
    
Gdata<-ObjectCv@ReduGdata
genes.names <- rownames(Gdata)

n.genes<-ObjectCv@n.genes
n.patients<-ncol(Gdata)


if (ObjectCv@n.genes<max(top)) stop("The max(top) should be less than or euqal to total number of genes")
            
         
HRPC<-array(NA,dim=c(ObjectCv@n.cv,3,length(top)))
HRPL<-array(NA,dim=c(ObjectCv@n.cv,3,length(top)))

#i=1
options( warn = -1)
for (j in 1: length(top)){
    for (i in 1:ObjectCv@n.cv){ 

            hrsetTE<-ObjectCv@HRTE[,,i]
            hrsetmT<-ObjectCv@HRT[,,i]
            
            Names.KGenesT<-rownames(Gdata)[order(hrsetTE[,1],decreasing=Decrease)]   
            index.Top.KGenesT <- order(hrsetTE[,1],decreasing=Decrease)
            index.Top.KGenes<-index.Top.KGenes[1:top[j]]
            
            #rxtract top K Genes
            #genei <- Gdata[Names.KGenesT[1:top[j]],]
            genei <-Gdata[intersect(rownames(Gdata),Names.KGenesT[1:top[j]]),]
            
            #PCA ----------------
            
            TrainTemp<-InterMediatePca(genei,ProgFact,SurvTime,Censor,ObjectCv@ind.train[i,])
            TestTemp<- InterMediatePca(genei,ProgFact,SurvTime,Censor,ObjectCv@ind.test[i,])   
                m1 <- TrainTemp$m0
                TrtandGene1<-summary(m1)[[7]][c("pc1"),1]
                rm(m1)
                p1.train     <- TrtandGene1[1]*TrainTemp$pc1
                p1.test     <-  TrtandGene1[1]*TestTemp$pc1
             TempGeneiTE <-EstimHR(p1.test,Sdata=data.frame(SurvTime=SurvTime[ObjectCv@ind.test[i,]],Censor=Censor[ObjectCv@ind.test[i,]]), 
                            ProgFact=ProgFact[ObjectCv@ind.test[i,],],Plots = FALSE, MedianCut = median(p1.train) )  

             HRPC[i,,j]<- (summary( TempGeneiTE$SurFit)[[8]][1,])[-2]
            
            #PLS -------------------
             TrainTemp1<-InterMediatePls(genei,ProgFact,SurvTime,Censor,ObjectCv@ind.train[i,])
             TestTemp1 <-InterMediatePls(genei,ProgFact,SurvTime,Censor,ObjectCv@ind.test[i,])   
                m1 <- TrainTemp1$m0
                TrtandGene2<-summary(m1)[[7]][c("pc1"),1]
                rm(m1)
                p1.train     <- TrtandGene2[1]*TrainTemp1$pc1
                p1.test     <- TrtandGene2[1]*TestTemp1$pc1
             TempGeneiTE <-EstimHR(p1.test,Sdata=data.frame(SurvTime=SurvTime[ObjectCv@ind.test[i,]],Censor=Censor[ObjectCv@ind.test[i,]]), 
                           ProgFact=ProgFact[ObjectCv@ind.test[i,],],Plots = FALSE, MedianCut = median(p1.train) )  

             HRPL[i,,j]<- (summary( TempGeneiTE$SurFit)[[8]][1,])[-2]
           
   #print("----------------------------------------")
   #message(i)           
   }
   #print("----------------------------------------")
   #message(j)         
}              




#par(mfrow=c(1,2))
#nn<-dim(HRPC)[3]
#PC.HRp<-HRPC[,1,1:nn]
#colnames(PC.HRp)<-top
#boxplot(PC.HRp,main="HR (PC1)",xlab="Top K Genes",ylab="HR",col=1:nn)
#
#
## summarize results ---------- based on PLS
##par(mfrow=c(1,2))
#PL.HRp<-HRPL[,1,1:nn]
#colnames(PL.HRp)<-top
#boxplot(PL.HRp,main="HR (PLS1)",xlab="Top K Genes",ylab="HR",col=1:nn)
#
#
#
##median HR
#mean.alpha<- sapply(1:nn,function(i) median(HRPC[,1,i],na.rm=T))
#se.alphal<- sapply(1:nn,function(i) quantile(HRPC[,1,i],na.rm=T,probs = c(0.025)))
#se.alphau<- sapply(1:nn,function(i) quantile(HRPC[,1,i],na.rm=T,probs = c(0.975)))
#mx1<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")
#
#
#mean.alpha<- sapply(1:nn,function(i) median(HRPL[,1,i],na.rm=T))
#se.alphal<- sapply(1:nn,function(i) quantile(HRPL[,1,i],na.rm=T,probs = c(0.025)))
#se.alphau<- sapply(1:nn,function(i) quantile(HRPL[,1,i],na.rm=T,probs = c(0.975)))
#mx3<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")
#
#
#HR<-data.frame(rbind(mx1,mx3))
#colnames(HR)<-top
#rownames(HR)<-c("HR  (PCA1)","HR  (PLS1)")
n.cv<-ObjectCv@n.cv
n.genes<-ObjectCv@n.genes
return(new("CVSeqInc",HRPC=HRPC,HRPL=HRPL,n.genes=n.genes,n.cv=n.cv,top=top))

}
