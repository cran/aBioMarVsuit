NullDistHR<-function(n.perm=100, case=2, Validation=c("PLSbased","PCAbased","L1based","MVbased"), SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL 
                              
){

options(warn=-1)
Validation <- match.arg(Validation)
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


HRlowPerm<-matrix(NA,nrow=n.perm,ncol=3)
HRlowObs<-as.vector(rep(NA,3))


#set.seed(123)
ind.s<-ind.gene<-ind.w<-matrix(NA,nrow=n.perm,ncol=n.patients)
            for (i in 1:n.perm) {
                ind.s[i,]<-1:n.patients
                ind.w[i,]<-1:n.patients
                ind.gene[i,]<-1:n.patients
                }

switch(case,
            
            {#case 2: permute survival
            for (i in 1:n.perm) {
                ind.s[i,]<-sample(c(1:n.patients),replace=F) 
                }
                        
            },
            
            {#case 3: permute survival, prognostic
            for (i in 1:n.perm) {
                ind.s[i,]<-sample(c(1:n.patients),replace=F) 
                }
            ind.w<-ind.s
            
            },
            
           
            {#case 4 b: permute survival, prognostic and permute Genes independently
            for (i in 1:n.perm) {
                ind.s[i,]<-sample(c(1:n.patients),replace=F) 
                }
                
                ind.w<-ind.s 
                for (i in 1:n.perm) {
                ind.gene[i,]<-sample(c(1:n.patients),replace=F) 
                }
                 
            },
            
            {#case 4 : permute  genes only
                for (i in 1:n.perm) {
                ind.gene[i,]<-sample(c(1:n.patients),replace=F) 
                }
            }
)

perProgFact<-NULL
   
   if (Validation=="PLSbased") {
            for (i in 1:n.perm) {
            if (!is.null(ProgFact)) perProgFact<-ProgFact[ind.w[i,],]
            
            Temp<-SurFitPlsClasif(SurvTime=SurvTime[ind.s[i,]],Gdata= ReduGdata[,ind.gene[i,]],Censor= Censor[ind.s[i,]], ReduceDim=FALSE,NuFeToBeSel=150, ProgFact=perProgFact, Plots = FALSE,  MedianCut = MedianCut )
            
            if (!is.null(ProgFact))  HRlowPerm[i,]<-(summary(Temp$SurFit)[[8]])[1,][-2]
            if ( is.null(ProgFact))  HRlowPerm[i,]<-(summary(Temp$SurFit)[[8]])[-2]
            }
   TempObs<-SurFitPlsClasif(SurvTime, ReduGdata, Censor, ReduceDim=FALSE,NuFeToBeSel=150, ProgFact=ProgFact, Plots = FALSE,  MedianCut = MedianCut )
            if (!is.null(ProgFact)) HRlowObs<-(summary(Temp$SurFit)[[8]])[1,][-2]
            if ( is.null(ProgFact)) HRlowObs<-(summary(Temp$SurFit)[[8]])[-2]
   }
   
   
   
   if (Validation=="PCAbased") {
            for (i in 1:n.perm) {
            if (!is.null(ProgFact)) perProgFact<-ProgFact[ind.w[i,],]
            Temp<-SurFitPcaClasif(SurvTime=SurvTime[ind.s[i,]], Gdata=ReduGdata[,ind.gene[i,]], Censor=Censor[ind.s[i,]], ReduceDim=FALSE,NuFeToBeSel=150, ProgFact=perProgFact, Plots = FALSE,  MedianCut = MedianCut )
          
            if (!is.null(ProgFact))  HRlowPerm[i,]<-(summary(Temp$SurFit)[[8]])[1,][-2]
            if ( is.null(ProgFact))  HRlowPerm[i,]<-(summary(Temp$SurFit)[[8]])[-2]
            
            }
   TempObs<-SurFitPcaClasif(SurvTime, ReduGdata, Censor, ReduceDim=FALSE,NuFeToBeSel=150, ProgFact=ProgFact, Plots = FALSE,  MedianCut = MedianCut )
  
            if (!is.null(ProgFact)) HRlowObs<-(summary(Temp$SurFit)[[8]])[1,][-2]
            if ( is.null(ProgFact)) HRlowObs<-(summary(Temp$SurFit)[[8]])[-2]
   }
   
   
   
   if (Validation=="MVbased") {
            for (i in 1:n.perm) {
            if (!is.null(ProgFact)) perProgFact<-ProgFact[ind.w[i,],]
            Ana1<-GeneSpecificCoxPh( SurvTime=SurvTime[ind.s[i,]], Gdata=ReduGdata[,ind.gene[i,]], Censor=Censor[ind.s[i,]], ReduceDim=FALSE,NuFeToBeSel=150, ProgFact=perProgFact,
                        MedianCut = MedianCut)
            
            Temp<-MajorityVotes(Ana1,ProgFact=ProgFact[ind.w[i,],], SurvTime=SurvTime[ind.s[i,]],Censor=Censor[ind.s[i,]],J=1)
            
            if (!is.null(ProgFact)) HRlowPerm[i,]<-(summary(Temp$m0)[[8]])[1,][-2]
            if ( is.null(ProgFact)) HRlowPerm[i,]<-(summary(Temp$m0)[[8]])[-2]
            }
            
           Ana2<-GeneSpecificCoxPh( SurvTime,
                        ReduGdata,
                        Censor,
                        ReduceDim=FALSE,
                        NuFeToBeSel=150,
                        ProgFact=ProgFact,
                        MedianCut = MedianCut)

   TempObs<-MajorityVotes(Ana2,ProgFact, SurvTime,Censor,J=1)

   if (!is.null(ProgFact)) HRlowObs<-(summary(TempObs$m0)[[8]])[1,][-2]
   if ( is.null(ProgFact)) HRlowObs<-(summary(TempObs$m0)[[8]])[-2]
   }
   
   
   
   if (Validation=="L1based") {
            for (i in 1:n.perm) {
            Temp<-NA
            if (!is.null(ProgFact)) perProgFact<-ProgFact[ind.w[i,],]
            try( Temp<-LassoElasticNetCoxPh(SurvTime=SurvTime[ind.s[i,]],Censor=Censor[ind.s[i,]],ReduGdata[,ind.gene[i,]], ProgFact=perProgFact, Plots = FALSE, 
                     MedianCut = MedianCut , GeneList=NULL, StZ=TRUE,  alpha=1)
             , silent = TRUE )        

           if ((!is.na(Temp))[1]) {
              if (!is.null(ProgFact)) HRlowPerm[i,]<-(summary(Temp$Results$SurFit)[[8]])[1,][-2]
              if ( is.null(ProgFact)) HRlowPerm[i,]<- summary(Temp$Results$SurFit)[[8]][-2]
           }
           if ((is.na(Temp))[1])  HRlowPerm[i,]<-NA
            message(i) 
            }
   TempObs<-LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata, ProgFact=ProgFact, Plots = FALSE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=1)         
                        
   if (!is.null(ProgFact)) HRlowObs<-(summary(TempObs$Results$SurFit)[[8]])[1,][-2]
   if ( is.null(ProgFact)) HRlowObs<-(summary(TempObs$Results$SurFit)[[8]])[-2]
   
   }
        
      



#
#if (!DymPlots) { 
#par(mfrow=c(2,2))
#hrDensSpec(HRGSplus.PCA,n=n.perm,vv=ObsHr[1],mymain="Density of HR GS+ \nPCA")
#hrDensSpec(HRGSminus.PCA,n=n.perm,vv=ObsHr[2] ,mymain="Density of HR GS- \nPCA")
#hrDensSpec(HRGSplus.PLS,n=n.perm,vv=ObsHr[3],mymain="Density of HR GS+ \nPLS")
#hrDensSpec(HRGSminus.PLS,n=n.perm,vv=ObsHr[4] ,mymain="Density of HR GS- \nPLS")
#}
#

return(new("Permutation",HRlowObs=HRlowObs,HRlowPerm=HRlowPerm,n.perm=n.perm,Validation=Validation))

}
