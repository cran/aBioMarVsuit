#-------
# Pushpike Thikarathne

#All methods and Classes

setGeneric("summary", function(object, ...) standardGeneric("summary"))


#=======================================================================================================


setClass("FurtherValidation",representation(RunTime="vector",fold="numeric", n.cv="numeric",n.innerCV="numeric",TopK="vector",
HRInner="array",HRTE="matrix",WeightMat="matrix"),
                           prototype=list(RunTime=c(NA),fold=3,n.cv=10, n.innerCV=100, TopK=c(NA),HRInner=array(NA,c(1,1,1)) , HRTE=matrix(1,1,1), WeightMat=matrix(1,1,1))
                           )

setMethod("show",signature="FurtherValidation"
, function(object){
          cat("Further Cross Valdiated Results for  \n Lasso and Elastic Net based Predictive Gene signature\n")
          cat("Number of Outer CV used: ", object@n.cv, "\n")
          cat("Number of Inner CV used: ", object@n.innerCV, "\n")
          cat("Number of Genes subject to shrinkage : ", length(object@TopK), "\n")
          })

setMethod("summary",signature="FurtherValidation", function(object){
        cat("Summary of Further Cross Validations\n")
        cat("Estimated Mean Weights for the Classifier\n")
        wmat<-colSums(object@WeightMat)/nrow(object@WeightMat)
        names(wmat)<-object@TopK
        print(wmat)
})



#=======================================================================================================


setClass("cvelasticnetcox",representation(coef.mat="matrix",Run.Time="vector",lambda="vector",n.g="vector",gene.mat="matrix",HRT="matrix",HRTE="matrix",pld="vector",ReduGdata="matrix"),
                           prototype=list(coef.mat=matrix(1,1,1),Run.Time=c(NA),lambda=c(NA), n.g=c(NA), gene.mat=matrix(1,1,1),HRT=matrix(1,1,1) , HRTE=matrix(1,1,1) , pld=c(NA),ReduGdata=matrix(1,1,1))
                           )

setMethod("show",signature="cvelasticnetcox"
, function(object){
          cat("Cross Valdiated Results for Lasso and Elastic Net based Predictive Gene signature\n")
          cat("Number of CV: ", length(object@lambda), "\n")
          })

setMethod("summary",signature="cvelasticnetcox", function(object){
        cat("Summary of Cross Validations\n")
                  
        cat("Estimated  quantiles of HR on test data\n")
        print(quantile(object@HRTE[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
        cat("\n")
        cat("Estimated quantiles of HR on train data\n")
        print(quantile(object@HRT[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
        cat("Mostly selected 30 features:\n")
        Freq=colSums(object@gene.mat)
        names(Freq)<-rownames(object@ReduGdata)
        sFreq<-sort(Freq,decreasing = TRUE)
        sFreq<-sFreq[sFreq>0]
        maxG<-length(sFreq)
        if (maxG>30) maxG<-30
   
   print(names(sFreq)[1:maxG])
})










#=======================================================================================================


setClass("Permutation",representation(HRlowObs="vector",HRlowPerm="matrix",n.perm="numeric",Validation="vector"),
                           prototype=list(HRlowObs=as.vector(rep(NA,3)),HRlowPerm=matrix(1,1,1),n.perm=100,Validation=c(NA))
                           )

setMethod("show",signature="Permutation"
, function(object){
          cat("Estimated Null Ditribution of the ",object@Validation,"\n")
          cat("Number of Permutations: ", object@n.perm, "\n")
          })

setMethod("summary",signature="Permutation", function(object){
          cat("Summary of Permutation Analysis\n")
          cat("validation scheme used :",object@Validation,"\n")
          cat("Number of Permutations: ", object@n.perm, "\n")
          cat("Estimated  quantiles of the null distribution of HR\n")
          print(quantile(object@HRlowPerm[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
          cat("\n")
          cat("Estimated HR on original data\n")
          ttt<-object@HRlowObs
          names(ttt)<-c("Estimate","lower95CI","Upper95CI")
          print(ttt)

})

#=======================================================================================================


setClass("CVMV",representation(HRp.train="matrix",HRp.test="matrix",nCV="numeric",ReduGdata="matrix",pFactors="vector"),
                           prototype=list(HRp.train=matrix(1,1,1),HRp.test=matrix(1,1,1),nCV=100,ReduGdata=matrix(1,1,1),pFactors=c(NA))
                           )

setMethod("show",signature="CVMV"
, function(object){
          cat("Cross Validated Majority Votes Based Classification Analysis\n")
          cat("Number of cross valdiations used: ", object@nCV, "\n")
          if (!is.na(object@pFactors)) cat("Prognostic factors used: ",object@pFactors,"\n")
          })



#=======================================================================================================


setClass("CVSeqInc",representation(HRPC="array",HRPL="array",n.genes="numeric",n.cv="numeric",top="numeric"),
                           prototype=list(HRPC=array(NA,dim=c(1,1,1)),HRPL=array(NA,dim=c(1,1,1)),n.genes=1,n.cv=3,top=seq(5,100,by=5))
                           )

setMethod("show",signature="CVSeqInc"
, function(object){
          cat("Cross Validated Top K1, K2, ..., Kn Analysis\n")
          cat("Number of Top K Genes used: ", object@top, "\n")
          cat("Number of cross valdiations used: ", object@n.cv, "\n")
          })


setMethod("summary",signature="CVSeqInc", function(object){
          cat("Results Based on Test Data\n")
          cat("Summary of Cross Validated Top K Genes  Analysis\n")
          cat("Estimated Median of the cross Validated HR for Top K genes \n")
nn<-dim(object@HRPC)[3]
mean.alpha<- sapply(1:nn,function(i) median(object@HRPC[,1,i],na.rm=T))
se.alphal<- sapply(1:nn,function(i) quantile(object@HRPC[,1,i],na.rm=T,probs = c(0.025)))
se.alphau<- sapply(1:nn,function(i) quantile(object@HRPC[,1,i],na.rm=T,probs = c(0.975)))
mx1<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


mean.alpha<- sapply(1:nn,function(i) median(object@HRPL[,1,i],na.rm=T))
se.alphal<- sapply(1:nn,function(i) quantile(object@HRPL[,1,i],na.rm=T,probs = c(0.025)))
se.alphau<- sapply(1:nn,function(i) quantile(object@HRPL[,1,i],na.rm=T,probs = c(0.975)))
mx3<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


HR<-data.frame(rbind(mx1,mx3))
colnames(HR)<-object@top
rownames(HR)<-c("HR  (PCA1)","HR  (PLS1)")
print(HR)
}
)


#=======================================================================================================


setClass("CVGbyG",representation(HRT="array",HRTE="array",ind.train="matrix",ind.test="matrix",n.genes="numeric",n.cv="numeric",ReduGdata="matrix"),
                           prototype=list(HRT=array(NA,dim=c(1,1,1)),HRTE=array(NA,dim=c(1,1,1)),ind.train=matrix(0,0,0),ind.test=matrix(0,0,0),n.genes=1,n.cv=3,ReduGdata=matrix(0,0,0))
                           )



setMethod("show",signature="CVGbyG"
, function(object){
          cat("Cross Validated Gene by Gene Analysis\n")
          cat("Number of Genes used: ", object@n.genes, "\n")
          cat("Number of cross valdiations used: ", object@n.cv, "\n")
          })

setMethod("summary",signature="CVGbyG", function(object,Which=1){
          cat("Summary of Cross Validated Gene by Gene Analysis\n")
          cat("Estimated Median of the cross Validated HR for Feature: ",Which,"\n")
  HRTEST  <- object@HRTE[Which,,][1,]
  HRTrain <- object@HRT[Which,,][1,]

  
  mean.alpha <- median(HRTrain,na.rm=T)
  se.alphal <- quantile(HRTrain,na.rm=T,probs = c(0.025))
  se.alphau <- quantile(HRTrain,na.rm=T,probs = c(0.975))
  cat("Estimated HR for Train Data \n")
  cat(paste(round(mean.alpha,4),"(",round(se.alphal,4),"-",
  round(se.alphau,4),")",sep=""))


  cat("\n")
  mean.alpha <- median(HRTEST,na.rm=T)
  se.alphal <- quantile(HRTEST,na.rm=T,probs = c(0.025))
  se.alphau <- quantile(HRTEST,na.rm=T,probs = c(0.975))
  cat("Estimated HR for Test Data \n")
  cat(paste(round(mean.alpha,4),"(",round(se.alphal,4),"-",
  round(se.alphau,4),")",sep=""))
  cat("\n")
            }
)







#=======================================================================================================

setClass("GeneSpecific",representation(res="list",HRp="matrix",gr="matrix",gNames="vector"),
                           prototype=list(res=list(1),HRp=matrix(0,0,0),gr=matrix(0,0,0),gNames=vector()))



setMethod("show",signature="GeneSpecific"
, function(object){
          cat("Gene by Gene CoxPh Model\n")
          cat("Number of Genes used: ", length(object@gNames), "\n")
          })




setMethod("summary",signature="GeneSpecific", function(object){
          cat("Summary of Gene by Gene CoxPh Models\n")
          cat("Number of Genes used: ", length(object@gNames), "\n")
          cat("Top 15 Genes out of ", length(object@gNames), "\n")
           cat("Estimated HR for low risk group\n")
          # top genes based on upper CI HR GS+
            Names.KGenes<-object@gNames  
            index.Top.KGenes <- order(object@HRp[,1],decreasing =FALSE)
            index.Top.KGenes<-index.Top.KGenes[1:15]
            Top.KGenes.GSplus<-data.frame(GeneNames=Names.KGenes[index.Top.KGenes],object@HRp[index.Top.KGenes,-2])
            
            colnames(Top.KGenes.GSplus)<-c("GeneNames","HR GS+","LCI","UCI")
            #Top.KGenes.GSplus
            
            # FDR corrected CI for top k  genes
            
            cilevel <- 1-0.05*15/nrow(object@HRp) 
            #qnorm(cilevel)
            HRpadj <- exp(log(object@HRp[index.Top.KGenes,1]) + log(object@HRp[index.Top.KGenes,c(3,4)])-log(object@HRp[index.Top.KGenes,1])*qnorm(cilevel)/1.96)  # 
            res.topkgenes<-data.frame(Top.KGenes.GSplus,HRpadj)
            
            colnames(res.topkgenes)<-c("GeneNames","HR","LCI","UCI","FDRLCI","FDRUCI")
           print(res.topkgenes)
          })

#=======================================================================================================



setClass("CVfordimMethods",representation(Results="data.frame",n="numeric",PCAorPLS="vector",cv.train="matrix",NuFeToBeSel="numeric"),
                           prototype=list(Results=data.frame(1),n=numeric(),PCAorPLS="PCA",cv.train=matrix(0,0,0),NuFeToBeSel=numeric()))




setMethod("show",signature="CVfordimMethods"
, function(object){
          cat("Cross Validated PLS and PCA based HR estimation methods\n")
          cat("Dimention Reduction Method: ", object@PCAorPLS,"\n", sep="")
          cat("Number of CVs : ", object@n,"\n", sep="")
          cat("Number of Genes: ", object@NuFeToBeSel, "\n")
          })



#=======================================================================================================
