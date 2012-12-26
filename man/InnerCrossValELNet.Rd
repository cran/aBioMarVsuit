\name{InnerCrossValELNet}
\alias{InnerCrossValELNet}
\title{Inner and Outer Cross Validations}


\description{
The function can be used to perform the further cross validations based on fixed gene list while classifier is evaluated on completely independent samples. 
The classifier is built on the weights obtain from the inner cross validations results and it is tested on out-of-bag data. 
These weights can be fixed or can be updated at each outer iteration. If weights are not fixed then patients are classified using majority votes.
 Otherwise, weights obtained from the inner cross validations are summarized by mean weights and used in the classifier. Inner cross validations are performed by calling to function \code{\link{CVLassoElasticNetCoxPh}}. Hazard ratio for low risk group is estimated using out-of-bag data.
}

\usage{
InnerCrossValELNet(fold=3,n.cv=50,n.innerCV=100,MixParAlpha=0.1,
                             GexprMat,TopGenes,WeightsFixed=FALSE,
                             SurvTime, Censor,ProgFact=NULL)
}


\arguments{
  \item{fold}{Number of folds and default 3} 
  \item{n.cv}{Number of outer cross validations and default 50.}
  \item{n.innerCV}{Number of inner cross validations.  Inner cross validations are performed using the function \code{CVLassoElasticNetCoxPh}. Default is 100.}
  \item{MixParAlpha}{Mixing parameter indicating Lasso or Elastic net penalties to be used in the inner cross validations.}
   \item{GexprMat}{Reduced or larger gene expression matrix. 
   However, function interanlly extract the gene expression of the \code{TopGenes} and use that reduced matrix for the rest of the analysis. } 
  \item{TopGenes}{Top list of genes. Usually this can be mostly selected genes by function \code{\link{CVLassoElasticNetCoxPh}}.}
 
  \item{WeightsFixed}{A logical flag indicates whether the fixed or non-fixed weights are being used during the classifier evaluations. Default is \code{FALSE}.} 
    \item{SurvTime}{A vector of survival time}
  \item{Censor}{Censoring indicator}
   \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
}



\value{
A class of \code{\link{FurtherValidation}} returns.
 }




\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{CVLassoElasticNetCoxPh}}, 
\code{\link{EstimHR}}, \code{\link{MajorityVotes}}, \code{\link{coxph}}, \code{\link{glmnet}}, \code{\link{LassoElasticNetCoxPh}}}


\examples{

\donttest{

data(exprLym)

GexprMatrix<-exprs(exprLym)
SurvData<-pData(exprLym)

#featureNames(exprLym)

PatietId<-rownames(SurvData[!is.na(SurvData[,c("IPI")]),])
Gdata<-GexprMatrix[,PatietId]
dim(Gdata)
#colnames(Gdata)
SurvTime<-SurvData[!is.na(SurvData[,c("IPI")]),c("FollowUpYears")]
Censor<-ifelse(SurvData[!is.na(SurvData[,c("IPI")]),c("StatusAtFollowUp")]=="Dead",1,0)

Gdata[is.na(Gdata)]<-mean(Gdata,na.rm=T)


nPatients<-ncol(Gdata)
ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))

#Select Top K genes (Mostly selected genes during the CV)
 TopGenes<-c("TPRT(23590)", "COL19A1(1310)","RPL7L1(285855)","DMD(1756)", "NR3C1(2908)",
 "ITGA6(3655)", "PRKD1(5587)","IFRD1(3475)","BMP6(654)","SFXN4(119559)","TOPBP1(11073)" 
 ,"LILRA4(23547)", "DNTTIP2(30836)" , "HBS1L(10767)"  ,  "MUS81(80198)")

#Example I-----------------------------------------
# without prognostic factors and Weights are updated

WOpWup<-InnerCrossValELNet(fold=3,n.cv=5,n.innerCV=3,MixParAlpha=1,
                             Gdata,TopGenes,WeightsFixed=FALSE,
                             SurvTime, Censor,ProgFact=NULL)
show(WOpWup)
summary(WOpWup)

#Example II-----------------------------------------                             
#with prognostic factors 
#and Weights are fixed
                             
WpWFtrue<-InnerCrossValELNet(fold=3,n.cv=20,n.innerCV=50,MixParAlpha=1,
                             Gdata,TopGenes,WeightsFixed=TRUE,
                             SurvTime, Censor,ProgFact)
show(WpWFtrue)
summary(WpWFtrue)   
plot(WpWFtrue,ptype=1)


#and Weights are NOT fixed
WpWFfalse<-InnerCrossValELNet(fold=3,n.cv=20,n.innerCV=50,MixParAlpha=1,
                             Gdata,TopGenes,WeightsFixed=FALSE,
                             SurvTime, Censor,ProgFact)
show(WpWFfalse)
summary(WpWFfalse)                            
plot(WpWFfalse,ptype=1)



# compare results based on classifier with fixed weights versus classifer with weights 
# changing in which  final classification based on majority votes

Res<-data.frame(HRTrue=slot(WpWFtrue, "HRTE")[,1],HRfalse=slot(WpWFfalse, "HRTE")[,1])
boxplot(Res,col="red",ylim=c(0,1),main="Estimated HR based on completely
independent Samples",names=c("Weights Fixed","Weights Updated"))

 
# estimated density of the HR 
plot(WpWFtrue,ptype=2)



}
}

\keyword{glmnet, Lasso, Elastic Net, CV}
