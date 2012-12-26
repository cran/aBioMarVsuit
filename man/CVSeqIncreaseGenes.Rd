\name{CVSeqIncreaseGenes}
\alias{CVSeqIncreaseGenes}
\title{Cross Valdiation for Top \code{k1, k2,...,kn} genes}


\description{
This function further processes the cross validated gene by gene analysis results.  And function sequentially considers top \code{k1, k2,...,kn} number of genes based on the estimated HR on the test data of the cross validated gene by gene analysis. 
For each top K genes function recompute first  PCA or PLS on train data and estimate risk scores on both test and train data only on the gene expression matrix with top k genes. 
Patients are then classified as having low or high risk 
based on the test data where the cutoff used is median of the risk score based on train data. Finally hazard ratio is estimated based on test data. The process is repeated for each top K gene sets.

The function can be interpreted as cross validated version of the function \code{\link{SeqIncreaseGenes}}. }

\usage{
CVSeqIncreaseGenes(ObjectCv,top=seq(5,100,by=5),SurvTime,Censor, ProgFact=NULL)
}



\arguments{
\item{ObjectCv}{Object of class \code{\link{CVGbyG}}.}
\item{top}{Top \code{k1, k2,...,kn} number of genes.}
  \item{SurvTime}{A vector of survival time}
  \item{Censor}{Censoring indicator}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
}


\value{
 A class object \code{\link{CVSeqInc}} returns. 
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{SeqIncreaseGenes}}, \code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, \code{\link{GeneSpecificCoxPh}}}


\examples{

\donttest{
ata(exprLym)

GexprMatrix<-exprs(exprLym)
SurvData<-pData(exprLym)

PatietId<-rownames(SurvData[!is.na(SurvData[,c("IPI")]),])
Gdata<-GexprMatrix[,PatietId]
dim(Gdata)
#colnames(Gdata)
SurvTime<-SurvData[!is.na(SurvData[,c("IPI")]),c("FollowUpYears")]
Censor<-ifelse(SurvData[!is.na(SurvData[,c("IPI")]),c("StatusAtFollowUp")]=="Dead",1,0)

Gdata[is.na(Gdata)]<-mean(Gdata,na.rm=T)

#-------------------
CVR1<-CVGeneSpecificCoxPh(fold=3,SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150,
                        ProgFact=NULL, MedianCut = NULL, n.cv=10)

plot(CVR1,ylim=c(0,2))

#-------------------
CVTopK<-CVSeqIncreaseGenes(CVR1,top=seq(5,100,by=5),SurvTime,Censor, ProgFact=NULL)

show(CVTopK)
plot(CVTopK)
summary(CVTopK)
}

}

\keyword{Survival, PCA, Risk, Hazard Ratio}
