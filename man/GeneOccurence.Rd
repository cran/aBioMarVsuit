\name{GeneOccurence}
\alias{GeneOccurence}
\title{Mostly selected genes by cross validated gene by gene analysis}


\description{
This function searches for mostly selected top genes during the cross validation of gene by gene analysis. Top genes are ranked based on estimated HR for low risk. Therefore top gene should have minimum HR estimate.
Number of top K genes need to be considered can be given in advanced. Finally it visualizes the genes that are selected at least 5\% times during the cross validations.
}

\usage{
GeneOccurence(ObjectCv,TopK=20,minFreq=5)
}



\arguments{
\item{ObjectCv}{Object returns from \code{\link{CVGeneSpecificCoxPh}}.}
\item{TopK}{Top K genes (20 by default) to be used in the sequential analysis.}
  \item{minFreq}{Minimum number of times a particular gene has been selected to be in the top gene list, 5 by default. 
  If the number specified is larger than the maximum selection frequency then \code{max(frequency)-1} will automatically be used as minFreq.} 

}


\value{
 A vector of frequencies of which genes are selected. 
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{CVGeneSpecificCoxPh}}, \code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, \code{\link{GeneSpecificCoxPh}}}


\examples{

\donttest{
data(exprLym)

GexprMatrix<-exprs(exprLym)
SurvData<-pData(exprLym)

PatietId<-rownames(SurvData[!is.na(SurvData[,c("IPI")]),])
Gdata<-GexprMatrix[,PatietId]
dim(Gdata)
#colnames(Gdata)
SurvTime<-SurvData[!is.na(SurvData[,c("IPI")]),c("FollowUpYears")]
Censor<-ifelse(SurvData[!is.na(SurvData[,c("IPI")]),c("StatusAtFollowUp")]=="Dead",1,0)

Gdata[is.na(Gdata)]<-mean(Gdata,na.rm=T)

#-- perform gene by gene cross validated analysis --
CVR1<-CVGeneSpecificCoxPh(fold=3,
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL,
                        n.cv=10)
show(CVR1)                         

#-- Mostly selected genes
GeneOccurence(CVR1,TopK=20,minFreq=2)


}

}

\keyword{Survival, PCA, Risk, Hazard Ratio}
