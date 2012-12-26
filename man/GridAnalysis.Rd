\name{GridAnalysis}
\alias{GridAnalysis}
\title{Classification, survival fit and HR estimation}


\description{
The function performs the sensitivity of the cut off values used in classification and percentiles of risk score obtained under \code{\link{SurFitPlsClasif}} or \code{\link{SurFitPcaClasif}} methods is used as fine grid points.
}

\usage{
GridAnalysis(SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150, ProgFact=NULL, 
Plots = FALSE,DimMethod=c("PLS","PCA"))
}



\arguments{
  \item{SurvTime}{A vector of survival time}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should be equal to the number of genes and number of columns should be equal to number of patients.
               } 
  \item{Censor}{Censoring indicator}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix and to be used in the classifier.
  } 
  \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the PCA or PLS based methods. This is valid only if ReduceDim=\code{TRUE}.}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{Plots}{A boolean parameter indicating the plot should be shown. Default is FALSE}
  \item{DimMethod}{Dimension reduction method PCA or PLS.} 
}


\value{
 A data frame with estimated HR along with 95\% CI at each cut off value. 
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}} }


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


GridAnalysis(SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150, 
ProgFact=NULL, Plots = TRUE,DimMethod="PLS")

GridAnalysis(SurvTime,Gdata,Censor,ReduceDim=TRUE, NuFeToBeSel=150, 
ProgFact=NULL, Plots = FALSE,DimMethod="PCA")

}

}

\keyword{Survival, PCA, Risk, Hazard Ratio}
