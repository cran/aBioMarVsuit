\name{SeqIncreaseGenes}
\alias{SeqIncreaseGenes}
\title{Sequentially increase the number of top K genes to be used in the classifier}


\description{
This function sequentially increases the number of top K genes to be used in the PCA or PLS methods in order to obtain the risk score. This function 
internally calls \code{\link{GeneSpecificCoxPh}} to 
rank the genes based on HR for each gene. Therefore genes can be ordered as HR for low risk is increasing. 
After that function takes few top K (15 is default) to be used in the sequential analysis.
}

\usage{
SeqIncreaseGenes(TopK=15,SurvTime,Gdata,Censor,ReduceDim=TRUE, 
NuFeToBeSel=150, ProgFact=NULL, Plots = FALSE,DimMethod=c("PLS","PCA"),...)
}



\arguments{
\item{TopK}{Top K genes (15 by default) to be used in the sequential analysis.}
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
  \item{\dots}{Further arguments for plotting and valid only if Plot=\code{TRUE}}
}


\value{
 A data frame with estimated HR along with 95\% CI at each cut off value for the sequential analysis. 
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{CVSeqIncreaseGenes}}, \code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, \code{\link{GeneSpecificCoxPh}}}


\references{
Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}, Springer-Verlag.  \cr

Vinzi, V.; Chin, W.W.; Henseler, J. et al., eds. (2010).  Handbook of Partial Least Squares. \emph{ISBN 978-3-540-32825-4}. \cr

Bair, E., Hastie, T., Debashis, P., and Tibshirani, R. (2006). Prediction
by supervised principal components. \emph{American Statistics Association,
101(473):119-137}.

  }

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


#- sequentially increase the number of top K genes -

SeqIncreaseGenes(TopK=15,SurvTime,Gdata,Censor,ReduceDim=TRUE, 
NuFeToBeSel=150, ProgFact=NULL, Plots = FALSE,DimMethod="PLS")
}

}

\keyword{Survival, PCA, Risk, Hazard Ratio}
