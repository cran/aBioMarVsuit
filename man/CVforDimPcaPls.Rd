\name{CVforDimPcaPls}
\alias{CVforDimPcaPls}
\title{Cross validations for PCA and PLS based methods }


\description{
Cross validations for the analysis performs by \code{\link{SurFitPlsClasif}} and \code{\link{SurFitPcaClasif}} functions where the dimension reduction methods are PCA and PLS.
}

\usage{
CVforDimPcaPls( fold=3,SurvTime, Gdata,Censor,ReduceDim=TRUE,NuFeToBeSel=150,
                ProgFact=NULL,Plots = FALSE,
                 n=5 ,  DR ="PCA",   mtitle="")
}


\arguments{
  \item{fold}{number of folds, 3 as default.} 
   \item{SurvTime}{A vector of survival time.}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should 
              be equal to the number of genes and number of columns should be equal to number of patients.} 
  \item{Censor}{Censoring indicator.}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix and to be used in the classifier.
  } 
  \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the PCA or PLS based methods. This is valid only if ReduceDim=\code{TRUE}.}
    \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{Plots}{A boolean parameter indicating the plots should be shown. Default is FALSE}
 \item{n}{number of validations, 5 as default.} 
   \item{DR}{A dimension reduction method. The values are \code{"PCA"} for Principle components analysis and \code{"PLS"} for Partial least squares.} 
  \item{mtitle}{A title to be used in the plots.} 
}


\value{
A object of class \code{\link{CVfordimMethods}} 
}

\author{Pushpike Thilakarathne and  Ziv Shkedy}
\seealso{\code{\link{plsr}}, \code{\link{princomp}}, \code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, \code{\link{GridAnalysis}} }

\references{
Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}, Springer-Verlag.  \cr

Vinzi, V.; Chin, W.W.; Henseler, J. et al., eds. (2010).  Handbook of Partial Least Squares. \emph{ISBN 978-3-540-32825-4}. \cr

Bair, E., Hastie, T., Debashis, P., and Tibshirani, R. (2006). Prediction
by supervised principal components. \emph{American Statistics Association,
101(473):119-137}.


}


\examples{

\donttest{
#load Survival data
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

nPatients<-ncol(Gdata)
ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))

          
R1<-CVforDimPcaPls(fold=3,
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=ProgFact,
                        Plots = FALSE,
                        n=100 ,
                         DR ="PCA", 
                         mtitle="")
plot(R1,ylim=c(0,3))
show(R1)                         
                         
                         
R2<-CVforDimPcaPls(fold=3,
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=ProgFact,
                        Plots = FALSE,
                        n=100 ,
                         DR ="PLS", 
                         mtitle="")
plot(R2,ylim=c(1,50))
show(R2)

}
}

\keyword{Cross Validations, PCA, PLS, CoxPh}
