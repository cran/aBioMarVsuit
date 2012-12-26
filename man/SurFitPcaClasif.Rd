\name{SurFitPcaClasif}
\alias{SurFitPcaClasif}
\title{PCA and Survival and Classification Function}


\description{
The function will perform the PCA on gene expression matrix and fit Cox proportional hazard model with covariates including first PCA as covariates. 
It can also be used to perform the grid analysis where the grid will be several cut off values and default is median cut-off. 
This function can handle single and multiple genes.  Function will reduce larger gene expression matrix to smaller version using supervised pca approach and this is by default done and can be control by 
argument ReduceDim.
Other prognostic factors can be included to the model.
}

\usage{
SurFitPcaClasif(SurvTime,
                Gdata,
                Censor,
                ReduceDim=TRUE,
                NuFeToBeSel=150,
                ProgFact=NULL,
                Plots = FALSE,
                MedianCut = NULL)
}


\arguments{
  \item{SurvTime}{A vector of survival time}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should be equal to the number of genes and number of columns should be equal to number of patients.
               } 
  \item{Censor}{Censoring indicator}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix and to be used in the classifier.
  } 
  \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the PCA based methods. This is valid only if ReduceDim=\code{TRUE}.}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{Plots}{A boolean parameter indicating the plots should be shown. Default is FALSE}
  \item{MedianCut}{A user defined cut-off value for the classifier and default is median cut-off.} 
}


\value{
 (Graphics) Three plots: Risk vs first PCA, Fitted survival curves for two groups (high risk and low risk) and targeting odds ratio of low risk, Kaplan-Meier curves and box plots showing distibution of survival time for low and high risk groups.
And a list with following items. \cr
  \item{SurFit}{A summary of survival fit.}
  \item{p1}{A vector of risk scores and size equal to no.of patents.}
  \item{gs}{A vector of labels of patients whether they have low (1) or high risk (-1).}
  \item{pc1}{First PCA scores based on reduced gene expression matrix.}
}




\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{princomp}}, \code{\link{SurFitPlsClasif}}, \code{\link{coxph}}, \code{\link{EstimHR}},  \code{\link{CVforDimPcaPls}}, \code{\link{GridAnalysis}}}

\references{
Venables, W. N. and B. D. Ripley (2002). \emph{Modern Applied Statistics with S}, Springer-Verlag.  \cr

Bair, E., Hastie, T., Debashis, P., and Tibshirani, R. (2006). Prediction
by supervised principal components. \emph{American Statistics Association,
101(473):119-137}.

  }






\examples{

\donttest{
#--------------------------------  example I-----------------------------

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



results1<-SurFitPcaClasif(
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=100,
                        ProgFact=NULL,
                        Plots = TRUE,
                        MedianCut = NULL 
                              
)
summary(results1$SurFit)

plot(survfit(results1$SurFit,newdata=data.frame(gid=as.factor(c(1,2)))),col=2:3,
main="Fitted Survival Curves [PCA]",xlab="Time",ylab="Probability")


#--------------------------------  example II -----------------------------


#generate some Survival data
SimData<-GenSynSurvData(nPatients=100,nGenes=150,Pi=0.5)

SurvTime<-SimData$SurvTime
Censor<-SimData$Censor
ProgFact<-SimData$ProgFact
Gdata<-SimData$Gdata
featurenames<-SimData$featurenames
          
          
results2<-SurFitPcaClasif(
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=FALSE,
                        NuFeToBeSel=150,
                        ProgFact=ProgFact,
                        Plots = TRUE,
                        MedianCut = NULL 
                              
)
summary(results2$SurFit)

windows()
plot(survfit(results2$SurFit,newdata=data.frame(gid=as.factor(c(1,2)),
Age=20,Stage=rep(0.5,4),sex=rep(0.5,2))),col=2:3,
main="Fitted Survival Curves",xlab="Time",ylab="Probability")

}

}

\keyword{PCA, CoxPh, Classification}
