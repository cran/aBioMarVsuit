\name{LassoElasticNetCoxPh}
\alias{LassoElasticNetCoxPh}
\title{Wapper function for glmnet}


\description{
This is a wrapper function for \code{glmnet} and fit model with all prognostic factors and genes. LASSO, Elastic net and Ridge regressions can be fitted. 
Optimum lambda will be used to select the non-zero shrinkage coefficients. The function can accommodates prognostic factors and they will be forced to be in the model if they are supplied. 
}

\usage{
LassoElasticNetCoxPh(SurvTime, Censor, ReduGdata,  ProgFact=NULL,  Plots = FALSE,
                        MedianCut = NULL, GeneList=NULL,  StZ=TRUE, alpha=1)
}


\arguments{
  \item{SurvTime}{A vector of survival time.}
  \item{Censor}{Censoring indicator.}
  \item{ReduGdata}{A reduced gene expression matrix by supervised PCA approach. } 
  \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{Plots}{A boolean parameter indicating the plots should be shown. Default is FALSE}
  \item{MedianCut}{A user defined cut-off value for the classifier and default is median of the risk scores.} 
  \item{GeneList}{A list of genes to be considered in the model usually smaller than the genes in the \code{ReduGdata} . Default is to use all genes available in \code{ReduGdata}.} 
  \item{StZ}{Logical flag for gene expression standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.}
  \item{alpha}{mixing parameter for glmnet (see \code{\link{glmnet}}) with \code{0<= alpha <= 1}.} 
}


\value{
  \item{RiskScore}{a vector of numeric values. That is estimated risk scores based on selected genes for each patients.}
  \item{Results}{Fitted CoxPh model object for low risk and high risk patients.}
  \item{betaS}{A vector of estimated non-zero shrinkage coefficients.}
  \item{gm}{list of genes selected.}
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{coxph}}, \code{\link{glmnet}}, \code{\link{princomp}}}

\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2009). glmnet: Lasso and
elastic-net regularized generalized linear models. \emph{R package version 1.8-5}. \cr


Hastie, T. and Efron, B. (2007). lars: Least Angle Regression, Lasso and
Forward Stagewise. \emph{R package version 0.9-7}. \cr

Trevor Hastie, Robert Tibshirani, and J. H. Friedman (2001). The elements of statistical learning: data mining, 
inference, and prediction: with 200 full-color illustrations. 
\emph{New York}: Springer-Verlag. \cr

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. \emph{Jour-
nal of the Royal Statistics Society}, 58:267-288. \cr

Zou, H. and Hastie, T. (2005). Regularization and variable selection via the
elastic net. \emph{Journal of the Royal Statistical Society Series B}, 67:301-320. \cr

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

#reduced gene expression matrix by superpc
NuFeToBeSel<-150

DataForReduction<-list(x=Gdata,y=SurvTime, censoring.status=Censor, 
featurenames=rownames(Gdata))
TentativeList<-names(sort(abs(superpc.train(DataForReduction, 
type="survival")$feature.scores),decreasing =TRUE))[1:NuFeToBeSel]
#TentativeList
               
ReduGdata<-Gdata[TentativeList,]

# application of Lasso with all genes
 set.seed(145)
LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata, ProgFact=NULL, Plots = FALSE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=1)
                        
# application of Lasso with some prognostic factors

#generate prognostic factors
nPatients<-ncol(ReduGdata)
ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))
          
LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata, ProgFact=ProgFact, Plots = TRUE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=1)
                        
# application of Lasso with 50 genes
LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata[1:50,], ProgFact=NULL, Plots = TRUE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=1)

# application of Ridge 
LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata, ProgFact=NULL, Plots = TRUE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=0.01)

# application of ElasticNet with all genes
LassoElasticNetCoxPh(SurvTime,Censor,ReduGdata, ProgFact=NULL, Plots = TRUE, 
                     MedianCut = NULL , GeneList=NULL,
                        StZ=TRUE,  alpha=0.3)

  }
}

\keyword{Lasso, ElasticNet, CoxPh}
