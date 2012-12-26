\name{CVGeneSpecificCoxPh}
\alias{CVGeneSpecificCoxPh}
\title{Cross validation for gene by gene Cox proportional hazard model}


\description{
This function performs the cross validation for gene by gene analysis. First data will be divided into train and test. Second, 
gene-specific model is fitted on train data and classifier is built.  And classifier is then evaluated on test data for that particular gene. Process is repeated for all genes. All these steps can be repeated many times (~1000).
}

\usage{
CVGeneSpecificCoxPh(fold=3, SurvTime, Gdata, Censor,  ReduceDim=TRUE,
                    NuFeToBeSel=150, ProgFact=NULL, MedianCut = NULL,  n.cv=3)
}


\arguments{
 \item{fold}{Number of folds. Default is 3.} 
  \item{SurvTime}{A vector of survival time}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should be equal to the number of genes and number of columns should be equal to number of patients.
               } 
  \item{Censor}{Censoring indicator}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix and to be used in the classifier.
  } 
    \item{NuFeToBeSel}{Number of top K genes (default is 150) to be selected from supervised PCA when there is a larger gene expression matrix. This is valid only if ReduceDim=\code{TRUE}.}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{MedianCut}{A user defined cut-off value for the classifier and default is median cut-off.} 
  \item{n.cv}{Number of cross validations. Default is 3 but it is recommended to have at least 100.} 
  
}


\value{ A object of Class \code{\link{CVGbyG}} returns.}



\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{\code{\link{CVSeqIncreaseGenes}}, 
\code{\link{GeneSpecificCoxPh}}, \code{\link{GridAnalysis}}}

\examples{
  
\donttest{
#generate some Survival data
SimData<-GenSynSurvData(nPatients=100,nGenes=1000,Pi=0.5)

SurvTime<-SimData$SurvTime
Censor<-SimData$Censor
ProgFact<-SimData$ProgFact
Gdata<-SimData$Gdata
featurenames<-SimData$featurenames


          
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
summary(CVR1)
summary(CVR1, Which=20)
plot(CVR1,ylim=c(0,5))
plot(CVR1, Which=20,ylim=c(0,5))                      
 }
}

\keyword{CV, CoxPh}
