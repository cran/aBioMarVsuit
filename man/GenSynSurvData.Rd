\name{GenSynSurvData}
\alias{GenSynSurvData}
\title{Generate Synthetic Survival Gene Expression Data }


\description{
The function generates the gene expression matrix where small set of genes (50) are informative and rest of them are set as noisy genes.  Next to that Survival time and Censoring information are generated 
based on first right singular vectors of \code{\link{svd}} of the gene expression matrix. It also generates other prognostic factors such as Age, Stage and sex which slightly correlated with survival time.
}

\usage{
GenSynSurvData(nPatients=100,nGenes=150,Pi=0.5)
}


\arguments{
  \item{nPatients}{Number of patients.} 
  \item{nGenes}{Number of Genes}
  \item{Pi}{A proportion of patients having low risk.} 
}


\value{
A list with following items \cr
  \item{Censor}{Event indicator.}
  \item{SurvTime}{Survival time.}
  \item{Gdata}{A gene expression matrix.}
   \item{featurenames}{A list of genes.}
   \item{ProgFact}{A data frame with prognostic factors.}
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{
\code{\link{coxph}}, 
\code{\link{SurFitPlsClasif}}, 
\code{\link{SurFitPcaClasif}},
\code{\link{GeneSpecificCoxPh}},
\code{\link{EstimHR}},
\code{\link{CVGeneSpecificCoxPh}},
\code{\link{CVforDimPcaPls}},
\code{\link{GridAnalysis}}
 }


\examples{
#generate some Survival data

SimData<-GenSynSurvData(nPatients=100,nGenes=150,Pi=0.5)

SurvTime<-SimData$SurvTime
Censor<-SimData$Censor
ProgFact<-SimData$ProgFact
Gdata<-SimData$Gdata
featurenames<-SimData$featurenames


}

\keyword{Survival, PCA, Risk, Hazard Ratio}
