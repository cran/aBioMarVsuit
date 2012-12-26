\name{GeneSpecificCoxPh}
\alias{GeneSpecificCoxPh}
\title{Gene by gene Cox proportional hazard model}


\description{
 This function fits gene by gene Cox proportional hazard model and perform the classification based on median risk score which has been estimated using a signle gene. Function is useful for majority 
 vote classification method and gene by gene analysis and also for top K genes. 
}

\usage{
GeneSpecificCoxPh(SurvTime,
                 Gdata,
                 Censor,
                 ReduceDim=TRUE,
                 NuFeToBeSel=150,
                 ProgFact=NULL,
                 MedianCut = NULL)
}


\arguments{
  \item{SurvTime}{A vector of survival time}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should be equal to the number of genes and number of columns should be equal to number of patients.
               } 
  \item{Censor}{Censoring indicator}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced 
  matrix and to be used in the classifier.
  } 
    \item{NuFeToBeSel}{Number of top K genes (default is 150) to be selected from supervised PCA. This is valid only if ReduceDim=\code{TRUE}.}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{MedianCut}{A user defined cut-off value for the classifier and default is median cut-off.} 
}


\value{
A object of class \code{\link{GeneSpecific}} 
}





\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{coxph}}, \code{\link{GridAnalysis}},  \code{\link{MajorityVotes}}}





\references{
Richard Simon, Roadmap for Developing and Validating
Therapeutically Relevant Genomic Classifiers. 
\emph{JOURNAL OF CLINICAL ONCOLOGY VOLUME 23 d NUMBER 29 d OCTOBER 10 2005}.\cr

Yasrebi H, Sperisen P, Praz V, Bucher P (2009).
 Can Survival Prediction Be Improved By Merging Gene Expression Data Sets? 
 \emph{PLoS ONE 4(10): e7431}. 
 \url{doi:10.1371/journal.pone.0007431}
  }




\examples{
  
\donttest{
#generate some Survival data

SimData<-GenSynSurvData(nPatients=100,nGenes=1000,Pi=0.5)

SurvTime<-SimData$SurvTime
Censor<-SimData$Censor
ProgFact<-SimData$ProgFact
Gdata<-SimData$Gdata



#example I------------       
Ana1<-GeneSpecificCoxPh( SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL)

show(Ana1)
plot(Ana1)

#example II------------
ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))
          
Ana2<-GeneSpecificCoxPh( SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=ProgFact,
                        MedianCut = NULL)
show(Ana2)
plot(Ana2)
summary(Ana2)


}
}

\keyword{CoxPh}
