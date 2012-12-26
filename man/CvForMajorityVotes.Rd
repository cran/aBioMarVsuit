\name{CvForMajorityVotes}
\alias{CvForMajorityVotes}
\title{Cross validation for majority votes}


\description{
  Performs cross validations for Majority votes based classification.
}

\usage{
CvForMajorityVotes(SurvTime,Censor, ProgFact=NULL, Gdata, 
                   ReduceDim=TRUE, NuFeToBeSel=150, fold=3, nCV=100)
}


\arguments{
\item{SurvTime}{A vector of survival time}
\item{Censor}{Censoring indicator}
 \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
 \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should 
              be equal to the number of genes and number of columns should be equal to number of patients.} 
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and reduced matrix will be used for cross validations.
  } 
  \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the cross validations. This is valid only if ReduceDim=\code{TRUE}.}
   \item{fold}{Number of folds. Defalut is 3. } 
  \item{nCV}{Number of cross validations. Default is 100.} 
}


\value{
A class of \code{\link{CVMV}} returns.
}




\author{Pushpike Thilakarathne and  Ziv Shkedy}


\seealso{\code{\link{MajorityVotes}}}


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

#without prognostic factors
cvmvres<-CvForMajorityVotes(SurvTime,Censor, ProgFact=NULL, 
Gdata, ReduceDim=TRUE, NuFeToBeSel=150, fold=3, nCV=50)
show(cvmvres)
plot(cvmvres)

#generate some prognotic factors
set.seed(1258)
ProgFact<-data.frame(Age=floor(abs(rnorm(ncol(Gdata),40,10))),
                     Stage=sample(1:4,ncol(Gdata),replace=T))

cvmvres1<-CvForMajorityVotes(SurvTime,Censor, ProgFact=ProgFact,
Gdata, ReduceDim=TRUE, NuFeToBeSel=150, fold=3, nCV=50)
show(cvmvres1)
plot(cvmvres1)
}

}

\keyword{CV, Majority votes}
