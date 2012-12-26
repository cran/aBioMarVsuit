\name{MajorityVotes}
\alias{MajorityVotes}
\title{Majority votes classification}


\description{
The function performs the classification based on majority votes and estimate the HR for low risk group. It also visualizes the gene-specific 
classifications for a set of (20) patients. Here 20 patients are considered for visualization purposes only. 
}

\usage{
MajorityVotes(gSpecificRes,ProgFact, SurvTime,Censor,J=1)
}


\arguments{
  \item{gSpecificRes}{object returned from function \code{\link{GeneSpecificCoxPh}}.}
  \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
 \item{SurvTime}{A vector of survival time}
  \item{Censor}{Censoring indicator}
 \item{J}{\code{j}th Set of patients (20) and default is \code{J=1}.  \code{J} should be less than \code{(no.of patients/20)}.}
}


\value{
A list with 

  \item{m0}{CoxPh Model fit object.}
  \item{ggr}{Majority vote classification for each patient.}
  \item{VoteMat}{Vote matrix returns from  \code{\link{GeneSpecificCoxPh}}.}
}




\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{GeneSpecificCoxPh}}, \code{\link{CvForMajorityVotes}}}

\references{
Trevor Hastie, Robert Tibshirani, and J. H. Friedman. The elements of statistical learning: 
data mining, inference, and prediction: with 200 full-color illustrations. 
\emph{New York}: Springer-Verlag, 2001. \cr

Bair, E., Hastie, T., Debashis, P., and Tibshirani, R. (2006). Prediction
by supervised principal components. \emph{American Statistics Association,
101(473):119-137}.

  }



\examples{
 data(exprLym)

GexprMatrix<-exprs(exprLym)
SurvData<-pData(exprLym)

PatietId<-rownames(SurvData[!is.na(SurvData[,c("IPI")]),])
Gdata<-GexprMatrix[,PatietId]

#colnames(Gdata)
SurvTime<-SurvData[!is.na(SurvData[,c("IPI")]),c("FollowUpYears")]
Censor<-ifelse(SurvData[!is.na(SurvData[,c("IPI")]),c("StatusAtFollowUp")]=="Dead",1,0)

Gdata[is.na(Gdata)]<-mean(Gdata,na.rm=TRUE)
          
Ana<-GeneSpecificCoxPh( SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL)
show(Ana)


MVres<-MajorityVotes(Ana,NULL, SurvTime,Censor,J=1)


}

\keyword{Majority votes}
