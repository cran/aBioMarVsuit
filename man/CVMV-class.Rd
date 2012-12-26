\name{CVMV-class}
\docType{class}
\alias{CVMV-class}
\alias{CVMV}
\alias{CVMV-method}
\alias{show,CVMV-method}
\alias{plot,CVMV,missing-method}

\title{The CVMV Class}

\description{
This class object contains cross validated  majority votes based classification and HR estimation results. It returns from \code{\link{CvForMajorityVotes}}. }

\section{Slots}{
 A object of Class \code{CVMV} returns \cr
\describe{

  \item{HRp.train}{A data frame with  estimated HR on train data for low risk group.}
  \item{HRp.test}{A data frame with  estimated HR on test data for low risk group.}
  \item{ReduGdata}{Reduced gene expression matrix.}
  \item{nCV}{No.of cross validations.}
  \item{pFactors}{Name(s) of the prognostic factor(s).} 
}
}

\usage{
\S4method{plot}{CVMV,missing}(x, y, ...)
\S4method{show}{CVMV}(object)
}


\arguments{
  \item{x}{A \code{CVMV} class object}
  \item{y}{missing}
  \item{object}{A \code{CVMV} class object}
  \item{\dots}{The usual extra arguments to generic plot functions.}
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CVMV-object)} for brief information.}

    \item{plot}{Use \code{plot(CVMV-object)}}              
     }
}


\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{
\code{\link{CVSeqIncreaseGenes}},
\code{\link{CVGeneSpecificCoxPh}},
\code{\link{CvForMajorityVotes}},
\code{\link{MajorityVotes}}
}



\keyword{CV, CoxPh}
