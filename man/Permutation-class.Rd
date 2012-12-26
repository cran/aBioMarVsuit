\name{Permutation-class}
\docType{class}
\alias{Permutation-class}
\alias{Permutation}
\alias{Permutation-method}
\alias{show,Permutation-method}
\alias{plot,Permutation,missing-method}
\alias{summary,Permutation-method}

\title{The Permutation Class}

\description{
This class object contains estimated null distribution of the HR for a given validation scheme in function \code{\link{NullDistHR}}. }

\section{Slots}{
 A object of Class \code{Permutation} returns \cr
\describe{

  \item{HRlowObs}{Estimated HR for low risk group on the original data}
  \item{HRlowPerm}{A data frame of estimated HR for low risk group on the permuted data.}
  \item{n.perm}{Number of permutations used.} 
  \item{Validation}{The validation scheme used.}
}
}

\usage{
\S4method{plot}{Permutation,missing}(x, y, ...)
\S4method{show}{Permutation}(object)
\S4method{summary}{Permutation}(object)
}


\arguments{
  \item{x}{A \code{Permutation} class object}
  \item{y}{missing}
  \item{object}{A \code{Permutation} class object}
  \item{\dots}{The usual extra arguments to generic plot functions.}
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(Permutation-object)} for brief information.}

    \item{plot}{Use \code{plot(Permutation-object)}}       
    
    \item{summary}{Use \code{summary(Permutation-object)}}        
     }
}


\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{\code{\link{NullDistHR}},
\code{\link{EstimHR}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, 
 \code{\link{MajorityVotes}},  \code{\link{LassoElasticNetCoxPh}}}




\keyword{Permutations, Null Distribution, CoxPh}
