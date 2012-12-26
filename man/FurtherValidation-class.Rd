\name{FurtherValidation-class}
\docType{class}
\alias{FurtherValidation-class}
\alias{FurtherValidation}
\alias{FurtherValidation-method}
\alias{show,FurtherValidation-method}
\alias{plot,FurtherValidation,missing-method}
\alias{summary,FurtherValidation-method}


\title{The FurtherValidation Class}

\description{
This class object contains estimated  HR on the out of bag data by function \code{\link{InnerCrossValELNet}}. }

\section{Slots}{
 A object of Class \code{FurtherValidation} returns \cr
\describe{

 \item{RunTime}{A vector of running times in seconds for each iteration.}
  \item{fold}{Number of folds used.}
  \item{n.cv}{Number of outer cross validations used.}
  \item{n.innerCV}{Number of inner cross validations used.}
  \item{TopK}{Top genes used.}
 
 \item{HRInner}{A 3-way array in which first, second, and third dimensions correspond to n.innerCV, 1, and n.cv respectively.  This contains estimated HR for low risk group on Out of Bag data.}
  \item{HRTE}{A matrix with three columns and first column represents the estimated HR  and last two gives lower and upper 95\% CI. And this HR is estimated based on Out Of Bag data.}
  \item{WeightMat}{a matrix with columns equals number of TopK genes and rows as many as n.cv. Note that Weights are estimated as colMeans of coefficients matrix return from the inner cross validations.}

}
}



\usage{
\S4method{plot}{FurtherValidation,missing}(x, y, ptype=1,  ...)
\S4method{show}{FurtherValidation}(object)
\S4method{summary}{FurtherValidation}(object)
}


\arguments{
  \item{x}{A \code{FurtherValidation} class object}
  \item{y}{missing}
  \item{ptype}{Plot type. 1 distribution of the inner cross validated HR under test data for each outer iterations and estimated HR on the Out of Bag data are superimposed. 
  2 Estimated HR Density for low Risk Group .}
  \item{object}{A \code{FurtherValidation} class object}
  \item{\dots}{The usual extra arguments to generic plot functions.}
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(FurtherValidation-object)} for brief information.}

    \item{plot}{Use \code{plot(FurtherValidation-object)}}       
    
    \item{summary}{Use \code{summary(FurtherValidation-object)}}        
     }
}




\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{\code{\link{CVLassoElasticNetCoxPh}}, \code{\link{InnerCrossValELNet}},
\code{\link{EstimHR}}, \code{\link{MajorityVotes}},  \code{\link{LassoElasticNetCoxPh}}}


\keyword{Out of Bag Data, Cross Validations, Classifier}
