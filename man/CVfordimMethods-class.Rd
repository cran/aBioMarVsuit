\name{CVfordimMethods-class}
\docType{class}
\alias{CVfordimMethods-class}
\alias{CVfordimMethods}
\alias{CVfordimMethods-method}
\alias{show,CVfordimMethods-method}
\alias{plot,CVfordimMethods,missing-method}
\alias{summary,CVfordimMethods-method}

\title{The CVfordimMethods Class}
\description{Object returned by the method \code{\link{CVforDimPcaPls}}.}
\section{Slots}{
     \describe{
   
  \item{Results}{A dataframe with two columns where the estimated HR of low risk group for train and test data are given.} 
  \item{n}{Number of cross valdiations  used.}
  \item{PCAorPLS}{Method of dimention reduction used.}
  \item{cv.train}{Train indeces matrix used for cross validations} 
  \item{NuFeToBeSel}{Number of genes have been considered for PCA or PLS method.} 
  }
}


\usage{
\S4method{plot}{CVfordimMethods,missing}(x, y, ...)
\S4method{summary}{CVfordimMethods}(object)
\S4method{show}{CVfordimMethods}(object)
}


\arguments{
  \item{x}{A \code{CVfordimMethods} class object}
  \item{y}{missing}
   \item{object}{A \code{CVfordimMethods} class object}
  \item{\dots}{The usual extra arguments to generic functions}
}




\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CVfordimMethods-object)} for brief information.}
    \item{summary}{Use \code{summary(CVfordimMethods-object)} to apply the
                   classic \code{summary()} function to the  \code{CVfordimMethods-object)} }
    \item{plot}{Use \code{plot(CVfordimMethods-object)}}
     }
}

\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{CVforDimPcaPls}}}

\keyword{Cross Validations, PCA, PLS, CoxPh}
