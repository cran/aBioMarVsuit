\name{GeneSpecific-class}
\docType{class}
\alias{GeneSpecific-class}
\alias{GeneSpecific}
\alias{GeneSpecific-method}
\alias{show,GeneSpecific-method}
\alias{plot,GeneSpecific,missing-method}
\alias{summary,GeneSpecific-method}

\title{The GeneSpecific Class}
\description{Object returned by function \code{\link{GeneSpecificCoxPh}}.}
\section{Slots}{
     \describe{
  \item{res}{A list of dataframes of each having output object of coxph for a certain gene.} 
  \item{HRp}{A dataframe with estimated gene-specific HR for low risk group and 95\% CI.}
  \item{gr}{A matrix of group labels (1: low risk and -1: high risk) with number of rows equals to number of genes used and number of columns equals to number of patients.}
  \item{gNames}{Genes names used.} 
  }
}


\usage{
\S4method{plot}{GeneSpecific,missing}(x, y, ...)
\S4method{summary}{GeneSpecific}(object)
\S4method{show}{GeneSpecific}(object)
}


\arguments{
  \item{x}{A \code{GeneSpecific} class object}
  \item{y}{missing}
   \item{object}{A \code{GeneSpecific} class object}
  \item{\dots}{The usual extra arguments to generic functions --- see 
               \code{\link[graphics]{plot}}, 
               \code{\link[graphics]{plot.default}}, 
               \code{\link[stats]{plot.stepfun}} }
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(GeneSpecific-object)} for brief information.}
    \item{summary}{Use \code{summary(GeneSpecific-object)} to apply the
                   classic \code{summary()} function to the  \code{GeneSpecific-object)} }
    \item{plot}{Use \code{plot(GeneSpecific-object)}}              
     }
}





\details{
\describe{
\item{plot}{\code{signature(x = "GeneSpecific", y = "missing")}: plots  
   of class GeneSpecific for cross validated PCA and PLS based methods to use classify patients. }
\item{plot}{\code{signature(x = "GeneSpecific", y = "missing")}: Plots for Gene-specific analysis results. }

}

Any parameters of \code{plot.default} may be passed on to this particular
\code{plot} method. 

}









\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{GeneSpecificCoxPh}}}

\keyword{Gene by Gene, Cross Validations, PCA, PLS, CoxPh}
