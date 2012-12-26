\name{CVGbyG-class}
\docType{class}
\alias{CVGbyG-class}
\alias{CVGbyG}
\alias{CVGbyG-method}
\alias{show,CVGbyG-method}
\alias{plot,CVGbyG,missing-method}
\alias{summary,CVGbyG-method}

\title{The CVGbyG Class}

\description{
This class object contains cross validated  gene by gene analysis results. It returns from \code{\link{CVGeneSpecificCoxPh}} }

\section{Slots}{
A object of Class \code{CVGbyG} returns.
\describe{ 

  \item{HRT}{3-way array in which first, second, and third dimensions correspond to n.gene, 4, and n.cv respectively. This contains the estimated HR for low risk on train data.}
  \item{HRTE}{3-way array in which first, second, and third dimensions correspond to n.gene, 4, and n.cv respectively. This contains the estimated HR for low risk on test data.}
  \item{ind.train}{A matrix which represents train indices.}
  \item{ind.test}{A matrix which represents test indices.}
  \item{n.genes}{Number of genes in the reduced matrix.}
  \item{n.cv}{No.of cross validations.}
  \item{ReduGdata}{A reduced gene expression matrix used in the analysis.} 
}
}

\usage{
\S4method{plot}{CVGbyG,missing}(x, y, Which=1, ...)
\S4method{summary}{CVGbyG}(object, Which=1)
\S4method{show}{CVGbyG}(object)
}


\arguments{
  \item{x}{A \code{CVGbyG} class object}
  \item{y}{missing}
  \item{Which}{Which gene for which estimated HR to be visualized. By default results of the first gene is used.}
   \item{object}{A \code{CVGbyG} class object}
  \item{\dots}{The usual extra arguments to generic functions --- see 
               \code{\link[graphics]{plot}}, 
               \code{\link[graphics]{plot.default}}, 
               \code{\link[stats]{plot.stepfun}} }
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CVGbyG-object)} for brief information.}
    \item{summary}{Use \code{summary(CVGbyG-object)} to apply the
                   classic \code{summary()} function to the  \code{CVGbyG-object)} }
    \item{plot}{Use \code{plot(CVGbyG-object)}}              
     }
}


\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{\code{\link{CVGeneSpecificCoxPh}}}



\keyword{CV, CoxPh}
