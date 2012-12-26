\name{CVSeqInc-class}
\docType{class}
\alias{CVSeqInc-class}
\alias{CVSeqInc}
\alias{CVSeqInc-method}
\alias{show,CVSeqInc-method}
\alias{plot,CVSeqInc,missing-method}
\alias{summary,CVSeqInc-method}

\title{The CVSeqInc Class}

\description{
This class object contains cross validated results of sequentially increasing of Top K genes analysis. It returns from \code{\link{CVSeqIncreaseGenes}}. }

\section{Slots}{
A object of Class \code{CVSeqInc} returns.
\describe{ 

  \item{HRPC}{3-way array in which first, second, and third dimensions correspond to n.gene, 3, and n.cv respectively. This contains the estimated HR on test data and dimension reduction method is PCA.}
  \item{HRPL}{3-way array in which first, second, and third dimensions correspond to n.gene, 3, and n.cv respectively. This contains the estimated HR on test data and dimension reduction method is PLS.}
  \item{n.genes}{Number of genes in the reduced matrix.}
  \item{n.cv}{No.of cross validations.}
  \item{top}{A sequence of top k genes considered. Default is \code{top=seq(5,100,by=5)}} 
}
}

\usage{
\S4method{plot}{CVSeqInc,missing}(x, y, ...)
\S4method{summary}{CVSeqInc}(object)
\S4method{show}{CVSeqInc}(object)
}


\arguments{
  \item{x}{A \code{CVSeqInc} class object}
  \item{y}{missing}
  \item{object}{A \code{CVSeqInc} class object}
  \item{\dots}{The usual extra arguments to generic functions --- see 
               \code{\link[graphics]{plot}}, 
               \code{\link[graphics]{plot.default}} }
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(CVSeqInc-object)} for brief information.}
    \item{summary}{Use \code{summary(CVSeqInc-object)} to apply the
                   classic \code{summary()} function to the  \code{CVSeqInc-object)} }
    \item{plot}{Use \code{plot(CVSeqInc-object)}}              
     }
}


\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{
\code{\link{CVSeqIncreaseGenes}},
\code{\link{CVGeneSpecificCoxPh}}
}



\keyword{CV, CoxPh}
