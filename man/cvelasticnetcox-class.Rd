\name{cvelasticnetcox-class}
\docType{class}
\alias{cvelasticnetcox-class}
\alias{cvelasticnetcox}
\alias{cvelasticnetcox-method}
\alias{show,cvelasticnetcox-method}
\alias{plot,cvelasticnetcox,missing-method}
\alias{summary,cvelasticnetcox-method}

\title{The cvelasticnetcox Class}

\description{
This class object contains estimated HR under testing and train data using \code{\link{CVLassoElasticNetCoxPh}}. }

\section{Slots}{
 A object of Class \code{cvelasticnetcox} returns \cr

\describe{
 \item{coef.mat}{A matrix with rows equals to number of cross validations and columns equals to number of genes.}
  \item{Run.Time}{A vector. Run time for each iteration and is given in seconds.}
  \item{lambda}{A vector. estimated optimum lambda for each iterations.}
  \item{n.g}{A vector. Number of non zero shrinkage coefficients at each iteration.}
  \item{gene.mat}{A matrix with 0 and 1. Number of rows equals to number of iterations and number of columns equals to number of genes.  1 indicates a particular gene had non zero coefficient and otherwise it is zero.}
  \item{HRT}{A matrix with three columns and first column represents the estimated HR  and last two gives lower and upper 95\% CI.  And this HR is estimated based on TRAINING data.}
  \item{HRTE}{A matrix with three columns and first column represents the estimated HR  and last two gives lower and upper 95\% CI.  And this HR is estimated based on TEST data.}
   \item{pld}{Partial likelihood deviance at each cross validations.}
   \item{ReduGdata}{Reduced gene expression matrix by supervised PCA. And this is the matrix used for cross validations.}
}
}

\usage{
\S4method{show}{cvelasticnetcox}(object)
\S4method{plot}{cvelasticnetcox,missing}(x, y, ptype=1, ...)
\S4method{summary}{cvelasticnetcox}(object)
}


\arguments{
  \item{x}{A \code{cvelasticnetcox} class object}
  \item{y}{missing}
   \item{ptype}{Plot type. 1 distribution of the HR under test and train data. 2 HR vs number genes selected. 3. Top mostly selected 30 genes.}
  \item{object}{A \code{cvelasticnetcox} class object}
  \item{\dots}{The usual extra arguments to generic plot functions.}
}


\section{Methods}{
  \describe{
    \item{show}{Use \code{show(cvelasticnetcox-object)} for brief information.}

    \item{plot}{Use \code{plot(cvelasticnetcox-object)}}       
    
    \item{summary}{Use \code{summary(cvelasticnetcox-object)}}        
     }
}


\author{Pushpike Thilakarathne and  Ziv Shkedy }


\seealso{\code{\link{CVLassoElasticNetCoxPh}}, \code{\link{coxph}}, \code{\link{glmnet}}, \code{\link{princomp}}}




\keyword{Permutations, Null Distribution, CoxPh}
