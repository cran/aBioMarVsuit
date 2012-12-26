\name{NullDistHR}
\alias{NullDistHR}
\title{Null Distribution of the Estimated HR}


\description{
This function generates the null distribution of the HR by permuting the observations.  Several ways of permutation setting are implemented. That is, function can be used to generate null distributions for four different  validation schemes, PLS based, PCA based, Majority votes based and Lasso based.
 Note this function internally calls function \code{\link{EstimHR}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, 
 \code{\link{MajorityVotes}}, and \code{\link{LassoElasticNetCoxPh}}.  
 It can handel larger gene expression matrices and they are reduced by supervised pca approach. Reduced matrix will be used for all the permutations.
 
}

\usage{
NullDistHR(n.perm=100, case=2, Validation=c("PLSbased","PCAbased","L1based","MVbased"), 
          SurvTime,  Gdata,  Censor, ReduceDim=TRUE, NuFeToBeSel=150, ProgFact=NULL,
            MedianCut = NULL)
}


\arguments{
  \item{n.perm}{Number of permutations to be used and default 100.}  
  \item{case}{There are seven different ways how to call this function:
    \enumerate{
     \item permute survival only.
    \item permute survival and rows of data frame of the prognostic factors.
    \item permute survival,  rows of data frame of the prognostic factors, columns of  genes expression matrix independently.
    \item permute genes expression matrix only.
    }
    case \code{2} is set as default one. 
}

\item{Validation}{There are four different validation schemes where the null distribution can be estimated. That is \code{c("PLSbased","PCAbased","L1based","MVbased")}. }

 \item{SurvTime}{A vector of survival time.}
  \item{Gdata}{A large or small gene expression matrix. Number columns resembles to number of patients.  A matrix with gene expressions where the number of rows should 
              be equal to the number of genes and number of columns should be equal to number of patients.} 
  \item{Censor}{Censoring indicator.}
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix 
  and to be used in the classifier. Note that use can directly used the reduced gene expression matrix by bypassing the supervised pca using \code{ReduceDim=FALSE}. } 
  \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the PCA or PLS based methods. This is valid only if ReduceDim=\code{TRUE}.}
    \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
    \item{MedianCut}{Default is \code{NULL} and median cut off will be used. And user define cutoff value can be provided.}
 
    
}

\value{
A object of Class \code{Permutation} returns together with estimated HR for low risk group on the permuted data and original data.
}




\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{EstimHR}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}}, 
 \code{\link{MajorityVotes}},  \code{\link{LassoElasticNetCoxPh}}}


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


#------- permuations Sampling --------------------------------------

#permute survival only and Validation scheme: PLS based
set.seed(123)
Perm<-NullDistHR(n.perm=100,case=2, Validation="PLSbased", SurvTime, Gdata,
                        Censor,  ReduceDim=TRUE,  NuFeToBeSel=150,
                        ProgFact=NULL,  MedianCut = NULL )

show(Perm)
plot(Perm)



#permute survival only and Validation scheme: PCA based
set.seed(123)
NullDistHR(n.perm=100,case=2, Validation="PCAbased", SurvTime, Gdata,
                        Censor,  ReduceDim=TRUE,  NuFeToBeSel=150,
                        ProgFact=NULL,  MedianCut = NULL )
                        
                        

#permute survival only and Validation scheme: Majority votes based
set.seed(123)
NullDistHR(n.perm=100,case=2, Validation="MVbased", SurvTime, Gdata,
                        Censor,  ReduceDim=TRUE,  NuFeToBeSel=150,
                        ProgFact=NULL,  MedianCut = NULL )



#permute survival only and Validation scheme: Lasso based
set.seed(123)
NullDistHR(n.perm=100,case=2, Validation="L1based", SurvTime, Gdata,
                        Censor,  ReduceDim=TRUE,  NuFeToBeSel=150,
                        ProgFact=NULL,  MedianCut = NULL )

    }
}

\keyword{permutation, PCA, PLS, CoxPh}
