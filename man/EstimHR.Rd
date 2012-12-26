\name{EstimHR}
\alias{EstimHR}
\title{Classification, survival fit and HR estimation}


\description{
The function does the classification based on the risk scores provided by the user and visualize survival fit along with HR estimates.
}

\usage{
EstimHR(RiskScore, Sdata, ProgFact = NULL, Plots = FALSE, MedianCut = NULL)
}


\arguments{
  \item{RiskScore}{Risk scores: a vector of size resembles to number of patients.} 
  \item{Sdata}{A data frame with survival time and censoring indications.}
  \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 
  \item{Plots}{A boolean parameter indicating the plots should be shown. Default is FALSE}
  \item{MedianCut}{Default is \code{NULL} and median cut off will be used. And user define cutoff value can be provided.}
}


\value{
  \item{SurFit}{Summary of the model fit with estimated HR.}
  \item{gs}{A vector with 1 (low risk)  and -1 (high risk) values and length is equal to number of patients.}
}



\author{Pushpike Thilakarathne and Ziv Shkedy}

\seealso{\code{\link{coxph}}, \code{\link{SurFitPlsClasif}}, \code{\link{SurFitPcaClasif}} }









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


#------------------------------------------------------------------
#load existing data files

NuFeToBeSel=150

DataForReduction<-list(x=Gdata,y=SurvTime, censoring.status=Censor, featurenames=rownames(Gdata))
TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),
decreasing =TRUE))[1:NuFeToBeSel]


ReduGdata<-Gdata[TentativeList,]

#generate some prognotic factors
nPatients<-ncol(ReduGdata)

ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))


#define some functions
  f.rvar <- function(x){apply(x,1,var)}

  f.rmean <- function(x, n = ncol(x)){c(x \%*\% rep(1, n))/n}

f.pca=function (x)
{
    ca <- match.call()
    if (ncol(x) > nrow(x)) {
        u = princomp(t(x))
        u$call = ca
        return(u)
    }
    xb <- x - (mu <- f.rmean(x))
    xb.svd <- svd(xb)
    u <- t(xb) \%*\% xb.svd$u
    dimnames(u)[[2]] <- paste("PC", 1:ncol(u), sep = "")
    l <- xb.svd$u
    dimnames(l) <- list(paste("V", 1:nrow(l), sep = ""), paste("Comp.",
        1:ncol(l), sep = ""))
    class(l) <- "loadings"
    sd = xb.svd$d/sqrt(ncol(x))
    names(sd) <- paste("Comp.", 1:length(sd), sep = "")
    u <- list(sdev = sd, loadings = l, center = mu, scale = rep(1,
        length(mu)), n.obs = ncol(x), scores = u, call = ca)
    class(u) <- "princomp"
    return(u)
}

#estimate risk score without prognotics factors
 
pc1 <- f.pca(as.matrix(ReduGdata))[[6]][,1]
cdata <- data.frame(SurvTime, Censor,pc1=pc1)
SurFit <- coxph(Surv(SurvTime, Censor==1) ~ pc1,data=cdata)
          
  TrtandPC1<-summary(SurFit)[[7]][c("pc1"),1]
  RiskScore1 <- TrtandPC1[1]*pc1

Sdata<-data.frame(SurvTime=SurvTime,Censor=Censor)

EstimHR(RiskScore1, Sdata,ProgFact=NULL,Plots = TRUE,MedianCut = NULL)
 
#estimate risk score WITH prognotic factors

cdata <- data.frame(SurvTime, Censor,pc1=pc1,ProgFact)
SurFit <- coxph(Surv(SurvTime, Censor==1) ~ pc1+Age+Stage+sex,data=cdata)
  TrtandPC1<-summary(SurFit)[[7]][c("pc1"),1]
  RiskScore2 <- TrtandPC1[1]*pc1

 EstimHR(RiskScore2, Sdata,ProgFact=ProgFact,Plots = TRUE,MedianCut = NULL)
}

}

\keyword{Survival, PCA, Risk, Hazard Ratio}
