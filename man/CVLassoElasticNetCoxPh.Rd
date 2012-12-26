\name{CVLassoElasticNetCoxPh}
\alias{CVLassoElasticNetCoxPh}
\title{Cross Validations for Lasso Elastic Net based predictive models}


\description{
The function performs the cross validations for LASSO and Elastic net models for Cox proportional hazard model. 
Top genes selected are being updated at each iteration and use in the classifier.  That means predictive gene signature is varied iteration to iteration.  The underline idea is to investigate the HR
under train and test data as well as mostly selected genes during this process. 
}

\usage{
CVLassoElasticNetCoxPh(n.cv=20,fold=3,SurvTime, Censor, Gdata,ReduceDim=TRUE, 
NuFeToBeSel=150,  ProgFact=NULL, GeneList=NULL,  StZ=TRUE, alpha=1)
}


\arguments{
  \item{n.cv}{Number of validations and default is 20.}
  \item{fold}{Number of folds to be used and default is 3.} 
    \item{SurvTime}{A vector of survival time.}
      \item{Censor}{Censoring indicator.}
  \item{Gdata}{A large or small gene expression matrix. A matrix with gene expressions where the number of rows should 
              be equal to the number of genes and number of columns should be equal to number of patients.} 
  \item{ReduceDim}{Default is \code{TRUE} and larger gene expression matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix and to be used in the classifier.
  } 
   \item{NuFeToBeSel}{Number of top K genes (default is 150) to be used in the PCA or PLS based methods. This is valid only if ReduceDim=\code{TRUE}.}
    \item{ProgFact}{A dataframe where the prognostic factors can be given to be used in the model.} 

  \item{GeneList}{A list of genes to be considered in the model usually smaller than the genes in the \code{Gdata}. If NULL all genes available in reduced gene expression matrix will be used. Otherwise only those genes are used in the predictive model.} 
  \item{StZ}{Logical flag for gene expression standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.}
  \item{alpha}{Mixing parameter for glmnet (see \code{\link{glmnet}}) with \code{0<= alpha <= 1}.} 
}


\value{
A class of \code{\link{cvelasticnetcox}} returns. }
 




\author{Pushpike Thilakarathne and  Ziv Shkedy}

\seealso{\code{\link{coxph}}, \code{\link{glmnet}}, \code{\link{princomp}}, \code{\link{InnerCrossValELNet}}, \code{\link{LassoElasticNetCoxPh}}}


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

# Example 1 :
#------------------------------------------------------------------
#Without Prognostic Factors
set.seed(123) 
CVres<-CVLassoElasticNetCoxPh(n.cv=50,fold=3,SurvTime, Censor, 
                              Gdata, NuFeToBeSel=100, ProgFact=NULL, 
                              ReduceDim=TRUE, GeneList=NULL,  StZ=TRUE, alpha=1)
                         

show(CVres)
summary(CVres)

#distribution of the HR
plot(CVres,ptype=1)

#estimated HR vs number of genes selected
plot(CVres,ptype=2)

#Mostly selected 30 genes
plot(CVres,ptype=3)


#With Prognostic Factors
set.seed(123) 
nPatients<-ncol(Gdata)
ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))


CVresP<-CVLassoElasticNetCoxPh(n.cv=10,fold=3,SurvTime, Censor, Gdata, 
                         NuFeToBeSel=150, ProgFact=ProgFact, 
                         ReduceDim=TRUE, GeneList=NULL,  StZ=TRUE, alpha=1)
summary(CVresP)


# Example 2 :
#------------------------------------------------------------------
#----------cross validations for Elastic Net ----------------------

# cross validations for Elastic Net with Alpha = 0.2
set.seed(123) 
ElNet<-CVLassoElasticNetCoxPh(n.cv=20,fold=3,SurvTime, 
  Censor, Gdata, NuFeToBeSel=150, ProgFact=NULL, 
  ReduceDim=TRUE, GeneList=NULL,  StZ=TRUE, alpha=0.2)

#Example 3
#------------------------------------------------------------------
# Apply Elastic Net on grid of Alpha values and cross validations

#This would take more than 5 hrs....
set.seed(123) 
grid.alpha<-seq(0.1,0.9,by=0.1)
CvElNetResults<-sapply(grid.alpha, 
               function(alpha) CVLassoElasticNetCoxPh(n.cv=50,fold=3,SurvTime, 
                                                    Censor, Gdata, NuFeToBeSel=100, ProgFact=NULL, 
                                                     ReduceDim=TRUE, GeneList=NULL,  StZ=TRUE, 
                                                     alpha=alpha))

Runtime<-sapply(1:length(grid.alpha),function(k) slot(CvElNetResults[[k]], "Run.Time"))
lambda <-sapply(1:length(grid.alpha),function(k) slot(CvElNetResults[[k]], "lambda"))
n.g    <-sapply(1:length(grid.alpha),function(k) slot(CvElNetResults[[k]], "n.g"))
HRT    <-sapply(1:length(grid.alpha), function(k) slot(CvElNetResults[[k]], "HRT")[,1] )
HRTE   <-sapply(1:length(grid.alpha), function(k) slot(CvElNetResults[[k]], "HRTE")[,1] )
HRTm   <-t(sapply(1:length(grid.alpha),function(k) quantile(slot(CvElNetResults[[k]], "HRT")[,1],
         na.rm=T,probs = c(0.025,0.5,0.975))))
HRTEm  <-t(sapply(1:length(grid.alpha),function(k) quantile(slot(CvElNetResults[[k]], "HRTE")[,1],
         na.rm=T,probs = c(0.025,0.5,0.975)))) 
freq   <-sapply(1:length(grid.alpha), function(k) colSums(slot(CvElNetResults[[k]], "gene.mat")))


# run time versus lambda

par(mfrow=c(3,3))
for (i in 1:length(grid.alpha)) {
Run.Time<-Runtime[,i]
lambda1<-lambda[,i]
plot(c(0,50),c(0,max(lambda1)),main= mtext(bquote(alpha == .(grid.alpha[i]))),
xlab="Run time (seconds)",
ylab=expression(lambda),cex.main=1.6,cex.lab=1.3,type="n")
points(Run.Time,lambda1,col=i+1)
}


#--Gene selection frequency across different alpha
par(mfrow=c(3,3))
for (i in 1:length(grid.alpha)) {
Run.Time<-Runtime[,i]
n.g1<-n.g[,i]
plot(c(0,length(n.g1)),c(0,max(n.g)),main= mtext(bquote(alpha == .(grid.alpha[i]))),
xlab="cv index",ylab="no.of Gene",cex.main=1.6,cex.lab=1.3,type="n")
lines(n.g1,col=i+1,type='s')
}


#----------HR 
par(mfrow=c(1,2))
colnames(HRT)<-grid.alpha
boxplot(HRT,ylim=c(0,max(HRT)),names=grid.alpha[1:length(grid.alpha)],
main="HR on Training Set \n 100 runs",col=2:(length(grid.alpha)+1),ylab="HR",
xlab=expression(alpha),cex.main=1.5,cex.lab=1.3)

colnames(HRTE)<-grid.alpha
boxplot(HRTE,ylim=c(0,max(HRTE)),names=grid.alpha[1:length(grid.alpha)],
main="HR on Testing Set \n 100 runs",col=2:(length(grid.alpha)+1),ylab="HR",
xlab=expression(alpha),cex.main=1.5,cex.lab=1.3)


# Estimated median HR on train and test data
#train
rownames(HRTm)<-grid.alpha
HRTm

#test
rownames(HRTEm)<-grid.alpha
HRTEm



#-- genes that are selected at least 50% times

par(mfrow=c(3,3),mar=c(1.75,1,1.5,0))    
  for (i in 1:length(grid.alpha)){
   Run.Time<-Runtime[,i]
  freq1<-freq[,i]

  freq1<-freq1[freq1>30]
  barplot(freq1,las=2,ylim=c(0,max(freq)),  ylab="",col=rainbow(length(freq1)),cex.names=0.35,
  main=mtext(bquote(alpha == .(grid.alpha[i]))),cex.lab=1,cex.main=1.5  )
  }


}
}


\keyword{CV, Lasso, ElasticNet, CoxPh}
