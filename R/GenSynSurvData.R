####**********************************************************************
####**********************************************************************
####
####  cross validated gene by gene analysis
####
####  Copyright 2012, CenStat, Uhasselt
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Pushpike J Thilakarathne, PhD.
####    Interuniversity Institute for Biostatistics and statistical Bioinformatics (I-BioStat) 
####    University of Hasselt and Catholic university of Leuven
####    3000 Leuven
####    Belgium
####
####    email:  pushpike@gmail.com
####    URL:    http://www.uhasselt.be/fiche?voornaam=Pushpike&naam=Thilakarathne
####
####    
####**********************************************************************
####**********************************************************************

########################################################################

GenSynSurvData<-function(nPatients=100,nGenes=150,Pi=0.5){

set.seed(425)
Gdata<-matrix(rnorm(nGenes*nPatients),ncol=nPatients)
Gdata[1:10,]<-matrix(rnorm(10*nPatients,sd=1,mean=3),ncol=nPatients)
Gdata[11:50,]<-matrix(rnorm(40*nPatients,sd=1,mean=2),ncol=nPatients)

SurvTime<-c(abs(10+svd(Gdata[1:floor(nPatients*Pi),])$v[1:floor(nPatients*Pi),1]+ 
                2*rnorm(floor(nPatients*Pi))),
                abs(8 +svd(Gdata[1:floor(nPatients*Pi),])$v[(floor(nPatients*Pi)+1):nPatients,1]+ 
                3*rnorm(nPatients-floor(nPatients*Pi))))
                
CenTime<-abs(10+svd(Gdata[1:floor(nPatients),])$v[1:floor(nPatients),1]+ 
                1*rnorm(floor(nPatients)))                
                
Censor<- ifelse(SurvTime>CenTime,0,1)

featurenames <- paste("feature",as.character(1:nGenes),sep="")
rownames(Gdata)<-featurenames 

ProgFact<-data.frame(Age=floor(SurvTime*0.68+rnorm(nPatients,30,10)),
          Stage=sample(1:4,nPatients,replace=T),sex=rbinom(nPatients, 1, 0.5))


return(list(Censor=Censor,SurvTime=SurvTime,featurenames=featurenames,Gdata=Gdata,ProgFact=ProgFact))
}
