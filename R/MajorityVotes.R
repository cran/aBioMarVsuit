####**********************************************************************
####**********************************************************************
####
#### 
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








MajorityVotes<-function(gSpecificRes,ProgFact, SurvTime,Censor,J=1){
        
        
if (class(gSpecificRes)!="GeneSpecific") stop("Invalid class object.")
if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
if (missing(Censor)) stop("Argument 'Censor' is missing...")


        ggr <- per.R<-per.NR <- NULL
        for (i in 1:length(SurvTime)){
        per.R[i]<-sum(gSpecificRes@gr[,i]==1)
        ggr[i]<-ifelse((nrow(gSpecificRes@gr)-per.R[i])>per.R[i],"high","low")
        }
        
        ggr<-as.factor(ggr)

        if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime,Censor,ggr)
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ ggr,data=cdata)
        }
        
        if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime,Censor,ggr,ProgFact)
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ ggr",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
             } else {
             
             stop(" Argument 'ProgFact' is NOT a data frame ")
             }

        }

VoteMat<-gSpecificRes@gr
ng<-nrow(VoteMat)
np<-ncol(VoteMat)
Jmax<-floor(np/20)
if (J<Jmax) {
            slist<-(1+(J-1)*20):(J*20)
            
            plot(0,0, xlim=c(-2,ng),ylim=c(1,25),xlab="genes",ylab="Patient Index",axes=FALSE,type="n",
            main="Gene-Specific Classification\n of 20  selected patients",cex.main=0.9)
            for(i in 1:ng){
            points(rep(i, 20)[VoteMat[i,slist]==1],which(VoteMat[i,slist]==1),col=3,pch=15)
            points(rep(i, 20)[VoteMat[i,slist]==-1],which(VoteMat[i,slist]==-1),col=2,pch=15)
            }
            axis(1);axis(2,at=1:20,slist);box()
            legend("topleft",c("Low Risk","High Risk"),pch=c(15,15),col=c(3,2))
            } else {stop("J should be less than (no.of patients / 20) ")
}

return(list(m0=m0,ggr=ggr,VoteMat=VoteMat))
                       
}
