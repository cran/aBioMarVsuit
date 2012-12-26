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



# function do the classification, surv fit and HR estimation and plots ------------------------------------


EstimHR<-function(   RiskScore,
                     Sdata,
                     ProgFact=NULL,
                    Plots = FALSE,
                    MedianCut = NULL ){
options(warn=-1)
if (missing(Sdata)) stop("Argument 'Sdata' is missing...")
if (missing(RiskScore)) stop("Argument 'RiskScore' is missing...")

SurvTime<-Sdata[,1]

    p1<-RiskScore
    n.patients<-length(RiskScore)
         if (is.null(MedianCut)) { MedianCut <- median(p1) }  
        pos <- c(1:length(p1))[p1 <= MedianCut] 
        neg <- c(1:length(p1))[p1 > MedianCut] 
                
        gs <- NULL
        gs[pos] <- 1  
        gs[neg] <- -1   
        #c1data <- data.frame(cdata,gs=as.factor(gs))
          if (Plots) {
            windows()
            plot(p1, p1,ylab="Score",main="Low vs High Risk",ask=F)
            points(p1[gs==-1],p1[gs==-1],col="red")
            points(p1[gs==1],p1[gs==1],col="green")
            abline(MedianCut,0)
            # click on the plot to give a location for legend
            legend("bottomleft",col=c("red","green"),c("High Risk","Low Risk"),pch=1)
        }

#-- surviaval fit -
        gid<- NULL
        gid[c(1:n.patients)[gs==1]]<- "low"
        gid[c(1:n.patients)[gs==-1]]<- "high"
        c2data <- data.frame(Sdata,gid=as.factor(gid))
        
        
          if (Plots) {
          windows() 
                    par(mfrow=c(1,2),mar=c(5,4,4,1))
                    mfit <- survfit(Surv(SurvTime, Censor == 1)~gid, data = c2data)
                    plot(mfit,ylab="Survival",xlab="Time",col=c("red","green"),main="Kaplan-Meier curves ") 
                    legend(1,0.4,col=c("green","red"),c("Low Risk","High Risk"),lty=rep(1,2),bty="n")       
                    
                    boxplot(SurvTime~gid,ylim=range(SurvTime),ylab="Survival Time",main="Low vs High Risk groups",col=c("red","green"))
        }
        
       if (is.null(ProgFact)) { 
        SurFit<- coxph(Surv(SurvTime, Censor==1) ~ gid,data=c2data)
       } 
       
       if (is.data.frame(ProgFact)) {
       nPrgFac<-ncol(ProgFact)
       c2data <- data.frame(Sdata,gid=as.factor(gid),ProgFact)
       NameProg<-colnames(ProgFact)
       eval(parse(text=paste( "SurFit <-coxph(Surv(SurvTime, Censor==1) ~ gid",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=c2data)" ,sep="")))  
       }
       
       
return(list(SurFit=SurFit,gs=gs))      
}

# END OF function do the classification, surv fit and HR estimation and plots ------------------------------------
#----------------------------------------------------------------------------------------------------------
