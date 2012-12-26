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

#---------------------------------------------------------------------------------------------------
#------------------------------------ function for  gene specific cox ph model ---------------------

# higher the risk lower the survival

GeneSpecificCoxPh<-function( SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        MedianCut = NULL){

if (missing(SurvTime)) stop("Argument 'SurvTime' is missing...")
if (missing(Gdata)) stop("Argument 'Gdata' is missing...")
if (missing(Censor)) stop("Argument 'Censor' is missing...")

if (ReduceDim) {
    DataForReduction<-list(x=Gdata,y=SurvTime, censoring.status=Censor, featurenames=rownames(Gdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:NuFeToBeSel]
    TentativeList
    
    ReduGdata<-Gdata[TentativeList,]
    } else {
    ReduGdata<-Gdata
}


n.genes<-nrow(ReduGdata)
n.patients<-ncol(ReduGdata)
gNames<-rownames(ReduGdata)  

    HRp <- matrix(0,n.genes,4)
    gr <- matrix(0, n.genes, n.patients)
    
    res <-  vector("list", n.genes)
    
    for (i in 1:n.genes){  #---------------------------  STRAT FOR LOOP ------------------------

                genei <- ReduGdata[i,]
                
                if (is.null(ProgFact)) {
                        
                cdata <- data.frame(SurvTime,Censor,genei)
                m0 <- coxph(Surv(SurvTime, Censor==1) ~ genei,data=cdata)
                }
                if (!is.null(ProgFact)) {
                    if (is.data.frame(ProgFact)) {
                    nPrgFac<-ncol(ProgFact)
                    cdata <- data.frame(SurvTime,Censor,genei,ProgFact)
                    NameProg<-colnames(ProgFact)
                    eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ genei",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
                    } else {
                             
                    stop(" Argument 'ProgFact' is NOT a data frame ")
                    }
                
                }
                #risk Score       
                TrtandPC1<-summary(m0)[[7]][c("genei"),1]
                p1 <- TrtandPC1*genei
              
              TempGenei <-EstimHR(p1,Sdata=cdata, ProgFact=ProgFact,Plots = FALSE, MedianCut = MedianCut )  
              res[[i]]<-TempGenei$SurFit     
              HRp[i,]<-summary(TempGenei$SurFit)[[8]][1,]
               gr[i,]<-TempGenei$gs
              
    }#---------------------------  END OF  FOR LOOP ------------------------

return(new("GeneSpecific",res=res,HRp=HRp,gr=gr,gNames=gNames))
}

#---------------------------------------------------------------------------------------------------
#------------------------------------END OF  function for  gene specific cox ph model --------------
