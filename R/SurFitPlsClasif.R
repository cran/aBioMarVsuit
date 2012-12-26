#---------------------------------------------------------------------------------------------------
#----------------------------------- function for PCA and survival and classification ---------------



SurFitPlsClasif<-function(                                                        
                        SurvTime,
                        Gdata,
                        Censor,
                        ReduceDim=TRUE,
                        NuFeToBeSel=150,
                        ProgFact=NULL,
                        Plots = FALSE,
                        MedianCut = NULL 
                              
){

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


    if (is.matrix(ReduGdata)) {
        
                        PLSforGSK<-data.frame(1:n.patients)
                        PLSforGSK$g<-as.matrix(t(ReduGdata))
                        colnames(PLSforGSK)[1]<-c("SurvTime")
                        PLSforGSK[,1]<-SurvTime
                       # plsr.1 <- plsr(SurvTime ~ stdize(g), ncomp = 2, data = PLSforGSK, validation =  "CV")
                        plsr.1 <- plsr(SurvTime ~ g, method="simpls",ncomp =2, scale =TRUE,data = PLSforGSK, validation =  "CV")
                        pc1<-scores(plsr.1)[,1] # extract the first com
        } else {
        pc1<-ReduGdata
        }         
        
    if (is.null(ProgFact)) {
        
        cdata <- data.frame(SurvTime,Censor,pc1)
        m0 <- coxph(Surv(SurvTime, Censor==1) ~ pc1,data=cdata)
    }
    if (!is.null(ProgFact)) {
         if (is.data.frame(ProgFact)) { 
             nPrgFac<-ncol(ProgFact)
             cdata <- data.frame(SurvTime,Censor,pc1,ProgFact)
             NameProg<-colnames(ProgFact)
            eval(parse(text=paste( "m0 <-coxph(Surv(SurvTime, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))  
         } else {
         
         stop(" Argument 'ProgFact' is NOT a data frame ")
         }

    }
              
       #risk Score       
                TrtandPC1<-summary(m0)[[7]][c("pc1"),1]
                p1 <- TrtandPC1*pc1
     TempRes<-EstimHR(p1,Sdata=cdata, ProgFact=ProgFact,Plots = Plots, MedianCut = MedianCut )
        
      tempp<-list(SurFit=TempRes$SurFit,p1=p1,gs=TempRes$gs,pc1=pc1)
class(tempp)<-"SurFitPlsOut"
return(tempp)
}


#---------------------------------------------------------------------------------------------------
#--------------------------------- END OF function for PC and survival and classification ----------
