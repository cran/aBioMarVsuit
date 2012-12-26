GeneOccurence<-function(ObjectCv,TopK=20,minFreq=5){

Decrease=FALSE
if (class(ObjectCv)!="CVGbyG") stop("Invalid class object.")

    gene <- ObjectCv@ReduGdata
    gene.mat<-matrix(0,ObjectCv@n.cv,ObjectCv@n.genes)
    genes.names <- rownames(gene)
for (i in 1:ObjectCv@n.cv){ #i=1

#genes are ranked based on testing data
            hrset<-ObjectCv@HRTE[,,i] #test
            hrsetm<-ObjectCv@HRT[,,i] # train 
            
            Names.KGenes<-rownames(gene)
            #index.Top.KGenes <- which(hrset[,4]<1)
            
            #if ( !is.null(dim(index.Top.KGenes))) {
            #Top.KGenes.GSplus<-data.frame(rownames(gene)[index.Top.KGenes],hrset[index.Top.KGenes,-2])
            
            Top.KGenes.GSplus<-data.frame(rownames(gene),hrset[,-2])
            colnames(Top.KGenes.GSplus)<-c("Gene","HR","LCI","UCI")
            #Top.KGenes.GSplus
            
            
            Top.KGenes.GSminus<-data.frame(rownames(gene),hrsetm[,-2])
            
            colnames(Top.KGenes.GSminus)<-c("Gene","HR","LCI","UCI")
            #Top.KGenes.GSminus
            
            
            # FDR corrected CI for top k  genes
            
            cilevel <- 1-0.05*nrow(Top.KGenes.GSplus)/ObjectCv@n.genes  
            #qnorm(cilevel)
            #2.524987
            
            HRpadj <- exp(log(hrset[,1]) + log(hrset[,c(3,4)])-log(hrset[,1])*qnorm(cilevel)/1.96)  # 
            res.topkgenes<-data.frame(Top.KGenes.GSplus,HRpadj,hrsetm[,1])
            
            colnames(res.topkgenes)<-c("Gene","HRTest","LCI","UCI","FDRLCI","FDRUCI","HRTrain")
            #res.topkgenes
            #order(res.topkgenes[,c("FDRUCI")])
            sort.topK<-res.topkgenes[order(res.topkgenes[,c("HRTest")],decreasing=Decrease),]
            
            topg20<-sort.topK[1:TopK,c("Gene")]
                    #topg20
            gene.mat[i,is.element(genes.names,topg20)]<-1



}

fr<-colSums(gene.mat)
names(fr)<-genes.names

if (max(fr)<minFreq) minFreq<-max(fr)-1

top.fr<-fr[fr>minFreq]
barplot(top.fr,las=2,ylim=c(0,ObjectCv@n.cv),  ylab="",col=1:length(top.fr),cex.names=0.6,main="Mostly selected as top",cex.lab=1,cex.main=1.5  )

return(top.fr)

}
