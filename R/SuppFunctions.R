
f.rmean=function(x, n = ncol(x)){
c(x %*% rep(1, n))/n
 }
#-----------------------------------END OF function for ro mean ------------------------------------


hrDens<-function(HR,vv=1,mymain=""){

        r <- density(HR[,1],from=0,to=3) 
        plot(r,main=mymain,ylim=c(0,3),xlab="Estimated HR")

        abline(v=vv,col=2)
        #CI for permuated cases
        qq<-quantile(sort(HR[,1]),prob=c(0.05,0.5,0.95))
        abline(v=qq[1],col=3)
        abline(v=qq[2],col=3)
        abline(v=median(HR[,1]),col=3,lwd=5)
        
        pvalue<-sum(vv>HR[,1])/nrow(HR)
        return(list(pvalue,qq))
        
}


#---------------------------------------------------------------------------------------------------
#----------------------------------- function for PCA ----------------------------------------------
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
    u <- t(xb) %*% xb.svd$u
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


f.rvar <- function(x){apply(x,1,var)}
