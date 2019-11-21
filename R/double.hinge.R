# lb.quantile=.05; ub.quantile=.95; mid.x=NULL; lower.y=0; upper.y=100
double.hinge=function(x, y, mid.x=NULL, lb.quantile=.05, ub.quantile=.95, lower.y=0, upper.y=100, 
    var.type=c("none","bootstrap"), ci.bootstrap.size=1000, alpha=0.05, save.boot=TRUE, ncpus=1
) {

    var.type<-match.arg(var.type)    
    
    coef.hat=double.hinge.fit (x, y, mid.x, lb.quantile, ub.quantile, lower.y, upper.y)
    
    res=list()
    res$x=x
    res$y=y
    res$coefficients=coef.hat
    res$lower.y=lower.y
    res$upper.y=upper.y
    
    #### bootstrap
    if (var.type=="bootstrap") {
        # save rng state before set.seed in order to restore before exiting this function
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
        set.seed(1)
        
        boot.out=boot::boot(x, R=ci.bootstrap.size, sim = "ordinary", stype = "i", parallel = ifelse(ncpus==1,"no","multicore"), ncpus=ncpus, statistic=function(dat, ii){
            double.hinge.fit (dat[ii], y[ii], mid.x, lb.quantile, ub.quantile, lower.y, upper.y)               
        }) 
        boot.samples=boot.out$t
        #str(boot.out$boot.samples)
        
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)             
        
        # symmetric percentile CI as defined in Hansen (2017)
        ci.boot.symm=sapply (1:length(coef.hat), function(i) {
            q.1=quantile(abs(boot.samples[,i]-coef.hat[i]), 1-alpha, na.rm=TRUE)
            coef.hat[i]+c(-q.1, q.1)
        })                
        colnames(ci.boot.symm)<-names(coef.hat) 
        rownames(ci.boot.symm)<-c((alpha/2*100)%.%"%",(100-alpha/2*100)%.%"%") 
        res$ci.boot=ci.boot.symm
        
    }
    
    class(res)=c("double.hinge", class(res))
    return (res)
}


double.hinge.fit=function (x, y, mid.x, lb.quantile, ub.quantile, lower.y, upper.y, useC=TRUE) {
    
    # sort x and y
    order.=order(x)
    x=x[order.]
    y=y[order.]
    mid.quantile=.5
    if (!is.null(mid.x)) mid.quantile=mean(x<mid.quantile)
    chngpts.1 = get.chngpts (x, lb.quantile, mid.quantile, n.chngpts=Inf, stratified.by.sorted=NULL)
    chngpts.2 = get.chngpts (x, mid.quantile, ub.quantile, n.chngpts=Inf, stratified.by.sorted=NULL)
    
    if(useC) {
        coef = .Call("double_hinge_fit", 
            as.double(x), as.double(y), 
            as.double(chngpts.1), as.double(chngpts.2), 
            as.double(lower.y), as.double(upper.y),
            as.integer(0)
        )    
        names(coef)=c("lower.x","upper.x","slope") # this is for the .Call
        
    } else {
        mse = sapply(chngpts.1, function(e.1) {
              sapply(chngpts.2, function(e.2) {
                y.hat=ifelse(x<=e.1, lower.y, ifelse(x>=e.2, upper.y, lower.y+(upper.y-lower.y)/(e.2-e.1)*(x-e.1) ))
                sum((y-y.hat)**2)
        })
        })
        
        which(min(mse)==mse)
        ind=arrayInd(which.min(mse), .dim=dim(mse)); ind
        
        e.1=chngpts.1[ind[2]]; e.2=chngpts.2[ind[1]]; e.1; e.2
        coef=c(lower.x=e.1, upper.x=e.2, slope=(upper.y-lower.y)/(e.2-e.1))    
        
    }
    
    coef
    
}


plot.double.hinge=function(x, lcol=1, lwd=2, ...) {
    fit=x
    plot(fit$x, fit$y, xlab="", ylab="")
    lines(x=c(min(fit$x,na.rm=T), fit$coefficients["lower.x"], fit$coefficients["upper.x"], max(fit$x,na.rm=T)), y=c(fit$lower.y, fit$lower.y, fit$upper.y, fit$upper.y), col=lcol, lwd=lwd)
    
}
