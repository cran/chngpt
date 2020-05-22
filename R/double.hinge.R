# lb.quantile=.05; ub.quantile=.95; mid.x=NULL; lower.y=0; upper.y=1; mid.quantile.search.step=0.05
double.hinge=function(x, y, lower.y=NULL, upper.y=NULL, 
    var.type=c("none","bootstrap"), ci.bootstrap.size=1000, alpha=0.05, save.boot=TRUE, ncpus=1
) {
    
    not.missing=!is.na(x) & !is.na(y)
    x=x[not.missing]
    y=y[not.missing]
    
    var.type<-match.arg(var.type)    
    
    if(is.null(lower.y)) lower.y=min(y)
    if(is.null(upper.y)) upper.y=max(y)
    
    coef.hat=double.hinge.fit (x, y, lower.y, upper.y)
    
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
            double.hinge.fit (dat[ii], y[ii], lower.y, upper.y)               
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
        if(save.boot) res$boot.samples=boot.samples        
    }
    
    class(res)=c("double.hinge", class(res))
    return (res)
}


double.hinge.fit=function (x, y, lower.y, upper.y) {    
    # sort x and y
    order.=order(x)
    x=x[order.]
    y=y[order.]
    
    #chngpts = get.chngpts (x, lb.quantile, ub.quantile, n.chngpts=Inf, stratified.by.sorted=NULL)
    
    coef = .Call("double_hinge_fit_2", 
        as.double(x), as.double(y), 
        as.double(lower.y), as.double(upper.y)
    )    
    names(coef)=c("lower.x","upper.x","slope","mse") # this is for the .Call
    coef    
}


#which=NULL; xlim=NULL; lwd=2; lcol="red"; lty=1; add.points=TRUE; add.ci=TRUE; breaks=20; mark.chngpt=FALSE; xlab=NULL; ylab=NULL
plot.double.hinge=function(x, which=NULL, xlim=NULL, lwd=2, lcol="red", lty=1, add.points=TRUE, add.ci=TRUE, breaks=20, mark.chngpt=FALSE, xlab=NULL, ylab=NULL, ...) {
    
    has.boot.samples=FALSE
    if(!is.null(x$boot.samples)) has.boot.samples=TRUE 
    
    if(is.null(which)) {
        par(mfrow=c(1+has.boot.samples,1))
                              plot(x,which=1,xlim=xlim,lwd=lwd,lcol=lcol,mark.chngpt=mark.chngpt,...)
        if (has.boot.samples) plot(x,which=2,xlim=xlim,lwd=lwd,lcol=lcol,...)
        return(invisible())
    }    
    
    fit=x
    if(is.null(xlim)) xlim=range(fit$x)
    
    out=list()
    if(which==1) {
    # scatterplot with lines
        plot(fit$x, fit$y, xlab="", ylab="", xlim=xlim)
        lines(x=c(min(fit$x,na.rm=T), fit$coefficients["lower.x"], fit$coefficients["upper.x"], max(fit$x,na.rm=T)), y=c(fit$lower.y, fit$lower.y, fit$upper.y, fit$upper.y), col=lcol, lwd=lwd)    
        
    } else if (which==2) {
    # histograms of threshold estimates from bootstrap
        if (!has.boot.samples) stop("boot samples not saved")
        if(!hasArg("main")) {
            hist(fit$boot.samples[,1:2], xlim=xlim, breaks=breaks, xlab="bootstrap threshold estimate", main="", ...)
        } else {
            hist(fit$boot.samples[,1:2], xlim=xlim, breaks=breaks, xlab="bootstrap threshold estimate", ...)
        }
        ci=summary(fit)$coefficients[1:2,c("(lower","upper)")]
        if(add.ci) abline(v=ci, lty=c(2,2)) 
        
    } else stop("wrong which")
    
    invisible(out)
    
}

print.double.hinge=function(x, ...) {
    print(x$coefficients)
}


summary.double.hinge=function(object, verbose=FALSE, ...) {    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"        
    fit=object
    
    if (is.null(fit$ci.boot)) {
        cat("No variance estimate available.\n\n")
        print(fit)
        return (invisible())
    } else {    
        ci.boot= fit$ci.boot
    }
    
    lb=ci.boot[1,]
    ub=ci.boot[2,]
    sd=(ub-lb)/1.96/2
    pval=pnorm(abs(fit$coefficients/sd), lower.tail=F)*2
    
    out=cbind(fit$coefficients, sd, lb, ub, pval)
    out=out[-nrow(out),] # remove mse row
    out[c("lower.x","upper.x"),ncol(out)]=NA # p-val does not make sense for thresholds
    colnames(out) = c("Estimate", "Std. Error*", "(lower", "upper)", "Pr(>|t|)*")
    
    res=list()
    res$coefficients=out    
    class(res) <- "summary.double.hinge"
    res
}

print.summary.double.hinge=function(x,...) {
    cat("Double hinge model \n\n")
    cat("Coefficients:\n")
    print(x$coefficients)
}

vcov.double.hinge=function(object, var.type=NULL, ...) {
    fit=object
    if (is.null(fit$ci.boot)) {
        cat("No variance estimate available.\n\n")
        return (invisible())
    } else {    
        ci.boot= fit$ci.boot
    }

    lb=ci.boot[1,]
    ub=ci.boot[2,]
    sd=(ub-lb)/1.96/2
    
    if (!is.null(fit$boot.samples)){
        # estimate cor using bootstrap samples, but estimate var by scaling CI
        cor=cor(fit$boot.samples, use="pairwise.complete.obs")
        vcov=diag(sd)%*%cor%*%diag(sd)
    } else {
        vcov=matrix(NA,4,4)
        diag(vcov)=sd*sd
    }
    rownames(vcov)<-colnames(vcov)<-names(fit$coefficients)
    vcov=vcov[1:3,1:3]
    vcov
}

residuals.double.hinge=function(object, ...) {
    object$y-fitted(object)
}

fitted.double.hinge=function(object, ...) {
    with(object, 
        ifelse (x<coefficients["lower.x"], 
            lower.y, 
            ifelse (x>coefficients["upper.x"], 
                upper.y, 
                (x-coefficients["lower.x"])*coefficients["slope"]
                ))
    )    
        
}
