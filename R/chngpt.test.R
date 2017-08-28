# chngpt.test is a streamlined version of chngpt.test.2 and only supports lr.mc method for making inference for both main effect-only model and interaction model
# mc.n=5e4; chngpts=NULL; lb.quantile=.1; ub.quantile=.9; chngpts.cnt=50; verbose=TRUE; prob.weights=NULL; family="logistic"
chngpt.test = function(formula.null, formula.chngpt, family=c("binomial","gaussian"), data, 
    type=c("step","hinge","segmented","stegmented"),
    main.method=c("lr","score"),
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, 
    chngpts.cnt=50, # this is set to 25 if int is weighted.two.sided or weighted.one.sided
    single.weight=1,
    mc.n=5e4, 
    prob.weights=NULL,
    compute.p.value=TRUE,
    verbose=FALSE 
) {    
    
    DNAME = deparse(substitute(data))
    if (missing(type)) stop("type mssing")
    family<-match.arg(family)        
    type<-match.arg(type)        
    main.method<-match.arg(main.method)        
    
    # keep only records that have no missing data for the null model and the change point model
    subset.1 = complete.cases(model.frame(formula.null, data, na.action=na.pass))
    subset.2 = complete.cases(model.frame(formula.chngpt, data, na.action=na.pass))
    data=data[subset.1 & subset.2,,drop=FALSE]
    
    tmp=as.matrix(model.frame(formula.chngpt, data))  
    col.Z=colnames(model.matrix(formula.null, data))  
    chngpt.var.name=setdiff(colnames(tmp), col.Z)[1]
    z.1.name=intersect(colnames(tmp), col.Z)
    chngpt.var = tmp[,chngpt.var.name]
    has.itxn = length(z.1.name)>0   
    
    
    #####################################################################################
    # make formula 
    
    f.null=if(type %in% c("segmented","stegmented")) update(formula.null, as.formula("~.+"%+%chngpt.var.name)) else formula.null    
    f.alt=update(f.null, as.formula(get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)))
    
    y=model.frame(f.null, data)[,1]
    Z=model.matrix(f.null, data)
    z.1 = Z[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    n=nrow(Z)
    p.null=ncol(Z)
    p.alt=switch(type, step=1, hinge=1, segmented=1, stegmented=2)
    p.alt.2=p.alt*ifelse(has.itxn,2,1)
    do.score= p.alt.2==1 & main.method=="score"
    
    if (is.null(prob.weights)) prob.weights=rep(1,n) 
    data$prob.weights=prob.weights # try put it in data to be found by glm
    
    if (is.null(chngpts)) chngpts=quantile(chngpt.var, seq(lb.quantile,ub.quantile,length=chngpts.cnt))
    M <- length(chngpts)  
    
    if(has.itxn & type!="step") stop("interaction model for this type not implemented yet: "%+%type)
    if(verbose) {
        myprint(n, p.null, p.alt, chngpt.var.name)        
        cat("Null model: "); print(f.null)
        cat("Change point model: "); print(formula.chngpt)
        if (has.itxn) cat("interaction var: ", z.1.name, "\n")
    }    
        
    if(p.alt.2==1 & main.method=="score") {
        method="Maximum of Score Statistics"
    } else {
        method="Maximum of Likelihood Ratio Statistics"
    }
    
    
    #####################################################################################
    # fit null model
    
    fit.null=keepWarnings(glm(formula=f.null, data=data, family=family, weights=prob.weights)) # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit.null$warning)!=0) {
        if(verbose) print(fit.null$warning)
        # ignore warning, often startsWith(fit.null$warning[[1]]$message,"non-integer #successes in a binomial glm!")
        fit.null=fit.null$value
    } else {
        fit.null=fit.null$value
    }
    
    
    #####################################################################################
    # compute LR statistics
    
    if(!do.score) {        
        glm.warn=FALSE
        QQ = numeric(M)
        for (m in 1:M) {
            data=make.chngpt.var(chngpt.var, chngpts[m], type, data)
            fit.alt=keepWarnings(glm(f.alt, data=data, family=family, weights=prob.weights)) 
            if(length(fit.alt$warning)!=0) {
                # how warning is handled greatly affects finite sample performance!
                # there are two main types of warnings: 
                #    glm.fit: fitted probabilities numerically 0 or 1 occurred
                #    glm.fit: algorithm did not converge
                if(verbose) { cat("At the ",m,"out of ",M," change point:"); print(fit.alt$warning) } 
                glm.warn=TRUE
                #print(summary(fit.alt$value))
                #return (list(chngpt=NA, p.value=NA, glm.warn=glm.warn))# need these fields for sim_test_batch
                QQ[m] = fit.null$aic - fit.alt$value$aic + 2*p.alt.2
            } else {
                QQ[m] = fit.null$aic - fit.alt$value$aic + 2*p.alt.2
                #QQ[m] = fit.null$deviance - fit.alt$value$deviance # this may not give the actual deviance, but up to a constant, and is family-dependent
            }
        }
        Q.max=max(QQ,na.rm=TRUE)
        max.id=which.max(QQ)
        chngpt=chngpts[max.id]
        # the following is used by chngptm, which call test to get starting value
        if(!compute.p.value) {
            res=list()
            res$chngpt=chngpt
            res$QQ=QQ
            res$chngpts=chngpts
            res$parameter=chngpt
            res$method=method
            res$fit.null.dev=fit.null$deviance
            class(res)=c("chngpt.test","htest",class(res))
            return (res) 
        }
        # repeat the fit at the chosen change point
        data=make.chngpt.var(chngpt.var, chngpt, type, data)
        fit.alt=glm(f.alt, data=data, family=family, weights=prob.weights)        
    }
    
    
    #####################################################################################
    # Compute W.null and V.S.hat, which are needed for computing p value and for computing test statistic in score test
    
    if (verbose) myprint("compute W.null")
    
    linear.predictors.null = fit.null$linear.predictors    
    beta.h = coef(fit.null)
    mu.h = fit.null$fitted.values
    
    # D.h is the inverse of the variance of working response variable
    D.h = if (family=="binomial") diag(c(mu.h*(1-mu.h))) else if (family=="gaussian") diag(length(mu.h)) / summary(fit.null)$dispersion 
    V.beta.h = solve(t(Z) %*% diag(prob.weights * diag(D.h)) %*% Z)
    V.eta.h = Z %*% V.beta.h %*% t(Z)
    A.h = diag(n) - D.h %*% V.eta.h %*% diag(prob.weights)
    
    # V.y is the variance of response variable
    V.y = if (family=="binomial") diag(c(mu.h*(1-mu.h))) else if (family=="gaussian") diag(length(mu.h)) * summary(fit.null)$dispersion 
    ADA = A.h %*% V.y %*% t(A.h) 
    
    W.M = matrix(0, nrow=n, ncol=p.alt*M)
    for (m in 1:M) W.M[,1:p.alt+p.alt*(m-1)] = make.chngpt.var(chngpt.var, chngpts[m], type)
    
    if(do.score) {
        W.null<-W.M
        
    } else {
        W.null = matrix(0, nrow=n, ncol=p.alt.2*M)
        for (m in 1:M) {
            data  = make.chngpt.var(chngpt.var, chngpts[m], type, data)
            X = model.matrix(f.alt, data) # make sure column order is the same as the coefficient order
            # fisher information for the full model, estimated under null. D.h is the diagonal matrix of variance estimated under NULL
            I.a = t(X) %*% D.h %*% X 
            I.bb.a = I.a[p.null+1:p.alt.2, p.null+1:p.alt.2] - I.a[p.null+1:p.alt.2, 1:p.null] %*% solve(I.a[1:p.null, 1:p.null], I.a[1:p.null, p.null+1:p.alt.2]) # the order depends on formula.1, hardcoded here
            #if(p.alt==1) W.null[,m] = W.M[,m] * I.bb.a**(-1/2) 
            if (p.alt.2>1) eig = eigen(solve(I.bb.a))
            I.bb.a.inv.sqrt=if (p.alt.2==1) I.bb.a**(-1/2) else eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
            W.null[,1:p.alt.2+p.alt.2*(m-1)] = X[,p.null+1:p.alt.2] %*% I.bb.a.inv.sqrt
        }       
    }    
    p = ncol(W.null)

    V.S.hat = t(W.null) %*% ADA %*% W.null
#    qr(V.S.hat, tol = 1e-8)$rank
#    isSymmetric(A.h)
    
    
    #####################################################################################
    # compute test statistics for score tests, this can only be done after V.S.hat is computed
    
    if(do.score) {        
        # scale to sd 1
        TT.0=c((t(W.M) %*% (y - mu.h)) / sqrt(diag(V.S.hat)))   
        TT = abs(TT.0)
        T.max=max(TT)
        max.id=which.max(TT)
        chngpt=chngpts[max.id]        
        # the following is used by chngptm, which call test to get starting value
        if(!compute.p.value) return (chngpt)                 
    } 
    
    
    #####################################################################################
    # compute p value
    #####################################################################################
    
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
    # simulate samples from mvrnorm
    set.seed(1)    
    sam=mvrnorm (mc.n, rep(0,p), cov2cor(V.S.hat)) # mvtnorm::rmvnorm can fail 
    
    if(p.alt.2==1 & main.method=="score") {
        
        sam=abs(sam)            
        tmp = apply(sam, 2, function(aux) list(aux))
        tmp = do.call(c, tmp)
        x.max=do.call(pmax, tmp)            
        p.value = mean(x.max>T.max)                    
        
        QQ=TT; Q.max=T.max # for consistent return
        glm.warn=NA
        fit.alt=NULL# if NA, then difficult to test is.na(fit.alt) because it is a test in LR test
        method="Maximum of Score Statistics"
                
    } else {
        
        # get p value
        sam=sam**2
        x.max = apply(sam, 1, function(aux) max( colSums(matrix(aux,nrow=p.alt.2)) )  )
        p.value = mean(x.max>Q.max)      
            
        method="Maximum of Likelihood Ratio Statistics"
        
        if(verbose==2) {
            myprint(p, dim(V.S.hat))
            #cat("V.S.hat\n"); print(V.S.hat)   
            #print(TT.0); print(round(V.S.hat[1:10,1:10],6))
        }
        
    }
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv) 
            
        
    #################################################################
    # return results
    res=list(chngpt=chngpt, statistic=Q.max, parameter=chngpt)
    names(res$statistic)="Maximal statistic"
    names(res$parameter)="change point"
    res$chngpts=chngpts    
    res$QQ=QQ
    res$method=method
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$data.name=DNAME
    res$V.S.hat = V.S.hat
    res$p.value=p.value    
    res$glm.warn=glm.warn
    res$fit.null=fit.null
    if (!is.null(fit.alt)) res$fit.alt=fit.alt
    class(res)=c("chngpt.test","htest",class(res))
    res
    
}




# note that when there is both main and interaction, there are two sets of statistics. the code cannot handle that for now, but remnant of it is in fold
plot.chngpt.test <- function(x, by.percentile=TRUE, both=FALSE, main=NULL, ...) {

    if(all(is.na(x$QQ))) {
        warning("all Q are NA, quit plotting")
        return()
    }
    fold=1
    perc=as.numeric(strtrim(names(x$chngpts),nchar(names(x$chngpts))-1))  
    idx=seq(1, length(x$chngpts), length.out=5)
    if(by.percentile) {
        # the primary x axis is percentile
        plot(perc, x$QQ, xlab=ifelse(both,"","change point (%)"), ylab="statistic", type="b", main=ifelse(is.null(main), x$method, main), xaxt="n", ...)
        axis(side=1, at=perc[idx], labels=round(perc[idx]))
        abline(v=(0:(fold-1))*100+as.numeric(strtrim(names(x$chngpt),nchar(names(x$chngpt))-1))  , lty=2)
        if (both) {
            axis(side=1, at=perc[idx], labels=round(x$chngpts[idx]), line=1.5, tick=FALSE, lwd=0)
        }
    } else {
        # the primary x axis is change point
        plot(x$chngpts, x$QQ, xlab=ifelse(both,"","change point"), ylab="statistic", type="b", main=ifelse(is.null(main), x$method, main), xaxt="n", ...)
        axis(side=1, at=x$chngpts[idx], labels=signif(x$chngpts[idx],2))
        abline(v=x$parameter, lty=2)
        if (both) {
            axis(side=1, at=x$chngpts[idx], labels=round(perc[idx]), line=1.5, tick=FALSE, lwd=0)
        }
    }
}
