# chngpt.test is a streamlined version of chngpt.test.2 and only supports lr.mc method for making inference for both main effect-only model and interaction model
# prec.weights has not been verified in MC experiments b/c it is not often used
# test.statistic="lr"; p.val.method="MC"; mc.n=5e4; chngpts=NULL; lb.quantile=.1; ub.quantile=.9; chngpts.cnt=50; verbose=TRUE; prec.weights=NULL; robust=FALSE
chngpt.test = function(formula.null, formula.chngpt, family=c("binomial","gaussian"), data, type=c("step","hinge","segmented","stegmented"),
    test.statistic=c("lr","score"), # support for score is gradually descreasing
    chngpts=NULL, lb.quantile=.1, ub.quantile=.9, chngpts.cnt=50, # this is set to 25 if int is weighted.two.sided or weighted.one.sided
    prec.weights=NULL,
    p.val.method=c("MC","param.boot"), 
    mc.n=5e4, # 1e3 won't cut it, the p values estimated could be smaller than nominal
    boot.B=1e4,
    robust=FALSE,
    keep.fits=FALSE, verbose=FALSE
) {    
    
    DNAME = deparse(substitute(data))
    if (missing(type)) stop("type mssing")
    family<-match.arg(family)        
    type<-match.arg(type)        
    test.statistic<-match.arg(test.statistic)        
    p.val.method<-match.arg(p.val.method)        
    
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
    
    ## It is pre-mature to use chngptm to do this because the fastgrid mode does not yet support coarse grids
    # fastgrid and C version of grid is implemented only for the following scenarios
    fastgrid.ok = family=="gaussian" & type %in% c("hinge","segmented") & !has.itxn
    
    # make formula. Note that f.null includes chngptvar if type is segmented or stegmented
    f.null=if(type %in% c("segmented","stegmented")) update(formula.null, as.formula("~.+"%.%chngpt.var.name)) else formula.null    
    f.alt=update(f.null, get.f.alt(type, chngpt.var.name, modified.by=if(has.itxn) z.1.name else NULL))
    
    # extract data
    y=model.frame(f.null, data)[,1]
    Z=model.matrix(f.null, data)
    z.1 = Z[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    n=nrow(Z)
    p.null=ncol(Z)
    p.alt=switch(type, step=1, hinge=1, segmented=1, stegmented=2)
    p.alt.2=p.alt*ifelse(has.itxn,2,1)
    do.score= p.alt.2==1 & test.statistic=="score"
    
    # weights
    if (is.null(prec.weights)) prec.weights=rep(1,n) # it is how glm.fits handles weights by default
    data$prec.weights=prec.weights # try put it in data to be found by glm
    prec.w.all.one=all(prec.weights==1) 
    
    # sorted data
    chngpt.var.sorted=sort(chngpt.var)
    Z.sorted=Z[order(chngpt.var),,drop=FALSE]
    y.sorted=y[order(chngpt.var)]
    w.sorted=data$prec.weights[order(chngpt.var)]
    data.sorted=data[order(chngpt.var),]    
    
    # change point candidates
    if (is.null(chngpts)) chngpts=get.chngpts(chngpt.var.sorted,lb.quantile,ub.quantile,chngpts.cnt)
    M <- length(chngpts)  
    
    if(has.itxn & type!="step") stop("interaction model for this type not implemented yet: "%.%type)
    if(verbose) {
        cat("Null: "); print(f.null)
        cat("Alt:  "); print(f.alt)
        myprint(family, n, p.null, p.alt, chngpt.var.name, mc.n, fastgrid.ok)
        if (has.itxn) cat("interaction var: ", z.1.name, "\n")
    }    
        
    if(p.alt.2==1 & test.statistic=="score") {
        method="Maximum of Scores"
    } else {
        method="Maximum of Likelihood Ratios"
    }
    
    
    #####################################################################################
    # fit null model. if glm gives a warning, we may do something
    
    fit.null=keepWarnings(glm(formula=f.null, data=data, family=family, weights=prec.weights)) 
    if(length(fit.null$warning)!=0) {
        if(verbose) print(fit.null$warning)
        # ignore warning, often startsWith(fit.null$warning[[1]]$message,"non-integer #successes in a binomial glm!")
        fit.null=fit.null$value
    } else {
        fit.null=fit.null$value
    }
    
#    # Debugging. The lesson is that svyglm properly handles sampling probability weights in estimated variability, but not glm
#    library(survey)
#    summary(glm(formula=f.null, data=data, family=family))
#    summary(glm(formula=f.null, data=data, family=family, weights=prec.weights))
#    summary(svyglm(f.null, family=binomial(), design=svydesign(id=~1,strata=~y, weights=rep(1,nrow(data)), data=data)))
#    summary(svyglm(f.null, family=binomial(), design=svydesign(id=~1,strata=~y, weights=~prec.weights, data=data)))
    
    #####################################################################################
    # compute LR statistics (score statistic is computed in a later section)
    
    if(!do.score) {       
#        if(fastgrid.ok) {
#            if (verbose) print("chgnpt.test: do fastgrid search")
#            tmpfit=chngptm (formula.null, formula.chngpt, family, data=data, type=type, est.method="grid", var.type="none", grid.search.max=Inf, verbose=verbose>1, keep.best.fit=TRUE, weights=prec.weights) 
#            chngpt=tmpfit$chngpt
#            fit.alt=tmpfit$best.fit
#            Q.max=fit.null$aic - fit.alt$aic + 2*p.alt.2 # this is difference in deviance (2 x log lik) between models 
#            QQ=tmpfit$logliks # this is only proportional to logliks
#            glm.warn=tmpfit$glm.warn
#            #chngpts=tmpfit$chngpts
#        } else {
            if (verbose) print("chgnpt.test: do grid search")
            glm.warn=FALSE
            QQ = numeric(M)
            for (m in 1:M) {
                data=make.chngpt.var(chngpt.var, chngpts[m], type, data)
                fit.alt=keepWarnings(glm(f.alt, data=data, family=family, weights=prec.weights)) 
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
            #plot(chngpts, QQ)
            Q.max=max(QQ,na.rm=TRUE)
            max.id=which.max(QQ)
            chngpt=chngpts[max.id]
            # repeat the fit at the chosen change point
            data=make.chngpt.var(chngpt.var, chngpt, type, data)
            fit.alt=glm(f.alt, data=data, family=family, weights=prec.weights)        
#        }                        
    }
    
    
    #####################################################################################
    # Compute W.null and V.S.hat, which are needed for computing p value and for computing test statistic in score test
    
    if (verbose) myprint("compute W.null")
    
    linear.predictors.null = fit.null$linear.predictors    
    beta.h = coef(fit.null)
    mu.h = fit.null$fitted.values
    
    # D.h is the inverse of the variance of working response variable
    D.h = if (family=="binomial") {
        diag(c(mu.h*(1-mu.h))) 
    } else if (family=="gaussian") {
        diag(length(mu.h)) * summary(fit.null)$dispersion # dispersion is variance estimate under gaussian family
    }
    
    DW = if(prec.w.all.one) D.h else diag(diag(D.h) * prec.weights)
    
    # compute R := I - DZ(Z'DWZ)^{-1}Z'W is the residual operator
    V.eta.h = Z %*% solve(t(Z) %*% DW %*% Z) %*% t(Z)
    R.h = diag(n) - D.h %*% V.eta.h %*% diag(prec.weights) # if we inline V.eta.h, numerically we get different results
    
    if(robust) {
        V.y = diag(resid(fit.null, type="response")**2)
        RDR = R.h %*% V.y %*% t(R.h)     
    } else {
        if (prec.w.all.one) {
            RDR = R.h %*% D.h
            # R.h %*% D.h %*% t(R.h) should simplify to  R.h %*% D.h, but numerically we get very slightly different results:
            # print((c(R.h %*% V.y)-c(R.h %*% D.h %*% t(R.h))))
        } else {
            RDR = R.h %*% D.h %*% t(R.h)      
        }
    }
    
    if(do.score) {
        W.null = matrix(0, nrow=n, ncol=p.alt*M)
        for (m in 1:M) W.null[,1:p.alt+p.alt*(m-1)] = make.chngpt.var(chngpt.var, chngpts[m], type)
        
    } else {
        W.null = matrix(0, nrow=n, ncol=p.alt.2*M)
        for (m in 1:M) {
            data  = make.chngpt.var(chngpt.var, chngpts[m], type, data)
            X = model.matrix(f.alt, data) # make sure column order is the same as the coefficient order
            # fisher information for the full model, estimated under null. D.h is the diagonal matrix of variance estimated under NULL
            I.a = t(X) %*% DW %*% X 
            I.bb.a = I.a[p.null+1:p.alt.2, p.null+1:p.alt.2] - I.a[p.null+1:p.alt.2, 1:p.null] %*% solve(I.a[1:p.null, 1:p.null], I.a[1:p.null, p.null+1:p.alt.2]) # the order depends on formula.1, hardcoded here
            if (p.alt.2>1) eig = eigen(solve(I.bb.a))
            I.bb.a.inv.sqrt=if (p.alt.2==1) I.bb.a**(-1/2) else eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
            W.null[,1:p.alt.2+p.alt.2*(m-1)] = X[,p.null+1:p.alt.2] %*% I.bb.a.inv.sqrt
        }       
    }    
    p = ncol(W.null)
    B=W.null
    if(!prec.w.all.one) B=diag(prec.weights) %*% B # following the notation from "logistic regression.pdf"
    
    V.S.hat = t(B) %*% RDR %*% B
#    qr(V.S.hat, tol = 1e-8)$rank
#    isSymmetric(R.h)
    
    
    #####################################################################################
    # compute test statistics for score tests, this can only be done after V.S.hat is computed
    
    if(do.score) {        
        # scale to standard deviation 1
        TT.0=c((t(W.null) %*% (y - mu.h)) / sqrt(diag(V.S.hat)))   
        TT = abs(TT.0)
        T.max=max(TT)
        max.id=which.max(TT)
        chngpt=chngpts[max.id]        
    } 
    
    
    #####################################################################################
    # compute p value
    #####################################################################################
    
    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    set.seed(1)    
    
    if(p.val.method=="MC") {
        # simulate samples from MASS::mvrnorm
        sam=mvrnorm (mc.n, rep(0,p), cov2cor(V.S.hat)) # mvtnorm::rmvnorm can fail 
        
        if(p.alt.2==1 & test.statistic=="score") {            
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
            sam=sam**2
            x.max = apply(sam, 1, function(aux) max( colSums(matrix(aux,nrow=p.alt.2)) )  )
            #str(Q.max)
            #print(max(x.max))
            p.value = mean(x.max>Q.max)      
                
            method="Maximum of Likelihood Ratio Statistics"
            
            if(verbose==2) {
                myprint(p, dim(V.S.hat))
                #cat("V.S.hat\n"); print(V.S.hat)   
                #print(TT.0); print(round(V.S.hat[1:10,1:10],6))
            }
            
        }
    
    } else if (p.val.method=="param.boot") {    
        if (!fastgrid.ok | test.statistic=="score") stop("only /lr method and fastgrid.ok is supported for parametric bootstrap for now")
        sd.null=sqrt(summary(fit.null)$dispersion)
        linear.predictors.null.sorted=linear.predictors.null[order(chngpt.var)]
        Q.max.boot=sapply (1:boot.B, function(b) {
            # simulate from fit.null
            y.b=rnorm(n, linear.predictors.null.sorted, sd.null) 
            llik.null.b=-n/2*log(mean(lm.fit(Z.sorted, y.b)$residuals**2)) # #tmpfit=lm(y~Girth, data.frame(y=y.b, Z)); logLik(tmpfit) + n/2*(1+log(2*pi)) # same as llik.null.b
            f.name="fastgrid_" %.% family
            yhy = .Call(f.name, 
                            cbind(Z.sorted,chngpt.var.sorted), 
                            as.double(y.b), 
                            as.double(w.sorted), 
                            prec.w.all.one, 
                            attr(chngpts,"index"),
                            0,# bootstrap size
                            type=="upperhinge"
                    ) 
            2 * (max(-n/2*log((sum(y.b**2) - yhy)/n)) - llik.null.b)
        })
        p.value = mean(Q.max.boot>Q.max)      
        
    }
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv) 
            
        
    #################################################################
    # return results
    res=list(chngpt=chngpt, statistic=Q.max, parameter=chngpt)
    names(res$statistic)="Maximal statistic"
    names(res$parameter)="threshold"
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
    if(keep.fits) {
        res$fit.null=fit.null
        if (!is.null(fit.alt)) res$fit.alt=fit.alt
    }
    
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
    idx=seq(1, length(x$chngpts), length.out=5)
    
    # the primary x axis is change point
    plot(x$chngpts, x$QQ, xlab=ifelse(both,"","threshold"), ylab="statistic", type="b", main=ifelse(is.null(main), "Test by "%.%x$method, main), xaxt="n", ...)
    axis(side=1, at=x$chngpts[idx], labels=signif(x$chngpts[idx],2))
    abline(v=x$parameter, lty=2)
    
#    perc=as.numeric(strtrim(names(x$chngpts),nchar(names(x$chngpts))-1))  
#    if(by.percentile) {
#        # the primary x axis is percentile
#        plot(perc, x$QQ, xlab=ifelse(both,"","change point (%)"), ylab="statistic", type="b", main=ifelse(is.null(main), "Test by "%.%x$method, main), xaxt="n", ...)
#        axis(side=1, at=perc[idx], labels=round(perc[idx]))
#        abline(v=(0:(fold-1))*100+as.numeric(strtrim(names(x$chngpt),nchar(names(x$chngpt))-1))  , lty=2)
#        if (both) {
#            axis(side=1, at=perc[idx], labels=round(x$chngpts[idx]), line=1.5, tick=FALSE, lwd=0)
#        }
#    } 
}
