# tol=1e-4; maxit=1e2; verbose=TRUE; chngpt.init=NULL; search.all.thresholds=TRUE
chngptm = function(formula.null, formula.chngpt, data, search.all.thresholds=NULL, tol=1e-4, maxit=1e2, chngpt.init=NULL, 
    prob.weights=NULL,
    verbose=FALSE) {
    
    form.all = update(formula.null, formula.chngpt)
    subset. = complete.cases(model.frame(form.all, data, na.action=na.pass))
    data=data[subset.,,drop=FALSE]
    
    # create a new column chngpt.var.name%+%"_dich_at_chngpt"
    y=model.frame(formula.null, data)[,1]
    Z=model.matrix(formula.null, data)
    tmp=model.matrix(formula.chngpt, data)[,-1,drop=F]
    n=nrow(Z)
    p=ncol(Z)
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    
    b.=-30
    
    if (is.null(prob.weights)) prob.weights=rep(1,n) 
    data$prob.weights=prob.weights # try put it in data to be found by glm
    
    if (is.null(search.all.thresholds)) search.all.thresholds = n<500 
    
    if (search.all.thresholds) {
        
        sorted.chng.var=sort(chngpt.var)
        if (!has.itxn) {
            formula.new = update(formula.null, as.formula("~.+"%+%chngpt.var.name%+%"_dich_at_chngpt"))
        } else {
            formula.new = update(formula.null, as.formula("~.+"%+%z.1.name%+%"*"%+%chngpt.var.name%+%"_dich_at_chngpt"))
        }
        logliks=sapply (sorted.chng.var[-1], function(e) {
            data[[chngpt.var.name%+%"_dich_at_chngpt"]] = ifelse(chngpt.var>=e, 1, 0)
            fit = glm(formula.new, data, family="binomial", weights=prob.weights)        
            as.numeric(logLik(fit))
        } )
        e=sorted.chng.var[-1][which.max(logliks)]
        data[[chngpt.var.name%+%"_dich_at_chngpt"]] = ifelse(chngpt.var>=e, 1, 0)
        fit = glm(formula.new, data, family="binomial", weights=prob.weights) 
        if (verbose==2) print(summary(fit))
        coef.hat=c(coef(fit), e)       
        
    } else {
        
        # do a test to get init value for change point to be used in estimation
        if (is.null(chngpt.init)) {
            test = chngpt.test (formula.null, formula.chngpt, data)
            if (verbose) print(test)
            e.init=test$chngpt
        } else {
            test=NULL # test is returned either way
            e.init=chngpt.init
        }
        names(e.init)="e"
        
        #
        coef.hat=rep(0, ncol(Z)+ifelse(has.itxn,3,2))
        n.iter=0
        converged=TRUE        
        while(TRUE){    
        
            n.iter=n.iter+1
            if (verbose) cat("iter ", n.iter, "\n")
            if (n.iter>maxit) {converged=FALSE; break}
            
            # remake the binary change point variable in every iteration based on the change point estimate from last iteration
            data[[chngpt.var.name%+%"_dich_at_chngpt"]] = ifelse(chngpt.var>e.init, 1, 0)
            
            # the following is dividied into has interaction and has not
            # in both, we iterate between estimating effect given change point and estimating change point and associated effect given other effects
            if(!has.itxn) {
                # no interaction
                
                fit.0 = glm(update (formula.null, as.formula("~.+"%+%chngpt.var.name%+%"_dich_at_chngpt")), data, family="binomial", weights=prob.weights)        
                beta.init=coef(fit.0)[p+1]; names(beta.init)="beta"
                alpha.hat=coef(fit.0)[-(p+1)]
                alpha.z = c(Z %*% alpha.hat)
                
                if (verbose) {
                    myprint(beta.init, e.init, b.)
                }
                optim.out = optim(par=c(beta.init, e.init), 
                      fn = function(theta,...) sum(deviance.chng.deriv(theta[1],theta[2],...)), 
                      # if we use analytical gradient function by deriv3, we can get situations like exp(100), which will be Inf, and Inf/Inf will be NaN
                      #gr = function(theta,...) colSums(attr(deviance.chng.deriv(theta[1],theta[2],...), "gradient")), 
                      gr = NULL,
                      chngpt.var, y, b., alpha.z, 
                      lower = c(-10, quantile(chngpt.var, .1)), 
                      upper = c(10, quantile(chngpt.var, .9)), 
                      method="L-BFGS-B", control = list(), hessian = F)
                #sum(deviance.chng.deriv(beta.init, e.init, chngpt.var, y, b., alpha.z))
                
                e.init=optim.out$par["e"]; names(e.init)="e"
                coef.tmp=c(alpha.hat, optim.out$par)
                if (verbose) print(coef.tmp)
                if (max(abs(coef.tmp-coef.hat))<tol) {
                    coef.hat=coef.tmp
                    break
                } else {
                    coef.hat=coef.tmp
                }
        
            } else {
                # has interaction
                
                fit.0 = glm(update (formula.null, as.formula("~.+"%+%chngpt.var.name%+%"_dich_at_chngpt*"%+%z.1.name)), data, family="binomial", weights=prob.weights)        
                beta.init=coef(fit.0)[c(p+1,p+2)]; names(beta.init)=c("beta1","beta2")
                alpha.hat=coef(fit.0)[-c(p+1,p+2)]
                alpha.z = c(Z %*% alpha.hat)
                
                optim.out = optim(par=c(beta.init, e.init), 
                      fn = function(theta,...) sum(deviance.chng.itxn.deriv(theta[1],theta[2],theta[3],...)), 
                      gr = function(theta,...) colSums(attr(deviance.chng.itxn.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
                      chngpt.var, y, b., alpha.z, z.1,
                      lower = c(-10, -10, quantile(chngpt.var, .1)), 
                      upper = c(10, 10, quantile(chngpt.var, .9)), 
                      method="L-BFGS-B", control = list(), hessian = TRUE)
                
                e.init=optim.out$par["e"]; names(e.init)="e"
                coef.tmp=c(alpha.hat, optim.out$par)
                if (verbose) myprint(coef.hat)
                if (verbose) myprint(coef.tmp)
                if (max(abs(coef.tmp-coef.hat))<tol) {
                    coef.hat=coef.tmp
                    break
                } else {
                    coef.hat=coef.tmp
                }
                
            }
            
        } # end while 
    
    } # end if search.all.thresholds else
    
    names(coef.hat)[length((coef.hat))]="chngpt"
    names(coef.hat)[p+1]="("%+%chngpt.var.name%+%">chngpt)"
    if (has.itxn) names(coef.hat)[p+2]="("%+%chngpt.var.name%+%">chngpt):"%+%z.1.name
    
    # variance-covariance matrix
    # expressions for use with optim with deriv3
    alpha.z.s="("%+% concatList("alpha"%+%1:p%+%"*z"%+%1:p%+%"","+") %+%")"
    if (!has.itxn) {
        deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - beta/(1+exp(b*(x-e))) ) ) "
        params=c("alpha"%+%1:p, "beta", "e")
    } else {
        deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + (beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - (beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) "
        params=c("alpha"%+%1:p, "beta1", "beta2", "e")
    }
    params.long=c(params,"x","y","b","z"%+%1:p,"z.1")
    if (verbose) print(deviance.s)
    if (verbose) myprint(params)
    loss.f=deriv3(parse(text=deviance.s), params, params.long)    
    param.list = c(as.list(coef.hat), list(chngpt.var), list(y), list(b.), lapply(1:ncol(Z), function (i) Z[,i]), list(z.1))
    names(param.list)=params.long    
    tmp=do.call(loss.f, param.list)
    hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)                
    var.est = try(solve(hess)) # should keep change point in, and not do hess[-ncol(hess), -ncol(hess)], otherwise lead to over estimation of sd
    rownames(var.est) <- colnames(var.est) <- names(coef.hat)
    
    fit=list(
          coefficients=coef.hat
        , vcov=var.est
        , data=data
        , formula.null=formula.null
        , formula.chngpt=formula.chngpt
        , chngpt.var=chngpt.var.name
        , chngpt=coef.hat["chngpt"]
        , search.all.thresholds = search.all.thresholds     
    )
    if (!search.all.thresholds) fit=c(fit, list(converged=converged, iter=n.iter, test=test))
    names(fit$chngpt) = round(100*mean(chngpt.var<fit$chngpt),1) %+% "%"
    class(fit)=c("chngptm", class(fit))
    
    fit    
}

print.chngptm=function(x, ...) {
    if (!x$search.all.thresholds) {
        if (!x$converged) cat("Warning: not converged\n")
    }
    print(x$coefficients)
}

coef.chngptm=function(object, ...) {
    object$coefficients[-length(object$coefficients)] # not return the chngpoint estimate
}

vcov.chngptm=function(object, ...) {
    object$vcov[-length(object$coefficients),-length(object$coefficients)] # not return the chngpoint estimate
}

summary.chngptm=function(object, ...) {
    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"    
    
    fit=object
    p=length(fit$coefficients)
    n=nrow(fit$data)
    
    # assuming the last of coefficients is always the change point
    res=list()
    res$coefficients=mysapply(1:(p-1), function (i) {
        c(
              "OR"=exp(unname(fit$coefficients[i]))
            , "p.value" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
            , "(lower" = exp(unname(fit$coefficients[i] - sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
            , "upper)" = exp(unname(fit$coefficients[i] + sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
#              "Estimate"=unname(fit$coefficients[i])
#            , "Std. Error" = sqrt(fit$vcov[i,i])
#            , "t value" = unname(fit$coefficients[i] / sqrt(fit$vcov[i,i]))
#            , "Pr(>|t|)" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
        )
    })
    rownames(res$coefficients)=names(fit$coefficients)[-p]
    
    i=p
    res$chngpt=c(
              fit$chngpt
            , "p.value" = unname(pt(abs(fit$coefficients[i] / sqrt(fit$vcov[i,i])), df=n-p, lower.tail=FALSE))
            , "(lower" = (unname(fit$coefficients[i] - sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
            , "upper)" = (unname(fit$coefficients[i] + sqrt(fit$vcov[i,i]) * qt(0.975, df=n-p, lower.tail=TRUE)))
    )
    
    res
}


#deviance.3pl <- expression( (1-y) * ( c + (d-c)/(1+exp(b*(t-e))) )  +  log( 1 + exp( -c - (d-c)/(1+exp(b*(t-e))) ) ) )
#deviance.3pl.deriv=deriv3(deviance.3pl, c("c","d","e"), c("c","d","e","t","y","b"))

deviance.chng <- expression( (1-y) * ( alpha.z + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - beta/(1+exp(b*(x-e))) ) ) )
deviance.chng.deriv=deriv3(deviance.chng, c("beta","e"), c("beta","e","x","y","b","alpha.z"))

deviance.chng.itxn <- expression( (1-y) * ( alpha.z + (beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) )
deviance.chng.itxn.deriv=deriv3(deviance.chng.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))


#    c.init=ifelse(beta<0, (alpha+beta), (alpha)); names(c.init)="c"
#    d.init=ifelse(beta>0, (alpha+beta), (alpha)); names(d.init)="d"            
#    lb=c(logit(1/n),logit(1/n),quantile(dat$x,.01))
#    ub=c(logit((n-1)/n),logit((n-1)/n),quantile(dat$x,.99))
#        
#    ee=quantile(dat$x,seq(.15,.85,length=10)); ee=c(ee, e.)
#    fits=lapply(ee, function (e.init){
#        names(e.init)="e"
#        optim(par=c(c.init,d.init,e.init), 
#              fn = function(theta,...) sum(deviance.3pl.deriv(theta[1],theta[2],theta[3],...)), 
#              gr = function(theta,...) colSums(attr(deviance.3pl.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
#              dat$x, dat$y, b., method="L-BFGS-B", lower = lb, upper = ub, control = list(), hessian = F)
#    })            
#    fit=fits[[which.min(sapply(fits, function(x) x$value))]]
#    res=c(res, "fit"=chngpt.score.stat(e.=fit$par["e"], y~x, dat, b.=b.)["Z.stat"]) # similar performance to max.Z
#    
#    
#    plot(y~x, dat, log=""); abline(v=exp(e.))
#    
#    # likelihood evaluated at init
#    tmp=sum(deviance.3pl.deriv(c.init,d.init,e.init,dat$x,dat$y)); sum(tmp)
#    
#    # plotting likelihood
#    ee=exp(seq(3,7,length=100))
#    llik=sapply(ee, function (e.1){
#        lik=ThreePL.chng.t(t=dat$x,c=c.init,d=d.init,e=e.1)
#        lik[dat$y==0]=1-lik[dat$y==0]
#        llik=sum(log(lik))
#    })
#    plot(ee, llik, log="x", type="l")
#    
#    lik=ThreePL.chng.t(t=dat$x,c=coef(fit)[1],d=coef(fit)[2],e=coef(fit)[3])
#    lik[dat$y==0]=1-lik[dat$y==0]
#    sum(log(lik))
    
#    # compare variance estimate computed by deriv3 and by optim
#    tmp=deviance.3pl.deriv(coef(fit)[1],coef(fit)[2],coef(fit)[3],dat$x,dat$y)
#    sum(tmp, na.rm=T)
#    solve(-apply(attr(tmp,"h"), 2:3, mean))
#    
#    fit$fit$val
#    solve(fit$fit$hessian) # same as vcov(fit)
#    
#    
#    if (fit$convergence!=0) { 
#        did.not.fit=T
#    } else {
#        # compute variance estimate using deriv3
#        tmp=deviance.3pl.deriv(fit$par["c"], fit$par["d"], fit$par["e"], dat$x, dat$y)
#        hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)                
#        var.est = try(solve(hess)[1:2,1:2])
#        if (class(var.est)=="try-error") {
#            did.not.fit=T
#        } else {
#            did.not.fit=F
#        }
#    }
#    
#    if (did.not.fit) {
#        res=c(res, sigm=NA, sigm.cover=NA)
#    } else {
#        d.minus.c = fit$par["d"] - fit$par["c"]
#        d.minus.c.var = c(-1,1) %*% var.est %*% c(-1,1) 
#        res=c(res, sigm = pt(abs(d.minus.c)/sqrt(d.minus.c.var), df=nrow(dat)-3, lower.tail=F) * 2 )                                 
#        res=c(res, sigm.cover = d.minus.c - 1.96*sqrt(d.minus.c.var) < beta & beta < d.minus.c + 1.96*sqrt(d.minus.c.var) ) # coverage
#        res=c(res, fit$par["e"]) # e
#    }
