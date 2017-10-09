# useC=T; tol=1e-4; maxit=1e2; verbose=TRUE; est.method="grid"; var.type="bootstrap"; lb.quantile=.1; ub.quantile=.9; grid.search.max=500; weights=NULL; chngpt.init=NULL; alpha=0.05; b.transition=Inf; search.bound=10; ci.bootstrap.size=500; m.out.of.n=F; aux.fift=NULL
chngptm = function(formula.1, formula.2, family, data, #family can be coxph or any glm family, but variance estimate is only available for binomial and gaussian (only model-based for latter)
  type=c("step","hinge","segmented","segmented2","stegmented"), # segmented2 is the model studied in Cheng 2008
  est.method=c("default","smoothapprox","grid"), 
  var.type=c("none","robust","model","smooth","robusttruth","bootstrap","all"), aux.fit=NULL, 
  test.inv.ci=TRUE, boot.test.inv.ci=FALSE, # test.inv.ci is passed to local functions, boot.test.inv.ci is global within this function
  lb.quantile=.1, ub.quantile=.9, grid.search.max=5000, ci.bootstrap.size=500, alpha=0.05, save.boot=FALSE, m.out.of.n=FALSE, # grid.search.max is the maximum number of grid points used in grid search
  b.transition=Inf,# controls whether threshold model or smooth transition model
  tol=1e-4, maxit=1e2, chngpt.init=NULL, search.bound=10,
  weights=NULL, verbose=FALSE, useC=TRUE, fast=TRUE, ...) 
{
    
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)    
    var.type<-match.arg(var.type)    
    est.method<-match.arg(est.method)    
    
    stopifnot(b.transition>0)
    b.search=if(b.transition==Inf) 30 else b.transition
    if(type=="segmented2" & b.transition==Inf) stop("for segmented2, b.transition should not be Inf") # limited implementation for segmented2
    
    # not all variance estimation methods are implemented for all families    
    if(!family %in% c("gaussian","binomial")) {
        if(!var.type %in% c("bootstrap")) stop("no analytical variance estimates are provided for this family: "%+%family)
    }     
    
    # grid search is implemented for all families, but smoothapprox is only implemented for gaussian and binomial
    if (est.method=="default") {
        est.method = if(family%in%c("binomial","gaussian") & is.null(weights)) "smoothapprox" else "grid"
    } else if (est.method=="smoothapprox" & !family%in%c("binomial","gaussian")) stop ("smoothapprox not implemented for this family")
    
    if (var.type %in% c("robust","robusttruth","all") & is.null(aux.fit)) stop("need an aux.fit for robust variance estimate")
    
    # remove missing observations
    form.all = update(formula.1, formula.2)
    subset. = complete.cases(model.frame(form.all, data, na.action=na.pass))
    data=data[subset.,,drop=FALSE]
    
    # decide whether there is interaction
    y=model.frame(formula.1, data)[,1]
    Z=model.matrix(formula.1, data)
    tmp=model.matrix(formula.2, data)[,-1,drop=F]
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    if(is.na(chngpt.var.name)) stop("Something is wrong. Check the formulal. ")
    if(verbose>=2) myprint(chngpt.var.name)
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    chngpt.var.sorted=sort(chngpt.var)
    #print(chngpt.var.sorted); chngpt.var.sorted=chngpt.var[order(chngpt.var)]; print(chngpt.var.sorted)
    Z.sorted=Z[order(chngpt.var),,drop=FALSE]
    y.sorted=y[order(chngpt.var)]
    data.sorted=data[order(chngpt.var),]
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    
    n=nrow(Z)
    p.z=ncol(Z)
    p.2=switch(type, step=1, hinge=1, segmented=2, segmented2=2, stegmented=3)
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    p=p.z+p.2.itxn+1 #total number of paramters, including threshold
    
    if (is.null(weights)) weights=rep(1,n) 
    data$weights=weights # try put it in data to be found by glm
    
    # make formula that includes all parameters but threshold
    formula.new = if (type %in% c("segmented","segmented2","stegmented")) update(formula.1, as.formula("~.+"%+%chngpt.var.name)) else formula.1
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    formula.new=update(formula.new, as.formula(f.alt))
    
    if (verbose) {
        myprint(type, est.method, has.itxn)     
        print(formula.new)
        myprint(p.z, p.2.itxn, p)
    }
        
    # change point candidates
    nLower=round(nrow(data)*lb.quantile)+1; nUpper=round(nrow(data)*ub.quantile); #myprint(nLower, nUpper)    
    chngpts=chngpt.var.sorted[nLower:nUpper]
#    # take the mid points between data points as change point candidates. For step models, use any point in between does not change likelihood, but it does change likelihood for segmented model
#    chngpts=(chngpts[1:(length(chngpts)-1)]+chngpts[-1])/2
    if (length(chngpts)>grid.search.max) {
        chngpts=chngpts[round(seq(1,length(chngpts),length=grid.search.max))]
        if(est.method=="grid") stop("some work needs to be done to resolve this issue")
    }
    # another way to define nLower/nUpper is as follows. It is equivalent
    #nLower=sum(chngpt.var.sorted<quantile(chngpt.var.sorted, lb.quantile))+1; nUpper=sum(chngpt.var.sorted<=quantile(chngpt.var.sorted, ub.quantile));myprint(nLower, nUpper)
    
    
    # 
    grid.search=function(){
        if(verbose) cat("perform grid search\n")
        if(useC & family=="gaussian") {
            # putting as.double around a matrix changes the matrix into a vector, which affects 
            # as.double is needed around y.sorted b/c if y.sort is integer, this throws an error b/c the cpp function expects real
            # logliks is actually -rss here
            logliks=-.Call("grid_search", cbind(Z.sorted,chngpt.var.sorted, if(type=="segmented") chngpt.var.sorted), as.double(y.sorted), nLower, nUpper) # this is a problem, Z may not have chngptvar as the last 
        } else {
            logliks=sapply (chngpts, function(e) {
                data=make.chngpt.var(chngpt.var, e, type, data, b.transition)
                fit = do.regression (formula.new, data, weights, family)            
                if(length(fit$warning)!=0) {
                    if(verbose>=3) print(fit$warning)
                    #print(summary(fit$value))
                    as.numeric(logLik(fit$value))
                } else as.numeric(logLik(fit$value))            
                # debug use
                #sum(resid(fit$value)**2)
            } )
        }
        #print(logliks)
    
        e=chngpts[which.max(logliks)]
        #myprint(which.max(logliks))
        #myprint(e)
        attr(e, "glm.warn")=any(is.na(logliks))
        attr(e, "logliks")=logliks
        if(verbose>=2) {
            plot(chngpts, logliks, type="b", xlab="change point")
            abline(v=e)
        }
        e        
    }
    
    # find e.final
    if (est.method=="grid") {
    #### grid search
        e.final=grid.search()           
        glm.warn=attr(e.final,"glm.warn")  
        
    } else if (est.method=="smoothapprox") {
    #### Newton-Raphson
        
        if(verbose) cat("smoothapprox search\n")
        # do a test to get init value for change point to be used in estimation 
        e.init=chngpt.init
        if (is.null(e.init)) {
            if(verbose) cat("looking for init\n")
            chngpt.test.0 = chngpt.test (formula.1, formula.2, family=family, data, type=ifelse(type=="segmented2","segmented",type), compute.p.value=FALSE)# note main.method is lr by default, which is faster
            if (verbose==2) plot(chngpt.test.0, by.percentile=FALSE)
            e.init=chngpt.test.0$chngpt
        } 
        
        if(length(e.init)==0) {
            e.final=grid.search()
            glm.warn=attr(e.final,"glm.warn")
        } else {
            names(e.init)="e"
            if (verbose) cat("init e: ", e.init, "\n")
            glm.warn=FALSE
            
            coef.hat=rep(0, 1+ncol(Z)+p.2.itxn)
            n.iter=0
            converged=TRUE   
            e.effective.maxit=numeric(maxit)     
            while(TRUE){    
            
                n.iter=n.iter+1
                if (n.iter>maxit) {converged=FALSE; break}
                if (verbose) cat("iter ", n.iter, "\t")
                
                # remake the binary change point variable in every iteration based on the change point estimate from last iteration
                data = make.chngpt.var(chngpt.var, e.init, type, data, b.search) # b.search is used here to be consistent with optim which uses b.search as well
                fit.0 = do.regression(formula.new, data, weights, family)$value
                
                # update threshold and associated coefficients
                beta.init=coef(fit.0)[p.z+1:p.2.itxn]; #names(beta.init)=c("beta1","beta2")# we may need this for variance calculation to work
                alpha.hat=coef(fit.0)[1:p.z]
                stopifnot(all(!is.na(alpha.hat)))
                alpha.z = c(Z %*% alpha.hat)
                
                # search for better e and slopes associated with x and thresholded x
                optim.out = try(optim(par=c(beta.init, e.init), 
                      fn = get("dev."%+%type%+%"."%+%ifelse(has.itxn,"itxn.","")%+%"f"), 
                      gr = NULL,
                      #gr = get("dev."%+%type%+%"."%+%ifelse(has.itxn,"itxn.","")%+%"deriv.f"), 
    #                  # if we use analytical gradient function by deriv3 in optim, we can get situations like exp(100), which will be Inf, and Inf/Inf will be NaN
    #                  fn = function(theta,...) sum(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...)), 
    #                  gr = function(theta,...) colSums(attr(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
                      chngpt.var, y, b.search, alpha.z, z.1, family,
                      lower = c(rep(-search.bound,length(beta.init)), quantile(chngpt.var, lb.quantile)), 
                      upper = c(rep( search.bound,length(beta.init)), quantile(chngpt.var, ub.quantile)), 
                      method="L-BFGS-B", control = list(), hessian = TRUE))
                      
                if(class(optim.out)=="try-error") {
                    if (verbose) cat("error doing smoothapprox search, switch to grid search\n")
                    e.final=grid.search()
                    glm.warn=attr(e.final,"glm.warn")
                    break;
                }
                            
                e.init=optim.out$par["e"]
                # to avoid algorithm ocillation due to artifact
                e.effective.maxit[n.iter]=which(chngpt.var.sorted>=e.init)[1]
    #            # set the "effective" change point to the data point that is just to the right of the optim output. This leads to more clear indications when perfect separate happens
                e.init=mean(chngpt.var.sorted[e.effective.maxit[n.iter]])
                # set the "effective" change point to the mid point between the two data points that sandwich the optim output. 
    #            e.init=mean(chngpt.var.sorted[e.effective.maxit[n.iter]+(-1):0])
                names(e.init)="e" # this is needed for optim to work
                
                coef.tmp=c(alpha.hat, optim.out$par)
                coef.tmp[length(coef.tmp)]=e.init # set to effective change point
                if (verbose) cat(coef.tmp, "\n")
                if (max(abs(coef.tmp - coef.hat)) < tol) {
                    coef.hat = coef.tmp
                    e.final=last(coef.hat)
                    break
                }
                else {
                    coef.hat = coef.tmp
                    e.final=last(coef.hat)
                }
                
            } # end while 
        } # end if else
    
    } # end if grid/smoothapprox
    # fit glm using e.final
    data = make.chngpt.var(chngpt.var, e.final, type, data, b.transition) # note that b.transition is used here instead of b.search
    fit = do.regression (formula.new, data, weights, family)$value
    coef.hat=c(coef(fit), "chngpt"=e.final)
    best.fit=fit
    
    names(coef.hat)[length((coef.hat))]="chngpt"
    if (type=="stegmented") {
        replacement="I("%+%chngpt.var.name%+%">chngpt)"
        replacement.2="("%+%chngpt.var.name%+%"-chngpt)+"
        new.names=sub("x.gt.e.2", replacement.2, names(coef.hat))    
        new.names=sub("x.gt.e", replacement, new.names)    
    } else {
        if(type=="step") {
            replacement="I("%+%chngpt.var.name%+%">chngpt)"
        } else if (type %in% c("hinge","segmented","segmented2")) {
            replacement="("%+%chngpt.var.name%+%"-chngpt)+"
        } 
        new.names=sub("x.gt.e", replacement, names(coef.hat))    
    }
    if (verbose) cat(new.names,"\n")
    names(coef.hat)=new.names
    

    # keep an empty line here between estimation and var estimate
    
                
    ###############################################################################################
    # variance-covariance 
    
    # only implemented when there is no interaction and for linear and logistic regression
    if (var.type!="bootstrap" & (has.itxn | !family %in% c("binomial","gaussian"))) {
        var.est=NULL
    } else {
        if (type %in% c("step","stegmented")) {
            # discontinous models
            var.est=vcov(best.fit)
            
        } else if (type %in% c("hinge","segmented","segmented2")) { # continuous change point models
            
            var.est.smooth=function(){
                if(verbose) cat("in var.est.smooth\n")
    
                # expressions for use with optim with deriv3, they are different from the set of expressions, dev.step etc, b/c each alpha has to be separate
                alpha.z.s="("%+% concatList("alpha"%+%1:p.z%+%"*z"%+%1:p.z%+%"","+") %+%")"
                if (family=="binomial") {
                    if (type=="hinge") {
                        deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + (x-e)*beta/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - (x-e)*beta/(1+exp(-b*(x-e))) ) ) "
                        params=c("alpha"%+%1:p.z, "beta", "e")
                    } else if (type=="segmented") {
                        deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + beta1*x + (x-e)*beta2/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - beta1*x - (x-e)*beta2/(1+exp(-b*(x-e))) ) ) "
                        params=c("alpha"%+%1:p.z, "beta1", "beta2", "e")
                    } else if (type=="segmented2") {
                        deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + beta1*x + x*beta2/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - beta1*x - x*beta2/(1+exp(-b*(x-e))) ) ) "
                        params=c("alpha"%+%1:p.z, "beta1", "beta2", "e")
                    }
                } else if (family=="gaussian") {
                    if (type=="hinge") {
                        deviance.s <- " (y- ( "%+% alpha.z.s %+%" + (x-e)*beta/(1+exp(-b*(x-e))) ))**2"
                        params=c("alpha"%+%1:p.z, "beta", "e")
                    } else if (type=="segmented") {
                        deviance.s <- " (y- ( "%+% alpha.z.s %+%" + beta1*x + (x-e)*beta2/(1+exp(-b*(x-e))) ))**2"
                        params=c("alpha"%+%1:p.z, "beta1", "beta2", "e")
                    } else if (type=="segmented2") {
                        deviance.s <- " (y- ( "%+% alpha.z.s %+%" + beta1*x + x*beta2/(1+exp(-b*(x-e))) ))**2"
                        params=c("alpha"%+%1:p.z, "beta1", "beta2", "e")
                    }
                }
                params.long=c(params,"x","y","b","z"%+%1:p.z,"z.1")
                if (verbose) print(deviance.s)
                if (verbose) myprint(params)
                loss.f=deriv3(parse(text=deviance.s), params, params.long)    
                # b.search instead of b.transition needs to be used below, otherwise hess is singular
                param.list = c(as.list(coef.hat), list(chngpt.var), list(y), list(b.search), lapply(1:ncol(Z), function (i) Z[,i]), list(z.1))
                names(param.list)=params.long    
                tmp=do.call(loss.f, param.list)
                hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)       
                print(hess)
                print(eigen(hess))         
                var.est = try(solve(hess)) # should keep change point in, and not do hess[-ncol(hess), -ncol(hess)], otherwise lead to over estimation of sd
                
                if (family=="gaussian") {
                    var.est=var.est*2*sigma(best.fit)^2
                }
                
                rownames(var.est) <- colnames(var.est) <- names(coef.hat)
                var.est
            }
            
            # model-based
            var.est.model=function(test.inv.ci, robust=FALSE){
                if(verbose) cat("in var.est.model\n")
                                
                # set up some variables
                e=coef.hat["chngpt"]
                x.gt.e=chngpt.var>e
                beta.4=coef.hat["("%+%chngpt.var.name%+%"-chngpt)+"]
                
                # x.tilda
                if(type %in% c("hinge","segmented")) {
                    if (b.transition==Inf) {
                        x.tilda=cbind(Z, if(type=="segmented") chngpt.var, (chngpt.var-e)*x.gt.e, -beta.4*x.gt.e)
                    } else {
                        x.tilda=cbind(Z, if(type=="segmented") chngpt.var, (chngpt.var-e)*expit(b.transition*(chngpt.var-e)), -beta.4*expit(b.transition*(chngpt.var-e))*(1+b.transition*(chngpt.var-e)*(1-expit(b.transition*(chngpt.var-e)))))
                    }               
                } else if (type=="segmented2") {
                    # earlier we checked that b.transition should not be Inf
                    x.tilda=cbind(Z, chngpt.var, chngpt.var*expit(b.transition*(chngpt.var-e)),                               -beta.4*expit(b.transition*(chngpt.var-e))        *b.transition*chngpt.var*(1-expit(b.transition*(chngpt.var-e))) )
                }
                
                if (family=="binomial") {
                    p.est=drop(expit(x.tilda[,-p] %*% coef.hat[-p]))                                                  
                    if(!robust) {
                        V=-tXDX(x.tilda, p.est*(1-p.est))/n    
                        var.est=solve(-V)/n                        
                    } else {
                        # sandwich estimation
                        B.inv=solve(tXDX(x.tilda, p.est*(1-p.est))/n)
                        M=tXDX(x.tilda, resid(best.fit,type="response")**2)/n    
                        var.est=B.inv %*% M %*% B.inv/n                        
                    }
                } else if (family=="gaussian") {
                    if(!robust) {
                        V=-t(x.tilda) %*% (x.tilda)/n/sigma(best.fit)^2
                        var.est=solve(-V)/n                        
                    } else {
                        # sandwich estimation
                        B.inv=solve(t(x.tilda) %*% (x.tilda)/n)
                        M=tXDX(x.tilda, resid(best.fit,type="response")**2)/n    
                        var.est=B.inv %*% M %*% B.inv/n                        
                    }
                }
                #myprint(e); print(eigen(V)$values); print(x.tilda[order(chngpt.var),])
                
                rownames(var.est) <- colnames(var.est) <- names(coef.hat)                
                
                # profile likelihood ratio test inversion CI
                if(test.inv.ci) {
                    chngpt.ci=ci.test.inv(qchisq(1-alpha, df=1) * 1)
                    attr(var.est,"chngpt.ci")=chngpt.ci
                }
                
                var.est
            }
    
            var.est.robust=function(aux.fit, test.inv.ci, true.param=FALSE){
                if(verbose) cat("in var.est.robust\n")
                                
                # compute density of X at change point. have to use coef.hat["chngpt"] for change point here b/c e is defined later
                den=density(chngpt.var)
                pt1=which(den$x>=coef.hat["chngpt"])[1]
                f.x.e = interpolate(c(den$x[pt1],den$y[pt1]), c(den$x[pt1-1],den$y[pt1-1]), coef.hat["chngpt"])   # interploate to find the density at change point
                
                # debugging only: use population parameter values to calculate variance
                if (true.param) {
                    # quadratic model                    
                    coef.hat=c("(Intercept)"=-2.8297600157,   "z"=0.3378716621,   "x"=0.5526491719,   "(x-chngpt)+"=1.3409043525,   "chngpt"=3.7541306365) 
                    f.x.e=dnorm(3.75, 4.7, 1.4) #0.22
                    ## sigmoid2, in sim studies, we don't need robust variance estimate for sigmoid2
                    #coef.hat=c("(Intercept)"=-0.09428639728, "z"=log(1.4), "x"=-log(0.67), "(x-chngpt)+"=log(0.4), "chngpt"=4.5) 
                    #f.x.e=dnorm(4.5, 4.7, 1.6) 
                }
                
                # set up some variables
                e=coef.hat["chngpt"]
                x.gt.e=chngpt.var>e
                beta.4=coef.hat["("%+%chngpt.var.name%+%"-chngpt)+"]
                
                # x.tilda
                x.tilda=cbind(Z, if (type=="segmented") chngpt.var, (chngpt.var-e)*x.gt.e, -beta.4*x.gt.e)
                
                antilink=get(family)()$linkinv
                mu.est=drop(antilink(x.tilda[,-p] %*% coef.hat[-p]))        
                M=tXDX(x.tilda, resid(best.fit,"response")^2)/n # has to use the robust, sandwich estimation version. Note the type needs to be response for glm
                # non-robust version is not justified because it depends on the assumption that mean model is correct
                #M=t(x.tilda)%*%x.tilda *sigma(best.fit)^2 /n 
                if (family=="binomial") {
                    V.1=-tXDX(x.tilda, mu.est*(1-mu.est))/n    
                } else if (family=="gaussian") {
                    V.1=-t(x.tilda) %*% (x.tilda)/n                    
                }
                # V.2 is the same between binomial and gaussian families
                V.2=matrix(0,nrow=p,ncol=p)
                # there are two ways to compute the off-diagonal element. when use true coefficients, it works better to use predicted response by true model; when use estimated, using y works better
                #V.2[p,p-1]<-V.2[p-1,p]<- if(is.null(attr(aux.fit,"truemodel"))) -mean((y-mu.est)*x.gt.e) else -mean((predict(aux.fit, data, "response")-mu.est)*x.gt.e)
                V.2[p,p-1]<-V.2[p-1,p]<- -mean((y-mu.est)*x.gt.e) 
                # compute V.2[p, p]
                newdata=data; newdata[[chngpt.var.name]]=e
                m.0=mean(predict(aux.fit, newdata, "response"))
                p.e.z.0=mean(antilink(cbind(Z, if (type=="segmented") e) %*% coef.hat[1:(p-2)]))                
                V.2[p, p]= beta.4 * f.x.e * (m.0-p.e.z.0)
                V=V.1+V.2    
                V.inv=solve(V) 
                
                var.est=V.inv%*%M%*%V.inv/n       
                rownames(var.est) <- colnames(var.est) <- names(coef.hat)                
                
                if(verbose==2){
                    print(V.1)
                    print(V.2)
                    myprint(diag(V))
                    myprint(diag(V.inv))
                    cat("eigen values: \n")
                    myprint(round(eigen(V.1)$values, 5))
                    myprint(round(eigen(V.2)$values, 5))
                    myprint(round(eigen(V)$values, 5))
#                    myprint(beta.4, f.x.e, m.0, p.e.z.0, m.0-p.e.z.0)
#                    myprint(mean(y*x.gt.e), mean(p.est*x.gt.e))
#                    myprint(mean(predict(aux.fit, data, "response") * x.gt.e))
                }    
                
                # profile likelihood ratio test inversion CI
                if(test.inv.ci) {
                    scale.chisq= last(diag(var.est))/(-last(diag(V.inv))/n) 
                    if (verbose>=2) {
                        myprint(last(diag(var.est)), -last(diag(V.inv))/n, scale.chisq)
                    }
                    if (scale.chisq<0) {
                        warning("scale.chisq is negative")
                        scale.chisq=abs(scale.chisq)
                    }                     
                    if (family=="binomial") {
                        # no need to do anything
                    } else if (family=="gaussian") {
                        scale.chisq=scale.chisq/sigma(best.fit)^2
                    }
                    chngpt.ci=ci.test.inv(qchisq(1-alpha, df=1) * scale.chisq)
                    attr(var.est,"chngpt.ci")=chngpt.ci
                }
    
                var.est
            }
            
            # c.alpha is the critical value for diff in deviance
            ci.test.inv=function(c.alpha){
                if (verbose) cat("in ci.test.inv\n")
                
                if (est.method=="grid") {
                    c(NA,NA)# todo
                } else if (est.method=="smoothapprox") {
                    
                    idx.chngpt=which(chngpt.var.sorted>=coef.hat["chngpt"])[1]
                    data = make.chngpt.var(chngpt.var, chngpt.var.sorted[idx.chngpt], type, data, b.transition)
                    fit=do.regression (formula.new, data, weights, family)
                    lik.max = as.numeric(logLik(fit$value))
                    
                    profile.liks=rep(NA, length(chngpt.var.sorted))
                    
                    # go left
                    lik=lik.max
                    idx=idx.chngpt
                    #myprint(idx, lik, c.alpha)
                    while (2*(lik.max-lik)<c.alpha) {
                        profile.liks[idx]=lik
                        idx=idx-1
                        if(idx==0) {
                            #if (verbose) warning("no more data on the left")
                            lik=NA; break
                        }
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.transition)
                        fit = do.regression (formula.new, data, weights, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            lik = NA; break
                        } else lik = as.numeric(logLik(fit$value))                             
                    } 
                    #lb=if(is.na(lik)) NA else chngpt.var.sorted[idx+1] 
                    lb=chngpt.var.sorted[idx+1] 
                    
                    # go right
                    lik=lik.max
                    idx=idx.chngpt
                    while (2*(lik.max-lik)<c.alpha) {
                        profile.liks[idx]=lik
                        idx=idx+1
                        if(idx>n) {
                            #if (verbose) warning("no more data on the right")
                            lik=NA; break
                        }
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.transition)
                        fit = do.regression (formula.new, data, weights, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            lik = NA; break
                        } else lik = as.numeric(logLik(fit$value))     
                    } 
                    #ub=if(is.na(lik)) NA else chngpt.var.sorted[idx-1]
                    ub=chngpt.var.sorted[idx-1]
                    
                    if (verbose==2) {
                        plot(chngpt.var.sorted, chngpt.test.0$fit.null.dev + 2*profile.liks, xlab="change points", ylab="deviance relative to null model", main="Profile Likelihood")
                    }
                    
                    c(lb,ub)
                }
                
            }
    
            # use F statistics to find test-inversion CI
            ci.test.inv.F=function(c.alpha){
                if (verbose) cat("in ci.test.inv.F\n")
                
                if (est.method=="grid") {
                    c(NA,NA)# todo
                } else if (est.method=="smoothapprox") {
                    
                    idx.chngpt=which(chngpt.var.sorted>=coef.hat["chngpt"])[1]
                    data = make.chngpt.var(chngpt.var, chngpt.var.sorted[idx.chngpt], type, data, b.transition)# it is ok that data is assigned to b/c the function adds columns
                    fit=do.regression (formula.new, data, weights, family)
                    sigsq.max = (sigma(fit$value))**2; #myprint(sigsq.max)
                    #sigsq.max = (sigma(best.fit))**2; myprint(sigsq.max)
                    
                    profile.sigsqs=rep(NA, length(chngpt.var.sorted))
                    
                    # go left
                    sigsq=sigsq.max
                    idx=idx.chngpt
                    #myprint(idx, sigsq, c.alpha)
                    while (n*(sigsq/sigsq.max-1)<c.alpha) {
                        profile.sigsqs[idx]=sigsq
                        idx=idx-1
                        if(idx==0) {
                            if (verbose) warning("no more data on the left")
                            sigsq=NA; break
                        }
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.transition)
                        fit = do.regression (formula.new, data, weights, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            sigsq = NA; break
                        } else sigsq = (sigma(fit$value))**2
                        if(verbose>=2) print(sigsq)
                    } 
                    #lb=if(is.na(sigsq)) NA else chngpt.var.sorted[idx+1] 
                    lb=chngpt.var.sorted[idx+1] 
                    
                    # go right
                    sigsq=sigsq.max
                    idx=idx.chngpt
                    while (n*(sigsq/sigsq.max-1)<c.alpha) {
                        profile.sigsqs[idx]=sigsq
                        idx=idx+1
                        if(idx>n) {
                            if (verbose) warning("no more data on the right")
                            sigsq=NA; break
                        }
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.transition)
                        fit = do.regression (formula.new, data, weights, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            sigsq = NA; break
                        } else sigsq = (sigma(fit$value))**2
                    } 
                    #ub=if(is.na(sigsq)) NA else chngpt.var.sorted[idx-1]
                    ub=chngpt.var.sorted[idx-1]
                    
                    if (verbose==2) {
                        plot(chngpt.var.sorted, chngpt.test.0$fit.null.dev + 2*profile.sigsqs, xlab="change points", ylab="F-statistic", main="Test Inversion")
                    }
                    
                    c(lb,ub)
                }
                
            }
                
            # the index may be hardcoded in this function
            ci.bootstrap=function(){
                if(verbose) cat("in ci.bootstrap\n")
                
                # bootstrapping
                # save rng state before set.seed in order to restore before exiting this function
                save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
                set.seed(1)     # this is only added on 7/30/2017, well after biometrics paper is published. should not change the results qualitatively
                
                if (est.method=="grid" & useC & family=="gaussian") {                    
                    boot.out=list()
                    class(boot.out)="boot"
                    boot.out$t0=coef.hat
                    boot.out$R=ci.bootstrap.size
                    boot.out$call=c("boot")
                    boot.out$sim="ordinary"
                    boot.out$stype="i"
                    boot.out$fake=TRUE
                    
                    # weirdly, f.name cannot be put inline b/c rcmdcheck throws an error otherwise
                    f.name=ifelse(fast,"boot_grid_search_fast", "boot_grid_search")
                    boot.out$t=.Call(f.name, cbind(Z.sorted,chngpt.var.sorted, if(type=="segmented") chngpt.var.sorted), as.double(y.sorted), nLower, nUpper, ci.bootstrap.size) # this is a problem, Z may not have chngptvar as the last 
                    boot.out$t=t(matrix(boot.out$t, ncol=ci.bootstrap.size))
                    
                } else {
                    boot.out=boot(data.sorted, R=ci.bootstrap.size, sim = "ordinary", stype = "i", statistic=function(dat, ii){
                        #print(ii)
                        if (m.out.of.n) ii=ii[1:(4*sqrt(n))] # m out of n bootstrap
                        # use chngpt.init so that the mle will be less likely to be inferior to the model conditional on e.hat, but it does not always work
                        # for studentized bootstrap interval, need variance estimates as well, however, stud performs pretty badly
                        fit.ii=try(chngptm (formula.1, formula.2, family, dat[ii,], type, est.method=est.method, var.type="none", b.transition=b.transition, chngpt.init=coef.hat["chngpt"], useC=useC, verbose=verbose>=3))
                        tmp=length(coef.hat)*1+ifelse(!boot.test.inv.ci, 0, ifelse(family=="gaussian",2,1)) # need to adjust depending on the last else clause
                        out = if(inherits(fit.ii,"try-error")) {
                            rep(NA,tmp) 
                        } else if (any(is.na(fit.ii$coefficients))) {
                            rep(NA,tmp) 
#                        # the next if is controversial, but at least has no impact on qudratic/250/gaussian
#                        } else if (any(abs(fit.ii$coefficients[1:(length(fit.ii$coefficients)-1)])>search.bound)) {
#                            rep(NA,tmp)
                        } else {
                            if(!boot.test.inv.ci) {
                                fit.ii$coefficients
                            } else {
                                # compute profile likelihood ratio
                                # pl can be negative even when grid search is used in both point estimate and bootstrap
                                # this is because the estimated chngpt may not be in the {x} of the bootstrapped dataset, and a value in between can be better than {x}
                                dat.tmp=make.chngpt.var(chngpt.var[ii], coef.hat["chngpt"], type, data[ii,], b.transition) # don't forget to do chngpt.var[ii]!
                                fit.ii.ehat=do.regression (formula.new, data=dat.tmp, weights, family); if(length(fit.ii.ehat$warning)!=0 & verbose) print(fit.ii.ehat$warning)
                                pl= 2*(as.numeric(logLik(fit.ii$best.fit)) - as.numeric(logLik(fit.ii.ehat$value)))
                                Fs=n*(sigma(fit.ii.ehat$value)^2/sigma(fit.ii$best.fit)^2-1) # F statistic, an approximation of pl, but may work better for linear regression
                                if(pl<0 & verbose) {
                                    print(summary(fit.ii$best.fit))
                                    print("------------------------------------------------------------")
                                    print(summary(fit.ii.ehat$value))
                                    print("------------------------------------------------------------")
                                    print(summary(fit.ii))
                                    print(sort(dat[ii,"x"]))
                                    #stop("debugging")
                                }
                                c(fit.ii$coefficients, pl, if(family=="gaussian") Fs) # for Gaussian, return both pl and Fs
                            }  
                        }
                        out                        
                    }) # end boot call
                    #str(boot.out$call)
                }
                # restore rng state 
                assign(".Random.seed", save.seed, .GlobalEnv)     
                
                boot.samples=boot.out$t
                # use boot.ci to extract basic, percentile and abc confidence intervals
                outs=list(perc=NULL, basic=NULL, bc=NULL) # this fixes the order of the list components
                for (i in 1:length(coef.hat)){
                    # for stud, index needs to be of length 2 and the second element needs to be var. If the second element is not present, will have warnings: index out of bounds; minimum index only used
                    #ci.out=try(boot.ci(boot.out, type=c("perc","basic",if(ci.bootstrap.size>n) "bca","stud"), index=c(i, i+length(coef.hat))))
                    ci.out=suppressWarnings({
                        #try(boot.ci(boot.out, type=c("perc","basic",if(ci.bootstrap.size>n & is.null(boot.out$fake)) "bca"), index=i), silent=TRUE)# bca requires a true boot.out object
                        try(boot.ci(boot.out, type=c("perc","basic"), index=i), silent=TRUE)# bca requires a true boot.out object, changed to this for performance timing to be fair to true boot call
                    })
                    # suppressWarnings changes the return type
                    outs$perc= cbind(outs$perc,  if(!inherits(ci.out,"bootci")) c(NA,NA) else ci.out$percent[1,4:5])
                    outs$basic=cbind(outs$basic, if(!inherits(ci.out,"bootci")) c(NA,NA) else ci.out$basic[1,4:5])                    
                    #outs$abc=  cbind(outs$abc,   if(!inherits(ci.out,"bootci")) c(NA,NA) else if(ci.bootstrap.size>n) ci.out$bca[1,4:5] else c(NA,NA) ) # for bca to work, the number of bootstrap replicates needs to be greater than sample size
                    #outs$stud= cbind(outs$stud,  if(!inherits(ci.out,"bootci")) c(NA,NA) else ci.out$student[1,4:5])
                }
                colnames(outs$perc)<-colnames(outs$basic)<-names(coef.hat) # <-colnames(outs$abc)
                # find bc CI from boot.out$t as ci.out does not provide bc
                outs$bc=sapply (1:length(coef.hat), function(i) {
                    z.0=qnorm(mean(boot.samples[,i]<coef.hat[i], na.rm=TRUE))                    
                    quantile(boot.samples[,i], c(pnorm(2*z.0+qnorm(alpha/2)), pnorm(2*z.0+qnorm(1-alpha/2))), na.rm=TRUE)
                })
                colnames(outs$bc)<-names(coef.hat) 
                # symmetric percentile CI as defined in Hansen (2017)
                outs$symm=sapply (1:length(coef.hat), function(i) {
                    q.1=quantile(abs(boot.samples[,i]-coef.hat[i]), 1-alpha, na.rm=TRUE)
                    coef.hat[i]+c(-q.1, q.1)
                })                
                colnames(outs$symm)<-names(coef.hat) 
                
                # compute test inversion CI for threshold based on bootstrap critical value
                if(boot.test.inv.ci) {
                    outs$testinv=matrix(NA,nrow=2,ncol=length(coef.hat))
                    colnames(outs$testinv)<-names(coef.hat)
                    if(verbose>=2) print(summary(boot.samples[,1+length(coef.hat)]))
                    if (family=="gaussian") {
                        # print both c.alpha for comparison, but use Fs-based, because it seems to work better
                        c.alpha=quantile(boot.samples[,1+length(coef.hat)], 1-alpha, na.rm=TRUE); if(verbose) myprint(c.alpha)
                        #c.alpha=quantile(boot.samples[,2+length(coef.hat)], 1-alpha, na.rm=TRUE); if(verbose) myprint(c.alpha)
                        if (!is.na(c.alpha)) outs$testinv[,length(coef.hat)]= ci.test.inv.F(c.alpha) 
                    } else {
                        c.alpha=quantile(boot.samples[,1+length(coef.hat)], 1-alpha, na.rm=TRUE); if(verbose) myprint(c.alpha)
                        if (!is.na(c.alpha)) outs$testinv[,length(coef.hat)]= ci.test.inv(c.alpha)                     
                    }
                }
         
                # the following percentile CI is numerically different from boot.ci results because the latter use norm.inter to do interpolation on the normal quantile scale
                #ci.perc=apply(boot.out$t, 2, function (x) quantile(x, c(alpha/2, 1-alpha/2), na.rm=TRUE))   
                                
#                # this takes a long time and returns non-sensible result
#                ci.abc=abc.ci(data, statistic=function(dat, ww){
#                        fit.ii=try(chngptm (formula.1, formula.2, family, dat, type, est.method="smoothapprox", var.type="none"))
#                        # note that we cannot use | in the following line
#                        if(inherits(fit.ii,"try-error")) 
#                            rep(NA,length(coef.hat)) 
#                        else if (any(abs(fit.ii$coefficients[1:(length(fit.ii$coefficients)-1)])>search.bound)) 
#                            rep(NA,length(coef.hat)) 
#                        else fit.ii$coefficients
#                    }, 
#                    index=4, conf=1-alpha
#                )
#                myprint(ci.abc)
                
                if(save.boot) outs[["boot.samples"]]=boot.samples
                
                outs                
            }
            
            # if perfect segregation happens, cann't compute variance estimate
            # what to do with glm.warn?
            if (var.type=="none") {
                var.est=NULL
            } else if (var.type=="bootstrap") {
                var.est=ci.bootstrap()
            } else {                
                if (any(abs(coef.hat[-c(1,length(coef.hat))])>=search.bound)) {
                    cat("point estimate over search bound. No analytical variance estimate produced\n")
                    var.est=NULL
                } else {
                    var.est=switch(var.type, 
                        smooth=var.est.smooth(), 
                        model=var.est.model(test.inv.ci), 
                        robust=     var.est.robust(aux.fit, test.inv.ci), 
                        robusttruth=var.est.robust(aux.fit, test.inv.ci, true.param=TRUE), 
                        all=list("smooth"=var.est.smooth(), 
                                  "model"=var.est.model(test.inv.ci), 
                                 "robust"=var.est.robust(aux.fit, test.inv.ci), 
                               "sandwich"=var.est.model(test.inv.ci=FALSE,robust=TRUE) )
                    ) 
                }
            }
        } else stop("type incorrect")
        
    } # end if (has.itxn | !family %in% c("binomial","gaussian"))  else 
    
        
    res=list(
          best.fit=best.fit # this is placed first b/c it is too long
        , coefficients=coef.hat
        , vcov=var.est
        , formula.1=formula.1
        , formula.2=formula.2
        , chngpt.var=chngpt.var.name
        , chngpt=coef.hat["chngpt"]
        , est.method = est.method
        , b.transition=b.transition
        , type = type
        , glm.warn = glm.warn
        , family=family    
        , var.type=var.type
    )
    names(res$chngpt) = round(100*mean(chngpt.var<res$chngpt),1) %+% "%"
    
    if (est.method=="smoothapprox") {
        res=c(res, list(converged=converged, iter=n.iter))
    } else if (est.method=="grid") {
        res=c(res, list(chngpts=chngpts, logliks=attr(e.final,"logliks")))
#        # good.soln needs to take as input a chngptm object that has chngpts etc
#        class(res)=c("chngptm", class(res)) 
#        tmp=good.soln(res, plot=verbose==3)
#        res=c(res, list(good.soln=tmp) )
    }
    
    class(res)=c("chngptm", class(res))
    res    
}


## make.chngpt.var is needed in both testing and estimation
## in estimation by grid search, b.transition is set to null when this function is called
## in estimation by smooth approx, b.transition value is saved to object
#
# The difference between the next two versions is that when b.transition is not infinity, if x<e, x is 0 in the older version but only close to 0 in the newer version
## before May 13, 2017
#make.chngpt.var=function(x, e, type, data=NULL, b.transition=NULL) {    
#    if(type=="step") {
#        out=ifelse(x>=e, 1, 0)
#    } else if (type=="hinge") {
#        out=ifelse(x>=e, x-e, 0)
#    } else if (type=="segmented") {
#        out=ifelse(x>=e, x-e, 0) # x also becomes part of the null model
#    } else if (type=="stegmented") {
#        out=cbind(ifelse(x>=e, 1, 0), ifelse(x>=e, x-e, 0))  # x also becomes part of the null model
#    }    
#    if (!is.null(b.transition)) out=out * 1/(1+exp(b.transition*(x-e))) 
#
## New on May 13, 2017 this function is modified to fit smooth transition model
make.chngpt.var=function(x, e, type, data=NULL, b.transition=Inf) {    
    transition=expit(b.transition*(x-e))
    if(b.transition==Inf) transition[is.nan(transition)]=1 # if this is not here, there will be NaN and that messes up things
    if(type=="step") {
        out=transition
    } else if (type=="hinge") {
        out=transition*(x-e)
    } else if (type=="segmented") {
        out=transition*(x-e) # x also becomes part of the null model
    } else if (type=="segmented2") {
        out=transition*x # x by itself also becomes part of the null model
    } else if (type=="stegmented") {
        out=cbind(transition, transition*(x-e))  # x also becomes part of the null model
    }    
    
    #str(out)    
    # unchanged on May 13, 2017
    if (is.null(data)) {
        out    
    } else {
        if (type=="stegmented") {
            data$x.gt.e = out[,1]
            data$x.gt.e.2 = out[,2]
        } else {
            data$x.gt.e = out
        }
        data
    }
}
    
do.regression=function(formula.new, data, weights, family){
    if (family=="coxph") {
        fit = keepWarnings(survival::coxph(formula.new, data=data, weights=weights))
    } else {
#    str(formula.new); str(data); str(weights); str(family)
        fit =             keepWarnings(glm(formula.new, data=data, weights=weights, family=family) )
    }
    fit
}
get.f.alt=function(type, has.itxn, z.1.name, chngpt.var.name) {
    if (type=="stegmented") {
        f.alt=if(has.itxn) "~.+(x.gt.e+x.gt.e.2)*"%+%z.1.name%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name else "~.+x.gt.e+x.gt.e.2"
    } else if (type %in% c("segmented","segmented2")) {
        f.alt=if(has.itxn) "~."%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name%+%"+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    } else if (type %in% c("step","hinge")) {
        f.alt=if(has.itxn) "~.+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    }
    f.alt
}
        
#deviance.3pl <- expression( (1-y) * ( c + (d-c)/(1+exp(b*(t-e))) )  +  log( 1 + exp( -c - (d-c)/(1+exp(b*(t-e))) ) ) )
#deviance.3pl.deriv=deriv3(deviance.3pl, c("c","d","e"), c("c","d","e","t","y","b"))


predict.chngptm=function (object, newdata = NULL, type = c("link", "response", "terms"), ...){    
    newdata = make.chngpt.var(newdata[[object$chngpt.var]], object$chngpt, object$type, newdata, object$b.transition)            
    predict(object$best.fit, newdata, type, ...)
}
print.chngptm=function(x, ...) {
    if (x$est.method=="smoothapprox") {
        if (!x$converged) cat("Warning: not converged\n")
    }
    print(x$coefficients)
}
coef.chngptm=function(object, ...) {
    object$coefficients[-length(object$coefficients)] # not return the chngpoint estimate
}
vcov.chngptm=function(object, ...) {
    if(object$type %in% c("hinge","segmented")) {
        object$vcov[-length(object$coefficients),-length(object$coefficients)] # not return the chngpoint estimate
    } else  {
        object$vcov
    }
}
getFixedEf.chngptm=function(object, ...) {
    capture.output({
        res=summary(object)
    })
    rbind(res$coefficients, chngpt=c(res$chngpt[1], NA, res$chngpt[2:3]))
}

summary.chngptm=function(object, var.type=NULL, verbose=FALSE, ...) {    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"        
    fit=object
    p.z=length(fit$coefficients)
    n=nrow(fit$best.fit$data)
    type=fit$type
    
    if (is.null(fit$vcov)) {
        cat("No variance estimate available.\n\n")
        print(fit)
        return (invisible())
    } else {
    
        boot.conf=FALSE        
        if(!is.null(var.type)) {
            if (var.type %in% c("perc","basic","bc","bcabc","symm")) {
                boot.conf=TRUE
            }
            vcov=fit$vcov[[var.type]]

        } else {            
            if(is.list(fit$vcov)){
                if (!is.null(fit$vcov$robust)) {
                    vcov=fit$vcov$robust
                } else if (!is.null(fit$vcov$symm)){
                    vcov=fit$vcov$symm
                    boot.conf=TRUE
                    #if (all(is.na(vcov))) vcov=fit$vcov$bc # if the bootstrap samples have too many NA, only bc quantiles can be estimated. However, bc quantiles may be misleading
                }
                
            } else {
                vcov=fit$vcov
            }
        } 
        if(is.null(vcov)) {
            cat("No variance estimate available.\n\n")
            print(fit)
            return (invisible())
        }         
    
    }
    str(vcov)
    
    if(object$type %in% c("hinge","segmented") & !boot.conf) {
        vcov.t=vcov[p.z,p.z]
        tmp=vcov[-p.z,-p.z] # not return the chngpoint estimate
        # the last line lost the attr chngpt.ci, need to make a copy
        if(!is.null(attr(vcov,"chngpt.ci"))) attr(tmp,"chngpt.ci")=attr(vcov,"chngpt.ci")
        vcov=tmp
    } else  {
        vcov
    }
    #print(vcov)
    
    # assuming the last of coefficients is always the change point
    # deal with coefficients and change point separately
    
    res=list()
    
    # coefficients
    transf=if(fit$family=="binomial") exp else if (fit$family=="gaussian") identity
    if (boot.conf){
        lb=transf(vcov[1,])
        ub=transf(vcov[2,])
        pval=rep(NA,p.z)        
    } else  {
        lb=transf(unname(fit$coefficients[1:(p.z-1)] - sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        ub=transf(unname(fit$coefficients[1:(p.z-1)] + sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        pval=unname(pt(abs(fit$coefficients[1:(p.z-1)] / sqrt(diag(vcov))), df=n-p.z, lower.tail=FALSE)) *2 # *2 is added on 8/2/2016
    }
    res$coefficients=mysapply(1:(p.z-1), function (i) {
        c(
              transf(unname(fit$coefficients[i]))
            , "p.value" = pval[i]
            , "(lower" = lb[i]
            , "upper)" = ub[i]
#              "Estimate"=unname(fit$coefficients[i])
#            , "Std. Error" = sqrt(vcov[i,i])
#            , "t value" = unname(fit$coefficients[i] / sqrt(vcov[i,i]))
#            , "Pr(>|t|)" = unname(pt(abs(fit$coefficients[i] / sqrt(vcov[i,i])), df=n-p.z, lower.tail=FALSE))
        )
    })
    rownames(res$coefficients)=names(fit$coefficients)[-p.z]
    colnames(res$coefficients)[1]=if(fit$family=="binomial") "OR" else if (fit$family=="gaussian") "Est"
    
    # change point
    i=p.z
    if (type %in% c("hinge","segmented")) {
        if (boot.conf){
            lb=vcov[1,p.z]
            ub=vcov[2,p.z]
        } else if(!is.null(attr(vcov,"chngpt.ci"))) {
            if(verbose) print("get test inversion CI for threshold")
            lb=unname(attr(vcov,"chngpt.ci")[1])
            ub=unname(attr(vcov,"chngpt.ci")[2])
        } else {
            lb=(unname(fit$coefficients[i] - sqrt(vcov.t) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
            ub=(unname(fit$coefficients[i] + sqrt(vcov.t) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        }
        res$chngpt=c(
              fit$chngpt
#            , "p.value" = unname(pt(abs(fit$coefficients[i] / sqrt(vcov[i,i])), df=n-p.z, lower.tail=FALSE))# makes no sense to have a pvalue for chngpt
            , "(lower" = lb
            , "upper)" = ub
        )
    } else if (type %in% c("step","stegmented")) {
        res$chngpt=c(fit$chngpt)
    } else stop("incorrect type")
    
    res$type=fit$type
    
    class(res) <- "summary.chngptm"
    res
}
print.summary.chngptm=function(x,...) {
    cat("Change point model type: ",x$type,"\n\n")
    cat("Coefficients:\n")
    print(x$coefficients)
    cat("\nThreshold:\n")
    print(x$chngpt)    
}


# the following functions are used as objective functions in smoothapprox, b/c we don't supply derivative to optim anymore. 
.lik.f=function(linear, y, family) {tmp=if(family=="binomial") (1-y)*linear + log(1+exp(-linear)) else if(family=="gaussian") (linear-y)**2; sum(tmp)}
# step change point model
dev.step.f <- function(theta,x,y,b,alpha.z,z.1,family)  {
    beta=theta[1]; e=theta[2]
    eta=beta/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
}
dev.step.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(beta1+beta2*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
}
# hinge change point model
dev.hinge.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta=theta[1]; e=theta[2]
    eta=(x-e)*beta/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
dev.hinge.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(x-e)*(beta1+beta2*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
}
# segmented change point model
dev.segmented.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + (x-e)*beta2/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
dev.segmented2.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + x*beta2/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
dev.segmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*(x-e)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
dev.segmented2.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*x/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
# stegmented change point model
dev.stegmented.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[2]; beta3=theta[3]; e=theta[4]
    eta=beta1*x + ((x-e)*beta3 + beta2)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
# the order of parameters in the following needs to be fixed
dev.stegmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; beta5=theta[5]; beta6=theta[6]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + ((beta3+beta4*z.1)*(x-e) + beta4+beta5*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family)
} 
# dev.step.deriv.f is experimental, not really used, does not seem to bring improvement when compared to gr=NULL
dev.step.deriv.f <- function(theta,x,y,b,alpha.z,z.1,family)  {
    if(family!="binomial") stop("only binomial supported")
    beta=theta[1]; e=theta[2]
    approxi=expit(b*(x-e))
    dBde=b*approxi*(1-approxi)
    pi=expit(alpha.z+beta*approxi)
    colSums((y-pi)*cbind(x>e, beta*dBde))
}
# expressions, binomial only
# step change point model
dev.step <- expression( (1-y) * ( alpha.z + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - beta/(1+exp(-b*(x-e))) ) ) )
dev.step.deriv=deriv3(dev.step, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.step.itxn <- expression( (1-y) * ( alpha.z + (beta1+beta2*z.1)/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -alpha.z - (beta1+beta2*z.1)/(1+exp(-b*(x-e))) ) ) )
dev.step.itxn.deriv=deriv3(dev.step.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))
# hinge change point model
dev.hinge <- expression( (1-y) * ( alpha.z + (x-e)*beta/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*beta/(1+exp(-b*(x-e))) ) ) )
dev.hinge.deriv=deriv3(dev.hinge, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.hinge.itxn <- expression( (1-y) * ( alpha.z + (x-e)*(beta1+beta2*z.1)/(1+exp(-b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*(beta1+beta2*z.1)/(1+exp(-b*(x-e))) ) ) )
dev.hinge.itxn.deriv=deriv3(dev.hinge.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))

# check whether the MLE is considered good, return Boolean
good.soln=function(fit, df=7, plot=FALSE) {
    if (class(fit)[1]!="chngptm") {
        warning("fit has to be of class chngptm")
        return (NA)
    }
    
    # fit a smoothing spline
    complete=!is.na(fit$logliks)
    fm1 <- smooth.spline(x=fit$chngpts[complete], y=fit$logliks[complete], df=df)
    # compute locations where the first derivatives are 0
    deriv.fm1=predict(fm1,deriv=1)
    tmp=deriv.fm1$y>0
    k=length(tmp)
    tmp.2=xor(tmp[1:(k-1)], tmp[2:k])
    optimum=which(tmp.2)
    optimum.1=ifelse(abs(deriv.fm1$y[optimum])<abs(deriv.fm1$y[optimum+1]), optimum, optimum+1)
    
    if(plot) {
        par(mfrow=c(2,1))
        plot(fit$chngpts, fit$logliks)
        lines(fm1, col=2, type="l")        
        abline(v=fit$chngpts[optimum.1])
        plot(deriv.fm1$x, deriv.fm1$y); abline(h=0)
        abline(v=fit$chngpts[optimum.1])
    }    
    
    if(length(optimum.1)!=1) return (FALSE) else {
        if (abs(match(fit$chngpt,fit$chngpts)-optimum.1)<=length(fit$chngpts)/10) return (TRUE) else return (FALSE) # 20 is hardcoded
    }        
}


performance.unit.test=function(formula.1, formula.2, family, data, B, I){
    y=model.frame(formula.1, data)[,1]
    Z=model.matrix(formula.1, data)
    tmp=model.matrix(formula.2, data)[,-1,drop=F]
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    if(is.na(chngpt.var.name)) stop("Something is wrong. Check the formulal. ")
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    chngpt.var.sorted=sort(chngpt.var)
    #print(chngpt.var.sorted); chngpt.var.sorted=chngpt.var[order(chngpt.var)]; print(chngpt.var.sorted)
    Z.sorted=Z[order(chngpt.var),,drop=FALSE]
    y.sorted=y[order(chngpt.var)]
    data.sorted=data[order(chngpt.var),]
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    
    type="segmented"
    n=nrow(Z)
    p.z=ncol(Z)
    p.2=switch(type, step=1, hinge=1, segmented=2, segmented2=2, stegmented=3)
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    p=p.z+p.2.itxn+1 #total number of paramters, including threshold
    
    # make formula that includes all parameters but threshold
    formula.new = if (type %in% c("segmented","segmented2","stegmented")) update(formula.1, as.formula("~.+"%+%chngpt.var.name)) else formula.1
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    formula.new=update(formula.new, as.formula(f.alt))
    
    if (FALSE) {
        myprint(type, est.method, has.itxn)     
        print(formula.new)
        myprint(p.z, p.2.itxn, p)
    }
        
    -.Call("performance_unit_test", cbind(Z.sorted,chngpt.var.sorted, if(type=="segmented") chngpt.var.sorted), as.double(y.sorted), B, I)
}



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
#              dat$x, dat$y, b.transition, method="L-BFGS-B", lower = lb, upper = ub, control = list(), hessian = F)
#    })            
#    fit=fits[[which.min(sapply(fits, function(x) x$value))]]
#    res=c(res, "fit"=chngpt.score.stat(e.=fit$par["e"], y~x, dat, b.transition=b.transition)["Z.stat"]) # similar performance to max.Z
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
#
#
#                # an alternative parameterization
#                # x.tilda
#                beta.3=coef.hat["x"]
#                beta.4=coef.hat["(x-chngpt)+"]+beta.3
#                if (type=="segmented") {
#                    x.tilda=cbind(Z, (chngpt.var-e)*(1-x.gt.e), (chngpt.var-e)*x.gt.e, -beta.3*(1-x.gt.e)-beta.4*x.gt.e)
#                } else stop("type not supported yet")                                
#                # V.1
#                V.1=-tXDX(x.tilda, p.est*(1-p.est))/n                    
#                # M
#                M= if(robust) tXDX((y-p.est)*x.tilda, rep(1,n))/n else -V.1
#                # V.2
#                V.2=matrix(0,nrow=p,ncol=p)
#                if (robust) V.2[p,p-1]<-V.2[p-1,p]<- -mean((y-p.est)*x.gt.e)
#                if (robust) V.2[p,p-2]<-V.2[p-2,p]<- -mean((y-p.est)*(1-x.gt.e))
#                # last diagonal element of V.2
#                den=density(chngpt.var)
#                pt1=which(den$x>=e)[1]
#                f.x.e = interpolate(c(den$x[pt1],den$y[pt1]), c(den$x[pt1-1],den$y[pt1-1]), e)   
#                mean.yz=mean(y - expit(cbind(Z, x=e, 0) %*% coef.hat[-p])) # note this parameterization is based on coef.hat
#                V.2[p, p]= (beta.4-beta.3) * f.x.e * mean.yz
    
