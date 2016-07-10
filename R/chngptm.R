# tol=1e-4; maxit=1e2; verbose=TRUE; est.method="smoothapprox"; lb.quantile=.1; ub.quantile=.9; grid.size=500; weights=NULL; chngpt.init=NULL; alpha=0.05
chngptm = function(formula.1, formula.2, family, data,
  type=c("step","hinge","segmented","stegmented"), 
  est.method=c("default","smoothapprox","grid"), 
  var.type=c("none","robust","model","smooth","robusttruth","bootstrap","all"), aux.fit=NULL, test.inv.ci=TRUE,
  lb.quantile=.1, ub.quantile=.9, grid.size=500, ci.bootstrap.size=500, alpha=0.05, save.boot=FALSE, m.out.of.n=FALSE,
  b.=-30,
  tol=1e-4, maxit=1e2, chngpt.init=NULL, 
  weights=NULL, verbose=FALSE,
  ...) 
{
    
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)    
    var.type<-match.arg(var.type)    
    est.method<-match.arg(est.method)    
    
    if (var.type %in% c("robust","robusttruth","all") & is.null(aux.fit)) stop("need an aux.fit for robust variance estimate")
    
    # remove missing observations
    form.all = update(formula.1, formula.2)
    subset. = complete.cases(model.frame(form.all, data, na.action=na.pass))
    data=data[subset.,,drop=FALSE]
    
    if (est.method=="default") {
        est.method = if(family=="binomial" & is.null(weights)) "smoothapprox" else est.method = "grid"
    }
    
    # decide whether there is interaction
    y=model.frame(formula.1, data)[,1]
    Z=model.matrix(formula.1, data)
    tmp=model.matrix(formula.2, data)[,-1,drop=F]
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    chngpt.var.sorted=sort(chngpt.var)
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    
    n=nrow(Z)
    p.z=ncol(Z)
    p.2=switch(type, step=1, hinge=1, segmented=2, stegmented=3)
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    p=p.z+p.2.itxn+1 #total number of paramters, including threshold
    
    if (is.null(weights)) weights=rep(1,n) 
    data$weights=weights # try put it in data to be found by glm
    
    # make formula that includes all parameters but threshold
    formula.new = if (type %in% c("segmented","stegmented")) update(formula.1, as.formula("~.+"%+%chngpt.var.name)) else formula.1
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    formula.new=update(formula.new, as.formula(f.alt))
    
    if (verbose) {
        myprint(type, est.method, has.itxn)     
        print(formula.new)
        myprint(p.z, p.2.itxn, p)
    }
    
    search.bound=10
    if (est.method=="grid") {
    #### grid search

        # change point candidates
        chngpts=chngpt.var.sorted[chngpt.var.sorted<quantile(chngpt.var.sorted, ub.quantile) & chngpt.var.sorted>quantile(chngpt.var.sorted, lb.quantile)]
        if (length(chngpts)>grid.size) chngpts=chngpts[round(seq(1,length(chngpts),length=grid.size))]
                
        logliks=sapply (chngpts, function(e) {
            data=make.chngpt.var(chngpt.var, e, type, data)
            fit = do.regression (formula.new, data, weights, family)
            if(length(fit$warning)!=0) {
                if(verbose) print(fit$warning)
                NA
            } else as.numeric(logLik(fit$value))
        } )
        glm.warn=any(is.na(logliks))
        e=chngpts[which.max(logliks)]
        if(verbose==2) {
            plot(chngpts, logliks, type="b", xlab="change point")
            abline(v=e)
        }
        data = make.chngpt.var(chngpt.var, e, type, data)
        fit = do.regression (formula.new, data, weights, family)$value
        coef.hat=c(coef(fit), "chngpt"=e)  
        best.fit=fit
             
        
    } else if (est.method=="smoothapprox") {
    #### Newton-Raphson
        
        # do a test to get init value for change point to be used in estimation 
        if (is.null(chngpt.init)) {
            chngpt.test.0 = chngpt.test (formula.1, formula.2, family=family, data, type=type, compute.p.value=FALSE)# note main.method is lr by default, which is faster
            e.init=chngpt.test.0$chngpt
            if (verbose==2) plot(chngpt.test.0, xlim=range(chngpt.var))
        } else {
            e.init=chngpt.init
        }
        names(e.init)="e"
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
            data = make.chngpt.var(chngpt.var, e.init, type, data, b.)
            fit.0 = do.regression(formula.new, data, weights, family)$value
            #if (verbose==2) print(coef(fit.0))
            
            # update threshold and associated coefficients
            beta.init=coef(fit.0)[p.z+1:p.2.itxn]; #names(beta.init)=c("beta1","beta2")# we may need this for variance calculation to work
            alpha.hat=coef(fit.0)[1:p.z]
            alpha.z = c(Z %*% alpha.hat)
            #if (verbose) myprint(beta.init, e.init, b., digits=6)
            
            optim.out = optim(par=c(beta.init, e.init), 
                  fn = get("dev."%+%type%+%"."%+%ifelse(has.itxn,"itxn.","")%+%"f"), 
                  gr = NULL,
                  #gr = get("dev."%+%type%+%"."%+%ifelse(has.itxn,"itxn.","")%+%"deriv.f"), 
#                  # if we use analytical gradient function by deriv3 in optim, we can get situations like exp(100), which will be Inf, and Inf/Inf will be NaN
#                  fn = function(theta,...) sum(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...)), 
#                  gr = function(theta,...) colSums(attr(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
                  chngpt.var, y, b., alpha.z, z.1,
                  lower = c(rep(-search.bound,length(beta.init)), quantile(chngpt.var, lb.quantile)), 
                  upper = c(rep( search.bound,length(beta.init)), quantile(chngpt.var, ub.quantile)), 
                  method="L-BFGS-B", control = list(), hessian = TRUE)
            
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
            if (max(abs(coef.tmp-coef.hat))<tol) {
                coef.hat=coef.tmp
                break
            } else {
                coef.hat=coef.tmp
            }
            
        } # end while 
        data = make.chngpt.var(chngpt.var, last(coef.hat), type, data, b.)
        fit = do.regression (formula.new, data, weights, family)$value
        coef.hat=c(coef(fit), "chngpt"=last(coef.hat))  
        best.fit=fit     
    
    } # end if est.method    
    names(coef.hat)[length((coef.hat))]="chngpt"
    if (type=="stegmented") {
        replacement="I("%+%chngpt.var.name%+%">chngpt)"
        replacement.2="("%+%chngpt.var.name%+%"-chngpt)+"
        new.names=sub("x.gt.e.2", replacement.2, names(coef.hat))    
        new.names=sub("x.gt.e", replacement, new.names)    
    } else {
        if(type=="step") {
            replacement="I("%+%chngpt.var.name%+%">chngpt)"
        } else if (type=="hinge" | type=="segmented") {
            replacement="("%+%chngpt.var.name%+%"-chngpt)+"
        } 
        new.names=sub("x.gt.e", replacement, names(coef.hat))    
    }
    if (verbose) print(new.names)
    names(coef.hat)=new.names
    
    
    
    ###############################################################################################
    # variance-covariance 
    
    if (!has.itxn & family=="binomial") {
    
        if (type %in% c("hinge","segmented")) { # continuous change point models
            
            var.est.smooth=function(){
                # if perfect segregation happens, cann't compute variance estimate
                if (glm.warn | any(abs(coef.hat[-1])>=search.bound) ) {
                    var.est=matrix(NA, nrow=length(coef.hat), ncol=length(coef.hat))
                    rownames(var.est) <- colnames(var.est) <- names(coef.hat)
                    warning("cannot estimate covariance")
                    return (var.est)
                }
    
                # expressions for use with optim with deriv3, they are different from the set of expressions, dev.step etc, b/c each alpha has to be separate
                alpha.z.s="("%+% concatList("alpha"%+%1:p.z%+%"*z"%+%1:p.z%+%"","+") %+%")"
                if (type=="hinge") {
                    deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + (x-e)*beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - (x-e)*beta/(1+exp(b*(x-e))) ) ) "
                    params=c("alpha"%+%1:p.z, "beta", "e")
                } else if (type=="segmented") {
                    deviance.s <- " (1-y) * ( "%+% alpha.z.s %+%" + beta1*x + (x-e)*beta2/(1+exp(b*(x-e))) )  +  log( 1 + exp( -"%+% alpha.z.s %+%" - beta1*x - (x-e)*beta2/(1+exp(b*(x-e))) ) ) "
                    params=c("alpha"%+%1:p.z, "beta1", "beta2", "e")
                }
                params.long=c(params,"x","y","b","z"%+%1:p.z,"z.1")
                if (verbose) print(deviance.s)
                if (verbose) myprint(params)
                loss.f=deriv3(parse(text=deviance.s), params, params.long)    
                param.list = c(as.list(coef.hat), list(chngpt.var), list(y), list(b.), lapply(1:ncol(Z), function (i) Z[,i]), list(z.1))
                names(param.list)=params.long    
                tmp=do.call(loss.f, param.list)
                hess=apply(attr(tmp,"h"), 2:3, sum, na.rm=T)       
                #print("smooth"); print(eigen(hess))         
                var.est = try(solve(hess)) # should keep change point in, and not do hess[-ncol(hess), -ncol(hess)], otherwise lead to over estimation of sd
                
                rownames(var.est) <- colnames(var.est) <- names(coef.hat)
                var.est
            }
            
            # model-based
            var.est.model=function(test.inv.ci){
                if(verbose>=2) print("in var.est.model")
                # if perfect segregation happens, cann't compute variance estimate
                if (glm.warn | any(abs(coef.hat[-1])>=search.bound) ) {
                    var.est=matrix(NA, nrow=length(coef.hat), ncol=length(coef.hat))
                    rownames(var.est) <- colnames(var.est) <- names(coef.hat)
                    attr(var.est,"chngpt.ci")=c(NA,NA)
                    return (var.est)
                }
                                
                # set up some variables
                e=coef.hat["chngpt"]
                x.gt.e=chngpt.var>e
                beta.4=coef.hat["("%+%chngpt.var.name%+%"-chngpt)+"]
                
                # x.tilda
                x.tilda=cbind(Z, if (type=="segmented") chngpt.var, (chngpt.var-e)*x.gt.e, -beta.4*x.gt.e)
                p.est=drop(expit(x.tilda[,-p] %*% coef.hat[-p]))
                                              
                V.1=-tXDX(x.tilda, p.est*(1-p.est))/n    
                V=V.1
                #myprint(e)
                #print(eigen(V.1)$values)
                #print(x.tilda[order(chngpt.var),])
                var.est=solve(-V)/n                        
                rownames(var.est) <- colnames(var.est) <- names(coef.hat)                
                
                # profile likelihood ratio test inversion CI
                if(test.inv.ci) {
                    chngpt.ci=ci.test.inv(1)
                    attr(var.est,"chngpt.ci")=chngpt.ci
                }
                
                var.est
            }

            var.est.robust=function(aux.fit, test.inv.ci, true.param=FALSE){
                if(verbose>=2) print("in var.est.robust")
                # if perfect segregation happens, cann't compute variance estimate
                if (glm.warn | any(abs(coef.hat[-1])>=search.bound) ) {
                    var.est=matrix(NA, nrow=length(coef.hat), ncol=length(coef.hat))
                    rownames(var.est) <- colnames(var.est) <- names(coef.hat)
                    attr(var.est,"chngpt.ci")=c(NA,NA)
                    if(verbose>=2) print("return NA")
                    return (var.est)
                }
                                
                # compute density of X at change point. have to use coef.hat["chngpt"] for change point here b/c e is defined later
                den=density(chngpt.var)
                pt1=which(den$x>=coef.hat["chngpt"])[1]
                f.x.e = interpolate(c(den$x[pt1],den$y[pt1]), c(den$x[pt1-1],den$y[pt1-1]), coef.hat["chngpt"])   # interploate to find the density at change point
                
                # use population parameter values to calculate variance for debugging
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
                p.est=drop(expit(x.tilda[,-p] %*% coef.hat[-p]))
    
                V.1=-tXDX(x.tilda, p.est*(1-p.est))/n    
                # V.2
                V.2=matrix(0,nrow=p,ncol=p)
                # there are two ways to compute the off-diagonal element. when use true coefficients, it works better to use predicted response by true model; when use estimated, using y works better
                #V.2[p,p-1]<-V.2[p-1,p]<- if(is.null(attr(aux.fit,"truemodel"))) -mean((y-p.est)*x.gt.e) else -mean((predict(aux.fit, data, "response")-p.est)*x.gt.e)
                V.2[p,p-1]<-V.2[p-1,p]<- -mean((y-p.est)*x.gt.e) 
    
                newdata=data; newdata[[chngpt.var.name]]=e
                m.0=mean(predict(aux.fit, newdata, "response"))
                p.e.z.0=mean(expit(cbind(Z, if (type=="segmented") e) %*% coef.hat[1:(p-2)]))                
                V.2[p, p]= beta.4 * f.x.e * (m.0-p.e.z.0)
                V=V.1+V.2    
                V.inv=solve(V) 
                M= tXDX((y-p.est)*x.tilda, rep(1,n))/n 
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
                    if (verbose==2) {
                        myprint(last(diag(var.est)), -last(diag(V.inv))/n, scale.chisq)
                    }
                    if (scale.chisq<0) {
                        warning("scale.chisq is negative")
                        scale.chisq=abs(scale.chisq)
                    } 
                    chngpt.ci=ci.test.inv(scale.chisq)
                    attr(var.est,"chngpt.ci")=chngpt.ci
                }
    
                var.est
            }

            ci.test.inv=function(scale.chisq=1){
                # if perfect segregation happens, cann't compute variance estimate
                if (glm.warn | any(abs(coef.hat[-1])>=search.bound) ) {# -1 gets rid of intercept
                    warning("cannot estimate covariance")
                    return (c(NA,NA))
                }
                
                if (est.method=="grid") {
                    
                } else if (est.method=="smoothapprox") {
                    
                    idx.chngpt=which(chngpt.var.sorted>=coef.hat["chngpt"])[1]
                    data = make.chngpt.var(chngpt.var, chngpt.var.sorted[idx.chngpt], type, data, b.)
                    fit=do.regression (formula.new, data, weights, family)
                    lik.max = as.numeric(logLik(fit$value))
                    
                    c.alpha= qchisq(1-alpha, df=1) * scale.chisq
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.)
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], type, data, b.)
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

            # the index are hardcoded in this function
            ci.bootstrap=function(){
                if(verbose>=2) print("in ci.bootstrap")
                # if perfect segregation happens, cann't compute variance estimate
                if (glm.warn | any(abs(coef.hat[-1])>=search.bound) ) {
                    warning("cannot estimate covariance")
                    return (NULL)
                }
                
                # save rng state before set.seed in order to restore before exiting this function
                save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                if (class(save.seed)=="try-error") {        
                    set.seed(1)
                    save.seed <- get(".Random.seed", .GlobalEnv)
                }                        
                boot.out=boot(data, 
                    statistic=function(dat, ii){
                        if (m.out.of.n) ii=ii[1:(4*sqrt(n))] # m out of n bootstrap
                        fit.b=try(chngptm (formula.1, formula.2, family, dat[ii,], type, est.method="smoothapprox", var.type="none"))
                        #fit.b=try(chngptm (formula.1, formula.2, family, dat[ii,], type, est.method="smoothapprox", var.type="robust", aux.fit=aux.fit, test.inv.ci=FALSE))# stud performs pretty badly
                        # note that we cannot use | in the following line
                        out = if(inherits(fit.b,"try-error")) {
                            rep(NA,length(coef.hat)*1) # need to adjust depending on below
                        }
                        else if (any(abs(fit.b$coefficients[1:(length(fit.b$coefficients)-1)])>search.bound)) 
                            rep(NA,length(coef.hat)*1) # need to adjust depending on below
                        else {
                            c(fit.b$coefficients)#, diag(fit.b$vcov))
                        }
                        out
                        
                    }, 
                    R=ci.bootstrap.size, sim = "ordinary", stype = "i"
                )
                # restore rng state 
                assign(".Random.seed", save.seed, .GlobalEnv)                 
                
                # go through parameter list and get different types of confidence intervals
                outs=list(perc=NULL, basic=NULL, bc=NULL)
                for (i in 1:length(coef.hat)){
                    ci.out=try(boot.ci(boot.out, type=c("perc","basic",if(ci.bootstrap.size>n) "bca"), index=c(i, i+length(coef.hat))))
                    #ci.out=try(boot.ci(boot.out, type=c("perc","basic",if(ci.bootstrap.size>n) "bca","stud"), index=c(i, i+length(coef.hat))))
                    outs$perc= cbind(outs$perc,  if(inherits(ci.out,"try-error")) c(NA,NA) else ci.out$percent[1,4:5])
                    outs$basic=cbind(outs$basic, if(inherits(ci.out,"try-error")) c(NA,NA) else ci.out$basic[1,4:5])
                    #outs$stud= cbind(outs$stud,  if(inherits(ci.out,"try-error")) c(NA,NA) else ci.out$student[1,4:5])
                    # for bca to work, the number of bootstrap replicates needs to be greater than sample size
                    outs$abc=  cbind(outs$abc,   if(inherits(ci.out,"try-error")) c(NA,NA) else if(ci.bootstrap.size>n) ci.out$bca[1,4:5] else c(NA,NA) )
                }
                # find bc CI from boot.out$t as ci.out does not provide bc
                boot.samples=boot.out$t
                outs$bc=sapply (1:length(coef.hat), function(i) {
                    z.0=qnorm(mean(boot.samples[,i]<coef.hat[i], na.rm=TRUE))                    
                    quantile(boot.samples[,i], c(pnorm(2*z.0+qnorm(alpha/2)), pnorm(2*z.0+qnorm(1-alpha/2))), na.rm=TRUE)
                })
    
         
                # the following percentile CI is numerically different from boot.ci results because the latter use norm.inter to do interpolation on the normal quantile scale
                #ci.perc=apply(boot.out$t, 2, function (x) quantile(x, c(alpha/2, 1-alpha/2), na.rm=TRUE))   
                                
#                # this takes a long time and returns non-sensible result
#                ci.abc=abc.ci(data, statistic=function(dat, ww){
#                        fit.b=try(chngptm (formula.1, formula.2, family, dat, type, est.method="smoothapprox", var.type="none"))
#                        # note that we cannot use | in the following line
#                        if(inherits(fit.b,"try-error")) 
#                            rep(NA,length(coef.hat)) 
#                        else if (any(abs(fit.b$coefficients[1:(length(fit.b$coefficients)-1)])>search.bound)) 
#                            rep(NA,length(coef.hat)) 
#                        else fit.b$coefficients
#                    }, 
#                    index=4, conf=1-alpha
#                )
#                myprint(ci.abc)
                
                if(save.boot) outs[["boot.samples"]]=boot.samples
                
                outs                
            }

            var.est=switch(var.type, 
                none=NA,
                smooth=var.est.smooth(), 
                model=var.est.model(test.inv.ci), 
                robust=var.est.robust(aux.fit, test.inv.ci), 
                robusttruth=var.est.robust(aux.fit, test.inv.ci, true.param=TRUE), 
                bootstrap=ci.bootstrap(), 
                all=list("smooth"=var.est.smooth(), "model"=var.est.model(test.inv.ci), "robust"=var.est.robust(aux.fit, test.inv.ci) )
            )
            
        } else if (type %in% c("step","stegmented")) {
            # discontinous models
            var.est=vcov(best.fit)
            
        } else stop("type incorrect")
        
    } else {
        var.est=NULL
    }
    
    res=list(
          coefficients=coef.hat
        , vcov=var.est
        , formula.1=formula.1
        , formula.2=formula.2
        , chngpt.var=chngpt.var.name
        , chngpt=coef.hat["chngpt"]
        , est.method = est.method
        , b.=if(est.method=="grid") NULL else b.
        , type = type
        , glm.warn = glm.warn
        , best.fit=best.fit    
    )
    names(res$chngpt) = round(100*mean(chngpt.var<res$chngpt),1) %+% "%"
    
    if (est.method=="smoothapprox") {
        res=c(res, list(converged=converged, iter=n.iter))
    } else if (est.method=="grid") {
        res=c(res, list(chngpts=chngpts, logliks=logliks))
#        # good.soln needs to take as input a chngptm object that has chngpts etc
#        class(res)=c("chngptm", class(res)) 
#        tmp=good.soln(res, plot=verbose==3)
#        res=c(res, list(good.soln=tmp) )
    }
    
    class(res)=c("chngptm", class(res))
    res    
}


predict.chngptm=function (object, newdata = NULL, type = c("link", "response", "terms"), ...){    
    newdata = make.chngpt.var(newdata[[object$chngpt.var]], object$chngpt, object$type, newdata, object$b.)            
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
summary.chngptm=function(object, ...) {    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"        
    fit=object
    p.z=length(fit$coefficients)
    n=nrow(fit$best.fit$data)
    type=fit$type
    
    boot.conf=FALSE
    
    if (is.null(fit$vcov)) {
        cat("No variance estimate available.\n\n")
        print(fit)
        return (invisible())
    } else {
        if(is.list(fit$vcov)){
            if (!is.null(fit$vcov$robust)) {
                vcov=fit$vcov$robust
            } else if (!is.null(fit$vcov$abc)){
                boot.conf=TRUE
                vcov=fit$vcov$abc
                #if (all(is.na(vcov))) vcov=fit$vcov$bc # if the bootstrap samples have too many NA, only bc quantiles can be estimated. However, bc quantiles may be misleading
            }
            
        } else {
            vcov=fit$vcov
        }
    }
    
    if(object$type %in% c("hinge","segmented") & !boot.conf) {
        vcov.t=vcov[p.z,p.z]
        vcov=vcov[-p.z,-p.z] # not return the chngpoint estimate
    } else  {
        vcov
    }
    #print(vcov)
    
    # assuming the last of coefficients is always the change point
    # deal with coefficients and change point separately
    
    res=list()
    
    # coefficients
    if (boot.conf){
        lb=exp(vcov[1,])
        ub=exp(vcov[2,])
        pval=rep(NA,p.z)        
    } else  {
        lb=exp(unname(fit$coefficients[1:(p.z-1)] - sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        ub=exp(unname(fit$coefficients[1:(p.z-1)] + sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        pval=unname(pt(abs(fit$coefficients[1:(p.z-1)] / sqrt(diag(vcov))), df=n-p.z, lower.tail=FALSE))
    }
    res$coefficients=mysapply(1:(p.z-1), function (i) {
        c(
              "OR"=exp(unname(fit$coefficients[i]))
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
    
    # change point
    i=p.z
    if (type %in% c("hinge","segmented")) {
        if (boot.conf){
            lb=vcov[1,p.z]
            ub=vcov[2,p.z]
        } else if(!is.null(attr(vcov,"chngpt.ci"))) {
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


# make.chngpt.var is needed in both testing and estimation
# in estimation by grid search, b. is set to null when this function is called
# in estimation by smooth approx, b. value is saved to object
make.chngpt.var=function(x, e, type, data=NULL, b.=NULL) {    
    if(type=="step") {
        out=ifelse(x>=e, 1, 0)
    } else if (type=="hinge") {
        out=ifelse(x>=e, x-e, 0)
    } else if (type=="segmented") {
        out=ifelse(x>=e, x-e, 0) # x also becomes part of the null model
    } else if (type=="stegmented") {
        out=cbind(ifelse(x>=e, 1, 0), ifelse(x>=e, x-e, 0))  # x also becomes part of the null model
    }
    
#    print(1/(1+exp(b.*(x-e))))
#    print(out)
    if (!is.null(b.)) out=out * 1/(1+exp(b.*(x-e))) 
#    print(out)
    
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
    
get.f.alt=function(type, has.itxn, z.1.name, chngpt.var.name) {
    if (type=="stegmented") {
        f.alt=if(has.itxn) "~.+(x.gt.e+x.gt.e.2)*"%+%z.1.name%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name else "~.+x.gt.e+x.gt.e.2"
    } else if (type=="segmented") {
        #f.alt=if(has.itxn) "~.+x.gt.e*"%+%z.1.name%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name else "~.+x.gt.e"
        f.alt=if(has.itxn) "~."%+%"+"%+%chngpt.var.name%+%":"%+%z.1.name%+%"+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    } else if (type %in% c("step","hinge")) {
        f.alt=if(has.itxn) "~.+x.gt.e*"%+%z.1.name else "~.+x.gt.e"
    }
    f.alt
}
        
#deviance.3pl <- expression( (1-y) * ( c + (d-c)/(1+exp(b*(t-e))) )  +  log( 1 + exp( -c - (d-c)/(1+exp(b*(t-e))) ) ) )
#deviance.3pl.deriv=deriv3(deviance.3pl, c("c","d","e"), c("c","d","e","t","y","b"))
    
# 
# step change point model
dev.step <- expression( (1-y) * ( alpha.z + beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - beta/(1+exp(b*(x-e))) ) ) )
dev.step.deriv=deriv3(dev.step, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.step.itxn <- expression( (1-y) * ( alpha.z + (beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) )
dev.step.itxn.deriv=deriv3(dev.step.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))
# hinge change point model
dev.hinge <- expression( (1-y) * ( alpha.z + (x-e)*beta/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*beta/(1+exp(b*(x-e))) ) ) )
dev.hinge.deriv=deriv3(dev.hinge, c("beta","e"), c("beta","e","x","y","b","alpha.z"))
dev.hinge.itxn <- expression( (1-y) * ( alpha.z + (x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) )  +  log( 1 + exp( -alpha.z - (x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) ) ) )
dev.hinge.itxn.deriv=deriv3(dev.hinge.itxn, c("beta1","beta2","e"), c("beta1","beta2","e","x","y","b","alpha.z","z.1"))
# the following functions are used as objective functions in smoothapprox, b/c we don't supply derivative to optim anymore. 
# step change point model
dev.step.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta=theta[1]; e=theta[2]
    eta=beta/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
# dev.step.deriv.f is experimental, not really used, does not seem to bring improvement when compared to gr=NULL
dev.step.deriv.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta=theta[1]; e=theta[2]
    approxi=1/(1+exp(b*(x-e)))
    dBde=b*approxi*(1-approxi)
    pi=expit(alpha.z+beta*approxi)
    colSums((y-pi)*cbind(x>e, beta*dBde))
}
dev.step.itxn.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(beta1+beta2*z.1)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
# hinge change point model
dev.hinge.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta=theta[1]; e=theta[2]
    eta=(x-e)*beta/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
dev.hinge.itxn.f <- function(theta,x,y,b,alpha.z,z.1)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(x-e)*(beta1+beta2*z.1)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
}
# segmented change point model
dev.segmented.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + (x-e)*beta2/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
dev.segmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*(x-e)/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
# stegmented change point model
dev.stegmented.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[2]; beta3=theta[3]; e=theta[4]
    eta=beta1*x + ((x-e)*beta3 + beta2)/(1+exp(b*(x-e))) 
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
} 
# the order of parameters in the following needs to be fixed
dev.stegmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; beta5=theta[5]; beta6=theta[6]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + ((beta3+beta4*z.1)*(x-e) + beta4+beta5*z.1)/(1+exp(b*(x-e)))
    linear=alpha.z+eta
    sum( (1-y)*linear + log(1+exp(-linear)) )
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
    
