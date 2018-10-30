# tol=1e-4; maxit=1e2; verbose=TRUE; est.method="fastgrid"; var.type="bootstrap"; lb.quantile=.1; ub.quantile=.9; grid.search.max=500; weights=NULL; chngpt.init=NULL; alpha=0.05; b.transition=Inf; search.bound=10; ci.bootstrap.size=500; m.out.of.n=F; aux.fit=NULL        
#family can be coxph or any glm family, but variance estimate is only available for binomial and gaussian (only model-based for latter)
chngptm = function(formula.1, formula.2, family, data, type=c("step","hinge","upperhinge","segmented","segmented2","stegmented"), # segmented2 is the model studied in Cheng 2008
  formula.strat=NULL,
  weights=NULL, # this arg is kept here due to historical reasons
  est.method=c("default","smoothapprox","grid","fastgrid","fastgrid2","gridC"), # gridC is for development use only 
  var.type=c("none","robust","model","robusttruth","bootstrap","all"), aux.fit=NULL, 
  lb.quantile=.1, ub.quantile=.9, grid.search.max=5000, 
  test.inv.ci=TRUE, boot.test.inv.ci=FALSE, # test.inv.ci is passed to local functions, boot.test.inv.ci is global within this function
  ci.bootstrap.size=1000, alpha=0.05, save.boot=FALSE, m.out.of.n=FALSE, # grid.search.max is the maximum number of grid points used in grid search
  b.transition=Inf,# controls whether threshold model or smooth transition model
  tol=1e-4, maxit=1e2, chngpt.init=NULL, search.bound=10,
  keep.best.fit=TRUE, # best.fit is needed for making prediction and plotting
  verbose=FALSE, ...) 
{
    
    if (missing(type)) stop("type missing")
    threshold.type<-match.arg(type)  
    var.type<-match.arg(var.type)    
    est.method<-match.arg(est.method)    
        
    # remove missing observations
    subset. = complete.cases(model.frame(update(formula.1, formula.2), data, na.action=na.pass))
    data=data[subset.,,drop=FALSE]
    
    # set b.search
    stopifnot(b.transition>0)
    b.search=if(b.transition==Inf) 30 else b.transition
    if(threshold.type=="segmented2" & b.transition==Inf) stop("for segmented2, b.transition should not be Inf") # limited implementation for segmented2
    
    # not all variance estimation methods are implemented for all families    
    if(!family %in% c("gaussian","binomial")) {
        if(!var.type %in% c("bootstrap","none")) stop("no analytical variance estimates are provided for this family: "%.%family)
    }     
    if (var.type %in% c("robust","robusttruth","all") & is.null(aux.fit)) stop("need an aux.fit for robust variance estimate")
    
    # extract design matrix, x, and y
    y=model.frame(formula.1, data)[,1] #if lhs of formula.1 is cbind(,), y is a matrix
    Z=model.matrix(formula.1, data)
    tmpdat=model.matrix(formula.2, data)[,-1,drop=F]
    chngpt.var.name=setdiff(colnames(tmpdat), colnames(Z))[1];     
    if(verbose) cat("variable to be thresholded: ", chngpt.var.name, "\n")
    if(is.na(chngpt.var.name)) stop("Something is wrong. Check the formula.")
    chngpt.var = tmpdat[,chngpt.var.name]
    # note that in chngpt.test, Z may already include chngptvar if threshold.type is segmented
    
    # interaction
    z.1.name=intersect(colnames(tmpdat), colnames(Z))
    has.itxn = length(z.1.name)>0
    z.1 = tmpdat[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    
    # stratification, a.k.a. covariate-dependent threshold regression
    if (is.null(formula.strat)) {
        stratified.by=NULL
        stratified=FALSE
    } else {
        if(!threshold.type %in% c("hinge","segmented")) stop("stratification only implemented for hinge and segmented threshold effects")
        if(class(formula.strat)!="formula") stop("formula.strat needs to be a formula")
        strat.dat=model.matrix(formula.strat, data)[,-1,drop=F]
        if (ncol(strat.dat)!=1) {str(strat.dat); stop("something is wrong with formula.strat") }
        stratified.by=strat.dat[,1]
        if (!all(sort(unique(stratified.by))==c(0,1))) {stop("Only 0/1 or T/F stratification variables are supported for now.") }
        stratified=TRUE
        if(!est.method %in% c("default","fastgrid2","grid")) stop("est.method not supported for stratified")
        if(!family %in% c("gaussian")) stop("family not supported for stratified")
        n.0=sum(stratified.by==0)
        n.1=sum(stratified.by==1)
    }
    
    # dimensions
    n=nrow(Z)
    p.z=ncol(Z)
    p.2=switch(threshold.type, step=1, hinge=1, upperhinge=1, segmented=2, segmented2=2, stegmented=3)
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    p=p.z+p.2.itxn+1 #total number of paramters, including threshold
    if(stratified) {
        # if these dimensions are needed, e.g. in smoothapprox, they need to be updated
    }
        
    # formula.new is the model to fit if threshold is known
    formula.new = if (threshold.type %in% c("segmented","segmented2","stegmented")) update(formula.1, as.formula("~.+"%.%chngpt.var.name)) else formula.1
    formula.new=update(formula.new, get.f.alt(threshold.type, chngpt.var.name, modified.by=if(has.itxn) z.1.name else NULL, stratified=stratified))
    
    # weights
    # w.all.one is needed in the next block. 
    # chngpt.glm.weights goes with the data because glm() is incredible that if the weights arg does not have global scope, it needs to be a variable name in the data frame
    # chngpt.glm.weights and does not get updated in the next block while weights are updated
    if (is.null(weights)) weights=rep(1,nrow(data)) 
    data$chngpt.glm.weights=weights
    w.all.one=all(weights==1)     
    # create a copy of chngpt.glm.weights in the function scope so that do.regression works
    chngpt.glm.weights=data$chngpt.glm.weights  
    do.regression=function(formula.new, data, family){
        if (family=="coxph") {
            fit = keepWarnings(survival::coxph(formula.new, data=data, weights=chngpt.glm.weights)) # interestingly and I don't understand why weights=data$chngpt.glm.weights fails
        } else {
    #    str(formula.new); str(data); str(chngpt.glm.weights); str(family)
            fit =             keepWarnings(glm(formula.new, data=data, weights=chngpt.glm.weights, family=family) )
        }
        fit
    }
    
    # deal with cbind() in lhs in logistic regression 
    # convert cbind() to single column and update weights (but not chngpt.glm.weights, b/c formula is still cbind())
    if(is.matrix(y)) {
        n.tmp <- y[,1]+y[,2]
        y <- ifelse(n.tmp == 0, 0, y[,1]/n.tmp)        
        weights <- n.tmp*weights
        w.all.one=FALSE
    }
    
    # est.method
    # smoothapprox cannot work when lhs of formula.1 is cbind() 
    # for now, fastgrid and C version of grid is implemented only for the following scenarios
    fastgrid.ok = family=="gaussian" & threshold.type %in% c("hinge","upperhinge","segmented") & !has.itxn
    if ((!family%in%c("binomial","gaussian")) & est.method=="smoothapprox") stop ("smoothapprox only implemented for binomial and guassian families and when lhs of formula is not cbind()")
    if(est.method%in%c("fastgrid","fastgrid2") & !fastgrid.ok) {
        warning ("Fast grid search only implemented for guassian family and hinge/upperhinge/segmented models without interaction. Switching to grid search ...")
        est.method="grid"
    }
    if (est.method=="default") {
        if(fastgrid.ok) { est.method="fastgrid2" 
        } else if(family=="binomial") { est.method="smoothapprox"
        } else { est.method="grid"
        }        
    }
    if (fastgrid.ok & !est.method%in%c("fastgrid","fastgrid2") & var.type!="none") cat("Note that for linear models, fastgrid2 is the best est.method.\n") # var.type condition is added so that this is not printed during bootstrap
    
    # sort data
    if(!stratified) {
        order.=order(chngpt.var)
        chngpt.var.sorted=chngpt.var[order.]
        Z.sorted=Z[order.,,drop=FALSE]
        y.sorted=y[order.]
        w.sorted=weights[order.]
        data.sorted=data[order.,]
        stratified.by.sorted=NULL
    } else {
        order.=order(stratified.by, chngpt.var)
        chngpt.var.sorted=chngpt.var[order.]
        Z.sorted=Z[order.,,drop=FALSE]
        y.sorted=y[order.]
        w.sorted=weights[order.]
        data.sorted=data[order.,]
        stratified.by.sorted=stratified.by[order.]
    }
    
    if (verbose) {myprint(family, threshold.type, est.method, var.type); myprint(has.itxn, p.z, p.2.itxn, p); print(formula.new)}
        
    # threshold candidates
    chngpts=get.chngpts(chngpt.var.sorted,lb.quantile,ub.quantile,n.chngpts=grid.search.max, stratified.by.sorted=stratified.by.sorted)
    
    
    grid.search=function(){
        if(verbose) myprint(fastgrid.ok, est.method)
        
        if(fastgrid.ok & est.method%in%c("gridC", "fastgrid", "fastgrid2")) {
#            # gridC:     Only used for a fair comparison of speed with fastgrid. It is not as robust as grid, which is an R implementation that handles weights more easily
#            # for upperhinge fastgrid2, there is an R implementation for a simple case: linear, without weights. It is even slower than grid search, 
#            # but it implements a good algorithm and can be used to debug fastgrid2
#            logliks=get.liks.upperhinge (y.sorted,Z.sorted,chngpt.var.sorted,chngpts)   
        
#            # for fastgrid2 development
#            A <- solve(t(Z.sorted) %*% Z.sorted)
#            H <- Z.sorted %*% A %*% t(Z.sorted)
#            r <- as.numeric((diag(n)-H) %*% y.sorted)
#            # find sqrt of A
#            a.eig <- eigen(A)   
#            A.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% t(a.eig$vectors) # solve can be replaced by t b/c A is symmetrical
#            B = Z.sorted %*% A.sqrt
                                    
#                # for fastgrid2 development, the first two argments may be replaced by:
#                cbind(B, chngpt.var.sorted), #design matrix, the last column is to be used as (x-e)+ in the C program
#                as.double(r), 
            
            # computes Y' H_e Y
            if(!stratified) {
                f.name=est.method%.%"_search" # weirdly, f.name cannot be put inline b/c rcmdcheck throws an error otherwise
                logliks = .Call(f.name, 
                    cbind(Z.sorted, chngpt.var.sorted, if(threshold.type=="segmented") chngpt.var.sorted), #design matrix, the last column is to be used as (x-e)+ in the C program
                    as.double(y.sorted), 
                    as.double(w.sorted), 
                    w.all.one, 
                    as.integer(attr(chngpts,"index")), # potential thresholds
                    0,# bootstrap size
                    threshold.type=="upperhinge"
                )             
            } else {
                #cbind(chngpt.var.sorted*c(rep.int(0,n.0),rep.int(1,n.1)), chngpt.var.sorted*c(rep.int(1,n.0),rep.int(0,n.1))), # added col
                logliks = .Call("twoD_search", 
                    cbind(Z.sorted, if(threshold.type=="segmented") chngpt.var.sorted), # X
                    as.double(chngpt.var.sorted), # threshold variable
                    as.double(y.sorted), 
                    as.double(w.sorted), 
                    w.all.one, 
                    as.integer(stratified.by.sorted),
                    as.double(lb.quantile),
                    as.double(ub.quantile),
                    0,# bootstrap size
                    threshold.type=="upperhinge"
                )
                logliks=matrix(logliks, nrow=length(chngpts[[1]]), ncol=length(chngpts[[2]]))
            }
            # remember to to put as.double around y.sorted since if y.sort is integer, this throws an error b/c the cpp function expects real 
            # don't put as.double around a matrix since changes the matrix into a vector, which affects the ordering of elements
            
        } else {
        
            # call glm repeatedly
            logliks = if(!stratified) {
                sapply (chngpts, function(e) {
                    data=make.chngpt.var(chngpt.var, e, threshold.type, data, b.transition, stratified.by)
                    fit = do.regression (formula.new, data, family)            
                    if(length(fit$warning)!=0 & verbose>=3) {myprint(e); print(fit$warning); print(coef(fit$value)); print((logLik(fit$value)))}                 
                    #cat(sum(resid(fit$value)**2), " ") # to compare with fastgrid, uncomment the next line
                    if(family=="gaussian") sum(w.sorted*y.sorted**2)-sum(resid(fit$value)**2) else as.numeric(logLik(fit$value))            
                })
            } else {
                sapply (chngpts[[2]], simplify="array", function(f) {
                sapply (chngpts[[1]], simplify="array", function(e) {
                    data=make.chngpt.var(chngpt.var, c(e,f), threshold.type, data, b.transition, stratified.by)
#                    myprint(e,f)
#                    print(formula.new)
#                    print(chngpts)
#                    print(data[order.,c("z","x","x.mod.e","x.mod.f")])
#                    stop("rest")
                    fit = do.regression (formula.new, data, family)            
                    if(length(fit$warning)!=0 & verbose>=3) {myprint(e); print(fit$warning); print(coef(fit$value)); print((logLik(fit$value)))}                 
                    #cat(sum(resid(fit$value)**2), " ") # to compare with fastgrid, uncomment the next line
                    if(family=="gaussian") sum(w.sorted*y.sorted**2)-sum(resid(fit$value)**2) else as.numeric(logLik(fit$value))            
                })
                })
            }
            
        }

        if(!stratified) {
            e.final=chngpts[which.max(logliks)]
        } else {
            e.final=arrayInd(which.max(logliks), .dim=dim(logliks))
            e.final=c(chngpts[[1]][e.final[1]], chngpts[[2]][e.final[2]])
        }
        #print(chngpts); print(logliks); myprint(which.max(logliks)); myprint(e.final)
        #attr(e.final, "glm.warn")=any(is.na(logliks))
        if(verbose>=2) {
            if(!stratified) {
                plot(chngpts, logliks, type="b", xlab="threshold")
                abline(v=e.final)
            } else {
                # cannot use points3d to add a red point after plot
                col=0*logliks; col[which.max(logliks)]=1; col=col+1
                plot3d(cbind(expand.grid(chngpts[[1]],chngpts[[2]]), c(logliks)), xlab="f", ylab="e", zlab="loglik", type="p", col=c(col))
            }
        }
        
        # fit glm using e.final
        data = make.chngpt.var(chngpt.var, e.final, threshold.type, data, b.transition, stratified.by) # note that b.transition is used here instead of b.search
        fit = do.regression (formula.new, data, family)$value
        # return e.final and logliks through attr. If we return them as elements of fit, upon removing these elements, the type of fit changes to a list, not glm anymore
        attr(fit,"e.final")=e.final
        attr(fit,"logliks")=logliks
        fit
    }
    
    # find e.final
    if (est.method %in% c("grid","gridC","fastgrid","fastgrid2")) {
    #### grid search
        best.fit=grid.search()   
        e.final=attr(best.fit,"e.final")        
        logliks=attr(best.fit,"logliks")
        attr(best.fit,"logliks")<-attr(best.fit,"e.final")<-NULL
        coef.hat=c(coef(best.fit), e.final)
        glm.warn=attr(e.final,"glm.warn")  
        
    } else if (est.method=="smoothapprox") {
    #### Newton-Raphson
        
        if(verbose) cat("smoothapprox search\n")
        e.init=chngpt.init
        if (is.null(e.init)) {
            if(verbose) cat("initializing through coarse grid search\n")
            # est.method has to be grid in the following, otherwise it will create an infinite loop
            tmp=chngptm (formula.1, formula.2, family=family, data, type=threshold.type, weights=weights, est.method="grid", lb.quantile=lb.quantile, ub.quantile=ub.quantile, grid.search.max=50, 
                        b.transition=b.transition, search.bound=search.bound, keep.best.fit=FALSE, verbose=ifelse(verbose>3,verbose-3,0))# pass verbose-3 directly seems to fail to do the job 
            e.init=tmp$chngpt
        } 
        
        if(length(e.init)==0) {
            # if there are no initial values for some reason, just do grid search
            e.final=attr(grid.search(), "e.final")
            glm.warn=attr(e.final,"glm.warn")
            
        } else {
            # smooth approximation 
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
                data = make.chngpt.var(chngpt.var, e.init, threshold.type, data, b.search) # b.search is used here to be consistent with optim which uses b.search as well
                fit.0 = do.regression(formula.new, data, family)$value
                
                # update threshold and associated coefficients
                beta.init=coef(fit.0)[p.z+1:p.2.itxn]; #names(beta.init)=c("beta1","beta2")# we may need this for variance calculation to work
                alpha.hat=coef(fit.0)[1:p.z]
                stopifnot(all(!is.na(alpha.hat)))
                alpha.z = c(Z %*% alpha.hat)
    
                # search for better e and slopes associated with x and thresholded x
                optim.out = try(optim(par=c(beta.init, e.init), 
                      fn = get("dev."%.%threshold.type%.%"."%.%ifelse(has.itxn,"itxn.","")%.%"f"), 
                      gr = NULL,
                      #gr = get("dev."%.%threshold.type%.%"."%.%ifelse(has.itxn,"itxn.","")%.%"deriv.f"), 
    #                  # if we use analytical gradient function by deriv3 in optim, we can get situations like exp(100), which will be Inf, and Inf/Inf will be NaN
    #                  fn = function(theta,...) sum(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...)), 
    #                  gr = function(theta,...) colSums(attr(dev.step.itxn.deriv(theta[1],theta[2],theta[3],...), "gradient")), 
                      chngpt.var, y, b.search, alpha.z, z.1, family, weights,
                      lower = c(rep(-search.bound,length(beta.init)), quantile(chngpt.var, lb.quantile)), 
                      upper = c(rep( search.bound,length(beta.init)), quantile(chngpt.var, ub.quantile)), 
                      method="L-BFGS-B", control = list(), hessian = TRUE))
                      
                if(class(optim.out)=="try-error") {
                    if (verbose) cat("error doing smoothapprox search, switch to grid search\n")
                    e.final=attr(grid.search(), "e.final")
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
                } else {
                    coef.hat = coef.tmp
                    e.final=last(coef.hat)
                }
                
            } # end while 
        } # end if else
    
        # fit glm using e.final
        data = make.chngpt.var(chngpt.var, e.final, threshold.type, data, b.transition) # note that b.transition is used here instead of b.search
        best.fit = do.regression (formula.new, data, family)$value
        coef.hat=c(coef(best.fit), "chngpt"=e.final)
        
    } else stop("wrong est.method") # end if grid/smoothapprox
    
    if(!stratified) {
        names(coef.hat)[length(coef.hat)]="chngpt"
        if (threshold.type=="stegmented") {
            replacement="I("%.%chngpt.var.name%.%">chngpt)"
            replacement.2="("%.%chngpt.var.name%.%"-chngpt)+"
            new.names=sub("x.mod.e.2", replacement.2, names(coef.hat))    
            new.names=sub("x.mod.e", replacement, new.names)    
        } else {
            if(threshold.type=="step") {
                replacement="I("%.%chngpt.var.name%.%">chngpt)"
            } else if (threshold.type %in% c("hinge","segmented","segmented2")) {
                replacement="("%.%chngpt.var.name%.%"-chngpt)+"
            } else if (threshold.type %in% c("upperhinge")) {
                replacement="("%.%chngpt.var.name%.%"-chngpt)-"
            } 
            new.names=sub("x.mod.e", replacement, names(coef.hat))    
        }
    } else {
        names(coef.hat)[length(coef.hat)-1:0]=c("chngpt.0","chngpt.1")
        if (threshold.type %in% c("hinge","segmented","segmented2")) {
            new.names=sub("x.mod.e", "("%.%chngpt.var.name%.%".0-chngpt.0)+", names(coef.hat))    
            new.names=sub("x.mod.f", "("%.%chngpt.var.name%.%".1-chngpt.1)+", new.names)    
        } 
    }
    if (verbose) cat(new.names,"\n")
    names(coef.hat)=new.names
    

    # keep an empty line here between estimation and var estimate
    
                
    ###############################################################################################
    # variance-covariance 
    
    warn.var=NULL
    
    # only implemented when there is no interaction and for linear and logistic regression
    if (var.type!="bootstrap" & (has.itxn | !family %in% c("binomial","gaussian"))) {
        var.est=NULL
    } else {
        if (threshold.type %in% c("step","stegmented")) {
            # discontinous models
            var.est=NULL
            
        } else if (threshold.type %in% c("hinge","upperhinge","segmented","segmented2")) { # continuous change point models
            
            
            # model-based
            var.est.model=function(test.inv.ci, robust=FALSE){
                if(verbose) cat("in var.est.model\n")
    
                # set up some variables
                e=coef.hat["chngpt"]
                beta.4=coef.hat["("%.%chngpt.var.name%.%"-chngpt)"%.%ifelse(threshold.type=="upperhinge","-","+")]
    
                # x.tilda
                if(threshold.type %in% c("hinge","segmented")) {
                    x.gt.e=chngpt.var>e
                    if (b.transition==Inf) {
                        x.tilda=cbind(Z, if(threshold.type=="segmented") chngpt.var, (chngpt.var-e)*x.gt.e, -beta.4*x.gt.e)
                    } else {
                        x.tilda=cbind(Z, if(threshold.type=="segmented") chngpt.var, (chngpt.var-e)*expit(b.transition*(chngpt.var-e)), -beta.4*expit(b.transition*(chngpt.var-e))*(1+b.transition*(chngpt.var-e)*(1-expit(b.transition*(chngpt.var-e)))))
                    }               
                } else if(threshold.type %in% c("upperhinge")) {
                    x.lt.e=chngpt.var<e
                    if (b.transition==Inf) {
                        x.tilda=cbind(Z, (chngpt.var-e)*x.lt.e, -beta.4*x.lt.e)
                    } else {
                        stop("var.est.model: not yet implemented")
                        #x.tilda=cbind(Z, (chngpt.var-e)*expit(b.transition*(chngpt.var-e)), -beta.4*expit(b.transition*(chngpt.var-e))*(1+b.transition*(chngpt.var-e)*(1-expit(b.transition*(chngpt.var-e)))))
                    }               
                } else if (threshold.type=="segmented2") {
                    # earlier we checked that b.transition should not be Inf
                    x.tilda=cbind(Z, chngpt.var, chngpt.var*expit(b.transition*(chngpt.var-e)),                               -beta.4*expit(b.transition*(chngpt.var-e))        *b.transition*chngpt.var*(1-expit(b.transition*(chngpt.var-e))) )
                } else stop("threshold.type not recognized in var.est.model")
                
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
                if (threshold.type=="upperhinge") {
                    x.lt.e=chngpt.var<e
                    beta.4=coef.hat["("%.%chngpt.var.name%.%"-chngpt)-"]                    
                    # x.tilda
                    x.tilda=cbind(Z, (chngpt.var-e)*x.lt.e, -beta.4*x.lt.e)
                } else {
                    x.gt.e=chngpt.var>e
                    beta.4=coef.hat["("%.%chngpt.var.name%.%"-chngpt)+"]                    
                    # x.tilda
                    x.tilda=cbind(Z, if (threshold.type=="segmented") chngpt.var, (chngpt.var-e)*x.gt.e, -beta.4*x.gt.e)
                }
                
                antilink=get(family)()$linkinv
                mu.est=drop(antilink(x.tilda[,-p] %*% coef.hat[-p]))        
                M=tXDX(x.tilda, resid(best.fit,"response")^2)/n # has to use the robust, sandwich estimation version. Note the threshold.type needs to be response for glm
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
                V.2[p,p-1]<-V.2[p-1,p]<- -mean((y-mu.est)*if(threshold.type=="upperhinge") x.lt.e else x.gt.e) 
                # compute V.2[p, p]
                newdata=data; newdata[[chngpt.var.name]]=e
                m.0=mean(predict(aux.fit, newdata, "response"))
                p.e.z.0=mean(antilink(cbind(Z, if (threshold.type=="segmented") e) %*% coef.hat[1:(p-2)]))                
                V.2[p, p]= ifelse(threshold.type=="upperhinge",-1,1)* beta.4 * f.x.e * (m.0-p.e.z.0)
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
                        cat("scale.chisq is negative\n")
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
                
                # TODO: we do a better job when it is grid-based search since we have already computed the profile likelihood 
#                if (est.method %in% c("grid","fastgrid","fastgrid2")) {
#                    c(NA,NA)# todo
#                } else if (est.method=="smoothapprox") {
                    
                    idx.chngpt=which(chngpt.var.sorted>=coef.hat["chngpt"])[1]
                    data = make.chngpt.var(chngpt.var, chngpt.var.sorted[idx.chngpt], threshold.type, data, b.transition)
                    fit=do.regression (formula.new, data, family)
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], threshold.type, data, b.transition)
                        fit = do.regression (formula.new, data, family)
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], threshold.type, data, b.transition)
                        fit = do.regression (formula.new, data, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            lik = NA; break
                        } else lik = as.numeric(logLik(fit$value))     
                    } 
                    #ub=if(is.na(lik)) NA else chngpt.var.sorted[idx-1]
                    ub=chngpt.var.sorted[idx-1]
                    
                    if (verbose==2) {
                        plot(chngpt.var.sorted, 2*profile.liks, xlab="change points", ylab="deviance relative to null model", main="Test Inversion CI")
                    }
                    
                    c(lb,ub)
#                } else stop("wrong est.method")
                
            } 
    
            # use F statistics to find test-inversion CI
            ci.test.inv.F=function(c.alpha){
                if (verbose) cat("in ci.test.inv.F\n")
                
                if (est.method %in% c("grid","fastgrid","fastgrid2")) {
                    c(NA,NA)# todo
                } else if (est.method=="smoothapprox") {
                    
                    idx.chngpt=which(chngpt.var.sorted>=coef.hat["chngpt"])[1]
                    data = make.chngpt.var(chngpt.var, chngpt.var.sorted[idx.chngpt], threshold.type, data, b.transition)# it is ok that data is assigned to b/c the function adds columns
                    fit=do.regression (formula.new, data, family)
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], threshold.type, data, b.transition)
                        fit = do.regression (formula.new, data, family)
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
                        data=make.chngpt.var(chngpt.var, chngpt.var.sorted[idx], threshold.type, data, b.transition)
                        fit = do.regression (formula.new, data, family)
                        if(length(fit$warning)!=0) {
                            #if(verbose) print(fit$warning)
                            sigsq = NA; break
                        } else sigsq = (sigma(fit$value))**2
                    } 
                    #ub=if(is.na(sigsq)) NA else chngpt.var.sorted[idx-1]
                    ub=chngpt.var.sorted[idx-1]
                    
                    c(lb,ub)
                    
                } else stop("wrong est.method")
                
            } 
                
            # the index may be hardcoded in this function
            ci.bootstrap=function(){
                if(verbose) cat("in ci.bootstrap\n")
                
                # bootstrapping
                # save rng state before set.seed in order to restore before exiting this function
                save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
                if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }      
                set.seed(1)     # this is only added on 7/30/2017, well after biometrics paper is published. should not change the results qualitatively
                
                #if (fastgrid.ok & est.method%in%c("fastgrid","fastgrid2")) {                    
                # for the upperhinge2 manuscript, we add grid here
                if (fastgrid.ok & est.method%in%c("fastgrid","fastgrid2","gridC")) {      
                    boot.out=list()
                    class(boot.out)="boot"
                    boot.out$t0=coef.hat
                    boot.out$R=ci.bootstrap.size
                    boot.out$call=c("boot")
                    boot.out$sim="ordinary"
                    boot.out$stype="i"
                    boot.out$fake=TRUE
                    
                    if(!stratified) {
                        f.name=est.method%.%"_search" # weirdly, f.name cannot be put inline b/c rcmdcheck throws an error otherwise
                        boot.out$t = .Call(f.name, 
                            cbind(Z.sorted, chngpt.var.sorted, if(threshold.type=="segmented") chngpt.var.sorted), #design matrix, the last column is to be used as (x-e)+ in the C program
                            as.double(y.sorted), 
                            as.double(w.sorted), 
                            w.all.one, 
                            as.integer(attr(chngpts,"index")), # potential thresholds
                            ci.bootstrap.size,# bootstrap size
                            threshold.type=="upperhinge"
                        )             
                    } else {
                        boot.out$t = .Call("twoD_search", 
                            cbind(Z.sorted, if(threshold.type=="segmented") chngpt.var.sorted), # X
                            as.double(chngpt.var.sorted), # threshold variable
                            as.double(y.sorted), 
                            as.double(w.sorted), 
                            w.all.one, 
                            as.integer(stratified.by.sorted),
        #                    as.integer(attr(chngpts[[1]],"index")), # potential thresholds 
        #                    as.integer(attr(chngpts[[2]],"index")), # potential thresholds
                            as.double(lb.quantile),
                            as.double(ub.quantile),
                            ci.bootstrap.size,# bootstrap size
                            threshold.type=="upperhinge"
                        )
                    }
                    boot.out$t=t(matrix(boot.out$t, ncol=ci.bootstrap.size))
                    
                } else {
                    boot.out=boot::boot(data.sorted, R=ci.bootstrap.size, sim = "ordinary", stype = "i", statistic=function(dat, ii){
                        if (m.out.of.n) ii=ii[1:(4*sqrt(n))] # m out of n bootstrap
                        # use chngpt.init so that the mle will be less likely to be inferior to the model conditional on e.hat, but it does not always work
                        # for studentized bootstrap interval, need variance estimates as well, however, stud performs pretty badly
                        fit.ii=try(chngptm (formula.1, formula.2, family, dat[ii,], type=threshold.type, 
                            formula.strat=formula.strat,
                            weights=w.sorted[ii], threshold.type, est.method=est.method, var.type="none", 
                            b.transition=b.transition, chngpt.init=coef.hat["chngpt"], verbose=verbose>=3, keep.best.fit=TRUE, lb.quantile=lb.quantile, ub.quantile=ub.quantile, 
                            grid.search.max=grid.search.max))
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
                                dat.tmp=make.chngpt.var(chngpt.var[ii], coef.hat["chngpt"], threshold.type, data[ii,], b.transition) # don't forget to do chngpt.var[ii]!
                                fit.ii.ehat=do.regression (formula.new, data=dat.tmp, family); if(length(fit.ii.ehat$warning)!=0 & verbose) print(fit.ii.ehat$warning)
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
                rownames(outs$symm)<-c((alpha/2*100)%.%"%",(100-alpha/2*100)%.%"%") 
                
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
#                        fit.ii=try(chngptm (formula.1, formula.2, family, dat, threshold.type, est.method="smoothapprox", var.type="none"))
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
            if (family=="binomial" & any(abs(coef.hat[-c(1,length(coef.hat))])>=search.bound)) {
                warn.var="Point estimate over search bound. Not computing variance estimate."
                cat(warn.var, "\n")
                var.est=NULL
            } else {
                if (var.type=="none") {
                    var.est=NULL
                } else if (var.type=="bootstrap") {
                    var.est=ci.bootstrap()
                } else {                
                    var.est=switch(var.type, 
#                        smooth=var.est.smooth(), 
                        model=var.est.model(test.inv.ci), 
                        robust=     var.est.robust(aux.fit, test.inv.ci), 
                        robusttruth=var.est.robust(aux.fit, test.inv.ci, true.param=TRUE), 
                        all=list(#"smooth"=var.est.smooth(), 
                                  "model"=var.est.model(test.inv.ci), 
                                 "robust"=var.est.robust(aux.fit, test.inv.ci), 
                               "sandwich"=var.est.model(test.inv.ci=FALSE,robust=TRUE) )
                    ) 
                }
            }
        } else stop("threshold.type incorrect")
        
    } # end if (has.itxn | !family %in% c("binomial","gaussian"))  else 
    
        
    res=list(
          coefficients=coef.hat
        , vcov=var.est
#        , formula.1=formula.1
#        , formula.2=formula.2
        , formula.new=formula.new
        , chngpt.var=chngpt.var.name
        , chngpt=e.final
        , est.method = est.method
        , b.transition=b.transition
        , threshold.type = threshold.type
        , glm.warn = glm.warn
        , family=family    
        , var.type=var.type
        , n=n # after removing rows with missing values
    )
    if (est.method=="smoothapprox") {
        res=c(res, list(
            converged=converged, 
            iter=n.iter))
    } else if (est.method %in% c("grid","fastgrid","fastgrid2","gridC")) {
        res=c(res, list(
            chngpts=chngpts, 
            logliks=logliks))
#        # good.soln needs to take as input a chngptm object that has chngpts etc
#        class(res)=c("chngptm", class(res)) 
#        tmp=good.soln(res, plot=verbose==3)
#        res=c(res, list(good.soln=tmp) )
    }
    if (!is.null(warn.var)) res[["warning"]]=warn.var
    if(keep.best.fit) res$best.fit=best.fit 
    if (!stratified) names(res$chngpt) = round(100*mean(chngpt.var<res$chngpt),1) %.% "%" else res$chngpt=unname(res$chngpt)
    
    class(res)=c("chngptm", class(res))
    res    
}
# for linear only
chngptm.xy = function(x, y, type=c("step","hinge","segmented","segmented2","stegmented"),  ...)  {
    tmp.dat=data.frame(x=x, y=y)
    chngptm(y~1, ~x, family="gaussian", tmp.dat, type=type, ...)
}


get.chngpts=function (chngpt.var.sorted, lb.quantile, ub.quantile, n.chngpts, stratified.by.sorted=NULL) {    
    if (is.null(stratified.by.sorted)) {
        n=length(chngpt.var.sorted)
        nLower=round(n*lb.quantile)+1; nUpper=round(n*ub.quantile)
        chngpts=chngpt.var.sorted[nLower:nUpper]
        attr(chngpts, "index")=nLower:nUpper
        
        # subset chngpts if needed
        if (length(chngpts)>n.chngpts) {
            ind=round(seq(1,length(chngpts),length=n.chngpts))
            chngpts=chngpts[ind]
            attr(chngpts, "index")=(nLower:nUpper)[ind]
        } 
        
    } else {
        # stratified
        chngpts = lapply (0:1, function(v) {
            n=sum(stratified.by.sorted==v)
            nLower=round(n*lb.quantile)+1; nUpper=round(n*ub.quantile)
            chngpts=chngpt.var.sorted[stratified.by.sorted==v][nLower:nUpper]
            attr(chngpts, "index")=nLower:nUpper + if (v==1) sum(stratified.by.sorted==0) else 0 # when v is 1, add a base so that the index correspond to the chngpt.var.sorted
            chngpts
        })
        
        if (length(chngpts)>n.chngpts) {
            stop("subsetting chngpts not supported for now for stratified")
        } 
    }
    
    chngpts
}
# an alternative and equivalent way to define nLower/nUpper: nLower=sum(chngpt.var.sorted<quantile(chngpt.var.sorted, lb.quantile))+1; nUpper=sum(chngpt.var.sorted<=quantile(chngpt.var.sorted, ub.quantile));myprint(nLower, nUpper)
#    # old from chngptm, take the mid points between data points as change point candidates. For step models, use any point in between does not change likelihood, but it does change likelihood for segmented model
#    chngpts=(chngpts[1:(length(chngpts)-1)]+chngpts[-1])/2
#    # old from chngpt.test, chngpts are mostly between observed values
#    chngpts=quantile(chngpt.var, seq(lb.quantile,ub.quantile,length=chngpts.cnt)) 


## make.chngpt.var is needed in both testing and estimation
## in estimation by grid search, b.transition is set to null when this function is called
## in estimation by smooth approx, b.transition value is saved to object
## New on May 13, 2017 this function is modified to fit smooth transition model
make.chngpt.var=function(x, e, threshold.type, data=NULL, b.transition=Inf, stratified.by=NULL) {    
    if (is.null(stratified.by)) {
        transition=expit(b.transition*(x-e))
        if(b.transition==Inf) transition[is.nan(transition)]=1 # if this is not here, there will be NaN and that messes up things
        
        if(threshold.type=="step") {
            out=transition
        } else if (threshold.type %in% c("hinge","segmented")) {
            out=transition*(x-e)
             # for segmented, x also becomes part of the null model
        } else if (threshold.type=="upperhinge") {
            out=(1-transition)*(x-e)
        } else if (threshold.type=="segmented2") {
            out=transition*x # x by itself also becomes part of the null model
        } else if (threshold.type=="stegmented") {
            out=cbind(transition, transition*(x-e))  # x also becomes part of the null model
        } else stop("wrong threshold.type")
        
        #str(out)    
        # unchanged on May 13, 2017
        if (is.null(data)) {
            out    
        } else {
            if (threshold.type=="stegmented") {
                data$x.mod.e = out[,1]
                data$x.mod.e.2 = out[,2]
            } else {
                data$x.mod.e = out
            }
            data
        }
        
    } else {
        if (length(e)!=2) stop("make.chngpt.var needs a vector of length two for the argument e when stratified.by is not null")
        if (length(stratified.by)!=length(x)) stop("make.chngpt.var: stratified.by needs to be a vector of 0/1 of the same length as x")
        # create thresholded var within each stratum
        # first col of out is for stratified.by==0, second for stratified.by==1
        out = cbind(1-stratified.by, stratified.by)
        out[stratified.by==0,1] = make.chngpt.var (x[stratified.by==0], e[1], threshold.type, b.transition=b.transition) 
        out[stratified.by==1,2] = make.chngpt.var (x[stratified.by==1], e[2], threshold.type, b.transition=b.transition)
        if (is.null(data)) {
            colnames(out)=c("x.mod.e","x.mod.f")
            out    
        } else {
            # assuming only segmented/hinge are supported
            data$x.mod.e = out[,1]
            data$x.mod.f = out[,2]
            data
        }
        
    }    
    
}
#
# The difference between the next two versions is that when b.transition is not infinity, if x<e, x is 0 in the older version but only close to 0 in the newer version
## before May 13, 2017
#make.chngpt.var=function(x, e, threshold.type, data=NULL, b.transition=NULL) {    
#    if(threshold.type=="step") {
#        out=ifelse(x>=e, 1, 0)
#    } else if (threshold.type=="hinge") {
#        out=ifelse(x>=e, x-e, 0)
#    } else if (threshold.type=="segmented") {
#        out=ifelse(x>=e, x-e, 0) # x also becomes part of the null model
#    } else if (threshold.type=="stegmented") {
#        out=cbind(ifelse(x>=e, 1, 0), ifelse(x>=e, x-e, 0))  # x also becomes part of the null model
#    }    
#    if (!is.null(b.transition)) out=out * 1/(1+exp(b.transition*(x-e))) 
#


# the added terms will be named as, say x.mod.e. these variables will be created by make.chngpt.var
get.f.alt=function(threshold.type, chngpt.var.name, modified.by=NULL, stratified=FALSE) {
    if (is.null(modified.by) & !stratified) {
        f.alt = if (threshold.type=="stegmented") "~.+x.mod.e+x.mod.e.2" else "~.+x.mod.e"
        
    } else if (!is.null(modified.by) & !stratified) {
        f.alt = if (threshold.type=="stegmented") {
            "~.+(x.mod.e+x.mod.e.2)*"%.%modified.by%.%"+"%.%chngpt.var.name%.%":"%.%modified.by 
        } else if (threshold.type %in% c("segmented","segmented2")) {
            "~."%.%"+"%.%chngpt.var.name%.%":"%.%modified.by%.%"+x.mod.e*"%.%modified.by
        } else if (threshold.type %in% c("step","hinge","upperhinge")) {
            "~.+x.mod.e*"%.%modified.by
        } else stop("wrong threshold.type within get.f.alt")
         
    } else if (stratified & is.null(modified.by)) {
        f.alt = if (threshold.type %in% c("segmented","hinge")) {
            "~.+x.mod.e+x.mod.f"
        } else stop("wrong threshold.type within get.f.alt")
        
    } else stop("modified.by and stratified.by cannot both be there")
    as.formula(f.alt)
}
        
#deviance.3pl <- expression( (1-y) * ( c + (d-c)/(1+exp(b*(t-e))) )  +  log( 1 + exp( -c - (d-c)/(1+exp(b*(t-e))) ) ) )
#deviance.3pl.deriv=deriv3(deviance.3pl, c("c","d","e"), c("c","d","e","t","y","b"))


predict.chngptm=function (object, newdata = NULL, type = c("link", "response", "terms"), ...){    
    if (is.null(object$best.fit)) stop("To make predictions, chngptm fit needs to have keep.best.fit=TRUE in the option.")
    newdata = make.chngpt.var(newdata[[object$chngpt.var]], object$chngpt, object$threshold.type, newdata, object$b.transition)            
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
residuals.chngptm=function(object, ...) {
    residuals(object$best.fit, ...)
}
vcov.chngptm=function(object, var.type=NULL, ...) {
#    if(object$threshold.type %in% c("hinge","segmented")) {
#        object$vcov[-length(object$coefficients),-length(object$coefficients)] # not return the chngpoint estimate
#    } else  {
#        object$vcov
#    }
    boot.conf=FALSE        
    if(!is.null(var.type)) {
        if (var.type %in% c("perc","basic","bc","bcabc","symm")) {
            boot.conf=TRUE
        }
        vcov=object$vcov[[var.type]]
    
    } else {            
        if(is.list(object$vcov)){
            if (!is.null(object$vcov$robust)) {
                vcov=object$vcov$robust
            } else if (!is.null(object$vcov$symm)){
                vcov=object$vcov$symm
                boot.conf=TRUE
                #if (all(is.na(vcov))) vcov=object$vcov$bc # if the bootstrap samples have too many NA, only bc quantiles can be estimated. However, bc quantiles may be misleading
            }
            
        } else {
            vcov=object$vcov
        }
    } 
    attr(vcov, "boot.conf")=boot.conf
    vcov
}
getFixedEf.chngptm=function(object, exp=FALSE, show.slope.post.threshold=FALSE, ...) {
    capture.output({
        res=summary(object, expo=exp, show.slope.post.threshold=show.slope.post.threshold, ...)
    })
    rbind(res$coefficients, chngpt=c(res$chngpt[1], NA, res$chngpt[2:3]))
}
lincomb=function(object, comb, alpha=0.05){
    if(!is.matrix(comb)) comb=matrix(comb, ncol=1)
    est=c(coef(object)%*%comb)
    
    if (is.matrix(object$vcov)) {
        # analytical
        if(length(comb)==ncol(object$vcov)-1) {
            # vcov has chngpt
            comb=c(comb, 0)
        }
        se=sqrt(comb%*%object$vcov%*%comb)
        c(est, est-1.96*se, est+1.96*se)
    } else {
        # bootstrap
        samples=object$vcov$boot.samples[,1:length(comb)] %*% comb        
        # symmetrical bootstrap CI
        q.1=quantile(abs(samples-est), 1-alpha, na.rm=TRUE)
        ci=est+c(-q.1, q.1)
        c(est, ci)    
    }
}

# which=1: scatterplot with fitted line, only works for simple regression
plot.chngptm=function(x, which=NULL, xlim=NULL, lwd=2, lcol="red", add=FALSE, add.points=TRUE, add.ci=TRUE, breaks=20, ...) {
    
    has.boot.samples=FALSE
    if(is.list(x$vcov)) if(!is.null(x$vcov$boot.samples)) has.boot.samples=TRUE 
    has.submodel.lik=FALSE
    if(!is.null(x$logliks)) has.submodel.lik=TRUE 
    
    if(is.null(which)) {
        par(mfrow=c(1+has.submodel.lik+has.boot.samples,1))
        plot(x,1,xlim=xlim,lwd=lwd,lcol=lcol,...)
        if (has.submodel.lik) plot(x,2,xlim=xlim,lwd=lwd,lcol=lcol,...)
        if (has.boot.samples) plot(x,3,xlim=xlim,lwd=lwd,lcol=lcol,...)
        return(invisible())
    }    
    
    fit=x
    linkinv=get(fit$family)()$linkinv
    if(is.null(xlim)) xlim=range(fit$best.fit$data[[fit$chngpt.var]])
    
    out=list()
    if(which==1) {
    # scatterplot with lines
        if(!add) plot(fit$best.fit$data[[fit$chngpt.var]], fit$best.fit$y, xlim=xlim, xlab=fit$chngpt.var, ylab=names(fit$best.fit$model)[1], type="n", ...)
        chngpt.est=fit$chngpt        
        intercept=coef(fit)[1]
        
        if(x$threshold.type=="upperhinge") {
            pre.slope=coef(fit)["("%.%fit$chngpt.var%.%"-chngpt)-"]
            myprint(chngpt.est)
            
            xx=seq(xlim[1],chngpt.est,length=100)
            yy=linkinv(intercept+pre.slope*(xx-chngpt.est))
            out[[1]]=cbind(xx,yy)
            lines(xx, yy, lwd=lwd, col=lcol)
            
            xx1=seq(chngpt.est, xlim[2], length=100)
            yy=rep(intercept, length(xx1))
            out[[1]]=rbind(out[[1]], cbind(xx1,yy))
            lines(xx1, yy, lwd=lwd, col=lcol)
            
        } else {
            if (!is.na(coef(fit)[fit$chngpt.var]))  pre.slope=coef(fit)[fit$chngpt.var] else pre.slope=0
            delta.slope=coef(fit)["("%.%fit$chngpt.var%.%"-chngpt)+"]            
            
            xx=seq(xlim[1],chngpt.est,length=100)
            yy=linkinv(intercept+pre.slope*xx)
            out[[1]]=cbind(xx,yy)
            lines(xx, yy, lwd=lwd, col=lcol)
            
            xx1=seq(chngpt.est, xlim[2], length=100)
            yy=linkinv(intercept+pre.slope*xx1+delta.slope*(xx1-chngpt.est))
            out[[1]]=rbind(out[[1]], cbind(xx,yy))
            lines(xx1, yy, lwd=lwd, col=lcol)
        }
        
        
        # add points again on top
        if(add.points) points(fit$best.fit$data[[fit$chngpt.var]], fit$best.fit$y)
        
    } else if (which==2) {
    # loglik vs threshold
        if (!has.submodel.lik) stop("no submodel likelihoods")
        plot(fit$chngpts, fit$logliks, xlim=xlim, xlab="threshold", ylab="log likelihood", type="b")
        abline(v=summary(fit)$chngpt[1], lty=c(1)) # point est 
        
    } else if (which==3) {
    # histograms of threshold estimates from bootstrap
        if (!has.boot.samples) stop("boot samples not saved")
        hist(fit$vcov$boot.samples[,1+length(coef(fit))], xlim=xlim, breaks=breaks, main="", xlab="bootstrap threshold estimate", add=add, ...)
        if(add.ci) abline(v=summary(fit)$chngpt[2:3], lty=c(2,2)) #CI
        
    } else stop("wrong which")
    
    invisible(out)
}

summary.chngptm=function(object, var.type=NULL, expo=FALSE, show.slope.post.threshold=FALSE, verbose=FALSE, ...) {    
    # "Estimate" "Std. Error" "t value" "Pr(>|t|)"        
    fit=object
    p.z=length(fit$coefficients) # this count includes threshold parameter
    n=fit$n
    threshold.type=fit$threshold.type
    
    if (is.null(fit$vcov)) {
        cat("No variance estimate available.\n\n")
        print(fit)
        return (invisible())
    } else {    
        vcov=vcov(object, var.type)
        boot.conf=attr(vcov, "boot.conf")
        if(is.null(vcov)) {
            cat("No variance estimate available.\n\n")
            print(fit)
            return (invisible())
        }             
    }    
    if(object$threshold.type %in% c("hinge","segmented","upperhinge") & !boot.conf) {
        vcov.t=vcov[p.z,p.z]
        tmp=vcov[-p.z,-p.z] # not return the chngpoint estimate
        # the last line lost the attr chngpt.ci, need to make a copy
        if(!is.null(attr(vcov,"chngpt.ci"))) attr(tmp,"chngpt.ci")=attr(vcov,"chngpt.ci")
        vcov=tmp
    } 
    if(verbose) str(vcov)
    
    # assuming the last of coefficients is always the change point
    # deal with coefficients and change point separately
    
    res=list()
    
    # coefficients
    transf=if(fit$family=="binomial" & expo) exp else identity
    if (boot.conf){
        lb=transf(vcov[1,])
        ub=transf(vcov[2,])
        # an approximate p value, prompted by comment about it being NA from Katie Chartrand, James Cook University
        sd.star=(vcov[2,]-vcov[1,])/1.96/2
        pval=pnorm(abs(fit$coefficients/sd.star), lower.tail=F)*2
    } else  {
        lb=transf(unname(fit$coefficients[1:(p.z-1)] - sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        ub=transf(unname(fit$coefficients[1:(p.z-1)] + sqrt(diag(vcov)) * qt(0.975, df=n-p.z, lower.tail=TRUE)))
        pval=unname(pt(abs(fit$coefficients[1:(p.z-1)] / sqrt(diag(vcov))), df=n-p.z, lower.tail=FALSE)) *2 # *2 is added on 8/2/2016
    }
    res$coefficients=mysapply(1:(p.z-1), function (i) {
        c(
              transf(unname(fit$coefficients[i]))
            , "p.value" = unname(pval[i])
            , "(lower" = unname(lb[i])
            , "upper)" = unname(ub[i])
#              "Estimate"=unname(fit$coefficients[i])
#            , "Std. Error" = sqrt(vcov[i,i])
#            , "t value" = unname(fit$coefficients[i] / sqrt(vcov[i,i]))
#            , "Pr(>|t|)" = unname(pt(abs(fit$coefficients[i] / sqrt(vcov[i,i])), df=n-p.z, lower.tail=FALSE))
        )
    })
    rownames(res$coefficients)=names(fit$coefficients)[-p.z]
    colnames(res$coefficients)[1]=if(fit$family=="binomial" & expo) "OR" else "est"
    if(boot.conf) colnames(res$coefficients)[2]=colnames(res$coefficients)[2]%.%"*"
    
    # may need to find linear combination
    if(show.slope.post.threshold){
        comb=rep(0,length=p.z-1)
        comb[p.z-2:1]=1
        tmp=lincomb(object, comb, alpha=0.05)
        res$coefficients[p.z-1,]=c(tmp[1],NA,tmp[2],tmp[3])
        rownames(res$coefficients)[p.z-1]=sub("\\+",">",rownames(res$coefficients)[p.z-1])
    }
    
    # change point
    i=p.z
    if (threshold.type %in% c("hinge","segmented","upperhinge")) {
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
    } else if (threshold.type %in% c("step","stegmented")) {
        res$chngpt=c(fit$chngpt)
    } else stop("incorrect threshold.type")
    
    res$threshold.type=fit$threshold.type
    
    class(res) <- "summary.chngptm"
    res
}
print.summary.chngptm=function(x,...) {
    cat("Change point model threshold.type: ",x$threshold.type,"\n\n")
    cat("Coefficients:\n")
    print(x$coefficients)
    cat("\nThreshold:\n")
    names(x$chngpt)[1]="est" # Thanks to Katie Chartrand, James Cook University
    print(x$chngpt)    
}


# the following functions are used as objective functions in smoothapprox, b/c we don't supply derivative to optim anymore. 
.lik.f=function(linear, y, family, weights) {
    tmp=if(family=="binomial") {
        (1-y)*linear + log(1+exp(-linear)) 
    } else if(family=="gaussian") {
        (linear-y)**2
    }
    sum(tmp*weights)
}
# step change point model
dev.step.f <- function(theta,x,y,b,alpha.z,z.1,family,weights)  {
    beta=theta[1]; e=theta[2]
    eta=beta/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
}
dev.step.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family,weights)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(beta1+beta2*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
}
# hinge change point model
dev.hinge.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta=theta[1]; e=theta[2]
    eta=(x-e)*beta/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
# upperhinge model
dev.upperhinge.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta=theta[1]; e=theta[2]
    eta=(x-e)*beta*(1-1/(1+exp(-b*(x-e))))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
dev.hinge.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family,weights)  {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=(x-e)*(beta1+beta2*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
}
# segmented change point model
dev.segmented.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + (x-e)*beta2/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
dev.segmented2.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[2]; e=theta[3]
    eta=beta1*x + x*beta2/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
dev.segmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*(x-e)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
dev.segmented2.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + (beta3+beta4*z.1)*x/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
# stegmented change point model
dev.stegmented.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[2]; beta3=theta[3]; e=theta[4]
    eta=beta1*x + ((x-e)*beta3 + beta2)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
} 
# the order of parameters in the following needs to be fixed
dev.stegmented.itxn.f <- function(theta,x,y,b,alpha.z,z.1,family,weights) {
    beta1=theta[1]; beta2=theta[3]; beta3=theta[2]; beta4=theta[4]; beta5=theta[5]; beta6=theta[6]; e=theta[5]# note the order change between beta and theta
    eta=(beta1+beta2*z.1)*x + ((beta3+beta4*z.1)*(x-e) + beta4+beta5*z.1)/(1+exp(-b*(x-e)))
    linear=alpha.z+eta
    .lik.f(linear,y,family,weights)
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
# upper hinge model
dev.upperhinge <- expression( (1-y) * ( alpha.z + (x-e)*beta*(1-1/(1+exp(-b*(x-e)))) )  +  log( 1 + exp( -alpha.z - (x-e)*beta * (1-1/(1+exp(-b*(x-e)))) ) ) )


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


# Adam Elder
# assum input are sorted
# was useful for debugging C implementation
get.liks.upperhinge=function(y,Z,x,cand_e){    
    n=length(y)
    
    ## Solving for the quantities that only need to be found once
    A <- solve(t(Z) %*% Z)
    H <- Z %*% A %*% t(Z)
    U <- as.numeric(t(y) %*% (H - diag(n)))
    M <- diag(n) - H
    
    ## change to pass this info in
#    ## Selecting the possible cutoff values
#    cand_e <- x[floor(n * 0.05):ceiling(n*0.95)]
    sudo_liks = rep(NA,length(cand_e))
    
    ## Creating a variable that stores the best cuttoff point, as
    ## a proxy for the maximum likelihood conditional on this cuttoff
    best_e <- c(cand_e[1], 0)
    for(i in 1:length(cand_e)){
        e=cand_e[i]
        nu_e <- pmax(e - x, 0) #the last row vector of our design matrix for a given cuttoff value
        #if(sum(nu_e) != 0){
        sudo_liks[i] <- as.numeric((t(U) %*% nu_e))**2/as.numeric((t(nu_e) %*% M %*% nu_e))
      #}
    }
    
    sudo_liks
}

## R code for validating 2D search
#X=cbind(1, dat$z, dat$x)
#A <- solve(t(X) %*% X)
#H <- X %*% A %*% t(X)
#dat$y %*% H %*% dat$y
#Y=dat$y
#r=Y-H%*%Y
#
#chngpt.var=dat$x
#stratified.by=dat$z>0
#order.=order(stratified.by, chngpt.var)
#chngpt.var.sorted=chngpt.var[order.]
#r=r[order.]
#v=cbind(chngpt.var.sorted,chngpt.var.sorted)
#v[12:20,1]=0
#v[1:11,1]=v[1:11,1]-v[2,1]
#v[1:2,1]=0
#v[1:11,2]=0
#v[12:20,2]=v[12:20,2]-v[13,2]
#v[11+1,2]=0
#v
#
#X1=X[order.,]
#
#t(v)%*%r
#t(v)%*%v
#t(t(v)%*%r) %*% solve(t(v)%*%v-t(v)%*%X1 %*% solve(t(X1) %*% X1) %*% t(X1)%*%v) %*% t(v)%*%r + dat$y %*% H %*% dat$y
#
#
#Xe=cbind(X[order.,],v)
#Y[order.] %*% Xe %*% solve(t(Xe) %*% Xe) %*% t(Xe) %*% Y[order.]
