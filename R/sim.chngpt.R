expit.2pl=function(x,e,b) sapply(x, function(x) 1/(1+exp(-b*(x-e))))
sim.chngpt = function (
    mean.model=c("thresholded","thresholdedItxn","quadratic","quadratic2b","cubic2b","exp","flatHyperbolic","z2"), 
    threshold.type=c("NA","step","hinge","segmented","segmented2","stegmented"),# segmented2 differs from segmented in parameterization, it is the model studied in Cheng 2008
    b.transition=Inf,
    family=c("binomial","gaussian"), 
    x.distr=c("norm","norm3","norm6","imb","lin","mix","gam","zbinary","gam1","gam2", "fixnorm", "fixnorm3", "fixnorm6"), # gam1 is a hack to allow e. be different
    e.=NULL, mu.x=4.7, sd.x=NULL, sd=0.3, 
    alpha=NULL, alpha.candidate=NULL, coef.z=log(1.4), beta=NULL, beta.itxn=NULL, 
    n, seed, 
    weighted=FALSE, # sampling weights
    heteroscedastic=FALSE,
    verbose=FALSE) 
{
    
    if (!requireNamespace("mvtnorm")) {print("mvtnorm does not load successfully"); return (NULL) }
    if (!is.numeric(n)) stop("n is not numeric")
    
    if (missing(threshold.type) & startsWith(mean.model,"thresholded")) stop("threshold.type mssing")
    if (missing(family)) stop("family mssing")
    threshold.type<-match.arg(threshold.type)    
    mean.model<-match.arg(mean.model)    
    family<-match.arg(family)    
    x.distr<-match.arg(x.distr)    
    
    if(is.null(sd.x)) sd.x=if (mean.model=="quadratic") sd.x=1.4 else 1.6
    
    if(startsWith(x.distr,"fix")) {
        set.seed(1) # "seed" will be set before simulation of y
    } else {
        set.seed(seed)
    }
    
    #######################################################################################
    # generate covariates
    
    if(x.distr=="imb") { # imbalance
        x=c(rnorm(n-round(n/3), mu.x, sd.x), mu.x-abs(rnorm(round(n/3), 0, sd.x)))
        z=rep(1,n)
    } else if(x.distr=="lin") { # unif
        x=runif(n)*4*sd.x + mu.x-2*sd.x
        z=rnorm(n, mean=0, 1)
    } else if(x.distr=="mix") { # mixture
        x=c(rnorm(n*.6, mu.x, sd.x), rep(mu.x-2*sd.x, n*.4))
        z=rep(1,n)
    } else if(x.distr %in% c("gam","gam1","gam2")) { # gamma
        x=1.4*scale(rgamma(n=n, 2.5, 1))+mu.x/2
        z=rnorm(n, mean=0, 1)
        if(x.distr=="gam") e.=2.2 else if(x.distr=="gam1") e.=1.5 else if(x.distr=="gam2") e.=1 
        # for thresholded, override input
        if (mean.model=="thresholded") {
            if(x.distr=="gam") {
                alpha= if (threshold.type=="hinge") -0.5 else if(threshold.type=="segmented") -1.3 else stop("wrong threshold.type") # to have similar number of cases
            } else if (x.distr=="gam1") {
                alpha= if (threshold.type=="hinge") -0.2 else if(threshold.type=="segmented") -1 else stop("wrong threshold.type") # to have similar number of cases
            } else if (x.distr=="gam2") {
                alpha= if (threshold.type=="hinge") 0.2 else if(threshold.type=="segmented") -0.6 else stop("wrong threshold.type") # to have similar number of cases
            }            
        }            
        x.distr="gam" # the number trailing gam is only used to change e.
        
    } else if(startsWith(x.distr,"norm") | startsWith(x.distr,"fixnorm")) { # normal, fixnorm means design matrix is not random
        if (x.distr=="norm" | x.distr=="fixnorm") {
             rho=0
        } else if (x.distr=="norm3" | x.distr=="fixnorm3") {
            rho=0.3 
        } else if (x.distr=="norm6" | x.distr=="fixnorm6") {
            rho=0.6
        } else {
            stop("x.distr not supported: "%+%x.distr)
        }        
        tmp=mvtnorm::rmvnorm(n, mean = c(mu.x,0), sigma = matrix(c(sd.x^2,sd.x*rho,sd.x*rho,1),2)) # use mvtnorm
        x=tmp[,1]
        z=tmp[,2]    
    } else if(startsWith(x.distr,"zbinary")) { 
        x=rnorm(n, mu.x, sd.x)
        z=rbern(n, 1/2)-0.5
    } else stop("x.distr not supported: "%+%x.distr)    
    
    if (is.null(e.) | !startsWith(mean.model,"thresholded")) e.=4.7 # hard code e. for mean models other than thresholded
    if (verbose) {print(e.); print(mean(x<e.))}
    
    x.star = expit.2pl(x, e=e., b=b.transition)  
    
        
    #######################################################################################
    # make design matrix and coefficients
    
    if (startsWith(mean.model,"thresholded")) {        
        
        # when mean.model is not found, alpha remains a null, and do not throw an error
        # but if mean.model is NULL, an error is thrown
        if(!is.null(alpha.candidate)) alpha=alpha.candidate # used to determine sim.alphas
        #cat(e., beta, "\n")
        if(is.null(alpha)) alpha=try(chngpt::sim.alphas[[mean.model%+%"_"%+%sub("fix","",x.distr)]][e.%+%"", ifelse(mean.model=="thresholdedItxn",beta.itxn,beta)%+%""], silent=TRUE)
        if(is.null(alpha) | inherits(alpha, "try-error")) stop("alpha not found, please check beta or provide a null") 
        
        X=cbind(1,     z,        x,   x.star,   if(threshold.type=="segmented2") x.star*x else x.star*(x-e.),     z*x,   z*x.star,   z*x.star*(x-e.))
        coef.=c(alpha, z=coef.z, x=0, x.star=0, x.hinge=0,                                                        z.x=0, z.x.star=0, z.x.hinge=0)
        
        if (mean.model=="thresholded") { 
            if (threshold.type=="step") {
                coef.[1:5]=c(alpha, coef.z,          0,    beta,     0) 
            } else if (threshold.type=="hinge") {
                coef.[1:5]=c(alpha, coef.z,          0,       0,  beta) 
            } else if (threshold.type=="segmented") {
                coef.[1:5]=c(alpha, coef.z,  -log(.67),       0,  beta) 
            } else if (threshold.type=="segmented2") {
                coef.[1:5]=c(alpha, coef.z,  -log(.67),       0,  beta) 
            } else if (threshold.type=="stegmented") {
                coef.[1:5]=c(2,     coef.z,   log(.67), log(.67), beta) # all effects of x in the same direction, subject to perfect separation, though that does not seem to be the main problem
                #coef.[1:5]=c(0, coef.z, -log(.67), log(.67), beta) # effects of x and x.star in different direction
            }            
            
        } else if (mean.model == "thresholdedItxn") { 
        # intercept + main effect + interaction
        # used to be "sigmoid3","sigmoid4","sigmoid5"
    #        beta.var.name=switch(threshold.type,step="x.star",hinge="x.hinge",segmented="x.hinge",stegmented="x.hinge")
    #        coef.[beta.var.name]=switch(mean.model,sigmoid3=log(.67),sigmoid4=-log(.67),sigmoid5=0)
    #        coef.["z."%+%beta.var.name]=beta
    #        if (threshold.type=="segmented") {coef.["x"]=tmp; coef.["z.x"]=log(.67) } 
            beta.var.name=switch(threshold.type,step="x.star",hinge="x.hinge",segmented="x.hinge",stegmented="x.hinge")
            coef.[beta.var.name]=beta
            coef.["z."%+%beta.var.name]=beta.itxn
            if (threshold.type=="segmented") {coef.["x"]=tmp; coef.["z.x"]=log(.67) }    
        }
    #    else if (mean.model=="sigmoid1") { 
    #    # intercept only
    #        coef.["x.star"]=beta
    #        coef.["z"]=0
    #    } else if (mean.model=="sigmoid6") { 
    #    # special treatment, model misspecification
    #        coef.=c(alpha, coef.z, log(.67),  beta)
    #        X=cbind(1, z, x.star, x.star*z^3)    
        
        
    } else if (mean.model=="quadratic") { 
    # x+x^2 
        X=cbind(1,     z,        x,   x*x)
        coef.=c(alpha=-1, z=coef.z, x=-1 , x.quad=0.3)
    
    } else if (mean.model=="quadratic2b") { 
        X=cbind(1,     z,        x,   x*x)
        coef.=c(alpha=-1, z=coef.z, x=-2*mu.x , x.quad=1)
    
    } else if (mean.model=="cubic2b") { 
    # x+x^2+x^3
        X=cbind(1,     z,        x,   x*x,   x*x*x)
        coef.=c(alpha=-1, z=coef.z, x=-1 , x.quad=beta, x.cube=1)
    
    } else if (mean.model=="exp") { 
        if(x.distr=="norm") {
            X=cbind(1,     z,        exp((x-5)/2.5))
            coef.=c(alpha=-5, z=coef.z, expx=4)
        } else if(x.distr=="gam") {
            X=cbind(1,     z,        exp((x-5)/2.5))
            coef.=c(alpha=-3, z=coef.z, expx=4)
        } else stop("wrong x.distr")
    
    } else if (mean.model=="flatHyperbolic") { 
    # beta*(x-e)+beta*sqrt((x-e)^2+g^2)
        if(x.distr=="norm") {
            g=1
            X=cbind(1,     z,        (x-e.)+sqrt((x-e.)^2+g^2) )
            coef.=c(alpha=-5, z=coef.z, 2)
        } else if(x.distr=="gam") {
            g=1
            X=cbind(1,     z,        (x-4)+sqrt((x-4)^2+g^2) )
            coef.=c(alpha=-2, z=coef.z, 2)
        }
    
    } else if (mean.model=="z2") { 
    # z^2
        X=cbind(1,     z,        z*z)
        coef.=c(alpha=alpha, z=coef.z, z.quad=0.3)
    
    } else stop("mean.model not supported: "%+%mean.model)     
    if (verbose) myprint(coef., digits=10)
    
    linear.predictors=drop(X %*% coef.)
    # simulate y
    if(startsWith(x.distr,"fix")) set.seed(seed)
    y=if(family=="binomial") {
        rbern(n, expit(linear.predictors)) 
    } else if(family=="gaussian") {
        if (!heteroscedastic) {
            rnorm(n, linear.predictors, sd)
        } else {
            rnorm(n, linear.predictors, sd*abs(linear.predictors))
        }
    }
    
    dat=data.frame (
        y=y,
        z=z,         
        
        x=x,
        x.sq=x*x,
        x.star=x.star,
        x.hinge=x.star*(x-e.),
        x.bin.med=ifelse(x>median(x), 1, 0),
        x.tri = factor(ifelse(x>quantile(x,2/3),"High",ifelse(x>quantile(x,1/3),"Medium","Low")), levels=c("Low","Medium","High")),
        
        x.tr.1=ifelse(x>log(100), x, 0) ,
        x.tr.2=ifelse(x>log(100), x, log(100)) ,
        x.tr.3=ifelse(x>3.5, x, 0) ,
        x.tr.4=ifelse(x>3.5, x, 3.5), 
        x.tr.5=ifelse(x>3.5, x, 3.5/2), 
        x.tr.6=ifelse(x>5, x, 5) ,
        x.ind =ifelse(x>3.5, 0, 1), 
        
        x.bin.35=ifelse(x>3.5, 1, 0), 
        x.bin.6=ifelse(x>6, 1, 0) ,
        x.bin.log100=ifelse(x>log(100), 1, 0),
        
        eta=linear.predictors
    )        
    
    # sampling
    if(weighted) {
        if (family=="gaussian") {
            # create strata
            dat$strata=ifelse(dat$y>4, 1, 2)
            #table(dat$strata)
            dat$sampling.p=ifelse(dat$strata %in% c(1), 1, 0.5)
            dat=dat[rbern(n, dat$sampling.p)==1,]
            
        } else if (family=="binomial") {
            # create strata
            dat$strata=2-dat$y # 1: cases; 2: controls
            dat$strata[dat$strata==2 & dat$z>0]=3
            #table(dat$strata)
            dat$sampling.p[dat$strata==1]=1
            dat$sampling.p[dat$strata==2]=1
            dat$sampling.p[dat$strata==3]=0.5
            dat=dat[rbern(n, dat$sampling.p)==1,]
            
        } else stop("weighted not supported here yet")
    } else {
        dat$sampling.p=rep(1, nrow(dat))
    }
    
    dat
}
