test.chngptm.linear <- function() {


library("chngpt")
library("RUnit")
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(tolerance)
verbose = FALSE


est.methods=c("grid","fastgrid2")
for (type in c("segmented","hinge")) {
#type="hinge"    
    print(type)
    dat=sim.chngpt("quadratic", n=20, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")
    out=sapply(est.methods, function(est.method){
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type="segmented", formula.strat=~I(z>0), est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        )
    })    
    out
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# thinned thresholds, fastgrid will not be suppported
est.methods=c("grid","gridC","fastgrid2")# fastgrid not implemented weights yet
for (type in c("segmented","hinge","upperhinge")) {
#type="hinge"    
    print(type)
    dat=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")      
    out=sapply(est.methods, function(est.method){
        fit.0 = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=10, verbose=verbose)
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks[1:3]),
          fit.0$coefficients, 
          fit.0$vcov$boot.samples[1,]
        )
    })    
    out
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# vanilla
est.methods=c("grid","gridC","fastgrid","fastgrid2")
for (type in c("segmented","hinge","upperhinge")) {
    print(type)
    dat=sim.chngpt (mean.model="thresholded", threshold.type=type, family="gaussian", sd=0.3, mu.z=0, alpha=0, coef.z=log(1.4), beta=-1, n=100, seed=1)     
    out=sapply(est.methods, function(est.method){
        # ci.bootstrap.size is set to 1 because grid uses boot pkg, which generates different bootsamples if bootstrap sample size is >1
        fit.0=chngptm (y~z, ~x, "gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=0) 
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks[1:3]),
          fit.0$coefficients, 
          fit.0$vcov$boot.samples[1,]
        ) 
    })    
    out
    checkEqualsNumeric(out[,"grid"], out[,"gridC"], tolerance=tolerance)    
    checkEqualsNumeric(out[,"grid"], out[,"fastgrid"], tolerance=tolerance)    
    checkEqualsNumeric(out[,"grid"], out[,"fastgrid2"], tolerance=tolerance)    
}



# weights. fastgrid2 not supported yet
est.methods=c("grid","gridC","fastgrid")
for (type in c("segmented","hinge","upperhinge")) {
    print(type)
    dat=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")      
    out=sapply(est.methods, function(est.method){
        fit.0 = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, weights=rep(1:2,each=30), verbose=verbose)
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks[1:3]),
          fit.0$coefficients, 
          fit.0$vcov$boot.samples[1,]
        )
    })    
    out
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# weights + thinned thresholds
# tbd


#make.chngpt.var (x=1:10, e=c(2,8), threshold.type="segmented", data=NULL, b.transition=Inf, stratified.by=rep(0:1,5)) 



}
