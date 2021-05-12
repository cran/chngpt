test.chngptm.linear <- function() {


library("RUnit")
library("chngpt")
  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(paste0("tol: ", tolerance), quote=FALSE)
verbose = 0


# autoregressive errors
dat=sim.chngpt(mean.model="thresholded", threshold.type="M20", n=100, seed=1, mu.x=5, beta=c(10,1),x.distr="lin", e.=5, family="gaussian", alpha=0, sd=3, coef.z=log(1.4), heteroscedastic=FALSE, ar=.5)
fit=  chngptm(y~z, ~x, type="M20", data=dat, family="gaussian", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, bootstrap.type="sieve")
checkEqualsNumeric(fit$vcov$boot.samples[1,], c(1.882445, 0.09883493,     4.70625,    -0.3412223, 5.937374), tolerance=tolerance)    
fit=  chngptm(y~z, ~x, type="M20", data=dat, family="gaussian", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, bootstrap.type="wild")
checkEqualsNumeric(fit$vcov$boot.samples[1,], c(0.7984787,-0.3663074,7.2876752,0.4433051,5.6787879), tolerance=tolerance)    
fit=  chngptm(y~z, ~x, type="M20", data=dat, family="gaussian", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, bootstrap.type="wildsieve")
checkEqualsNumeric(fit$vcov$boot.samples[1,], c(1.0655015,1.6554599,8.0428826,0.5240296,5.5494949), tolerance=tolerance)    
fit=  chngptm(y~z, ~x, type="M20", data=dat, family="gaussian", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, bootstrap.type="awb")
checkEqualsNumeric(fit$vcov$boot.samples[1,], c(1.3466068,1.5963537,7.9991371,0.5712195,5.6141414), tolerance=tolerance)    


print("########  segmented")
type="segmented"    
dat = sim.chngpt(mean.model="thresholded", threshold.type=type, n=200, seed=1, beta=c(2,2),       x.distr="lin", e.=5, family="gaussian", alpha=0, sd=3, coef.z=1)
attr(dat,"coef")
#plot(y~x, dat)
fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method="fastgrid2", var.type="bootstrap", ci.bootstrap.size=2, verbose=verbose)
if(verbose) plot(fit); fit
checkEqualsNumeric(coef(fit), c(-0.6677877,0.8992538,2.2789307,1.9569163), tolerance=tolerance)    
est=lincomb(fit, comb=c(0,0,1,1), alpha=0.05); print(est)


## test mclapply support for fastgrid, cannot do it on windows
#dat = sim.threephase(n = 20, seed=10)
#fit.0=chngptm(y~z, ~x, type="M111", data=dat, family="gaussian", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=3, verbose=1, ncpus=3); fit.0


# M111
est.methods=c("fastgrid","grid")
dat = sim.threephase(n = 20, seed=10)
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(y~z, ~x, type="M111", data=dat, family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose); fit.0
    if(verbose) plot(fit.0)
    out=cbind(out, c(
      fit.0$coefficients
      , fit.0$vcov$boot.samples[1,]
    ))
}
colnames(out)=est.methods
if(verbose) print(out)
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    


print("########  random intercept")
dat=sim.twophase.ran.inte(threshold.type="segmented", n=50, seed=1)
fit = chngptm (formula.1=y~z+(1|id), formula.2=~x, family="gaussian", dat, type="segmented", est.method="grid", var.type="bootstrap", ci.bootstrap.size=1)
checkEqualsNumeric(fit$coefficients, c(2.7154145,0.3514853,1.7894006,2.5695986,5.1571429), tolerance=tolerance)    
checkEqualsNumeric(fit$vcov$boot.samples[1,], c(3.799664,0.3901292,1.5940243,2.2716004,5.1571429), tolerance=tolerance)    
plot(fit)

# also works for M20
#fit = chngptm (formula.1=y~z+(1|id), formula.2=~x, family="gaussian", dat, type="M20", est.method="grid", var.type="bootstrap", ci.bootstrap.size=1)


print("########  M12c")
type="M12c"
est.methods=c("fastgrid2","grid")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, 
        var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=c(rep(1,9),rep(5,10)))
    if(verbose) plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    ))
}
colnames(out)=est.methods
if(verbose) print(out)
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    


print("########  M21c")
type="M21c"
dat = sim.chngpt(mean.model="thresholded", threshold.type="M21", n=250, seed=1, beta=if(type=="M21c") c(1,2,1) else c(1,1,2),       x.distr="lin", e.=5, family="gaussian", alpha=0, sd=3, coef.z=1)
est.methods=c("fastgrid2","grid")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose)
    if(verbose) plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    ))
}
colnames(out)=est.methods
if(verbose) print(out)
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    


print("########  offset ")
for (type in c("hinge","M02")) {
    est.methods=c("fastgrid2","grid")
    out=sapply(est.methods, function(est.method){
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="none", offset=rep(1,nrow(pressure)), ci.bootstrap.size=1, verbose=verbose, weights=c(rep(1,9), rep(5,10)))
        if(verbose) plot(fit.0); fit.0
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        )
    })    
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  step model m out of n bootstrap")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=20, seed=1, beta=-log(.67), alpha=1)
est.method="fastgrid2"
fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", m.out.of.n=10, ci.bootstrap.size=10, 
    lb.quantile=.1, ub.quantile=.9, verbose=verbose)
if (verbose) plot(fit.0); fit.0
out=c(
  fit.0$logliks[1],
  diff(fit.0$logliks)[1:3],
  fit.0$coefficients,
  fit.0$vcov$boot.samples[1,]
)    
checkEqualsNumeric(out, c(35.52723714,-0.15529784,-0.18206414,0.02324697,0.99337520,0.38349876,0.48532422,4.67409558,1.3862616,0.1185507,0.5859363,6.4998895), tolerance=tolerance)


print("########  step model subsampling bootstrap")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=20, seed=1, beta=-log(.67), alpha=1)
est.method="fastgrid2"
fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", subsampling=10, ci.bootstrap.size=10, verbose=verbose, 
    lb.quantile=.1, ub.quantile=.9)
if (verbose) plot(fit.0); fit.0
out=c(
  fit.0$logliks[1],
  diff(fit.0$logliks)[1:3],
  fit.0$coefficients,
  fit.0$vcov$boot.samples[1,]
)    
checkEqualsNumeric(out, c(35.52723714,-0.15529784,-0.18206414,0.02324697,0.99337520,0.38349876,0.48532422,4.67409558,1.03713972,0.27900925,0.65087769,  5.32027458), tolerance=tolerance)


print("########  step")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=250, seed=1, beta=-log(.67), alpha=1)
est.methods=c("fastgrid2","grid")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, weights=rep(c(1,5),each=nrow(dat)/2), verbose=verbose)
    if (verbose) plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    ))
}    
colnames(out)=est.methods
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)


print("########  segmented hinge upperhinge") 
est.methods=c("fastgrid2","grid"); names(est.methods)=est.methods
for (type in c("segmented","hinge","upperhinge")) {
# type="segmented"
    if (verbose) print(type)
    dat=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")      
    fits=lapply(est.methods, function(est.method){
        fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, weights=rep(c(1,5),each=nrow(dat)/2), verbose=verbose)
    })
    out=sapply(fits, function(fit){
        c(
          fit$logliks[1],
          diff(fit$logliks[1:3]),
          fit$coefficients 
          ,fit$vcov$boot.samples[1,]
        )
    })    
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  M02 and M20")
for (type in c("M20","M02")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="M20") c(32,2) else c(10,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0)
    if (verbose) plot(y~x, dat)
    est.methods=c("fastgrid2","grid")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="none", save.boot=T, ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
        if(verbose) plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  M12")
for (type in c("M12")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=c(10,10,20), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0)
    if (verbose) plot(y~x, dat)
    est.methods=c("fastgrid2","grid")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="none", save.boot=T, ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
        if(verbose) plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  M21 and M12")
for (type in c("M21","M12")) {
    est.methods=c("fastgrid2","grid")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, 
            lb.quantile=.1, ub.quantile=.9, verbose=verbose, weights=c(rep(1,9),rep(5,10)))
        if(verbose) plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  M22")
dat=sim.chngpt(mean.model="thresholded", threshold.type="M22", n=100, seed=1, beta=c(32,2,10,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("fastgrid2","grid")
out=sapply(est.methods, function(est.method){
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
    if(verbose) plot(fit.0); fit.0
    c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    )
})    
if(verbose) print(out)
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    



print("########  M22c")
dat=sim.chngpt(mean.model="thresholded", threshold.type="M22c", n=100, seed=1, beta=c(32,2,32,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("fastgrid2","grid")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22c", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
    if(verbose) plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    ))
}    
colnames(out)=est.methods
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)




print("########  M30 and M03")
for (type in c("M30","M03")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="M30") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.methods=c("fastgrid2","grid")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
        if(verbose) plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


print("########  M33c")
dat=sim.chngpt(mean.model="thresholded", threshold.type="M33c", n=100, seed=1, beta=c(32,2,-5,5), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("fastgrid2","grid")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M33c", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
    if(verbose) plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1]
      ,diff(fit.0$logliks)[1:3]
      ,fit.0$coefficients
      ,fit.0$vcov$boot.samples[1,]
    ))
} 
colnames(out)=est.methods
if(verbose) print(out)
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)


print("########  M31 and M13")
for (type in c("M31","M13")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="M31") c(32,2,-5,5) else c(-10,-.25,1,-2), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.methods=c("fastgrid2","grid")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
        if(verbose) plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks)[1:3],
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    if(verbose) print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# grid search only
print("########  M04, M40")
for(t in c("M04", "M40")){
    dat=sim.chngpt(mean.model="thresholded", threshold.type=ifelse(t=="M04","M03","M30"), n=250, seed=1, beta=if(t=="M40") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.method="grid"
    fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=t, est.method=est.method, var.type="none", save.boot=T, ci.bootstrap.size=1, verbose=verbose)
    if(verbose) plot(fit.0)
    if (t=="upperhinge") checkEqualsNumeric(fit.0$coefficients, c(0.03832806,0.36759453,31.75427906,1.60274645,-5.15254933,-0.02319485,5.00493816), tolerance=tolerance)
    if (t=="hinge") checkEqualsNumeric(fit.0$coefficients, c(-0.03536995,0.37418059,-9.81715060,-0.29244682,0.89944391,0.03742214,5.00493816), tolerance=tolerance)
}




print("########  thinned thresholds, only for grid")
est.method="grid"
type="segmented"    
dat=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")    
fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=10, verbose=verbose, 
    lb.quantile=.1, ub.quantile=.9)
out=c(
      fit$logliks[1],
      diff(fit$logliks[1:3]),
      fit$coefficients, 
      fit$vcov$boot.samples[1,]
    )  
if(verbose) print(out)
checkEqualsNumeric(out, c(544.7867707,1.3788716,2.1591613,-4.6760790,0.3340772,1.1283779,1.5775026,4.8043910,-4.6452678,0.3701558,1.1152278,1.7350673,4.9145547), tolerance=tolerance)    


# stratified
# segmented is not working because it is not clear how the model should be 
for (type in c("hinge","upperhinge")) {
#type="hinge"; est.method="grid"
    print(type)
    dat=sim.chngpt("quadratic", n=20, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")
    est.methods=c("fastgrid2","grid")    
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, formula.strat=~I(z>0), est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=verbose, 
    lb.quantile=.1, ub.quantile=.9)
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    #print(out)
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}



}# end of test function
