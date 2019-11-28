test.chngptm.linear <- function() {


library("RUnit")
library("chngpt")
  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(paste0("tol: ", tolerance), quote=FALSE)
verbose = FALSE


print("########  step")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=250, seed=1, beta=-log(.67), alpha=1)
est.methods=c("grid","fastgrid2")
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


print("########  step model m out of n subsampling")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=20, seed=1, beta=-log(.67), alpha=1)
est.method="fastgrid2"
fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", m.out.of.n=10, ci.bootstrap.size=10, verbose=verbose)
if (verbose) plot(fit.0); fit.0
out=c(
  fit.0$logliks[1],
  diff(fit.0$logliks)[1:3],
  fit.0$coefficients,
  fit.0$vcov$boot.samples[1,]
)    
checkEqualsNumeric(out, c(35.52723714,-0.15529784,-0.18206414,0.02324697,0.99337520,0.38349876,0.48532422,4.67409558,1.03713972,0.27900925,0.65087769,  5.32027458), tolerance=tolerance)




print("########  segmented hinge upperhinge") 
est.methods=c("grid","fastgrid2"); names(est.methods)=est.methods
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
for (type in c("quadupperhinge","quadhinge")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="quadupperhinge") c(32,2) else c(10,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0)
    if (verbose) plot(y~x, dat)
    est.methods=c("grid","fastgrid2")
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
    est.methods=c("grid","fastgrid2")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=c(rep(1,9),rep(5,10)))
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
est.methods=c("grid","fastgrid2")
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
est.methods=c("grid","fastgrid2")
out=sapply(est.methods, function(est.method){
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22c", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose, weights=rep(c(1,5),each=nrow(dat)/2))
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



print("########  M30 and M03")
for (type in c("cubicupperhinge","cubichinge")) {
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="cubicupperhinge") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.methods=c("grid","fastgrid2")
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
est.methods=c("grid","fastgrid2")
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
    est.methods=c("grid","fastgrid2")
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
print("########  x4hinge, x4upperhinge")
for(t in c("hinge", "upperhinge")){
    dat=sim.chngpt(mean.model="thresholded", threshold.type=paste0("cubic",t), n=250, seed=1, beta=if(t=="upperhinge") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.method="grid"
    fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=paste0("x4",t), est.method=est.method, var.type="none", save.boot=T, ci.bootstrap.size=1, verbose=verbose)
    if(verbose) plot(fit.0)
    if (t=="upperhinge") checkEqualsNumeric(fit.0$coefficients, c(0.03832806,0.36759453,31.75427906,1.60274645,-5.15254933,-0.02319485,5.00493816), tolerance=tolerance)
    if (t=="hinge") checkEqualsNumeric(fit.0$coefficients, c(-0.03536995,0.37418059,-9.81715060,-0.29244682,0.89944391,0.03742214,5.00493816), tolerance=tolerance)
}



print("########  offset ")
for (type in c("hinge","quadhinge")) {
    est.methods=c("grid","fastgrid2")
    out=sapply(est.methods, function(est.method){
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", offset=rep(1,nrow(pressure)), ci.bootstrap.size=1, verbose=verbose, weights=c(rep(1,9), rep(5,10)))
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


print("########  thinned thresholds, fastgrid will not be suppported")
est.method="grid"
type="segmented"    
dat=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")    
fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat,  type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=10, verbose=verbose)
out=c(
      fit$logliks[1],
      diff(fit$logliks[1:3]),
      fit$coefficients, 
      fit$vcov$boot.samples[1,]
    )  
if(verbose) print(out)
checkEqualsNumeric(out, c(544.7867707,1.3788716,2.1591613,-4.6760790,0.3340772,1.1283779,1.5775026,4.8043910,-4.6452678,0.3701558,1.1152278,1.7350673,4.9145547), tolerance=tolerance)    


## stratified
#for (type in c("segmented","hinge")) {
##type="upperhinge"
#    print(type)
#    dat=sim.chngpt("quadratic", n=20, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")
#    est.methods=c("grid","fastgrid2")
#    out=sapply(est.methods, function(est.method){
#        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, formula.strat=~I(z>0), est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=verbose)
#        c(
#          fit.0$logliks[1],
#          diff(fit.0$logliks),
#          fit.0$coefficients,
#          fit.0$vcov$boot.samples[1,]
#        )
#    })    
#    out
#    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
#}



}
