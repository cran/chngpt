test.chngptm.linear <- function() {

library("chngpt")
library("RUnit")
  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(tolerance)
verbose = FALSE


# step model m out of n subsampling
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=20, seed=1, beta=-log(.67), alpha=1)
est.method="fastgrid2"
fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", m.out.of.n=10, ci.bootstrap.size=10, verbose=verbose)
plot(fit.0); fit.0
out=c(
  fit.0$logliks[1],
  diff(fit.0$logliks)[1:3],
  fit.0$coefficients,
  fit.0$vcov$boot.samples[1,]
)    
checkEqualsNumeric(out, c(35.52723714,-0.15529784,-0.18206414,0.02324697,0.99337520,0.38349876,0.48532422,4.67409558,1.03713972,0.27900925,
0.65087769,  5.32027458), tolerance=tolerance)



# step model
print("step")
dat=sim.chngpt("thresholded", threshold.type="step", family="gaussian", n=250, seed=1, beta=-log(.67), alpha=1)
est.methods=c("grid","fastgrid2")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, family="gaussian", dat, type="step", est.method=est.method, var.type="bootstrap", m.out.of.n=FALSE, ci.bootstrap.size=1, verbose=verbose)
    plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1],
      diff(fit.0$logliks)[1:3],
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    ))
}    
colnames(out)=est.methods
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)



# x4hinge, x4upperhinge
for(t in c("hinge", "upperhinge")){
    dat=sim.chngpt(mean.model="thresholded", threshold.type=paste0("cubic",t), n=250, seed=1, beta=if(t=="upperhinge") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.method="grid"
    fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=paste0("x4",t), est.method=est.method, var.type="none", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
    plot(fit.0)
    if (t=="upperhinge") checkEqualsNumeric(fit.0$coefficients, c(0.03832806,0.36759453,31.75427906,1.60274645,-5.15254933,-0.02319485,5.00493816), tolerance=tolerance)
    if (t=="hinge") checkEqualsNumeric(fit.0$coefficients, c(-0.03536995,0.37418059,-9.81715060,-0.29244682,0.89944391,0.03742214,5.00493816), tolerance=tolerance)
}



# M33c
dat=sim.chngpt(mean.model="thresholded", threshold.type="M33c", n=100, seed=1, beta=c(32,2,-5,5), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("grid","fastgrid2")
out=NULL
for (est.method in est.methods) {
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M33c", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=1)
    plot(fit.0); fit.0
    out=cbind(out, c(
      fit.0$logliks[1]
      ,diff(fit.0$logliks)
      ,fit.0$coefficients
      ,fit.0$vcov$boot.samples[1,]
    ))
} 
colnames(out)=est.methods
print(head(out))
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)


# M31 and M13
for (type in c("M31","M13")) {
    print(type)
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="M31") c(32,2,-5,5) else c(-10,-.25,1,-2), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.methods=c("grid","fastgrid2")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
        plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    print(head(out))
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# cubicupperhinge and cubichinge
for (type in c("cubicupperhinge","cubichinge")) {
    print(type)
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="cubicupperhinge") c(32,2,-5) else c(-10,-.25,1), x.distr="lin", e.=5, b.transition=Inf, family="gaussian", alpha=0); plot(y~x, dat)
    est.methods=c("grid","fastgrid2")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
        plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    print(head(out))
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# M22c
dat=sim.chngpt(mean.model="thresholded", threshold.type="M22c", n=100, seed=1, beta=c(32,2,32,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("grid","fastgrid2")
out=sapply(est.methods, function(est.method){
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22c", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=1)
    plot(fit.0); fit.0
    c(
      fit.0$logliks[1],
      diff(fit.0$logliks),
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    )
})    
print(head(out))
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)



# M22
dat=sim.chngpt(mean.model="thresholded", threshold.type="M22", n=100, seed=1, beta=c(32,2,10,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0, sd=0, coef.z=0)
est.methods=c("grid","fastgrid2")
out=sapply(est.methods, function(est.method){
    fit.0=chngptm(formula.1=y~z, formula.2=~x, dat, type="M22", family="gaussian", est.method=est.method, var.type="bootstrap", ci.bootstrap.size=1, verbose=1)
    plot(fit.0); fit.0
    c(
      fit.0$logliks[1],
      diff(fit.0$logliks),
      fit.0$coefficients,
      fit.0$vcov$boot.samples[1,]
    )
})    
print(head(out))
for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    



# M21 and M12
for (type in c("M21","M12")) {
    est.methods=c("grid","fastgrid2")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", offset=rep(1,nrow(pressure)), ci.bootstrap.size=1, verbose=1)
        plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    print(head(out))
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


#type="quadhinge"
for (type in c("quadupperhinge","quadhinge")) {
    print(type)
    dat=sim.chngpt(mean.model="thresholded", threshold.type=type, n=250, seed=1, beta=if(type=="quadupperhinge") c(32,2) else c(10,10), x.distr="norm", e.=6, b.transition=Inf, family="gaussian", alpha=0)
    plot(y~x, dat)
    est.methods=c("grid","fastgrid2")
    out=NULL
    for (est.method in est.methods) {
        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
        plot(fit.0); fit.0
        out=cbind(out, c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        ))
    }
    colnames(out)=est.methods
    print(head(out))
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


# offset 
for (type in c("hinge","quadhinge")) {
    est.methods=c("grid","fastgrid2")
    out=sapply(est.methods, function(est.method){
        fit.0=chngptm(formula.1=pressure~-1, formula.2=~temperature, pressure, type=type, family="gaussian", est.method=est.method, var.type="bootstrap", offset=rep(1,nrow(pressure)), ci.bootstrap.size=1, verbose=1)
        plot(fit.0); fit.0
        c(
          fit.0$logliks[1],
          diff(fit.0$logliks),
          fit.0$coefficients,
          fit.0$vcov$boot.samples[1,]
        )
    })    
    print(head(out))
    for (m in est.methods) checkEqualsNumeric(out[,"grid"], out[,m], tolerance=tolerance)    
}


## stratified
#for (type in c("segmented","hinge")) {
##type="upperhinge"
#    print(type)
#    dat=sim.chngpt("quadratic", n=20, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="gaussian")
#    est.methods=c("grid","fastgrid2")
#    out=sapply(est.methods, function(est.method){
#        fit.0=chngptm (formula.1=y~z, formula.2=~x, family="gaussian", dat, type=type, formula.strat=~I(z>0), est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, grid.search.max=1e5, verbose=verbose)
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


# thinned thresholds, fastgrid will not be suppported
est.methods=c("grid","fastgrid2")# fastgrid not implemented weights yet. gridC used to be supported, but not when we make all .C upperhinge
for (type in c("segmented","hinge","upperhinge")) {
#type="segmented"    
    if (verbose) print(type)
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
#type="segmented"
    if (verbose) print(type)
    dat=sim.chngpt (mean.model="thresholded", threshold.type=type, family="gaussian", sd=0.3, mu.z=0, alpha=0, coef.z=log(1.4), beta=-1, n=100, seed=1)     
    out=sapply(est.methods, function(est.method){
        # ci.bootstrap.size is set to 1 because grid uses boot pkg, which generates different bootsamples if bootstrap sample size is >1
        fit.0=chngptm (y~z, ~x, "gaussian", dat, type=type, est.method=est.method, var.type="bootstrap", save.boot=T, ci.bootstrap.size=1, verbose=verbose) 
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
    if (verbose) print(type)
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
