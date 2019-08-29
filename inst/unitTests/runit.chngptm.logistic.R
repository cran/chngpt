library("chngpt")
library("RUnit")
library("splines")
library("kyotil")

test.chngptm.logistic <- function() {

  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(tolerance)
verbose = FALSE


########################################
# test multiple interacting variables
n=1e3
set.seed(1)
X=data.frame(matrix(rnorm(3*n),ncol=3))
f=~I(X1>0)+X2*I(X1>0)+X3*I(X1>0)
# (Intercept)"      "I(X1 > 0)TRUE"    "X2"               "X3"               "I(X1 > 0)TRUE:X2" "I(X1 > 0)TRUE:X3"
lincomb=c(model.matrix(f, X) %*% c(0,1,0,0,1,1))
y=rbern(n,expit(lincomb))
dat=data.frame(y,X)
# use a binomial to test
fit=chngptm(formula.1=y~X2+X3, formula.2=~X2*X1+X3*X1, dat, type="step", family="binomial", grid.search.max=5, var.type="none")
# (Intercept)              X2              X3    I(X1>chngpt) X2:I(X1>chngpt) X3:I(X1>chngpt)          chngpt 
checkEqualsNumeric(c(coef(fit),fit$chngpt), c(0.20860034,-0.02784660,-0.08843246,0.65396890,1.02101488,1.24566539,0.02800216), tolerance=tolerance)    



########################################
# hinge model

data=sim.chngpt("thresholded", threshold.type="hinge", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="binomial")  
data$xx=data$x; data$zz=data$z; data$yy=data$y

## commented on 8/28/2019 due to a check error after submission to CRAN 
## fastgrid
#fit.2 = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="smoothapprox", var.type="none", verbose=verbose); fit.2
#fit.1 = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="grid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose); fit.1$vcov$boot.samples
#fit.3 = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=1, verbose=verbose); fit.3$vcov$boot.samples
#checkEqualsNumeric(coef(fit.1), coef(fit.3), tolerance=tolerance)    
#checkEqualsNumeric(fit.1$vcov$boot.samples, fit.3$vcov$boot.samples, tolerance=tolerance)    


fit.0=glm(yy~zz+ns(xx,df=3), data, family="binomial")
fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="smoothapprox", var.type="all", verbose=verbose, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
#    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.399209, tolerance=tolerance)    
checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 0.614832, tolerance=tolerance)    
checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 2.311749, tolerance=tolerance)        

#bootstrap
fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data[1:120,], type="hinge", est.method="smoothapprox", var.type="bootstrap", verbose=verbose, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, ci.bootstrap.size=10, boot.test.inv.ci=TRUE)
checkEqualsNumeric(fit$vcov$symm[1,], c(-0.4888763,0.6802883,-6.0128577,4.4632592), tolerance=tolerance)
checkEqualsNumeric(fit$vcov$testinv[,4], c(1.043623,6.013954), tolerance=tolerance)



########################################
# upperhinge model
dat=sim.chngpt (mean.model="thresholded", threshold.type="upperhinge", family="binomial", sd=0.3, mu.z=0, alpha=0, coef.z=log(1.4), beta=-1, n=100, seed=1)     
fit.0 = chngptm (formula.1=y~z, formula.2=~x, family="binomial", dat,  type="upperhinge", est.method="smoothapprox", var.type="model", save.boot=T, ci.bootstrap.size=1, verbose=verbose)
fit.1 = chngptm (formula.1=y~z, formula.2=~x, family="binomial", dat,  type="upperhinge", est.method="smoothapprox", var.type="robust", save.boot=T, ci.bootstrap.size=1, verbose=verbose, aux.fit=glm(y~z+ns(x,2), family="binomial", dat))
checkEqualsNumeric(fit.0$coef, c(0.2278046,   0.2595208,  -0.8485878, 5.4202994), tolerance=tolerance)    
checkEqualsNumeric(diag(fit.0$vcov), c(0.11440831,0.06818376,0.15998558,0.56502668), tolerance=tolerance)    
checkEqualsNumeric(diag(fit.1$vcov), c(0.16060624,0.05749946,0.24441672,1.06928181), tolerance=tolerance)        


data=sim.chngpt("quadratic", n=60, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="binomial")      
# weighted logistic models fit
fit.d = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data,  type="segmented", est.method="grid",         var.type="none", weights=rep(1:2,each=30), verbose=verbose)
fit.e = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data,  type="segmented", est.method="smoothapprox", var.type="none", weights=rep(1:2,each=30), verbose=verbose)
checkEqualsNumeric(coef(fit.d), c(-8.4120932,0.8365985,1.9086715,135.5846880), tolerance=tolerance)    
checkEqualsNumeric(coef(fit.e), c(-8.5621192,0.8431925,1.9465604,1713.3986973), tolerance=tolerance)    

# cbind() in formula
dat.2=sim.chngpt("thresholded", "step", n=200, seed=1, beta=1, alpha=-1, x.distr="norm", e.=4, family="binomial")
set.seed(1)
dat.2$success=rbinom(nrow(dat.2), 10, 1/(1 + exp(-dat.2$eta)))
dat.2$failure=10-dat.2$success
fit.2a=chngptm(formula.1=cbind(success,failure)~z, formula.2=~x, family="binomial", dat.2, type="step", est.method="grid", verbose=verbose)
fit.2b=chngptm(formula.1=cbind(success,failure)~z, formula.2=~x, family="binomial", dat.2, type="step", est.method="smoothapprox")
checkEqualsNumeric(fit.2a$coefficients, c(-0.8634819,0.3477191,0.9316376,4.0116612), tolerance=tolerance)    
checkEqualsNumeric(fit.2b$coefficients, c(-0.8637481,0.3477286,0.9255444,3.9907330), tolerance=tolerance)    
checkEqualsNumeric(diag(vcov(fit.2a$best.fit)), c(0.008068637,0.002473882,0.011023127), tolerance=tolerance)    
checkEqualsNumeric(diag(vcov(fit.2b$best.fit)), c(0.008207570,0.002474872,0.011138584), tolerance=tolerance)    
# compare cbind with single column outcome specification
n <- dat.2$success + dat.2$failure
dat.2$y.2 <- ifelse(n == 0, 0, dat.2$success/n)
dat.2$weights <- n
fit.2a1=chngptm(formula.1=y.2~z, formula.2=~x, family="binomial", dat.2, type="step", est.method="grid", weights=dat.2$weights)
fit.2b1=chngptm(formula.1=y.2~z, formula.2=~x, family="binomial", dat.2, type="step", est.method="smoothapprox", weights=dat.2$weights)
checkEqualsNumeric(fit.2a1$coefficients, c(-0.8634819,0.3477191,0.9316376,4.0116612), tolerance=tolerance)    
checkEqualsNumeric(fit.2b1$coefficients, c(-0.8637481,0.3477286,0.9255444,3.9907330), tolerance=tolerance)    
checkEqualsNumeric(diag(vcov(fit.2a1$best.fit)), c(0.008068637,0.002473882,0.011023127), tolerance=tolerance)    
checkEqualsNumeric(diag(vcov(fit.2b1$best.fit)), c(0.008207570,0.002474872,0.011138584), tolerance=tolerance)    

#    n <- dat.2$success + dat.2$failure
#    dat.2$y.2 <- ifelse(n == 0, 0, dat.2$success/n)
#    dat.2$weights <- n
#    stats::summary.glm(stats::glm(y.2~z+ I(z^2), family="binomial", dat.2, weights=dat.2$weights))
#    stats::summary.glm(stats::glm(cbind(success,failure)~z+ I(z^2), family="binomial", dat.2))


########################################
# segmented model

data=sim.chngpt("thresholded", threshold.type="segmented", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf, family="binomial")  
data$xx=data$x
data$zz=data$z
data$yy=data$y
fit.0=glm(yy~zz+xx+I(xx^2), data, family="binomial")
fit.0$coefficients=c(alpha=-1, z=log(1.4), x=-1 , x.quad=0.3)
fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="segmented", est.method="smoothapprox", var.type="all", verbose=verbose, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
#    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.9401789, tolerance=tolerance)    
checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 1.479249, tolerance=tolerance)    
checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 0.8593607, tolerance=tolerance)    


########################################
# step model

data=sim.chngpt("thresholdedItxn", threshold.type="step", family="binomial", n=250, seed=1, beta=-log(.67), beta.itxn=0, x.distr="norm", e.=3.4, b.transition=Inf, verbose=verbose)

# fit main effect

#    fit = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth", keep.best.fit=TRUE)
#    checkEqualsNumeric(mean(fit$coefficients), 0.9414254, tolerance=tolerance)
#    checkEqualsNumeric(mean(vcov(fit$best.fit)), 0.004475924, tolerance=tolerance)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth", keep.best.fit=TRUE)
#    checkEqualsNumeric(fit$coefficients, c( -0.8522692,0.3352071,0.5475595,3.7352043), tolerance=tolerance)
#    checkEqualsNumeric(mean(vcov(fit$best.fit)), 0.004475924, tolerance=tolerance)

#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth")    
#    checkEqualsNumeric(mean(fit$coefficients), 1.450894, tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit)), 0.02103244, tolerance=tolerance)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")    
#    checkEqualsNumeric(fit$coefficients, c(-0.5198031,0.3037050,0.2920897,5.7275861), tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit)), 0.02103244, tolerance=tolerance)

# this is an example that test leads to less noise
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=verbose, var.type="smooth")
#    checkEqualsNumeric(mean(fit$coefficients), 1.365853, tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit)), -0.006630736, tolerance=tolerance)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=verbose, var.type="smooth")
#    checkEqualsNumeric(fit$coefficients, c(-0.76090367,0.31253713,0.06159024,0.39690474,6.81840691), tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit)), -0.006539154, tolerance=tolerance)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=verbose, var.type="smooth", keep.best.fit=TRUE)
#    checkEqualsNumeric(mean(fit$coefficients), 0.8001427, tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit$best.fit)), 0.05447952, tolerance=tolerance)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=verbose, var.type="smooth", keep.best.fit=TRUE)
#    checkEqualsNumeric(fit$coefficients, c(-0.06705628,0.32633105,-0.28688468,0.76675233,0.32650919,3.73520431), tolerance=tolerance)    
#    checkEqualsNumeric(mean(vcov(fit$best.fit)), 0.05447952, tolerance=tolerance)
    
# fit interaction effect

fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox")
checkEqualsNumeric(mean(fit$coefficients), 1.356519, tolerance=tolerance)
fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0, ub.quantile=1, est.method="grid", verbose=verbose)
checkEqualsNumeric(fit$coefficients, c( -0.4815401,0.3042468,16.0476083,-0.3042468,8.3927654), tolerance=tolerance)

fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid")
checkEqualsNumeric(fit$coefficients, c(-0.5368853,0.2618009,0.2526734,0.1231553,5.3939234), tolerance=tolerance)
# Error in chngpt.test(formula.1, formula.2, data, type = type) : 
#  interaction model for this type not implemented yet: hinge
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox")
#    checkEqualsNumeric(fit$coefficients, c(-0.5247333,0.2731954,0.3058931,0.1291960,5.7165258), tolerance=tolerance)

fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=verbose)
checkEqualsNumeric(fit$coefficients, c(-1.27487926,2.11715992,0.20647617,-0.07942202,-0.74984522,0.81418233,2.65745247), tolerance=tolerance)
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=verbose)
#    checkEqualsNumeric(fit$coefficients, c(-0.78767807,0.71956572,0.05555813,0.16831871,-0.11267073,0.30441830,5.31461397), tolerance=tolerance)

fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=verbose)
checkEqualsNumeric(fit$coefficients, c(-0.77952531,0.89672361,0.05137674,0.41156747,-0.09905192,1.46771433,-0.30845593,-0.17022170,6.00352437), tolerance=tolerance)


########################################
# stegmented model

#    data=sim.chngpt("thresholded", threshold.type="stegmented", family="binomial", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.transition=Inf)
#    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")
#    checkEqualsNumeric(fit$coefficients, c(2.5067137,0.4783413,-0.5892674,1.3710912,-1.6098422,5.8467320), tolerance=tolerance)

}
