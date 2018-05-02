library("RUnit")
library("chngpt")

test.hinge.test <- function() {

    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-6
    # R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
    if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
    verbose=0
    
        
    #########################################################################

    # logistic regression

    # Tl-B Tl-FDB under null with multivariate z  
    dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.4),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
    
    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
    checkEqualsNumeric(test$p.value, 0.70, tolerance=tolerance) 
    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
    checkEqualsNumeric(test$p.value, 0.76, tolerance=tolerance) 
    
    
#    # Tl-B Tl-FDB under alternative e=0 with univariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'binomial',thres=0,X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.4),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.33, tolerance=tolerance) 
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.26, tolerance=tolerance)
#    
#    # Tl-B Tl-FDB under alternative e=0.67 with univariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'binomial',thres=0.67,X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.4),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.49, tolerance=tolerance) 
#    
#    # Tl-B Tl-FDB under alternative e=0 with multivariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'binomial',thres=0,X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.4),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.72, tolerance=tolerance) 
#    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.68, tolerance=tolerance)
#    
#    # Tl-B Tl-FDB under null with univariate z  
#    dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.4),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.77, tolerance=tolerance) 
#    
#    
#    # l-B l-FDB under null size for thres=0 with univariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.4),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.39, tolerance=tolerance) 
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.41, tolerance=tolerance) 
#
#    # l-B l-FDB under null size for thres=0.67 with univariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.4),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = 0.67,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.41, tolerance=tolerance) 
#    test=hinge.test(Y~z, "x", family="binomial", data=dat,  thres = 0.67,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.46, tolerance=tolerance) 
#
#    # l-B l-FDB under null size for thres=0 with multivariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'binomial',thres='NA',X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.4),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.67, tolerance=tolerance) 
#    test=hinge.test(Y~X1+X2, "x", family="binomial", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.64, tolerance=tolerance) 
#    
    
    #########################################################################
    # linear regression

    # TL-B TL-FDB under null with multivariate z
    dat=sim.hinge(threshold.type = 'NA',family = 'gaussian',thres='NA',X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.5),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
    
    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
    checkEqualsNumeric(test$p.value, .9,tolerance=tolerance)
    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
    checkEqualsNumeric(test$p.value, .97,tolerance=tolerance)    
    
#    # L-B L-DB L-FDB under null size for thres=0 with multivariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'gaussian',thres='NA',X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.5),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.85, tolerance=tolerance)
#    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='DB',boot.B=1e2,B2=50); test$p.value
#    checkEqualsNumeric(test$p.value, 0.76, tolerance=tolerance)
#    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .93,tolerance=tolerance)
#    
#    
#    # L-B L-DB L-FDB under null size for thres=0 with univariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'gaussian',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.5),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.52, tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='DB',boot.B=1e2,B2=50); test$p.value
#    checkEqualsNumeric(test$p.value, 0.5, tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .5,tolerance=tolerance)
#    
#    # L-B L-DB L-FDB under null size for thres=0.67
#    dat=sim.hinge(threshold.type = 'NA',family = 'gaussian',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.5),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0.67,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, 0.29, tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0.67,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='DB',boot.B=1e2,B2=50); test$p.value
#    checkEqualsNumeric(test$p.value, 0.24, tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = 0.67,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .22,tolerance=tolerance)
#
#    # TL-B TL-FDB under alternative e=0 with univariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'gaussian',thres=0,X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.5),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .07,tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .12,tolerance=tolerance)
#    
#    # TL-B TL-FDB under alternative e=0.67 with univariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'gaussian',thres=0.67,X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.5),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .05,tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .07,tolerance=tolerance)
#
#    # TL-B TL-FDB under alternative e=0 with multivariate z
#    dat=sim.hinge(threshold.type = 'hinge',family = 'gaussian',thres=0,X.ditr = 'norm',mu.X = c(0,0,0),coef.X = c(0,.5,.5,.5),cov.X = diag(3),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .78,tolerance=tolerance)
#    test=hinge.test(Y~X1+X2, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .83,tolerance=tolerance)
#    
#    # TL-B TL-FDB under null with univariate z
#    dat=sim.hinge(threshold.type = 'NA',family = 'gaussian',thres='NA',X.ditr = 'norm',mu.X = c(0,0),coef.X = c(0,.5,.5),cov.X = diag(2),eps.sd = 1,seed = 1,n=100)
#    
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='B',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .98,tolerance=tolerance)
#    test=hinge.test(Y~z, "x", family="gaussian", data=dat,  thres = NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,'method'='FDB',boot.B=1e2); test$p.value
#    checkEqualsNumeric(test$p.value, .98,tolerance=tolerance)

}
