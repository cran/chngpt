test.chngptm <- function() {

    library("chngpt")
    library("RUnit")
    library("splines")
    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-1
    # R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
    if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
    print(tolerance)
    
#    # performance unit testing
#    data=sim.chngpt("quadratic", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)      
#    system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 1)); system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 2))
#    data=sim.chngpt("quadratic", n=500, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)      
#    system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 1)); system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 2))
#    data=sim.chngpt("quadratic", n=1000, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)      
#    system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 1)); system.time(performance.unit.test (formula.1=y~z, formula.2=~x, family="gaussian", data, 200, 2))
    
    # ci.bootstrap.size being 1 or 5 affects the first bootstrap sample because of the way random numbers are generated and used. 
    data=sim.chngpt("quadratic", n=50, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)      

    # weighted linear models fit
    fit.a = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="grid", useC=FALSE, var.type="none", weights=rep(1:2,each=25))
    fit.b = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="fastgrid",         var.type="none", weights=rep(1:2,each=25), verbose=TRUE)
    fit.c = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="smoothapprox",     var.type="none", weights=rep(1:2,each=25), verbose=3)
    checkEqualsNumeric(coef(fit.a), coef(fit.b), tolerance=tolerance)    
    checkEqualsNumeric(coef(fit.b), coef(fit.c), tolerance=tolerance)    
    # weighted logistic models fit
    fit.d = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data,  type="segmented", est.method="grid", useC=FALSE, var.type="none", weights=rep(1:2,each=25), verbose=3)
    fit.e = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data,  type="segmented", est.method="smoothapprox",     var.type="none", weights=rep(1:2,each=25), verbose=3)
    checkEqualsNumeric(coef(fit.d), coef(fit.e), tolerance=tolerance)    
    
    fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="smoothapprox", var.type="bootstrap", ci.bootstrap.size=1, boot.test.inv.ci=TRUE, save.boot=TRUE)
    checkEqualsNumeric(fit$vcov$boot.samples[1,], c(-1.2739250,0.0393964,0.4140933,-0.4286910,5.4976075,26.1100387,34.2866744), tolerance=tolerance)        
    fit = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="smoothapprox", var.type="bootstrap", ci.bootstrap.size=5, boot.test.inv.ci=TRUE, save.boot=TRUE)
    checkEqualsNumeric(fit$vcov$boot.samples[1,], c(-1.12622283,0.09335172,0.39529622,-0.43627796,5.49760748,9.19697651,10.09715698), tolerance=tolerance)        
    checkEqualsNumeric(fit$vcov$symm[1,], c(-1.776669105,-0.009065366,0.104055843,-0.612082626,5.161310881), tolerance=tolerance)            
    # useC or not leads to different bootstrap results b/c bootstrap data are generated differently
    fit.1 = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="grid", var.type="bootstrap", ci.bootstrap.size=5, boot.test.inv.ci=TRUE, save.boot=TRUE, useC=FALSE)
    checkEqualsNumeric(fit.1$vcov$symm[1,], c(-2.2037845,-0.0121389,-0.0834454,-1.0730307,3.7849741), tolerance=tolerance)    
    fit.2 = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="grid", var.type="bootstrap", ci.bootstrap.size=5, boot.test.inv.ci=TRUE, save.boot=TRUE, useC=TRUE)
    checkEqualsNumeric(fit.2$vcov$symm[1,], c(-2.470325465,0.009685969,-0.073207089,-0.707083210,4.702745101), tolerance=tolerance)    
    fit.3 = chngptm (formula.1=y~z, formula.2=~x, family="gaussian", data,  type="segmented", est.method="fastgrid", var.type="bootstrap", ci.bootstrap.size=5, boot.test.inv.ci=TRUE, save.boot=TRUE, useC=TRUE)
    checkEqualsNumeric(fit.3$vcov$symm[1,], c(-2.470325465,0.009685969,-0.073207089,-0.707083210,4.702745101), tolerance=tolerance)    
    checkEqualsNumeric(coef(fit.2), coef(fit.3), tolerance=tolerance)    
    # Also note that sym CI also depends on est, thus est.method makes a difference. Although, now that est.method.boot is removed as an argument, this may not have any real impact

    #######################################
    # linear model
    fit = chngptm (formula.1=Volume~1, formula.2=~Girth, family="gaussian", trees,  type="segmented", est.method="grid", var.type="all", verbose=2, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, aux.fit=glm(Volume~ns(Girth,df=2),trees,family="gaussian"))
    checkEqualsNumeric(fit$coefficients, c(-24.614440,3.993966,4.266618,16.000000), tolerance=tolerance)    
    checkEqualsNumeric(diag(fit$vcov$model), c(18.2158018,0.1258957,1.1668183,0.4427713), tolerance=tolerance)    
    checkEqualsNumeric(diag(fit$vcov$sandwich), c(6.94919944,0.06118311,0.22564000,0.33254293), tolerance=tolerance)    
    checkEqualsNumeric(diag(fit$vcov$robust), c(779.922082,6.380663,8.091818,48.945107), tolerance=tolerance)    
    
    fit = chngptm (formula.1=Volume~1, formula.2=~Girth, family="gaussian", trees,  type="segmented", est.method="smoothapprox", var.type="model", verbose=2, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, aux.fit=1)
    checkEqualsNumeric(diag(fit$vcov), c(18.2158018,0.1258957,1.1668183,0.4427713), tolerance=tolerance)    
    checkEqualsNumeric(fit$coefficients, c(-24.614440,3.993966,4.266618,16.000000), tolerance=tolerance)    
    checkEqualsNumeric(attr(fit$vcov,"chngpt.ci"), c(13.3,16.3), tolerance=tolerance)        
    
    fit = chngptm (formula.1=Volume~1, formula.2=~Girth, family="gaussian", trees,  type="segmented", est.method="smoothapprox", var.type="model", verbose=2, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, b.transition=1)
    checkEqualsNumeric(diag(fit$vcov), c(18.7190518,0.1254883,1.1123345,0.5228118), tolerance=tolerance)    
    checkEqualsNumeric(fit$coefficients, c(-26.192998,4.158016,3.950498,16.000000), tolerance=tolerance)    
    fit = chngptm (formula.1=Volume~1, formula.2=~Girth, family="gaussian", trees,  type="segmented2", est.method="smoothapprox", var.type="model", verbose=2, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, b.transition=1)
    checkEqualsNumeric(diag(fit$vcov), c(22.1371998,0.1588287,0.0604084,0.3344677), tolerance=tolerance)    
    checkEqualsNumeric(fit$coefficients, c(-22.722485,3.816606,1.075842,17.900000), tolerance=tolerance)    
    

    ########################################
    # hinge model
    
    data=sim.chngpt("thresholded", threshold.type="hinge", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)  
    data$xx=data$x; data$zz=data$z; data$yy=data$y
    
    fit.0=glm(yy~zz+ns(xx,df=3), data, family="binomial")
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="smoothapprox", var.type="all", verbose=2, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.399209, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 0.614832, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 2.311749, tolerance=tolerance)        
    
    #bootstrap
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data[1:100,], type="hinge", est.method="smoothapprox", var.type="bootstrap", verbose=0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, ci.bootstrap.size=10, boot.test.inv.ci=TRUE)
    checkEqualsNumeric(fit$vcov$symm[1,], c(-0.7218213,0.5248393,-242.0967555,5.9602234), tolerance=tolerance)
    checkEqualsNumeric(fit$vcov$testinv[,4], c(2.484917,6.170364), tolerance=tolerance)
    
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="gaussian", data[1:100,], type="hinge", est.method="grid", var.type="bootstrap", verbose=0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, ci.bootstrap.size=100, save.boot=TRUE)
    checkEqualsNumeric(fit$vcov$boot.samples[2,], c(0.4613017,0.18780891,-0.1779347,4.416634,16.134385), tolerance=tolerance)
    
    
    ########################################
    # segmented model
    
    data=sim.chngpt("thresholded", threshold.type="segmented", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.transition=Inf)  
    data$xx=data$x
    data$zz=data$z
    data$yy=data$y
    fit.0=glm(yy~zz+xx+I(xx^2), data, family="binomial")
    fit.0$coefficients=c(alpha=-1, z=log(1.4), x=-1 , x.quad=0.3)
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="segmented", est.method="smoothapprox", var.type="all", verbose=2, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.9401789, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 1.479249, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 0.8593607, tolerance=tolerance)    


    ########################################
    # step model
    
    data=sim.chngpt("thresholdedItxn", threshold.type="step", family="binomial", n=250, seed=1, beta=-log(.67), beta.itxn=0, x.distr="norm", e.=3.4, b.transition=Inf, verbose=TRUE)
    
    # fit main effect
    
    fit = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 0.9414254, tolerance=tolerance)
    checkEqualsNumeric(mean(vcov(fit)), 0.004475924, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c( -0.8522692,0.3352071,0.5475595,3.7352043), tolerance=tolerance)
    checkEqualsNumeric(mean(vcov(fit)), 0.004475924, tolerance=tolerance)
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth")    
    checkEqualsNumeric(mean(fit$coefficients), 1.450894, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.02103244, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")    
    checkEqualsNumeric(fit$coefficients, c(-0.5198031,0.3037050,0.2920897,5.7275861), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.02103244, tolerance=tolerance)
    
    # this is an example that test leads to less noise
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 1.365853, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), -0.006630736, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(-0.76090367,0.31253713,0.06159024,0.39690474,6.81840691), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), -0.006539154, tolerance=tolerance)
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 1.034438, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.02996687, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(-0.06705628,0.32633105,-0.28688468,0.76675233,0.32650919,3.73520431), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.05447952, tolerance=tolerance)
        
    # fit interaction effect
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox")
    checkEqualsNumeric(mean(fit$coefficients), 1.356519, tolerance=tolerance)
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0, ub.quantile=1, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c( -0.4815401,0.3042468,16.0476083,-0.3042468,8.3927654), tolerance=tolerance)
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid")
    checkEqualsNumeric(fit$coefficients, c(-0.5368853,0.2618009,0.2526734,0.1231553,5.3939234), tolerance=tolerance)
# Error in chngpt.test(formula.1, formula.2, data, type = type) : 
#  interaction model for this type not implemented yet: hinge
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox")
#    checkEqualsNumeric(fit$coefficients, c(-0.5247333,0.2731954,0.3058931,0.1291960,5.7165258), tolerance=tolerance)
    
    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-1.27487926,2.11715992,0.20647617,-0.07942202,-0.74984522,0.81418233,2.65745247), tolerance=tolerance)
#    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=FALSE)
#    checkEqualsNumeric(fit$coefficients, c(-0.78767807,0.71956572,0.05555813,0.16831871,-0.11267073,0.30441830,5.31461397), tolerance=tolerance)

    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE)
    checkEqualsNumeric(fit$coefficients, c(-0.77952531,0.89672361,0.05137674,0.41156747,-0.09905192,1.46771433,-0.30845593,-0.17022170,6.00352437), tolerance=tolerance)


    ########################################
    # stegmented model
    
    data=sim.chngpt("thresholded", threshold.type="stegmented", family="binomial", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.transition=Inf)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(2.5067137,0.4783413,-0.5892674,1.3710912,-1.6098422,5.8467320), tolerance=tolerance)
        



}
