library("RUnit")
library("chngpt")


test.chngptm <- function() {

    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-1
    if(file.exists("D:/gDrive/3software/_checkReproducibility") & R.Version()$system %in% c("x86_64, mingw32")) tolerance=1e-6 # 32 bit system gives different results from 64 bit system
    print(tolerance)

    
    #######################################
    # linear model
    
    fit = chngptm (formula.1=Volume~1, formula.2=~Girth, family="gaussian", trees,  type="segmented", est.method="grid", var.type="none", verbose=2, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
    


    ########################################
    # hinge model
    
    data=sim.chngpt("sigmoid2", type="hinge", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.=-Inf)  
    data$xx=data$x
    data$zz=data$z
    data$yy=data$y
    fit.0=glm(yy~zz+ns(xx,df=3), data, family="binomial")
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="hinge", est.method="smoothapprox", var.type="all", verbose=2, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.3995812, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 0.6152094, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 2.329935, tolerance=tolerance)        
    #bootstrap
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data[1:100,], type="hinge", est.method="smoothapprox", var.type="bootstrap", verbose=2, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3, ci.bootstrap.size=10)
    checkEqualsNumeric(sum(diag(fit$vcov[["perc"]])), 1.22741, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["basic"]])), 0.9132618, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["bc"]])), 1.353463, tolerance=tolerance)    


    ########################################
    # segmented model
    
    data=sim.chngpt("sigmoid2", type="segmented", n=250, seed=1, beta=log(0.4), x.distr="norm", e.=4.1, b.=-Inf)  
    data$xx=data$x
    data$zz=data$z
    data$yy=data$y
    fit.0=glm(yy~zz+xx+I(xx^2), data, family="binomial")
    fit.0$coefficients=c(alpha=-1, z=log(1.4), x=-1 , x.quad=0.3)
    fit = chngptm (formula.1=yy~zz, formula.2=~xx, family="binomial", data, type="segmented", est.method="smoothapprox", var.type="all", verbose=2, aux.fit=fit.0, lb.quantile=0.1, ub.quantile=0.9, tol=1e-4, maxit=1e3)
    checkEqualsNumeric(sum(diag(fit$vcov[["smooth"]])), 0.9398962, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["model"]])), 1.479077, tolerance=tolerance)    
    checkEqualsNumeric(sum(diag(fit$vcov[["robust"]])), 0.8593664, tolerance=tolerance)    


    ########################################
    # step model, main effect
    
    data=sim.chngpt("sigmoid4", type="step", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)

    fit = chngptm (formula.1=y~z, formula.2=~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 0.9373754, tolerance=tolerance)
    checkEqualsNumeric(mean(vcov(fit)), 0.004481922, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c( -0.8522692,0.3352071,0.5475595,3.7352043), tolerance=tolerance)
    checkEqualsNumeric(mean(vcov(fit)), 0.004475924, tolerance=tolerance)
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", var.type="smooth")    
    checkEqualsNumeric(mean(fit$coefficients), 1.450928, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.007047492, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="hinge",     lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")    
    checkEqualsNumeric(fit$coefficients, c(-0.5198031,0.3037050,0.2920897,5.7275861), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.007046526, tolerance=tolerance)
    
    # this is an example that test leads to less noise
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 1.365648, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.02452156, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="segmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(-0.76090367,0.31253713,0.06159024,0.39690474,6.81840691), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.02455961, tolerance=tolerance)
    
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(mean(fit$coefficients), 1.035061, tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.0301887, tolerance=tolerance)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented", lb.quantile=0.1, ub.quantile=0.9, est.method="grid", verbose=FALSE, var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(-0.06705628,0.32633105,-0.28688468,0.76675233,0.32650919,3.73520431), tolerance=tolerance)    
    checkEqualsNumeric(mean(vcov(fit)), 0.05447952, tolerance=tolerance)
    
    
    
    ########################################
    # step model, interaction effect
    
    data=sim.chngpt("sigmoid4", type="step", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)

    fit = chngptm (y~z, ~x*z, family="binomial", data, tol=1e-4, maxit=1e3, type="step",      lb.quantile=0.1, ub.quantile=0.9, est.method="smoothapprox")
    checkEqualsNumeric(mean(fit$coefficients), 1.355475, tolerance=tolerance)
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
    
    data=sim.chngpt("sigmoid2", type="stegmented", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)
    fit = chngptm (y~z, ~x, family="binomial", data, tol=1e-4, maxit=1e3, type="stegmented",      lb.quantile=0.1, ub.quantile=0.9, est.method="grid", var.type="smooth")
    checkEqualsNumeric(fit$coefficients, c(2.5067137,0.4783413,-0.5892674,1.3710912,-1.6098422,5.8467320), tolerance=tolerance)
        



}
