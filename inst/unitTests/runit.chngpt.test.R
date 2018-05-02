library("RUnit")
library("chngpt")

test.chngpt.test <- function() {

    RNGkind("Mersenne-Twister", "Inversion")    
    tolerance=1e-1
    # R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
    if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
    verbose=0
    
    #################################
    # linear regression

    # if cnt is 10, there is only one peak, if cnt is 100, a second, higher peak appears
    test = chngpt.test (formula.null=Volume~1, formula.chngpt=~Girth, family="gaussian", data=trees, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr"); # test; plot(test)
    checkEqualsNumeric(test$p.value, c(1e-04), tolerance=tolerance) 
    test = chngpt.test (formula.null=Volume~1, formula.chngpt=~Girth, family="gaussian", data=trees, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="score")
    checkEqualsNumeric(test$p.value, c(0.0017), tolerance=tolerance) 
    test = chngpt.test (formula.null=Volume~1, formula.chngpt=~Girth, family="gaussian", data=trees, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="score", robust=TRUE)
    checkEqualsNumeric(test$p.value, c(0.0107), tolerance=tolerance) 
    # param.boot    
    #formula.null=Volume~1; formula.chngpt=~Girth; family="gaussian"; data=trees; type="segmented"; mc.n=1e4; verbose=verbose; chngpts.cnt=100; test.statistic="lr"; p.val.method="param.boot"
    test = chngpt.test (formula.null=Volume~1, formula.chngpt=~Girth, family="gaussian", data=trees, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr", p.val.method="param.boot")
    
    data=sim.chngpt(mean.model="thresholded", family="gaussian", threshold.type="step", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=Inf, alpha=0)
    test = chngpt.test (formula.null=y~1, formula.chngpt=~z, family="gaussian", data=data, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr")
    checkEqualsNumeric(test$p.value, 0.1745, tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z+I(z^2), formula.chngpt=~x, family="gaussian", data=data, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr")
    checkEqualsNumeric(test$p.value, 0.5628, tolerance=tolerance) 
    test = chngpt.test (formula.null=y~1, formula.chngpt=~z, family="gaussian", data=data, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr", p.val.method="param.boot")
    checkEqualsNumeric(test$p.value, 0.1764, tolerance=tolerance) 
#    test = chngpt.test (formula.null=y~1, formula.chngpt=~z, family="gaussian", data=data, type="hinge", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr")
#    checkEqualsNumeric(test$p.value, 0, tolerance=tolerance) 
#    test = chngpt.test (formula.null=y~1, formula.chngpt=~z, family="gaussian", data=data, type="hinge", mc.n=1e4, verbose=verbose, chngpts.cnt=100, test.statistic="lr", p.val.method="param.boot")
#    checkEqualsNumeric(test$p.value, 0, tolerance=tolerance) 
    

    #################################
    # logistic regression
    data=sim.chngpt(mean.model="thresholded", family="binomial", threshold.type="step", n=250, seed=1, beta=-log(.67), x.distr="norm", e.=3.4, b.=Inf, alpha=chngpt::sim.alphas$thresholded_gam[1,1])
    
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="hinge", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="lr"); test
    checkEqualsNumeric(test$p.value, c(.3474), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="step", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="lr"); test
    checkEqualsNumeric(test$p.value, c(.689), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="lr"); test
    checkEqualsNumeric(test$p.value, c(0.7706), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="stegmented", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="lr"); test
    checkEqualsNumeric(test$p.value, c(0.8799), tolerance=tolerance) 
    
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="hinge", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="score"); test
    checkEqualsNumeric(test$p.value, c(.3355), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="step", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="score"); test
    checkEqualsNumeric(test$p.value, c(.6921), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, family="binomial", data, type="segmented", mc.n=1e4, verbose=verbose, chngpts.cnt=10, test.statistic="score"); test
    checkEqualsNumeric(test$p.value, c(0.7765), tolerance=tolerance) 


#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x, data, type="step", mc.n=1e3, interaction.method="lr.mc", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.618), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x, data, type="hinge", mc.n=1e3, interaction.method="lr.mc", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.36), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e4, interaction.method="lr.mc", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.6196), tolerance=tolerance) 
#    test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e4, verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.6422), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="lr.pastor", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(0.8552631), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="weighted.two.sided", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.646), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="weighted.one.sided", verbose=verbose, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.557), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="main.itxn", verbose=verbose, chngpts.cnt=50); test
#    checkEqualsNumeric(test$p.value, c(.679), tolerance=tolerance) 




}
