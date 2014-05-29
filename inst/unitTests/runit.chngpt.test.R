library("RUnit")
library("chngpt")

test.chngpt.test <- function() {


tolerance=1e-1
RNGkind("Mersenne-Twister", "Inversion")

data=sim.sigmoid("sigmoid4", n=250, seed=1, alpha=sim.alphas[["sigmoid4_norm"]]["3.4","0"], beta=0, x.distr="norm", e.=3.4, b.=-Inf)
if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance=1e-6

test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="lr.mc", verbose=0, chngpts.cnt=10); test
checkEqualsNumeric(test$p.value, c(.638), tolerance=tolerance) 
test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="lr.pastor", verbose=0, chngpts.cnt=10); test
checkEqualsNumeric(test$p.value, c(0.8552631), tolerance=tolerance) 
test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="weighted.two.sided", verbose=0, chngpts.cnt=10); test
checkEqualsNumeric(test$p.value, c(.646), tolerance=tolerance) 
test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="weighted.one.sided", verbose=0, chngpts.cnt=10); test
checkEqualsNumeric(test$p.value, c(.557), tolerance=tolerance) 
#test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="score", verbose=0, chngpts.cnt=10); test
#checkEqualsNumeric(test$p.value, c(.642), tolerance=tolerance) 
#test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="score.norm", verbose=0, chngpts.cnt=10); test
#checkEqualsNumeric(test$p.value, c(.672), tolerance=tolerance) 
#test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="score.power", verbose=0, chngpts.cnt=10); test
#checkEqualsNumeric(test$p.value, c(.494), tolerance=tolerance) 
#test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="score.power.norm", verbose=0, chngpts.cnt=10); test
#checkEqualsNumeric(test$p.value, c(.494), tolerance=tolerance) 
test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="main.itxn", verbose=0, chngpts.cnt=50); test
checkEqualsNumeric(test$p.value, c(.679), tolerance=tolerance) 




}
