library("RUnit")
library("chngpt")
    
test.chngpt.test <- function() {

tolerance=1e-3
# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance=1e-6
} else if (substr(R.Version()$system, 1, 10)=="sparc-sun-" | substr(R.Version()$system, 1, 4)=="i686") { 
    tolerance=1e-1
}
RNGkind("Mersenne-Twister", "Inversion")

data=sim.sigmoid("sigmoid4", n=250, seed=1, alpha=sim.alphas[["sigmoid4_norm"]]["3.4","0"], beta=0, x.distr="norm", e.=3.4, b.=-Inf)

set.seed(1)
x.pre=runif(1)
set.seed(1)
test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="main.itxn", chngpts.cnt=50)
x.aft=runif(1)
checkEqualsNumeric(test$p.value, c(.679), tolerance=tolerance) # can depend on a number things, including 0/1 or logistic weights
checkEqualsNumeric(x.pre, x.aft, tolerance=tolerance) # check if calling chngpt.test changes the rng state


test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, mc.n=1e3, interaction.method="lr", verbose=0, chngpts.cnt=10); test

# check if chngpt.test is invariant to affine transformation of z
data$z.1=100000*data$z+2000
test.1 = chngpt.test (formula.null=y~z.1, formula.chngpt=~x*z.1, data, mc.n=1e3, interaction.method="lr", verbose=0, chngpts.cnt=10); test.1


}
