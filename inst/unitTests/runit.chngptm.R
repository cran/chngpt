### --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("chngpt")
}

test.chngptm <- function() {

tolerance=1e-3

# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance=1e-6
}
 
RNGkind("Mersenne-Twister", "Inversion")


data=sim.sigmoid("sigmoid4", n=250, seed=1, alpha=sim.alphas[["sigmoid4_norm"]]["3.4","0"], beta=0, x.distr="norm", e.=3.4, b.=-Inf)
formula.null=y~z
formula.chngpt=~x*z

fit.1 = chngptm (formula.null, formula.chngpt, data, tol=1e-4, maxit=1e3)
checkEqualsNumeric(fit.1$coefficients, c( -0.4815401,    0.3042468,   16.0476083,   -0.3042468,    8.3927654), tolerance=tolerance)

fit.1 = chngptm (formula.null, formula.chngpt, data, tol=1e-4, maxit=1e3, search.all.thresholds=FALSE)
checkEqualsNumeric(fit.1$coefficients, c(-0.5391406,0.2158219,0.4210706,0.5553366,6.0611921), tolerance=tolerance)



}
