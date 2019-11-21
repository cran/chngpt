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


dat=sim.pastor(seed=1)[1:200,]
fit=double.hinge(x=-dat$x.star, y=dat$x.star.expit, lower.y=0, upper.y=1, var.type="bootstrap", ci.bootstrap.size=10)
plot.double.hinge(fit, lcol="red")
checkEqualsNumeric(fit$coefficients, c(-0.5554024,-0.4524178,9.7101859), tolerance=tolerance)    
checkEqualsNumeric(c(fit$ci.boot), c(-0.5724834, -0.5383214, -0.4636889, -0.4411466,  7.3517072, 12.0686646), tolerance=tolerance)    


# decreasing slope
dat=sim.pastor(seed=1)[1:200,]
fit=double.hinge(x=dat$x.star, y=dat$x.star.expit, lower.y=1, upper.y=0, var.type="bootstrap", ci.bootstrap.size=10)
checkEqualsNumeric(fit$coefficients, c(0.4524178,0.5554024,-9.7101859), tolerance=tolerance)    
plot.double.hinge(fit, lcol="red")



}
