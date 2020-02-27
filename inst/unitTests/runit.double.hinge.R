test.double.hinge <- function() {

library("chngpt")
library("RUnit")
  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")    
tolerance=1e-1
# R.Version()$system is needed b/c 32 bit system gives different results from 64 bit system
if((file.exists("D:/gDrive/3software/_checkReproducibility") | file.exists("~/_checkReproducibility")) & R.Version()$system %in% c("x86_64, mingw32","x86_64, linux-gnu")) tolerance=1e-6 
print(paste0("tol: ", tolerance))
verbose = FALSE


dat=sim.pastor(seed=1)[1:200,]
fit=double.hinge(x=-dat$x.star, y=dat$x.star.expit, var.type="bootstrap", ci.bootstrap.size=1e1)
plot.double.hinge(fit, lcol="red")
checkEqualsNumeric(fit$coefficients, c(-0.5413842, -0.4572128, 11.8805124,   0.06287312), tolerance=tolerance)    
checkEqualsNumeric(c(fit$ci.boot), c(-0.5475406,-0.5352278,-0.4615826,-0.4528430,10.4503459,13.3106790,0.02341853, 0.10232771), tolerance=tolerance)    

plot(fit$x, residuals(fit))
plot(fit$x, fitted(fit))
diag(vcov(fit))

# decreasing slope, need to specify lower.y and upper.y
dat=sim.pastor(seed=1)[1:200,]
fit=double.hinge(x=dat$x.star, y=dat$x.star.expit, lower.y=1, upper.y=0, var.type="bootstrap", ci.bootstrap.size=10)
checkEqualsNumeric(fit$coefficients, c(0.4572128,   0.5413842, -11.8805124, 0.06287312), tolerance=tolerance)    
plot.double.hinge(fit, lcol="red")



}
