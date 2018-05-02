sim.hinge=function(threshold.type=c("NA","hinge"),family=c("binomial","gaussian"),thres='NA',X.ditr='norm',mu.X,coef.X,cov.X,eps.sd,seed,n){
  threshold.type=match.arg(threshold.type)
  family=match.arg(family)
  npar=length(mu.X)
  B=coef.X
  set.seed(seed)
  if(X.ditr=='norm'){
    Z=mvrnorm(n, mu=mu.X, Sigma=cov.X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    z <- Z[,-npar]
    x <- Z[,npar]
    if(threshold.type=='hinge'){
      e <- thres
      Xe <- cbind(1,z,pmax(0,x-e))
      if(family=='binomial'){
        pr = 1/(1+exp(-Xe%*%B))
        Y = rbinom(n,1,pr)
      } else if(family=='gaussian'){
        Y=Xe%*%B+rnorm(n,0,eps.sd)
      }
 
    } else if(threshold.type=='NA'){
      X <- cbind(1,z,x)
      if(family=='binomial'){
        pr = 1/(1+exp(-X%*%B))
        Y = rbinom(n,1,pr)
      } else if(family=='gaussian'){
        Y=X%*%B+rnorm(n,0,eps.sd)
      }
    }
    dat=data.frame(Y,z,x)
    return(dat)
  }
}
