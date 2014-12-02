expit.2pl=function(t,e,b) sapply(t, function(t) 1/(1+exp(b*(t-e))))

sim.sigmoid = function (label, n, seed, alpha, beta, coef.z=log(1.4), x.distr="norm", e.=NULL, b.=NULL) {
    
    set.seed(seed)
    require(mvtnorm)
    
    mu=4.7 
    sd.x=1.6
    
    # distribution of x
    if(x.distr=="imb") { # imbalance
        x=c(rnorm(n-round(n/3), mu, sd=sd.x), mu-abs(rnorm(round(n/3), 0, sd=sd.x)))
    } else if(x.distr=="lin") { # unif
        x=runif(n)*4*sd.x + mu-2*sd.x
    } else if(x.distr=="mix") { # mixture
        x=c(rnorm(n*.6, mu, sd=sd.x), rep(mu-2*sd.x, n*.4))
    } else if(x.distr=="gam") { # gamma
        x=sd.x*scale(rgamma(n=n, 2.5, 1))+mu
        z=scale(rgamma(n=n, 2.5, 1))
        
    } else if(startsWith(x.distr,"norm")) { # normal
    
        if (x.distr=="norm") {
            rho=0
        } else if (x.distr=="norm3") {
            rho=0.3 
        } else if (x.distr=="norm6") {
            rho=0.6
        } else {
            stop("x.distr not supported: "%+%x.distr)
        }
        
        tmp=mvtnorm::rmvnorm(n, mean = c(mu,0), sigma = matrix(c(sd.x^2,sd.x*rho,sd.x*rho,1),2)) # use mvtnorm
        x=tmp[,1]
        z=tmp[,2]
    
    } else stop("x.distr not supported: "%+%x.distr)
    
    x.star = expit.2pl(x, e=e., b=b.)        
    
    if (label=="sigmoid1") {
        
        # null model is intercept only
        X=cbind(1, x.star)
        eta = X %*% c(alpha, beta)
        z=rep(1,nrow(X))# just so that every dataset has a z
    
    } else {
        
        if (label=="sigmoid2") {
            coef.=c(alpha, coef.z,     beta, 0)
        } else if (label=="sigmoid3" | label=="sigmoid6") {        
            coef.=c(alpha, coef.z, log(.67),  beta)
        } else if (label=="sigmoid4") {        
            coef.=c(alpha, coef.z, -log(.67), beta)         
        } else if (label=="sigmoid5") {        
            coef.=c(alpha, coef.z, 0,         beta)         
        } else {
            stop("label not supported: "%+%label)    
        }
        
        if (label=="sigmoid6") {
            #X=cbind(1, z, x.star, x.star*(z+3)^3)
            #X=cbind(1, z, x.star, x.star*z^2)
            X=cbind(1, z, x.star, x.star*z^3)
        } else {
            X=cbind(1, z, x.star, x.star*z)
        }
        
        eta = X %*% coef.
    
    }
    
    y=rbern(n, expit(eta))
    
    data.frame (
        y=y,
        z=z,         
        
        x=x,
        x.star=x.star,
        x.bin.med=ifelse(x>median(x), 1, 0),
        x.tri = factor(ifelse(x>quantile(x,2/3),"High",ifelse(x>quantile(x,1/3),"Medium","Low")), levels=c("Low","Medium","High")),
        
        x.tr.1=ifelse(x>log(100), x, 0) ,
        x.tr.2=ifelse(x>log(100), x, log(100)) ,
        x.tr.3=ifelse(x>3.5, x, 0) ,
        x.tr.4=ifelse(x>3.5, x, 3.5), 
        x.tr.5=ifelse(x>3.5, x, 3.5/2), 
        x.tr.6=ifelse(x>5, x, 5) ,
        x.ind =ifelse(x>3.5, 0, 1), 
        
        x.bin.35=ifelse(x>3.5, 1, 0), 
        x.bin.6=ifelse(x>6, 1, 0) ,
        x.bin.log100=ifelse(x>log(100), 1, 0)        
    )        
}
