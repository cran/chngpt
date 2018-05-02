hinge.test <- function(formula, cov.interest, family=c("binomial","gaussian"), data, 
    thres=NA,lb.quantile=.1,ub.quantile=.9,chngpts.cnt=10,method=c("FDB","B","DB"),boot.B=1e4,B2=NA, verbose=FALSE){
    
    DNAME = deparse(substitute(data))
    family<-match.arg(family)        
    method<-match.arg(method)        
    
    # keep only records that have no missing data for the null model and the change point model
    subset.1 = complete.cases(model.frame(formula, data, na.action=na.pass))
    subset.2 = complete.cases(model.frame(as.formula("~"%+%cov.interest), data, na.action=na.pass))
    if(!all(subset.1 & subset.2)) warning("Missing data are being exluded.")
    data=data[subset.1 & subset.2,,drop=FALSE]
    
    # Y is the response, 
    # z is the variables without the variable supposed to have a threshold, 
    # x is the variable we want to test the existence of a threshold, 
    Y=model.frame(formula, data)[,1]
    z=model.matrix(formula, data)[,-1]
    x=data[[cov.interest]]
    n=length(Y)
        
    M <- chngpts.cnt

    # save rng state before set.seed in order to restore before exiting this function
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    set.seed(1)    
    
    # parametric bootstrap and fast double bootstrap
    if(method!='DB'){
    
        ###################################
        # maximize over candidate thresholds
        if(is.na(thres)){
                
            em <- quantile(x, probs = seq(lb.quantile,ub.quantile,length.out = M))
            X <- cbind(1,z,x)
            xem <- matrix(0, nrow = nrow(X), ncol = M)
            for(m in 1:M){
            xem[,m] <- pmax(0,x-em[m])
            }
            
            fit <- glm(Y~z+x, family = family)
            fitem <- 0
            for(m in 1:M){
                fitem[m] <- logLik(glm(Y~z+xem[,m], family = family))
            }
                
            TeM <- 1/n*(max(fitem)-logLik(fit))
            
            TeMb <- 0
            mean.b=matrix(0,nrow=n,ncol=boot.B)
            sd.b=rep(0,boot.B)
            for(b in 1:boot.B){
                if(family=='binomial'){
                    Yb <- rbinom(n, 1, fit$fitted.values)
                } else if(family=='gaussian'){
                    Yb <- rnorm(n, mean=fit$fitted.values, sd=mean(fit$residuals^2))
                }
                fitb <- glm(Yb~z+x, family = family)
                fitemb <- 0
                for(m in 1:M){
                    fitemb[m] <- logLik(glm(Yb~z+xem[,m], family = family))
                }
                
                TeMb[b] <- 1/n*(max(fitemb)-logLik(fitb))
                mean.b[,b] <- fitb$fitted.values
                sd.b[b]=mean(fitb$residuals^2)
            }
            p <- mean(TeMb > TeM)
            p.value <- p
            
            statistic=TeM
        
            # fast double bootstrap
            if(method=='FDB'){
                TeMfdb <- 0
                for(b in 1:boot.B){
                if(family=='binomial'){
                    Yfdb <- rbinom(n, 1, mean.b[,b])
                } else if(family=='gaussian'){
                    Yfdb <- rnorm(n,mean=mean.b[,b],sd=sd.b[b])
                }
                fitfdb <- glm(Yfdb~z+x, family = family)
                fitemfdb <- 0
                for(m in 1:M){
                    fitemfdb[m] <- logLik(glm(Yfdb~z+xem[,m], family = family))
                }
                
                TeMfdb[b] <- 1/n*(max(fitemfdb)-logLik(fitfdb))
                }
                # fast double bootstrap p-value
                adj.p <- mean(TeMb > quantile(TeMfdb, 1-p))
                p.value <- adj.p
            }
            
    
        ###################################
        # candidate threshold provided
        } else if(!is.na(thres)){        
        
            e=thres
            xe <- pmax(0,x-e)
            X <- cbind(1,z,x)
            
            fit <- glm(Y~z+x, family = family)
            fite <- glm(Y~z+xe, family = family) 
            
            Te <- 1/n*(logLik(fite)-logLik(fit))
            statistic=Te
            
            Teb <- 0
            mean.b=matrix(0,nrow=n,ncol=boot.B)
            sd.b=rep(0,boot.B)
            for(b in 1:boot.B){
                if(family=='binomial'){
                    Yb <- rbinom(n, 1, fit$fitted.values)
                } else if(family=='gaussian'){
                    Yb <- rnorm(n, mean=fit$fitted.values, sd=mean(fit$residuals^2))
                }
                fitb <- glm(Yb~z+x, family = family)
                fiteb <- glm(Yb~z+xe, family = family) 
                Teb[b] <- 1/n*(logLik(fiteb)-logLik(fitb))
                mean.b[,b] <- fitb$fitted.values
                sd.b[b]=mean(fitb$residuals^2)
            }  
            p <- mean(Teb > Te)
            p.value <- p
            
            if(method=='FDB'){
                # fast double bootstrap
                Tefdb <- 0
                for(b in 1:boot.B){
                    if(family=='binomial'){
                      Yfdb <- rbinom(n, 1, mean.b[,b])
                    } else if(family=='gaussian'){
                      Yfdb <- rnorm(n,mean=mean.b[,b],sd=sd.b[b])
                    }
                                    
                    fitfdb <- glm(Yfdb~z+x, family = family)
                    fitefdb <- glm(Yfdb~z+xe, family = family) 
                    
                    Tefdb[b] <- 1/n*(logLik(fitefdb)-logLik(fitfdb))
                }
                # fast double bootstrap p-value
                adj.p <- mean(Teb > quantile(Tefdb, 1-p))
                p.value <- adj.p
            }
        } #end if (!is.na(thres))
      

    ###################################
    # double bootstrap
    } else if(method=='DB'){
        if(family!="gaussian" | is.na(thres)) stop("Only Gaussian implemented for method DB and for a single threshold")
        
        e=thres
        xe <- pmax(0,x-e)
        X <- cbind(1,z,x)
        Xe <- cbind(1,z,xe)
        
        H <- X%*%solve(t(X)%*%X)%*%t(X)
        He <- Xe%*%solve(t(Xe)%*%Xe)%*%t(Xe)
        P <- diag(n)-H
        Pe <- diag(n)-He
        
        sig2.hat <- drop(1/n*(t(Y)%*%P%*%Y))
        sige2.hat <- drop(1/n*(t(Y)%*%Pe%*%Y))
        Te <- 1/2*(log(sig2.hat)-log(sige2.hat))
        statistic=Te
        
        Yb <- mvrnorm(boot.B, mu=H%*%Y, Sigma=diag(sig2.hat,n), tol = 1e-6, empirical = F, EISPACK = FALSE)
        sig2.hat.b <- 1/n*rowSums((Yb%*%P)*Yb)
        sige2.hat.b <- 1/n*rowSums((Yb%*%Pe)*Yb)
        Teb <- 1/2*(log(sig2.hat.b)-log(sige2.hat.b))
        p <- mean(Teb > Te)
        # double bootstrap
        pb <- 0
        for(b in 1:B2){
            Ymb <- mvrnorm(boot.B, mu=H%*%Yb[b,], Sigma=diag(sig2.hat.b[b],n), tol = 1e-6, empirical = F, EISPACK = FALSE)
            sig2.hat.mb <- 1/n*rowSums((Ymb%*%P)*Ymb)
            sige2.hat.mb <- 1/n*rowSums((Ymb%*%Pe)*Ymb)
            Temb <- 1/2*(log(sig2.hat.mb)-log(sige2.hat.mb))
            pb[b] <- mean(Temb > Teb[b])
        }
        
        p.value <- mean(pb < p)        
    } 
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv) 
    
    res=list(family=family,method="Maximal likelihood ratio test by "%+%method,p.value=p.value, statistic=statistic, data.name=DNAME)
    names(res$statistic)="Maximal likelihood ratio statistic"
    
    class(res)=c('htest',class(res))
    return(res) 
} 
