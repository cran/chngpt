# chngpts=NULL; lb.quantile=.1; ub.quantile=.9; chngpts.cnt=50; b.=-Inf; mc.n=5e4; ret.p.val=TRUE
chngpt.score.test = function(formula.null, formula.chngpt, data, 
    interaction.only=TRUE, chngpts=NULL, lb.quantile=.1, ub.quantile=.9, chngpts.cnt=50, b.=-30, mc.n=5e4, ret.p.val=TRUE, verbose=TRUE 
) {
    
    DNAME = deparse(substitute(data))
    
    y=model.frame(formula.null, data)[,1]
    Z=model.matrix(formula.null, data)
    tmp=as.matrix(model.frame(formula.chngpt, data,na.action=na.pass))
    if (nrow(Z)!=nrow(tmp)) stop("number of records do not match between formula.null and formula.chngpt")
    
    # tmp may have NA, remove those rows from all data
    y=y[complete.cases(tmp)]
    Z=Z[complete.cases(tmp),,drop=FALSE]
    data=data[complete.cases(tmp),,drop=FALSE]
    tmp=tmp[complete.cases(tmp),,drop=FALSE]
    n=nrow(Z)

    chng.var = tmp[,setdiff(colnames(tmp), colnames(Z))[1]]
    z.1 = tmp[,intersect(colnames(tmp), colnames(Z))] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itx = length(intersect(colnames(tmp), colnames(Z)))>0        
        
    if (is.null(chngpts)) {
        chngpts=quantile(chng.var, seq(lb.quantile,ub.quantile,length=chngpts.cnt))
    }
    p=length(chngpts)  
    if (has.itx & !interaction.only) p = p*2  
    
    fit=keepWarnings(glm(formula.null,  data=data, family="binomial"))    # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit$warning)!=0) {
        return (NA)
    } else {
        fit=fit$value
    }
    
    beta.h = coef(fit)
    mu.h = drop(expit(Z %*% beta.h))    
    D.h = diag(c(mu.h*(1-mu.h)))
    V.beta.h = solve(t(Z) %*% D.h %*% Z)
    V.eta.h = Z %*% V.beta.h %*% t(Z)
    A.h = diag(n) - D.h %*% V.eta.h
    V.tmp = A.h %*% D.h %*% t(A.h)
    
    # when beta.0 is not null, compute V.tmp using true values of beta
#    if (!is.null(beta.0)) {
#        mu = drop(expit(Z %*% beta.0))    
#        D = diag(c(mu*(1-mu)))
#        V.beta = solve(t(Z) %*% D %*% Z)
#        V.eta = Z %*% V.beta %*% t(Z)
#        A = diag(n) - D %*% V.eta
#        V.tmp = A %*% D %*% t(A)
#    }
    
    W=sapply(chngpts, function (e.){
        u = exp(b.*(chng.var-e.))
        w = (1/(1+u))/sum(1/(1+u))
    })
    if (has.itx) {
        if (interaction.only) {
            W = W * z.1
        } else {
            W = cbind(W, W * z.1)
        }
    }
    
    V.S.hat = t(W) %*% V.tmp %*% W
    
#    qr(V.S.hat, tol = 1e-8)$rank 
#    isSymmetric(A.h)
    
    TT.0=c(t(W) %*% (y - mu.h) / sqrt(diag(V.S.hat)))
    TT = abs(TT.0)
    
    T.max=max(TT)
    names(T.max)="Maximum score statistics"
    
    res=list()
    res$chngpts=chngpts
    res$TT=TT
    res$statistic=T.max
    res$chngpt = chngpts[which.max(TT)]
    res$parameter=NULL
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$method="Max-Score Change Point Test"
    res$data.name=DNAME
    res$V.S.hat = V.S.hat
    res$p.value=NA
    
    if (ret.p.val) {                
        #### method 1 for findind p-value
        # one can also find the quantile of max of multivariate normal by the following, but it actually takes 3 times as long
        #qmvnorm(.975, interval = c(-10, 10), tail = c("lower.tail"), mean = 0, corr = cov2cor(Sigma), sigma = NULL, maxpts = 25000, abseps = 0.001, releps = 0)        
        require(MASS)
        set.seed(1)
        sam=mvrnorm (n=mc.n, mu=rep(0,p), Sigma=cov2cor(V.S.hat))
        sam=abs(sam)            
        # there are several programming methods to get the max
        #x.max=rowMaxs(sam) #from matrixStats is slowest
        #x.max=pmax(sam[,1], sam[,2], sam[,3], sam[,4], sam[,5], sam[,6], sam[,7], sam[,8], sam[,9], sam[,10]) # faster than apply, but hard coded
        # the following is a little slower than doing pmax as above, but acceptable for n=1e5
        tmp = apply(sam, 2, function(aux) list(aux))
        tmp = do.call(c, tmp)
        x.max=do.call(pmax, tmp)            
        p.value = mean(x.max>T.max)
        res$p.value = p.value
        
#        if (ret.p.val.2) {
#            #### method 2 for finding p-value, not very good
#            sqrt.inv = solve(chol(cov2cor(V.S.hat))) # t(sqrt.inv) %*% cov2cor(V.S.hat) %*% sqrt.inv = identity, i.e. t(sqrt.inv) %*% TT has identity covariance matrix
#            TT.std = t(sqrt.inv) %*% TT.0
#            res$p.value.2 = pchisq(sum(TT.std**2), df=p, lower.tail = FALSE)                
#        }
    } 
    
    class(res)=c("chngpt.score.test","htest",class(res))
    res
    
}

plot.chngpt.score.test <- function(x, ...) {
    fold=length(x$TT)/length(x$chngpts)
    for (i in 1:fold) {
        plot(x$chngpts, x$TT[1:length(x$chngpts)+(i-1)*length(x$chngpts)], xlab="change point", ylab="T", type="b", ...)
        #perc=as.numeric(strtrim(names(x$chngpts),nchar(names(x$chngpts))-1))  
    }
}

antoka.test=function(formula, data, chngpt.var, plot.=FALSE) {
    
    DNAME = deparse(substitute(data))
    
    dat.1=data[order(data[[chngpt.var]]),]
    
    tmp=model.frame(formula, dat.1)
    y=tmp[[1]]
    
    fit=keepWarnings(glm(formula,  data=dat.1, family="binomial"))    # if glm gives a warning, use sandwich estimator to get variance est
    if(length(fit$warning)!=0) {
        return (NA)
    } else {
        fit=fit$value
    }
    
    n = nrow(dat.1)
    mu.hat = expit(predict(fit))
    
    T. = sapply (1:(n-1), function(k) {
        S.k.0 = sum((y-mu.hat)[1:k])
        V.hat = mean(mu.hat*(1-mu.hat)) * k * (n-k) / n # this V.hat is a simplified version of the actual variance if the formula has more than an intercept as predictors
        abs(S.k.0)/sqrt(V.hat)
    })    
    if(plot.) plot(T., type="b")
    
    # compute p-value
    T.max=max(T.)
    names(T.max)="Maximum statisticsZ"
    loglogn=log(log(n))
    T1=sqrt(2*loglogn)*T.max - 2*loglogn - 1/2*log(loglogn) + 1/2*log(pi)
    p.value=exp(-2*exp(-T1))
    
    res=list()
    res$statistic=T.max
    res$p.value=p.value
    res$k=which.max(T.)
    res$parameter=NULL
    res$conf.int=NULL
    res$estimate=NULL
    res$null.value=NULL
    res$alternative="two-sided"
    res$method="Antoka Change Point Test"
    res$data.name=DNAME
        
    class(res)=c("htest",class(res))
    res
}
