# source ("http://sites.google.com/site/youyifong/misc/youtil.R")

# author: Youyi Fong, if not specified otherwise

# most useful functions:
#   mypostscript/mypdf
#   mytex
#   %+% # this operator add up two strings, e.g. "a" %+% "b"
#   getFormattedSummary

# easy to forget useful R functions
#   expand.grid


#########################################################
### string functions
#########################################################

# paste two strings together
# e.g. "a" %+% "b"
"%+%" <- function (a, b) {
    if (is.numeric(a) & is.numeric(b)) out=sum(a,b)
    else out=paste(a,b,sep="")
    out
}


# many code in the packages expect the default behavior, thus it is dangerous to do this
#"[" <- function (x, ...) {
#    base::"["(x, ...,drop=FALSE)
#}

firstIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
            break;
        }
    }
    ret
}

lastIndex =function (s1, s2) {
    k=nchar (s2)
    ret=-1;
    for (i in 1:(nchar(s1)-k+1) ) {
        if (substr(s1, i, i+k-1)==s2) {
            ret=i;
        }
    }
    ret
}

# return TRUE if s1 starts with s2? 
startsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, 1, nchar(s2)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}

# return TRUE if s1 ends with s2, s1 can be a vector
endsWith=function(s1, s2){
    sapply (s1, function (s) {
        if ( substring (s, nchar(s)-nchar(s2)+1, nchar(s)) == s2 ) {
            return (TRUE);
        } else {
            return (FALSE);
        }
    })
}

# return TRUE if s1 contains s2
contain =function (s1, s2) {
    sapply (s1, function (s) {
        k=nchar (s2)
        matched=FALSE
        for (i in 1:(nchar(s)-k+1) ) {
            if (substr(s, i, i+k-1)==s2) {
                matched=TRUE
            }
        }
        matched
    })
}

####################################################
# Matrix functions
####################################################

# serial covariance matrix
AR1 = function (p, w) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            m [i,j]=w**abs(i-j)
        }
    }
    m
}

# exchangeable covariance matrix
EXCH = function (p, rho) {
    m = matrix(1, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            if (i!=j) m [i,j]=rho
        }
    }
    m
}

# matrix trace
tr=function (m) {
    s=0
    for (i in 1:length(m[1,])) {
        s=s+m[i,i]  
    }
    s
}

getUpperRight = function (matri, func=NULL) {
    n=nrow (matri)
    out= numeric ( (n-1)*n/2 )
    index=0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            index=index+1
            out[index]=matri[i,j]
        }
    }
    if (is.null(func)) {
        out
    } else {
        func(out)
    }
}


#repeat a matrix in a block diagonal fashion
rep.matrix.block = function (.matrix, times=2) {
    orig.dim = nrow (.matrix)     
    m = matrix (0, orig.dim * times, orig.dim * times)
    for (i in 1: times) {
        m[(1+(i-1)*orig.dim):(i*orig.dim), (1+(i-1)*orig.dim):(i*orig.dim)]=.matrix
    }
    m
}


#it does not work on data.frame
rep.matrix = function (.matrix, times=1, each=1, by.row=TRUE) {
    if (by.row) {
        colnames.=colnames(.matrix)
        new.matrix=matrix(0, nrow(.matrix)*each*times, ncol(.matrix) )
        for (i in 1:nrow(.matrix)) {
            for (j in 1:each) {
                new.matrix[(i-1)*each + j,] = .matrix[i,]
            }
        }
        
        if(times>1) {
            for (i in 2:times) {
                new.matrix[ ((i-1)*nrow(.matrix)*each+1) : (i*nrow(.matrix)*each), ] = new.matrix[1:(nrow(.matrix)*each), ]
            }
        }

        dimnames(new.matrix)[[2]]=colnames.
        new.matrix
    }
    else {
        t ( rep.matrix(t(.matrix), times, each, by.row=TRUE) )
    }    
}

# rep.data.frame(chi[21,], 2)
rep.data.frame = function (.data.frame, times=1){
    out = .data.frame
    if (time==1) return (out)
    for (i in 2:times) {
        out = rbind (out, .data.frame)
    }
    out
}


##############################################################
# misc functions
##############################################################

#binomial coefficient
binom.coef=function(n,m) { prod((n-m+1):n)/prod(1:m) }

sign.test = function (diff) {
    binom.test( sum( diff > 0 , na.rm=TRUE) , sum( diff != 0 , na.rm=TRUE), p = 0.5,
           alternative = c("two.sided", "less", "greater"),
           conf.level = 0.95)   
}

# DONT CHANGE THIS, USED By RUMINEX!
# if ret.mat is set to true, always return a matrix
# in the output, each row corresponds to one element of X, instead of each column
mysapply=function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, ret.mat=TRUE) 
{
    if (is.null(names(X)) & is.numeric(X)) names(X)=X%+%""
    FUN <- match.fun(FUN)
    answer <- lapply(X, FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1) {
         if (common.len >= 1) 
            if (common.len == 1 & !ret.mat) 
                unlist(answer, recursive = FALSE)
            else 
                t(array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                    length(X)), dimnames = if (!(is.null(n1 <- names(answer[[1]])) & 
                    is.null(n2 <- names(answer)))) 
                    list(n1, n2)))
        else t(answer)
    }
    else t(answer)
}
## test mysapply
#sapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )
#mysapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )


# This function process all columns of x together instead of processing them one at a time
# FUN can return an array or a list. It does not have to return a scalar. This
#    saves from having to redo grouping for every col that has to be returned
#    also this eliminates the necessity to process a column of x at a time
myaggregate = function (x, by, FUN, new.col.name="aggregate.value", ...) 
{
    if (!is.data.frame(x)) 
        x <- as.data.frame(x)
    if (!is.list(by)) 
        stop(paste(sQuote("by"), "must be a list"))
    if (is.null(names(by))) 
        names(by) <- paste("Group", seq(along = by), sep = ".")
    else {
        nam <- names(by)
        ind <- which(nchar(nam) == 0)
        names(by)[ind] <- paste("Group", ind, sep = ".")
    }
    
    #original
    #y <- lapply(x, tapply, by, FUN, ..., simplify = FALSE)
    #modified
    z=mytapply(x, by, FUN, ...)
    
    #original
    #if (any(sapply(unlist(y, recursive = FALSE), length) > 1)) 
    #    stop(paste(sQuote("FUN"), "must always return a scalar"))
    #z <- y[[1]]
    
    d <- dim(z)
    w <- NULL
    for (i in seq(along = d)) {
        j <- rep.int(rep.int(seq(1:d[i]), prod(d[seq(length = i - 
            1)]) * rep.int(1, d[i])), prod(d[seq(from = i + 1, 
            length = length(d) - i)]))
        w <- cbind(w, dimnames(z)[[i]][j])
    }
    w <- w[which(!unlist(lapply(z, is.null))), , drop = FALSE]
        
    #original
    #y <- data.frame(w, lapply(y, unlist, use.names = FALSE))
    #modified
    if(nrow(w)!=nrow(matrix(unlist(z), nrow=nrow(w), byrow=TRUE))) stop("SOMETHING WRONG IN myaggregate")
    y <- data.frame(w, matrix(unlist(z), nrow=nrow(w), byrow=TRUE))
    #original
    #names(y) <- c(names(by), names(x))
    #modified
    names(y) <- c(names(by), new.col.name)
    y
}

# This function can handle X being matrix instead of just a vector
mytapply = function (X, INDEX, FUN = NULL, ..., simplify = TRUE) 
{
    FUN <- if (!is.null(FUN)) 
        match.fun(FUN)
    if (!is.list(INDEX)) 
        INDEX <- list(INDEX)
    nI <- length(INDEX)
    namelist <- vector("list", nI)
    names(namelist) <- names(INDEX)
    extent <- integer(nI)
    
    #original
    #nx <- length(X)
    #modified
    nx = ifelse(!(is.data.frame(X) | is.matrix(X)), length(X), length(X[,1]) )
    
    one <- as.integer(1)   
    group <- rep.int(one, nx)
    ngroup <- one
    for (i in seq(INDEX)) {
        index <- as.factor(INDEX[[i]])
        if (length(index) != nx) 
            stop("arguments must have same length")
        namelist[[i]] <- levels(index)
        extent[i] <- nlevels(index)
        group <- group + ngroup * (as.integer(index) - one)
        ngroup <- ngroup * nlevels(index)
    }
    if (is.null(FUN)) 
        return(group)
    ans <- lapply(split(X, group), FUN, ...)
    index <- as.numeric(names(ans))
    if (simplify && all(unlist(lapply(ans, length)) == 1)) {
        ansmat <- array(dim = extent, dimnames = namelist)
        ans <- unlist(ans, recursive = FALSE)
    }
    else {
        ansmat <- array(vector("list", prod(extent)), dim = extent, 
            dimnames = namelist)
    }
    names(ans) <- NULL
    ansmat[index] <- ans
    ansmat
}

expit=function (x) {1/(1+exp(-x))}
logit=function (x) {log(x/(1-x))}


# commented out b/c as.date has no global visible binding
## set the classes for columns of a dataframe
#SetColType = function (df, colClasses) {
#    colNum = length(names(df))
#    for (i in 1:colNum) {
#        if(colClasses[i]=="factor") 
#            df[,i]=as.factor(df[,i])
#        if(colClasses[i]=="Date") 
#            df[,i]=as.date(as.character(df[,i]))
#        if(colClasses[i]=="integer") 
#            df[,i]=as.integer(df[,i])
#        if(colClasses[i]=="numeric") 
#            df[,i]=as.numeric(df[,i])
#        if(colClasses[i]=="character") 
#            df[,i]=as.character(df[,i])
#        if(colClasses[i]=="logical") 
#            df[,i]=as.logical(df[,i])
#    }
#    df
#}
#
quick.t.test = function (x, y, var.equal=FALSE) {
    mean.x = mean(x)
    mean.y = mean(y)
    m=length(x)
    n = length(y)
    if (var.equal) {
        (mean.x-mean.y)/sqrt( ( sum((x-mean.x)**2) + sum((y-mean.y)**2) ) * (1/m + 1/n) /(m+n-2) )
    }    else {
        (mean.x-mean.y)/sqrt( sum((x-mean.x)**2)/(m-1)/m + sum((y-mean.y)**2)/(n-1)/n )         
    }
}

vector.t.test = function (mean.x, mean.y, var.x, var.y, n) {
    new.var = (var.x + var.y) /n
    t.stat = abs(mean.x-mean.y)/sqrt(new.var)
    names(t.stat)=NULL
    t.stat
}

tukey.mtest = function (mu, ms, n) {
    #mu = c(45, 58, 46, 45, 56 )
    #ms = 5
    #n=3
    m=length(mu)
    cutoff = qtukey(p=.01, nmeans=m, df=(n-1)*(m-1), nranges = 1, lower.tail = FALSE, log.p = FALSE)/sqrt(2)
    
    t=matrix(0, m,m)
    for (i in 1:m) {
        for (j in i:m) {
            t[i,j] =  ( (mu[i]-mu[j])/ sqrt( 2/n * ms ) )
        }
    }
    
    cat ("The t statistics for Tukey method are calculated below:\n")
    print (signif(t,3))
    cat ("\n")
    
    cat ("By comparing the t values with ", signif (cutoff,3), ", Tukey method declares that the following pairs are significantly different: ")
    for (i in 1:m) {
        for (j in i:m) {
            if (abs(t[i,j]) > cutoff) {
                if (t[i,j]>0) cat (i, "&", j)
                else cat (j, "&", i)
                if (i==m & j==n) cat (". ")    
                else cat (", ")   
            }
        }
    }
    cat ("In other words, the following pairs are not significantly different: ")
    for (i in 1:m) {
        for (j in i:m) {
            if (abs(t[i,j]) <= cutoff & i!=j) {
                cat (i, "&", j)
                if (i==m & j==n) cat (". ")    
                else cat (", ")   
            }
        }
    }
    cat("\n")

}

as.binary <- function(n, base=2 , r=FALSE)
{
   out <- NULL
   while(n > 0) {
     if(r) {
       out <- c(out , n%%base)
     } else {
       out <- c(n%%base , out)
     }   
     n <- n %/% base
   }
   return(out)
}


############################################################
# output functions
############################################################

# print a matrix/table or a list of them to a latex file as xtable
# note file.name can not have space in it
# e.g. mytex(matrix(0,2,2));
# e.g. mytex(matrix(0,2,2), digits=4);
# e.g. mytex(list(matrix(0,2,2), c(1,1))); 
# default arguments: file.name="temp"; digits=NULL; display=NULL; align="r"; append=FALSE; preamble=""; keep.row.names=TRUE
mytex.begin=function(file.name,preamble=""){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%+%".tex"
    cat ("\\documentclass{article}\n", file=file.name, append=FALSE)
    cat (preamble, file=file.name, append=TRUE)
    cat("\n\\usepackage{geometry}\n", file=file.name, append=TRUE)    
    cat("\n\\begin{document}\n", file=file.name, append=TRUE)    
}
mytex.end=function(file.name){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%+%".tex"
    cat ("\n\\end{document}", file=file.name, append=TRUE)
}
mytex=function(dat=NULL, file.name="temp", digits=NULL, display=NULL, align="r", append=FALSE, preamble="", keep.row.names=TRUE, floating=FALSE, lines=TRUE, ...) {
    
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    
    require(xtable)
    if (!endsWith(file.name,".tex")) file.name=file.name%+%".tex"
    
    if(is.data.frame(dat)) dat=list(dat)
    if (!is.list(dat)) dat=list(dat)
    if (!append) {
        mytex.begin(file.name, preamble)
    } 
#    else {
#        fil1 = file(file.name, "r")
#        # copy content to a new file until hit \end{document}, assuming that is the last line
#        n=as.numeric(strsplit(system ("wc -l "%+%file.name, intern=TRUE), " ")[[1]][1])
#        print(n)
#        tmp = readLines(fil1, n-1)
#        close(fil1)
#        cat (concatList(tmp, "\n"), file=file.name, append=FALSE)
#    }
    
    if (length(dat)>0) {
        names(dat)=gsub("_"," ",names(dat))
        for (i in 1:length(dat)) {
            cat ("\\leavevmode\\newline\\leavevmode\\newline\n", file=file.name, append=TRUE)
            dat1 = dat[[i]]        
            .ncol=ncol(dat1)
            if (is.null(.ncol)) {
                if (is.null(nrow(dat1))) .ncol=1
                else .ncol=nrow(dat1)
            }
            
            if (!is.matrix(dat1) & is.character(dat1)) {
                cat (dat1%+%"\n\n\n", file=file.name, append=TRUE)
            } else {        
                if (is.vector(dat1)) dat1=as.matrix(dat1)
                
                cat (names(dat)[i]%+%"\n\n", file=file.name, append=TRUE)
                if (!is.null(dat1)) {
                    if (!is.null(attr(dat1,"caption"))) caption=attr(dat1,"caption") else caption=NULL
                    
                    if (lines) hline.after=c(-1,0,nrow(dat1)) else hline.after=c()
                    if (length(align)==1) align=rep(align,.ncol+1)
                    if (!keep.row.names) rownames(dat1)=1:nrow(dat1) # there is no way to not print rownames by xtable, here we make it not so distracting
                    print(..., xtable(dat1, 
                        digits=(if(is.null(digits)) rep(3, .ncol+1) else digits), # cannot use ifelse here!!!
                        display=(if(is.null(display)) rep("f", .ncol+1) else display), # or here
                        align=align, caption=caption, ...), 
                        hline.after=hline.after,
                            type = "latex", file = file.name, append = TRUE, floating = floating )
                }
                cat ("\n", file=file.name, append=TRUE)
            }
        
        }
    }
    
    if(!append) mytex.end(file.name)
    print ("data saved to "%+%getwd()%+%"/"%+%file.name)
}
#x=matrix(0,2,2)
#attr(x,"caption")="cap"
#mytex(x, floating=TRUE)


# write a table that contains mean and sd to temp.tex in the current working directory, getwd()
# models can be a list of models, or a single model
make.latex.coef.table = function (models, model.names=NULL, row.major=FALSE, round.digits=NULL) {
# e.g.: models=list(gam1, gam2); round.digits= c(3,3,3,3,3); model.names=c("gam1", "gam2");  row.major=TRUE   
    if (! ("list" %in% class (models) ) ) {models=list(models)}
    
    numParams = nrow (getFixedEf(models[[1]]))
    numModels = length (models)
    
    if (is.null (model.names)) {model.names=rep("",numModels)}
    if (is.null(round.digits)) round.digits=rep(3,numParams)    
    
    coef.table = mysapply (1:numModels, function (i.model) {
        temp = getFixedEf(models[[i.model]]) [,1:2,drop=FALSE]
        for (i.param in 1:numParams) {
            temp[i.param,] = round (temp[i.param,], round.digits[i.param])
        }
        temp2 = paste (format(temp[,1]), "(", format(temp[,2]), ")")
        names (temp2) = dimnames(temp)[[1]]
        temp2
    })
    dimnames (coef.table)[[1]] = model.names
    
    if (row.major) mytex ( coef.table, align="r" ) 
    else mytex (t(coef.table), align="r") 
}

# convert a factor to integer using its value, e.g. 1 to 1, 2 to 2
ftoi = function (f) {
    as.integer (as.character (f) )
}

##returns day of year from a date
#day.of.year = function (date1) {
#   temp = date.mdy (date1)
#    date.new.year = as.date(paste("1/1/",as.character(temp$year)))
#    date1-date.new.year + 1
#}

between = function (x, a, b){
    x>=a & x<=b
}

# convert temp from f to c
f2c = function (f) {
    (f-32)*5/9
}

# like lag, move vector to the right/left by given number of steps
# x is a vector
shift.right = function (x, k=1) {
    c(rep(NA,k), x[1:(length(x)-k)])
}
shift.left = function (x, k=1) {
    c(x[(k+1):length(x)], rep(NA,k))
}

# return a subset of data that is 1 row every thin.factor rows
ThinRows = function (dat, thin.factor=10) {
    NumRows = nrow(dat)
    dat[1:(NumRows/thin.factor)*thin.factor,]
}
thin.rows=ThinRows

#mix two arrays in an interlacing way
mix = function (a, b) {
    if (length(a)!=length(b)) print ("Length of two arguments to mix function not equal.")
    out = rep (a, each=2)
    for (i in 1:length(a)) {
        out[2*i]=b[i]
    }
    out
}

# this is written before I know about matplot {graphics}
#note that the first column must be the time axis
#plot several series in one plot
#dat is a matrix, first column is x, second is y1, third is y2, ...
#it is hard to change color to blue and red
# e.g. plotSeries (cbind(1:10, 2:11, 3:12))
plotSeries = function(dat, main="", legend=NULL,ylim=NULL,col=NULL,lty=NULL,ylab=NULL,type="b"
    ,legend.cex=1,legend.x="topleft",legend.inset=0,legend.bty="n",pch=NULL,     ...) {
    if (is.null (ylim)) ylim=range(dat[,2:ncol(dat)], na.rm=TRUE)
    if (is.null (col)) col=1:(ncol(dat)-1)
    if (is.null (lty)) lty=rep(1,(ncol(dat)-1))
    if (is.null (ylab)) ylab=""
    
    plot (dat[,1],rep(0,nrow(dat)),type="n",ylim=ylim, ylab=ylab, main=main, ...)
    for (i in 2:ncol(dat)) {
        lines (dat[,1], dat[,i], type=type, col=col[i-1], lty=lty[i-1], pch=pch, ...)
    }
    
    if (!is.null(legend)) 
        if (is.null(pch)) 
            legend(x=legend.x, col=col, legend=legend, lty=lty, cex=legend.cex,inset=legend.inset,bty=legend.bty)
        else 
            legend(x=legend.x, col=col, legend=legend, pch=pch, cex=legend.cex,inset=legend.inset,bty=legend.bty, pt.cex=2)
}


#generating rejective sample
rejective.sampling = function (N, n, pik) {
    s = sample(N, n, replace = TRUE, prob = pik)
    while (length(unique(s)) != n) {
        s = sample(N, n, replace = TRUE, prob = pik)
    } 
    s   
}

#mysystem can call any exe file
mysystem = function (cmd, ...) {
    system ( paste(Sys.getenv("COMSPEC")," /c ",cmd) , invisible =TRUE, intern=FALSE, ...)
}

getModeFromLocfit=function (fit) {
    tmp=preplot(fit)
    tmp$xev[[1]][rank (tmp$fit) == length(tmp$fit)]
}

endWith = function (s1, c1) {
    substr(s1, nchar(s1), nchar(s1) )==c1
}

getMfrow=function (len) {
    ret=NULL
    if (len==1) { 
        ret=c(1,1)
    } else if (len==2) { 
        ret=c(1,2)
    } else if (len==3) { 
        ret=c(1,3)
    } else if (len<=4) { 
        ret=c(2,2)
    } else if (len<=6) { 
        ret=c(2,3)
    } else if (len<=9) { 
        ret=c(3,3)
    } else if (len<=12) { 
        ret=c(3,4)
    } else if (len<=16) { 
        ret=c(4,4)
    } else if (len<=20) { 
        ret=c(4,5)
    } else if (len<=25) { 
        ret=c(5,5)
    }
    ret
}

concatList = function (lis, sep=""){
    out=lis[[1]]
    i=2
    while (i<=length(lis)){
        out=out%+%sep%+%lis[[i]]
        i=i+1
    }
    out
}

escapeUnderline=function (name) {
    gsub("_", "\\_", name)
}

#Usage eg: sendEmail("ying@gmail.com", "re: simulation", "body 1")
sendEmail=function (to, subject, body) {
    if (!unix())        
        system ("C:/1Youyisstuff/3Software/C++/SendEmail/Release/SendEmail.exe "
            %+% to %+%" \""
            %+% subject %+%"\" \""
            %+% body %+% "\"",
        wait=FALSE, show.output.on.console=FALSE, intern=TRUE)
    else {
        cmd="echo '"%+%body%+%"' | /usr/bin/mail -s '"%+%subject%+%"' " %+% to
        system (cmd)
    }
}

# it is best to resize a figure in latex, because the proportions won't look right, if we change the wide and height
# is.prezo controls whether or not it is used for a presentation
mypostscript=function (file="temp", mfrow=c(1,1), mfcol=NULL, width=NULL, height=NULL, ext="eps", is.prezo=FALSE, oma=NULL, mar=NULL,main.outer=FALSE, save2file=TRUE, ...) {    
    
    print(paste(getwd(),"/",file,sep=""))
    
    if (!is.null(mfcol)) {
        nrow=mfcol[1]; ncol=mfcol[2]        
    } else {
        nrow=mfrow[1]; ncol=mfrow[2]
    }
    
    #if (nrow>4) warning ("nrow > 4 will not fit a page without making the figures hard to see")
    
    # sca controls how much to scale down for use in a paper
    if(is.null(width) | is.null(height))  {
        if (nrow==1 & ncol==1) {width=6.7; height=6.7
        } else if (nrow==1 & ncol==2) {width=9.7; height=5.2
        } else if (nrow==1 & ncol==3) {width=9.7; height=3.4
        } else if (nrow==1 & ncol==4) {width=14; height=3.4
        } else if (nrow==2 & ncol==3) {width=9.7; height=6.7
        } else if (nrow==4 & ncol==6) {width=15; height=10
        } else if (nrow==2 & ncol==4) {width=13; height=6.7
        } else if (nrow==3 & ncol==6) {width=17.5; height=9
        } else if (nrow==4 & ncol==8) {width=17.5; height=9
        } else if (nrow==4 & ncol==9) {width=20; height=9
        } else if (nrow==3 & ncol==5) {width=15; height=9.6
        } else if (nrow==3 & ncol==4) {width=12; height=9.6
        } else if (nrow==4 & ncol==5) {width=15; height=12.5
        } else if (nrow==5 & ncol==6) {width=9; height=8.3
        } else if (nrow==2 & ncol==2) {width=8; height=8.5
        } else if (nrow==3 & ncol==3) {width=9.7; height=10.3
        } else if (nrow==4 & ncol==4) {width=9.7; height=10.3
        } else if (nrow==6 & ncol==5) {width=18; height=17
        } else if (nrow==5 & ncol==5) {width=15; height=15
        } else if (nrow==5 & ncol==3) {width=9; height=15
        } else if (nrow==4 & ncol==2) {width=6; height=13
        } else if (nrow==6 & ncol==3) {width=9; height=19
        } else if (nrow==7 & ncol==3) {width=9; height=22
        } else if (nrow==8 & ncol==5) {width=10; height=16
        } else if (nrow==6 & ncol==4) {width=12; height=19
        } else if (nrow==7 & ncol==5) {width=18; height=19
        } else if (nrow==5 & ncol==4) {width=12; height=15
        } else if (nrow==2 & ncol==1) {width=6.7; height=9.7
        } else if (nrow==5 & ncol==1) {width=5; height=13
        } else if (nrow==3 & ncol==2) {width=6.7; height=10.3
        } else if (nrow==4 & ncol==3) {width=9; height=12
        } else stop ("nrow x ncol not supported: "%+%nrow%+%" x "%+%ncol)
    }
    
    if(save2file){      
        if (ext=="pdf") pdf (paper="special", file=file%+%"."%+%ext, width=width, height=height, ...)
        else postscript (paper="special", horizontal=FALSE, file=file%+%"."%+%ext, width=width, height=height, ...)
    } else {
        print("not saving to file")
    }
    
    if (!is.null(mfcol)) par(mfcol=mfcol)
    else par(mfrow=mfrow)    
    
    if (!is.null(oma)) par(oma=oma)
    if (!is.null(mar)) {
        par(mar=mar)
    }
    
    if (main.outer) {
        tmp=par()$oma
        tmp[3]=tmp[3]+1
        par(oma=tmp)
    }
    
}
mypdf=function (...) {mypostscript(ext="pdf",...)}
##test
#mypdf(mfrow=c(1,1),file="test1x1");plot(1:10,main="LUMX",xlab="t",ylab="y");dev.off()
#mypdf(mfrow=c(1,2),file="test1x2");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);dev.off()
#mypdf(mfrow=c(2,2),file="test2x2");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(1,3),file="test1x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(2,3),file="test2x3");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);dev.off()
#mypdf(mfrow=c(4,4),file="test4x4");plot(1:10,main="LUMX",xlab="t",ylab="y");plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10);plot(1:10,main="Luminex");dev.off()

# default values for the control parameters have been tuned so that the eps figure is suitable for a latex file
# example: mypostscript(); plot(1:10); dev.off()
mypostscript.old=function (file="temp", mfrow=c(1,1), mfcol=NULL, width=NULL, height=NULL, my.oma=c(0,0,0,0),
    my.mar=c(3,3,2.5,1.5), type="eps", ...) {    
#    if(exists("figurePath") && file.exists(figurePath)) {
#        file=figurePath%+%"/"%+%file
#    } else {
#        file=file
#    }
    
    print(paste(getwd(),"/",file,sep=""))
    
    if (!is.null(mfcol)) {
        nrow=mfcol[1]; ncol=mfcol[2]        
    } else {
        nrow=mfrow[1]; ncol=mfrow[2]
    }
    
    if(is.null(width))  {
        if (ncol>=3) width=6.5
        else if (ncol==2) width=5
        else width=4
    }
    if(is.null(height)) {
        if (ncol>=4) height=width/ncol*nrow * 1
        else height=width/ncol*nrow       
    }
    
    if (type=="pdf") pdf (paper="special", file=file%+%"."%+%type, width=width, height=height, ...)
    else postscript (paper="special", horizontal=FALSE, file=file%+%"."%+%type, width=width, height=height, ...)
    
    if (!is.null(mfcol)) par(mfcol=mfcol)
    else par(mfrow=mfrow)    
    
    my.mgp=c(1.6,.5,0); 
    if (ncol>=4)      par (oma=my.oma, mar=my.mar, mgp=my.mgp, cex=.5, cex.main=1, cex.axis=1, cex.lab=1) 
    else if (ncol==3) par (oma=my.oma, mar=my.mar, mgp=my.mgp, cex=.75, cex.main=.75, cex.axis=.75, cex.lab=.85) 
    else if (ncol==2) par (oma=my.oma, mar=my.mar, mgp=my.mgp, cex=.75, cex.main=.75, cex.axis=.75, cex.lab=.85)
    else              par (oma=my.oma, mar=my.mar, mgp=my.mgp, cex=1, cex.main=.75, cex.axis=.75, cex.lab=.85)
}

normalize=function (a) {
    a/sum(a)
}

# returns binary representation of an integer
binary<-function(i) if (i) paste(binary(i %/% 2), i %% 2, sep="") else "" 

# returns binary representatin of an integer with leading 0, the length of string is n
binary2<-function(i, n) {
    a<-2^((n-1):0)
    b<-2*a
    sapply(i,function(x) paste(as.integer((x %% b)>=a),collapse=""))
} 

##-- Stirling numbers of the 2nd kind
##-- (Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)

##> S^{(m)}_n = number of ways of partitioning a set of $n$ elements into $m$
##> non-empty subsets

Stirling2 <- function(n,m)
{
    ## Purpose:  Stirling Numbers of the 2-nd kind
    ##      S^{(m)}_n = number of ways of partitioning a set of
    ##                      $n$ elements into $m$ non-empty subsets
    ## Author: Martin Maechler, Date:  May 28 1992, 23:42
    ## ----------------------------------------------------------------
    ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
    ## Closed Form : p.824 "C."
    ## ----------------------------------------------------------------
    ## maechler@_stat.math.ethz.ch

    if (0 > m || m > n) stop("'m' must be in 0..n !")     
    k <- 0:m
    sig <- rep(c(1,-1)*(-1)^m, length= m+1)
    # 1 for m=0; -1 1 (m=1)     
    ## The following gives rounding errors for (25,5) :     
    ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )     
    ga <- gamma(k+1)
    round(sum( sig * k^n /(ga * rev(ga)))) 
}

myprint <- function(object, ...) UseMethod("myprint") 

setSeedAsInR=function() {
# choose a rng so that results can be compared with that from C program
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = NULL)
    set.seed(1, kind = "Marsaglia-Multicarry")
    myprint(runif(1))
    myprint("the number above should be 0.006153224")
}

lotto=function(max, min=1) {
    cat ("    ")
    cat (floor(runif(1, min, max+1)))
    cat("    : a random number from "%+%min%+%" to "%+%max)
    cat ("\n")
}

# the input is a list, typically the output from a sapply call that should be matrix, but have different length
fill.jagged.array=function(a) {
    # don't check is.matrix, because for some reason, it will return true
    max.len=max(sapply(a, length))
    mysapply(a, function (e) {
        c(e, rep(NA, max.len-length(e)))
    })    
}

# don't have to transpose x
mywrite=function(x, ...){
    if (is.list(x)) x=fill.jagged.array(x)
    if (is.null(ncol(x))) i=length(x)
    else i=ncol(x)
    write (t(x), ncolumns=i, ...)
}

empty.plot=function () {
    plot(1,1,type="n",xlab="",ylab="",xaxt="n", yaxt="n", bty="n")
}

logSumExp=function (logx){
    logMeanExp(logx, 1)
}
logSumExpFor2=function (logx, logy){
    c=max(logx, logy)
    dif=abs(logx-logy)
    if (dif>300) return (c)
    else {
        log(sum(exp(logx-c), exp(logy-c)))+c
    }
}
# log( exp(logx1)-exp(logx2) )
logDiffExp=function (logx1,logx2){
    c=logx1
    c+ log(1-exp(logx2-logx1))
}

logMeanExp=function (logx,B=NULL){
# mean function for small numbers
# logx is a vector of large negative values
# return log (sum(exp(logx))/B)
# calculate log of the mean of a vector which contains 0 and very small real numbers (logged)
# return log of the mean
    if (is.null(B)) B=length(logx)
    c=max(logx)
    log(sum(exp(logx-c))/B)+c
}
# logMeanExp(log(1:5), 5) # test, should return log(3)

logDiffExp=function (logx1, logx2){
# diff function for small numbers
# logx1 and logx2 are typically large negative values, logx1>logx2
# return log (exp(logx1)-exp(logx2))
    if (logx1<logx2) {cat("\nlWarning [logDiffExp]: first argument smaller than second, return NaN.\n\n"); return (NaN);}
    c=max(logx1, logx2)
    log (exp(logx1-c)-exp(logx2-c))+c
}
# logDiffExp(log(2), log(1)) # test, should return 0

# from combinat package
permn=function (x, fun = NULL, ...) 
{
    if (is.numeric(x) && length(x) == 1 && x > 0 && trunc(x) == 
        x) 
        x <- seq(x)
    n <- length(x)
    nofun <- is.null(fun)
    out <- vector("list", gamma(n + 1))
    p <- ip <- seqn <- 1:n
    d <- rep(-1, n)
    d[1] <- 0
    m <- n + 1
    p <- c(m, p, m)
    i <- 1
    use <- -c(1, n + 2)
    while (m != 1) {
        out[[i]] <- if (nofun) 
            x[p[use]]
        else fun(x[p[use]], ...)
        i <- i + 1
        m <- n
        chk <- (p[ip + d + 1] > seqn)
        m <- max(seqn[!chk])
        if (m < n) 
            d[(m + 1):n] <- -d[(m + 1):n]
        index1 <- ip[m] + 1
        index2 <- p[index1] <- p[index1 + d[m]]
        p[index1 + d[m]] <- m
        tmp <- ip[index2]
        ip[index2] <- ip[m]
        ip[m] <- tmp
    }
    out
}

last = function (x, n=1, ...) {
    if (is.character(x)) tail (readLines(x), n=n, ...) # read file
    else if (is.vector(x)) x[length(x)]
    else if (is.array(x)) x[length(x)]
    else if (is.list(x)) x[[length(x)]]
    else stop ("last(): x not supported")
}

# dat is all positive and we want density to have positive support as well
get.density.boundary.corrected=function(dat) {
    # first log transform dat to real line
    dat1=log(dat)
    tmp=density(dat1)
}


unix=function (){
    substr(Sys.getenv("R_HOME") , 1,1)=="/"
}


#############################################################################
# regression functions or functions written for 570, 571
#############################################################################

#mycoef <- function(object, ...) UseMethod("mycoef") 
#
#"mycoef.drc" <-
#function(object, ...)
#{
#    if (!is.null(object$"coefficients"))
#    {
#        out=object$"coefficients"
#        names(out)=substr(names(out), 1, 1)
#        return(out)
#    } else {
#        retVec <- object$fit$par
#        names(retVec) <- object$parNames[[2]]
#        return(retVec)
#    }
#}
#
#
# return a matrix, first column coef, second column se, 
getFixedEf <- function(object, ...) UseMethod("getFixedEf") 

getFixedEf.MIresult=function(mir) {
    cbind(coef(mir), sqrt(diag(vcov(mir))))
}

getFixedEf.ltm=function (fit) {
    fit[-1,1:2]
}

#add CI to summary
mysummary <- function(object, ...) UseMethod("mysummary") 

#add prediction interval to summary
mypredict <- function(object, ...) UseMethod("mypredict") 

# returns sandwich estimator of variance matrix
# from Thomas Lumley
infjack.glm<-function(glm.obj,groups){
    umat<-estfun.glm(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    modelv<-summary(glm.obj)$cov.unscaled
    modelv%*%(t(usum)%*%usum)%*%modelv
}

# from Thomas Lumley
jack.glm<-function(glm.obj,groups){
    umat<-jackvalues(glm.obj)
    usum<-rowsum(umat,groups,reorder=FALSE)
    t(usum)%*%usum*(nrow(umat)-1)/nrow(umat)
}

# from Thomas Lumley
jackvalues<-function(glm.obj){
    db<-lm.influence(glm.obj)$coef
    t(t(db)-apply(db,2,mean))
}   

# from Thomas Lumley
estfun.glm<-function(glm.obj){
    if (is.matrix(glm.obj$x)) 
        xmat<-glm.obj$x
    else {
        mf<-model.frame(glm.obj)
        xmat<-model.matrix(terms(glm.obj),mf)       
    }
    residuals(glm.obj,"working")*glm.obj$weights*xmat
}

coef.geese = function  (geese1, ...) {
    tmp = summary(geese1)$mean[,1]
    names (tmp)=names (geese1$beta)
    tmp
}
vcov.geese = function  (geese1, ...) {
    tmp = geese1$vbeta
    dimnames (tmp)=list (names (geese1$beta), names (geese1$beta))
    tmp
}
residuals.geese = function (geese1, y, x) {
    y - x %*% geese1$beta   
}
predict.geese = function (geese1, x) {
    x %*% geese1$beta 
}

##get estimates, variances, sd from lmer fit
#getFixedEf.mer = function (lmerfit) {
#    Vcov <- lme4::vcov(lmerfit, useScale = FALSE) 
#    betas <- lme4::fixef(lmerfit) 
#    se <- sqrt(diag(Vcov)) 
#    zval <- betas / se 
#    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
#    cbind("Estimate"=betas, se, zval, pval) 
#}

#get estimates, variances, sd from lme fit
getFixedEf.lme = function (lme1) {
    betas <- lme1$coef$fixed
    se <- sqrt (diag (lme1$varFix))
    zval <- betas / se 
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE) 
    cbind(betas, se, zval, pval) 
}

getFixedEf.geese = function (geese1) {
    summary(geese1)$mean
}
    
getFixedEf.glm = function (glm1, exp=FALSE) {
    out=summary(glm1)$coef
    if(exp) out[,1]=exp(out[,1])
    out
}

getFixedEf.logistf = function (fit, exp=FALSE) {
    temp = summary(fit)
    out = cbind (coef=temp$coef, se=NA, p.value=temp$prob)
    if(exp) out[,1]=exp(out[,1])
    out
}
vcov.logistf=function(fit){
    fit$var
}

getFixedEf.gam = function (gam1) {
    temp = summary(gam1)
    cbind (temp$p.coef, temp$se[1:length(temp$p.coef)])
}

getFixedEf.lm = function (lm1) {
    summary(lm1)$coef
}

getFixedEf.gee = function (gee1) {
    summary(gee1)$coef
}

getFixedEf.inla = function (inlafit) {
    tmp = summary(inlafit)$fixed
    n=nrow(tmp)
    tmp.name = row.names(tmp)[n]
    # move intercept to the top
    if (tmp.name=="intercept") {
        tmp = rbind (tmp[n,],tmp)[1:n,,drop=FALSE]
        dimnames (tmp)[[1]][1] = tmp.name
    }
    # rename first column
    dimnames (tmp)[[2]][1] = "Estimate"
    tmp
}

getFixedEf.coxph=function (fit){
    sum.fit<-summary(fit)
    sum.fit$coef[,c(1,3)]
#    round(sqrt(diag(attr(fit$var,"phases")$phase1)),3)
#    round(sqrt(diag(attr(fit$var,"phases")$phase2)),3)    
}

# used to get mean and sd from a jags or winbugs sample 
# each column of samples is a variable
getFixedEf.matrix=function (samples) {
    t(apply(samples, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x))))
}

# add CI
getFixedEf2 = function (object, ...) {
    temp=getFixedEf (object, ...)
    temp = cbind(temp, "lower bound"=temp[,1]-1.96*temp[,2])
    temp = cbind(temp, "upper bound"=temp[,1]+1.96*temp[,2])
    temp
}

mysummary.lm = function (object, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
    z <- object
    p <- z$rank
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/(n - p)
        ans <- z[c("call", "terms")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0, 6)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate", 
            "Std. Error", "lower CI", "higher CI", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no terms nor qr component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if (rdf != z$df.residual) 
        warning("inconsistent residual degrees of freedom. -- please report!")
    p1 <- 1:p
    r <- z$residuals
    f <- z$fitted
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    #youyi
    ci.l = est - qt(.975, rdf, lower.tail=TRUE)*se
    ci.r = est + qt(.975, rdf, lower.tail=TRUE)*se
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval), 
        rdf, lower.tail = FALSE), ci.l, ci.r)
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "CI left", "CI right"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept")) 
            1
        else 0
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, 
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.lm"
    ans
}

#mysummary.ols=function (x, digits = 4, long = FALSE, ...) 
#{
#    oldopt <- options(digits = digits)
#    on.exit(options(oldopt))
#    cat("\n")
#    cat("Linear Regression Model\n\n")
#    dput(x$call)
#    cat("\n")
#    if (!is.null(z <- x$na.action)) 
#        naprint(z)
#    stats <- x$stats
#    if (lst <- length(stats)) {
#        if (.R.) {
#            cstats <- character(lst)
#            names(cstats) <- names(stats)
#            for (i in 1:lst) cstats[i] <- format(stats[i])
#            print(cstats, quote = FALSE)
#        }
#        else print(x$stats)
#        cat("\n")
#    }
#    pen <- length(x$penalty.matrix) > 0
#    resid <- x$residuals
#    n <- length(resid)
#    p <- length(x$coef) - (names(x$coef)[1] == "Intercept")
#    if (length(x$stats) == 0) 
#        cat("n=", n, "   p=", p, "\n\n", sep = "")
#    ndf <- x$stats["d.f."]
#    df <- c(ndf, n - ndf - 1, ndf)
#    r2 <- x$stats["R2"]
#    sigma <- x$stats["Sigma"]
#    rdf <- df[2]
#    if (rdf > 5) {
#        cat("Residuals:\n")
#        if (length(dim(resid)) == 2) {
#            rq <- apply(t(resid), 1, quantile)
#            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", 
#                "Max"), dimnames(resid)[[2]])
#        }
#        else {
#            rq <- quantile(resid)
#            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
#        }
#        print(rq, digits = digits, ...)
#    }
#    else if (rdf > 0) {
#        cat("Residuals:\n")
#        print(resid, digits = digits, ...)
#    }
#    if (nsingular <- df[3] - df[1]) 
#        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", 
#            sep = "")
#    else cat("\nCoefficients:\n")
#    se <- sqrt(diag(x$var))
#    z <- x$coefficients/se
#    P <- 2 * (1 - pt(abs(z), rdf))
#    ci.l = x$coefficients - qt(.975, rdf, lower.tail=TRUE)*se
#    ci.r = x$coefficients + qt(.975, rdf, lower.tail=TRUE)*se
#    co <- cbind(x$coefficients, se, z, P, ci.l, ci.r)
#    dimnames(co) <- list(names(x$coefficients), c("Value", "Std. Error", 
#        "t", "Pr(>|t|)", "CI left", "CI right"))
#    print(co)
#    if (pen) 
#        cat("\n")
#    else cat("\nResidual standard error:", format(signif(sigma, 
#        digits)), "on", rdf, "degrees of freedom\n")
#    rsqa <- 1 - (1 - r2) * (n - 1)/rdf
#    if (length(x$stats) == 0) 
#        cat("Multiple R-Squared:", format(signif(r2, digits)), 
#            " ")
#    cat("Adjusted R-Squared:", format(signif(rsqa, digits)), 
#        "\n")
#    if (!pen) {
#        if (long && p > 0) {
#            correl <- diag(1/se) %*% x$var %*% diag(1/se)
#            dimnames(correl) <- dimnames(x$var)
#            cat("\nCorrelation of Coefficients:\n")
#            ll <- lower.tri(correl)
#            correl[ll] <- format(round(correl[ll], digits), ...)
#            correl[!ll] <- ""
#            print(correl[-1, -(p + 1), drop = FALSE], quote = FALSE, 
#                digits = digits, ...)
#        }
#    }
#    cat("\n")
#    invisible()
#}
#
#summary.ols = function (ols1) {
#    print(ols1)
#}

mysummary.glm = function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion)) 
        dispersion <- if (any(object$family$family == c("poisson", 
            "binomial"))) 
            1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0)) 
                warning("observations with zero weight ", "not used for calculating dispersion")
            sum(object$weights * object$residuals^2)/df.r
        }
        else Inf
    p <- object$rank
    if (p > 0) {
        p1 <- 1:p
        Qr <- object$qr
        aliased <- is.na(coef(object))
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf   <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        
        #youyi
        ci.l = coef.p - qt(.975, df.r, lower.tail=TRUE)*s.err
        ci.r = coef.p + qt(.975, df.r, lower.tail=TRUE)*s.err
    
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue, ci.l, ci.r)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "z value", "Pr(>|z|)", "lower CI", "upper CI"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue, ci.l, ci.r)
            dimnames(coef.table) <- list(names(coef.p), c(dn, 
                "t value", "Pr(>|t|)", "lower CI", "upper CI"))
        }
        else {
            coef.table <- cbind(coef.p, Inf)
            dimnames(coef.table) <- list(names(coef.p), dn)
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0, 4)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0, 0)
        aliased <- is.na(coef(object))
        df.f <- length(aliased)
    }
    ans <- c(object[c("call", "terms", "family", "deviance", 
        "aic", "contrasts", "df.residual", "null.deviance", "df.null", 
        "iter")], list(deviance.resid = residuals(object, type = "deviance"), 
        coefficients = coef.table, aliased = aliased, dispersion = dispersion, 
        df = c(object$rank, df.r, df.f), cov.unscaled = covmat.unscaled, 
        cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
}

mypredict.glm = function (object, newdata = NULL, type = c("link", "response", 
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL, 
    na.action = na.pass, ...) 
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors, 
                response = object$fitted, terms = predict.lm(object, 
                  se.fit = se.fit, scale = 1, type = "terms", 
                  terms = terms))
            if (!is.null(na.act)) 
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                type = ifelse(type == "link", "response", type), 
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg")) 
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0) 
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict.lm(object, newdata, se.fit, scale = residual.scale, 
            type = ifelse(type == "link", "response", type), 
            terms = terms, na.action = na.action, interval="confidence")#youyi adds interval here
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}

mypredict.lm=function (object, newdata, se.fit = FALSE, scale = NULL, df = Inf, 
    interval = c("none", "confidence", "prediction"), level = 0.95, 
    type = c("response", "terms"), terms = NULL, na.action = na.pass, 
    ...) 
{
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    }
    else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, 
            xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
            eval(attr(tt, "variables")[[off.num + 1]], newdata)
        else if (!is.null(object$offset)) 
            eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq(len = p)
    piv <- object$qr$pivot[p1]
    if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) 
        warning("prediction from a rank-deficient fit may be misleading")
    beta <- object$coefficients
    predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
    if (!is.null(offset)) 
        predictor <- predictor + offset
    interval <- match.arg(interval)
    type <- match.arg(type)
    if (se.fit || interval != "none") {
        res.var <- if (is.null(scale)) {
            r <- object$residuals
            w <- object$weights
            rss <- sum(if (is.null(w)) 
                r^2
            else r^2 * w)
            df <- n - p
            rss/df
        }
        else scale^2
        if (type != "terms") {
            if (p > 0) {
                XRinv <- if (missing(newdata) && is.null(w)) 
                  qr.Q(object$qr)[, p1, drop = FALSE]
                else X[, piv] %*% qr.solve(qr.R(object$qr)[p1, 
                  p1])
                ip <- drop(XRinv^2 %*% rep(res.var, p))
            }
            else ip <- rep(0, n)
        }
    }
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        if (attr(tt, "intercept") > 0) 
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        hasintercept <- attr(tt, "intercept") > 0
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if (!mmDone) {
                mm <- model.matrix(object)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx[piv] * beta[piv])
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (se.fit || interval != "none") {
                ip <- matrix(ncol = nterms, nrow = NROW(X))
                dimnames(ip) <- list(rownames(X), names(asgn))
                Rinv <- qr.solve(qr.R(object$qr)[p1, p1])
            }
            if (hasintercept) 
                X <- sweep(X, 2, avx)
            unpiv <- rep.int(0, NCOL(X))
            unpiv[piv] <- p1
            for (i in seq(1, nterms, length = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0] <- 0
                predictor[, i] <- if (any(iipiv) > 0) 
                  X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
                if (se.fit || interval != "none") 
                  ip[, i] <- if (any(iipiv) > 0) 
                    as.matrix(X[, iipiv, drop = FALSE] %*% Rinv[ii, 
                      , drop = FALSE])^2 %*% rep.int(res.var, 
                      p)
                  else 0
            }
            if (!is.null(terms)) {
                predictor <- predictor[, terms, drop = FALSE]
                if (se.fit) 
                  ip <- ip[, terms, drop = FALSE]
            }
        }
        else {
            predictor <- ip <- matrix(0, n, 0)
        }
        attr(predictor, "constant") <- if (hasintercept) 
            termsconst
        else 0
    }
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip), 
            prediction = sqrt(ip + res.var))
        if (type != "terms") {
            predictor <- cbind(predictor, predictor + hwid %o% 
                c(1, -1))
            colnames(predictor) <- c("fit", "lwr", "upr")
        }
        else {
            lwr <- predictor + hwid
            upr <- predictor - hwid
        }
    }
    if (se.fit || interval != "none") 
        se <- sqrt(ip)
    if (missing(newdata) && !is.null(na.act <- object$na.action)) {
        predictor <- napredict(na.act, predictor)
        if (se.fit) 
            se <- napredict(na.act, se)
    }
    if (type == "terms" && interval != "none") {
        if (missing(newdata) && !is.null(na.act)) {
            lwr <- napredict(na.act, lwr)
            upr <- napredict(na.act, upr)
        }
        list(fit = predictor, se.fit = se, lwr = lwr, upr = upr, 
            df = df, residual.scale = sqrt(res.var))
    }
    else if (se.fit) 
        list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var))
    else predictor
}

# predict.gam returns a list of fit and se.fit. This function transforms it to a matrix: 1st col is fit, 2nd col is LI, 3rd col is UI, 4th col is se.fit
mypredict.gam = function (...) {
    pred = predict(...)
    x=matrix (0, length (pred$fit), 4)
    x[,1]=pred$fit
    x[,2]=pred$fit - 2 * pred$se.fit
    x[,3]=pred$fit + 2 * pred$se.fit
    x[,4]=pred$se.fit
    dimnames (x) = list (NULL, list("prediction", "LI", "UI", "se"))
    x
}

getMidPoints=function(x){
    ((c(0,x)+c(x,0))/2)[2:length(x)] 
}

# from a list of fits, say lmer, inla fits, return formatted summary controlled by "type"
# for a matrix, return Monte Carlo variance

# random=TRUE returns variance components
# type=1: est
# type=2: est (se)
# type=3: est (2.5%, 97.5%)
# type=4: est   se
getFormattedSummary=function(fits, type=1, est.digits=2, se.digits=2, random=FALSE){
    
    # should not use mysapply here 
    t(mysapply(fits, function (fit) {
        if (random) {
            tmp = getVarComponent (fit)
            if (class(fit)=="mer" & type!=1) {
                warning ("only point estimate is available for variance components from lmer fit, forcing type to 1")
                type=1
            }
        } else {
            tmp = getFixedEf (fit)
        }
        
        if (type==1)
            # this is the best way to format: first round, then nsmall
            format(round(tmp[,1], est.digits), nsmall=est.digits, scientific=FALSE) 
        else if (type==2)
            format(tmp[,1], digits=1, nsmall=est.digits, scientific=FALSE) %+% " (" %+% 
                format(tmp[,2], digits=1, nsmall=se.digits, scientific=FALSE) %+% ")"
        else if (type==3)
            format(tmp[,1], digits=1, nsmall=est.digits, scientific=FALSE) %+% " (" %+% 
                format(tmp[,3], digits=1, nsmall=est.digits, scientific=FALSE) %+% ", " %+% 
                    format(tmp[,4], digits=1, nsmall=est.digits, scientific=FALSE) %+% ")" 
        else if (type==4)
            # a space is inserted between est and se, they could be post-processed in swp
            format(tmp[,1], digits=1, nsmall=est.digits, scientific=FALSE) %+% " " %+% 
                format(tmp[,2], digits=1, nsmall=se.digits, scientific=FALSE)
        else 
            stop ("getFormattedSummaries(). type not supported: "%+%type)
    }))
}

getVarComponent <- function(object, ...) UseMethod("getVarComponent")

#getVarComponent.mer = function (lmer1) {
#    tmp=lme4::VarCorr(lmer1)
#    mysapply(tmp, function (comp) attr(comp, "stddev") )
#}

# comment out b/c building other packages, we will need to include dependency on some other packages that define VarCorr
#getVarComponent.lme = function (lme.fit) {
#    VarCorr(lme.fit)
#}


# used to get mean and sd from a jags or winbugs sample, getVarComponent.matrix and getFixedEf.matrix do the same thing
# each column of samples is a variable
getVarComponent.matrix = function (samples) {
    t(apply(samples, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}

getFixedEf.matrix = function (samples) {
    t(apply(samples, 2, function (x) c("Estimate"=mean(x), "sd"=sd(x), "2.5%"=quantile(x,.025), "97.5"=quantile(x,.975))))
}

# returns estimate of standard deviation and the estimated sd of that estimate
getVarComponent.hyperpar.inla = function (hyper1, transformation=NULL) {
    marginals = hyper1$marginals
    out = mysapply(1:length(marginals), function (i) {  
        # this is a little precarious, but hey
        if (startsWith(names(marginals)[i],"Prec")) {
            if (is.null (transformation)) {      
                getMeanSd(marginals[[i]],"inversesqrt")
            } else {
                getMeanSd(marginals[[i]],transformation)
            }
        } else if (startsWith(names(marginals)[i],"Rho")) {
            hyper1$summary[i, c(1,2,3,5)]
        } else {
            stop ("don't know what to do with this names(marginals)[i]: "%+% names(marginals)[i] )
        }
    })
    dimnames (out)[[1]]="sigma.or.rho."%+%dimnames (out)[[1]]
    out
}


###########################################
# for use with INLA
###########################################
# dat is all positive and we want density to have positive support as well
# remember to transform the density!
get.density.boundary.corrected=function(dat) {
    # first log transform dat to real line
    if (min(dat)<0) return (density(dat))
    else {
        dat1=log(dat)
        tmp=density(dat1)
        list (x=exp(tmp$x), y=tmp$y/exp(tmp$x))
    }
}

# return the mean, sd, CI of the transformed variable
getMeanSd=function (marginal, f="identity") {
    
    # interpolations suggested by Havard: do it on the original scale
    logtao=log(marginal[,1]); p.logtao=marginal[,2]*marginal[,1]
    fun = splinefun(logtao, log(p.logtao)) 
    h=0.001
    x = seq(min(logtao),max(logtao),by=h) 
    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)/2,max(logtao)+sd(logtao)/2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao),max(logtao)+sd(logtao),by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
#    x = seq(min(logtao)-sd(logtao)*2,max(logtao)+sd(logtao)*2,by=h) 
#    pmf = exp(fun(x))*h
#    sum (pmf) # Prob
    
    lower.boundary = rle(cumsum(pmf)>.025)$lengths[1]+1
    upper.boundary = rle(cumsum(pmf)>.975)$lengths[1]+1
#    if (pmf[lower.boundary]>.04) {
#        #stop ("getMeanSd(): pmf too large at lower boundary: "%+%pmf[lower.boundary])
#        return (rep(NA, 4))
#    }
#    if (pmf[upper.boundary]>.04) {
#        stop ("getMeanSd(): pmf too large at upper boundary"%+%pmf[upper.boundary])
#        return (rep(NA, 4))
#    }
    
    if (f=="identity") {
        func=function(x) { exp(x) }
    } else if (f=="inverse") {
        func=function(x) { exp(x)**-1 }
    } else if (f=="inversesqrt") {
        func=function(x) { exp(x)**-.5 }
    } else 
        stop ("getMeanSd(): function not supported "%+%f)
    
    mu = sum( pmf * func(x))
    stdev = (sum( pmf * func(x)**2 ) - mu**2) ** .5
    out = c("mean"=mu, "stddev"=stdev, 0,0 ) # we may run into problems with the boundary
    #out = c("mean"=mu, "stddev"=stdev, sort(func(x)[c(lower.boundary, upper.boundary)]) ); 
    names(out)[3]="2.5%"; names(out)[4]="97.5%";
    out
}

get.inla.f.param=function(arg) {
    if (arg$distr=="Gamma")
        inla.f.param=c(arg$shape, arg$rate)
    else if (arg$distr=="2DWishart")
        inla.f.param=c(arg$r, arg$R[1,1], arg$R[2,2], arg$R[1,2])
    else if (arg$distr=="3DWishart")
        inla.f.param=c(arg$r, arg$R[1,1], arg$R[2,2], arg$R[3,3], arg$R[1,2], arg$R[1,3], arg$R[2,3])
    else 
        stop ("get.inla.f.param(). distr not supported: "%+%arg$distr) 
    inla.f.param
}

getSamplesFileName=function (data.model,seed,chain,prior){
    seedstr=ifelse(is.null(seed) | is.na(seed),"","_seed"%+%seed)
    data.model%+%"/"%+%data.model%+%"_prior"%+%prior%+%"_chain"%+%chain%+%seedstr%+%".jags"
}

# plot marginal distribution from jags and inla
plot.marginals=function (.data.model, prior, seed=NA) {
    # first combine samples
    samples = NULL
    for (i in 1:.data.model@gibbs.chains) samples=rbind(.data.model@gibbs.samples[[i]])
    
    seedstr=""
    if (.data.model@is.sim) seedstr="_seed"%+%seed
    mypostscript(file=.data.model@name%+%"_prior"%+%prior%+%seedstr%+%".eps" )
    for (i in 1:length(.data.model@gibbs.watch.list)) {
        # this is a little precarious, but hey
        if (startsWith(.data.model@gibbs.watch.list[i],"sigma")) {
            jags.marginal = get.density.boundary.corrected ( 1/samples[,i]**2 ) 
        } else if (startsWith(.data.model@gibbs.watch.list[i],"rho")) {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else if (startsWith(.data.model@gibbs.watch.list[i],"log")) {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else if (.data.model@gibbs.watch.list[i]=="z") {
            jags.marginal = get.density.boundary.corrected ( samples[,i] ) 
        } else {
            next; # only do variance component parameters
        }
        
        hyper.marginal=.data.model@hyper.fit$marginals[[i]]
        # sometimes it is necessary to do transformation
        if (startsWith(.data.model@gibbs.watch.list[i],"log")) {
            hyper.marginal = cbind(log(hyper.marginal[,1]), hyper.marginal[,2]*hyper.marginal[,1])
        } else if (.data.model@gibbs.watch.list[i]=="z") {
            hyper.marginal = cbind(logit(.5+hyper.marginal[,1]/2), hyper.marginal[,2]*(.5-hyper.marginal[,1]**2/2))
        } else {
            hyper.marginal = hyper.marginal 
        }
        
        plot(jags.marginal$x, jags.marginal$y, type="l", xlab="", ylab="density", ylim=c(0,max(jags.marginal$y,hyper.marginal[,2])))
            #, xlim=range(c(marginal2[,1])) )
        lines(hyper.marginal[,1], hyper.marginal[,2], col=3)
#        plot(hyper.marginal[,1], hyper.marginal[,2], col=3, type="l", xlab="", ylab="density")
#        lines(jags.marginal$x, jags.marginal$y )
        legend(legend=c("mcmc","hyperpar"), col=c(1,3),x="topright", inset=0, bty="o", lty=1, cex=.6)
        title(main=.data.model@gibbs.watch.list[i])
    }
    dev.off()
}

load.jags.samples=function (.data.model, jags.samples.tmp, prior,seed=NA) {
    for (i in 1:.data.model@gibbs.chains) {
        load(getSamplesFileName(.data.model@name,seed,i,prior))
        .data.model@gibbs.samples[[i]] = jags.samples.tmp[,.data.model@gibbs.watch.list,drop=FALSE]
        print(getVarComponent(.data.model@gibbs.samples[[i]]))
    }
    .data.model
}

# calculate entropy
# p can be count vector or probability vector, but not a vector of membership indicator
H=function (p) { 
    if (sum(p)!=1) p=p/sum(p) # if p is count, transform to probability
    p=p[p!=0] # remove zero entries
    sum(-p*log2(p)) 
}

# each row is one observation in euclidean space
sum.of.square<-withinss<-function(x) {
    mean1 = colMeans(x)
    ss1 = sum( apply(x, 1, function(x) sum((x-mean1)**2) ) )
    ss1
}

# remove file extension from file name
rmExt = function(name){
    substr(name, 1, lastIndex(name, ".")-1 )
}
getExt = function(name){
    substr(name, lastIndex(name, ".")+1, nchar(name) )
}




# generate random number from (generalized) Bernoulli dist
# if generalized, output is 1,2,3,...; otherwise, output is 1 or 0
# n is number of random numbers to generate, prob can be one number or a vector
rbern=function(n, prob, generalized=FALSE) {
    if (!generalized){
        # binary outcome
        if (length(prob)==1) {
            rbinom(n,1,prob)
        } else {
            prob=cbind(1:n,prob)[,2]
            rbinom(n,prob,size=1)
#            sapply(prob, function(p) rbinom(1,1,p))
        }
    } else {
        x=rmultinom(n, 1, prob)
        out = apply(x, 2, function(y) which(y==1))
        out
    }
}

dbern=function(x,prob,log=FALSE){
    out=ifelse(x==1, prob, 1-prob)
    ifelse(log, log(out), out)
}

# correlated Bernoulli
dcorbern=function(x, p, a, log=FALSE){
    out=dbern(x[1],p)
    if (length(x)>1) {
        for (i in 2:length(x)) {
            out = out * dbern(x[i], p+a*(x[i-1]-p))
        }
    }
    ifelse(log, log(out), out)
}


roundup=function (value, digits) {
    format(round(value, digits), nsmall=digits, scientific=FALSE) 
}
format.int=function (value, digits, fill="0") {
    formatC(value, format="d", flag=fill, width=digits) 
}

# from Cai Tianxi
## count how many YY's are smaller or equal to yy
N.L.E <- function(yy, YY)  ## sum I(YY <= yy[i])
{
   rank(c(yy+1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.L <- function(yy, YY)  ## sum I(YY < yy[i])
{
   rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
}
N.G.E <- function(yy, YY)  ## sum I(YY >= yy[i])
{
   length(YY)-(rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy))
}


        
## next function copied from DiagnosisMed, had to modify it
#ROC=function (gold, test, CL = 0.95, Cost = 1, Prevalence = 0, Plot = TRUE, 
#    Plot.point = "Min.ROC.Dist", p.cex = 1, Full = FALSE, Print = TRUE) 
#{
#    if (any(is.na(test) | is.na(gold))) {
#        stop("It seems there are NAs either in the index test or in the reference test. Consider imputing or removing NAs!")
#    }
#    test.table <- table(test, gold)
#    if (dim(test.table)[2] != 2) {
#        stop("It seems that your gold standard has more than 2 categories")
#    }
#    CL <- CL
#    cost <- Cost
#    sample.size <- sum(test.table)
#    sample.prevalence <- (sum(test.table[, 2])/sample.size)
#    if (Prevalence == 0) {
#        pop.prevalence <- sample.prevalence
#    }
#    if (Prevalence > 0) {
#        (pop.prevalence <- Prevalence)
#    }
#    if (is.numeric(gold) == TRUE) {
#        X <- sort(test[gold == 0])
#        Y <- sort(test[gold == 1])
#        AUC <- ((as.double(length(test[gold == 0]))) * (as.double(length(test[gold == 
#            1]))) + ((as.double(length(test[gold == 0]))) * ((as.double(length(test[gold == 
#            0]))) + 1))/2 - sum(rank(test, ties.method = "average")[gold == 
#            0]))/((as.double(length(test[gold == 0]))) * (as.double(length(test[gold == 
#            1]))))
#        AUC[AUC < 0.5] <- 1 - AUC
#    }
#    if (is.factor(gold) == TRUE) {
#        X <- sort(test[gold == "negative"])
#        Y <- sort(test[gold == "positive"])
#        AUC <- ((as.double(length(test[gold == "negative"]))) * 
#            (as.double(length(test[gold == "positive"]))) + ((as.double(length(test[gold == 
#            "negative"]))) * ((as.double(length(test[gold == 
#            "negative"]))) + 1))/2 - sum(rank(test, ties.method = "average")[gold == 
#            "negative"]))/((as.double(length(test[gold == "negative"]))) * 
#            (as.double(length(test[gold == "positive"]))))
#        AUC[AUC < 0.5] <- 1 - AUC
#    }
#    m <- as.double(length(X))
#    n <- as.double(length(Y))
#    test.summary <- round(c(summary(test), sd(test)), digits = 5)
#    test.summary <- rbind(test.summary, round(c(summary(X), sd(X)), 
#        digits = 5))
#    test.summary <- rbind(test.summary, round(c(summary(Y), sd(Y)), 
#        digits = 5))
#    colnames(test.summary) <- c("Min.", "1st Qu.", "Median", 
#        "Mean", "3rd Qu.", "Max.", "SD")
#    rownames(test.summary) <- c("Overall summary", "Without disease", 
#        "With disease")
#    D10X <- function(Xi) {
#        (1/n) * sum(Y >= Xi[1])
#    }
#    D01Y <- function(Yi) {
#        (1/m) * sum(Yi[1] >= X)
#    }
#    VAR.AUC <- sum((tapply(X, X, "D10X") - AUC)^2)/(m * (m - 
#        1)) + sum((tapply(Y, Y, "D01Y") - AUC)^2)/(n * (n - 1))
#    SD.AUC <- sqrt(VAR.AUC)
#    alpha <- 1 - CL
#    AUC.summary <- c(AUC - qnorm(1 - alpha/2) * SD.AUC, AUC, 
#        AUC + qnorm(1 - alpha/2) * SD.AUC)
#    D <- sum(test.table[, 2])
#    ND <- sum(test.table[, 1])
#    test.values <- (as.numeric(rownames(unclass(test.table))))
#    test.diag.table <- as.data.frame(test.values)
#    for (i in 1:nrow(test.diag.table)) {
#        test.diag.table$TP[i] <- sum(test.table[i:nrow(test.table), 
#            2])
#        test.diag.table$FN[i] <- sum(test.table[1:i - 1, 2])
#        test.diag.table$FP[i] <- sum(test.table[i:nrow(test.table), 
#            1])
#        test.diag.table$TN[i] <- sum(test.table[1:i - 1, 1])
#    }
#    test.diag.table$Sensitivity <- round(test.diag.table$TP/D, 
#        digits = 4)
#    test.diag.table$Se.inf.cl <- round(binom.wilson(test.diag.table$TP, 
#        D, conf.level = CL)[4]$lower, digits = 4)
#    test.diag.table$Se.sup.cl <- round(binom.wilson(test.diag.table$TP, 
#        D, conf.level = CL)[5]$upper, digits = 4)
#    test.diag.table$Specificity <- round(test.diag.table$TN/ND, 
#        digits = 4)
#    test.diag.table$Sp.inf.cl <- round(binom.wilson(test.diag.table$TN, 
#        ND, conf.level = CL)[4]$lower, digits = 4)
#    test.diag.table$Sp.sup.cl <- round(binom.wilson(test.diag.table$TN, 
#        ND, conf.level = CL)[5]$upper, digits = 4)
#    test.diag.table$PPV <- round(test.diag.table$TP/(test.diag.table$TP + 
#        test.diag.table$FP), digits = 4)
#    test.diag.table$PPV.inf.cl <- round(binom.wilson(test.diag.table$TP, 
#        (test.diag.table$TP + test.diag.table$TP), conf.level = CL)[4]$lower, 
#        digits = 4)
#    test.diag.table$PPV.sup.cl <- round(binom.wilson(test.diag.table$TP, 
#        (test.diag.table$TP + test.diag.table$FN), conf.level = CL)[5]$upper, 
#        digits = 4)
#    test.diag.table$NPV <- round(test.diag.table$TN/(test.diag.table$TN + 
#        test.diag.table$FN), digits = 4)
#    test.diag.table$NPV.inf.cl <- round(binom.wilson(test.diag.table$TN, 
#        (test.diag.table$TN + test.diag.table$FN), conf.level = CL)[4]$lower, 
#        digits = 4)
#    test.diag.table$NPV.sup.cl <- round(binom.wilson(test.diag.table$TN, 
#        (test.diag.table$TN + test.diag.table$FN), conf.level = CL)[5]$upper, 
#        digits = 4)
#    test.diag.table$PLR <- round(test.diag.table$Sensitivity/(1 - 
#        test.diag.table$Specificity), digits = 2)
#    test.diag.table$PLR.inf.cl <- round(exp(log(test.diag.table$PLR) - 
#        (qnorm(1 - ((1 - CL)/2), mean = 0, sd = 1)) * sqrt((1 - 
#            test.diag.table$Sensitivity)/((D) * test.diag.table$Specificity) + 
#            (test.diag.table$Specificity)/((ND) * (1 - test.diag.table$Specificity)))), 
#        digits = 2)
#    test.diag.table$PLR.sup.cl <- round(exp(log(test.diag.table$PLR) + 
#        (qnorm(1 - ((1 - CL)/2), mean = 0, sd = 1)) * sqrt((1 - 
#            test.diag.table$Sensitivity)/((D) * test.diag.table$Specificity) + 
#            (test.diag.table$Specificity)/((ND) * (1 - test.diag.table$Specificity)))), 
#        digits = 2)
#    test.diag.table$NLR <- round((1 - test.diag.table$Sensitivity)/test.diag.table$Specificity, 
#        digits = 2)
#    test.diag.table$NLR.inf.cl <- round(exp(log(test.diag.table$NLR) - 
#        (qnorm(1 - ((1 - CL)/2), mean = 0, sd = 1)) * sqrt((test.diag.table$Sensitivity)/((D) * 
#            (1 - test.diag.table$Sensitivity)) + (1 - test.diag.table$Specificity)/((ND) * 
#            (test.diag.table$Specificity)))), digits = 2)
#    test.diag.table$NLR.sup.cl <- round(exp(log(test.diag.table$NLR) + 
#        (qnorm(1 - ((1 - CL)/2), mean = 0, sd = 1)) * sqrt((test.diag.table$Sensitivity)/((D) * 
#            (1 - test.diag.table$Sensitivity)) + (1 - test.diag.table$Specificity)/((ND) * 
#            (test.diag.table$Specificity)))), digits = 2)
#    test.diag.table$Accuracy <- (test.diag.table$TN + test.diag.table$TP)/sample.size
#    test.diag.table$DOR <- ((test.diag.table$TN) * (test.diag.table$TP))/((test.diag.table$FP) * 
#        (test.diag.table$FN))
#    test.diag.table$DOR <- ifelse(test.diag.table$DOR == Inf, 
#        NA, test.diag.table$DOR)
#    test.diag.table$Error.rate <- ((test.diag.table$FP) + (test.diag.table$FN))/sample.size
#    test.diag.table$Accuracy.area <- ((test.diag.table$TP) * 
#        (test.diag.table$TN))/(D * ND)
#    test.diag.table$Max.Se.Sp <- test.diag.table$Sensitivity + 
#        test.diag.table$Specificity
#    test.diag.table$Youden <- test.diag.table$Sensitivity + test.diag.table$Specificity - 
#        1
#    test.diag.table$Se.equals.Sp <- abs(test.diag.table$Specificity - 
#        test.diag.table$Sensitivity)
#    test.diag.table$MinRocDist <- (test.diag.table$Specificity - 
#        1)^2 + (1 - test.diag.table$Sensitivity)^2
#    test.diag.table$Efficiency <- (test.diag.table$Sensitivity * 
#        (pop.prevalence)) + ((1 - (pop.prevalence)) * test.diag.table$Specificity)
#    test.diag.table$MCT <- (1 - (pop.prevalence)) * (1 - test.diag.table$Specificity) + 
#        (cost * (pop.prevalence)) * (1 - test.diag.table$Sensitivity)
#    test.best.cutoff <- as.data.frame(rbind(test.diag.table[which.max(test.diag.table$Accuracy), 
#        c(1, 6:11, 18:20)], test.diag.table[which.max(test.diag.table$DOR), 
#        c(1, 6:11, 18:20)], test.diag.table[which.min(test.diag.table$Error.rate), 
#        c(1, 6:11, 18:20)], test.diag.table[which.max(test.diag.table$Accuracy.area), 
#        c(1, 6:11, 18:20)], test.diag.table[which.max(test.diag.table$Max.Se.Sp), 
#        c(1, 6:11, 18:20)], test.diag.table[which.max(test.diag.table$Youden), 
#        c(1, 6:11, 18:20)], test.diag.table[which.min(test.diag.table$Se.equals.Sp), 
#        c(1, 6:11, 18:20)], test.diag.table[which.min(test.diag.table$MinRocDist), 
#        c(1, 6:11, 18:20)], test.diag.table[which.max(test.diag.table$Efficiency), 
#        c(1, 6:11, 18:20)], test.diag.table[which.min(test.diag.table$MCT), 
#        c(1, 6:11, 18:20)]))
##    rownames(test.best.cutoff) <- c("Max. Accuracy", "Max. DOR", 
##        "Min. Error rate", "Max. Accuracy area", "Max. Sens+Spec", 
##        "Max. Youden", "Se=Sp", "Min. ROC distance", "Max. Efficiency", 
##        "Min. MCT")
#    reteval <- list(pop.prevalence = pop.prevalence, sample.size = sample.size, 
#        sample.prevalence = sample.prevalence, test.summary = test.summary, 
#        AUC.summary = AUC.summary, test.table = test.table, test.best.cutoff = test.best.cutoff, 
#        test.diag.table = test.diag.table, CL = CL, cost = cost)
#    class(reteval) <- "ROC"
#    if (Print == TRUE) {
#        if (Full == TRUE) {
#            print(reteval, Full = TRUE)
#        }
#        else {
#            print(reteval)
#        }
#    }
#    if (Plot == TRUE) {
#        plot(reteval, Plot.point = Plot.point, p.cex = p.cex)
#    }
#    invisible(reteval)
#}


# both dat must have two columns, each row is dat from one subject
# x.ori=0; xaxislabels=rep("",2); cex.axis=1; add=FALSE; xlab=""; ylab=""; pcol=NULL; lcol=NULL
my.interaction.plot=function(dat, x.ori=0, xaxislabels=rep("",2), cex.axis=1, add=FALSE, xlab="", ylab="", pcol=NULL, lcol=NULL, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,2),ylim=range(dat), ylab=ylab, xlab=xlab, xaxt="n", ...)
    cex=.25; pch=19
    if (is.null(lcol)) lcol=ifelse(dat[,1]>dat[,2],"red","black")
    for (i in 1:nrow(dat)) {
        points (1+x.ori, dat[i,1], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,1]))
        points (2+x.ori, dat[i,2], cex=cex, pch=pch, col=ifelse(is.null(pcol), 1, pcol[i,2]))
        lines (1:2+x.ori, dat[i,], lwd=.25, col=lcol[i])
    }
    axis(side=1, at=1:2+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}

myboxplot <- function(object, ...) UseMethod("myboxplot") 

# myboxplot.formula and myboxplot.list make a boxplot with data points and do inferences for two group comparions. 
# cex=.5; ylab=""; xlab=""; main=""; box=FALSE; highlight.list=NULL; at=NULL;pch=1;col=1;
myboxplot.formula=function(formula, data, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test=c("t","w"), ...){
    
    if (box) {
        boxplot(formula, data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL,...)
    } else {
        boxplot(formula, data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL, border="white", ...)
    }
        
    dat.tmp=model.frame(formula, data)
    if(is.null(at)){
        xx=jitter(as.numeric(as.factor(dat.tmp[,2])))
    } else{
        xx=jitter(at[as.factor(dat.tmp[,2])])
    }    
    points(xx, dat.tmp[[1]], cex=cex,pch=pch,col=col)
    
    # inference
    x.unique=unique(dat.tmp[[2]])
    if (length(test)>0 & length(x.unique)==2) {
        y.1=dat.tmp[[1]][dat.tmp[[2]]==x.unique[1]]
        y.2=dat.tmp[[1]][dat.tmp[[2]]==x.unique[2]]
        sub="p-value: "
        if ("t" %in% test) sub=sub%+%signif(t.test(y.1, y.2)$p.value,2) 
        if (length(test)==2) sub=sub%+%"|"
        if ("w" %in% test) sub=sub%+%signif(wilcox.test(y.1, y.2)$p.value,2)
        title(sub=sub)
    }
    
}
myboxplot.list=function(data, cex=.5, ylab="", xlab="", main="", box=TRUE, at=NULL, pch=1, col=1, test=c("t","w"), ...){
    
    if (box) {
        boxplot(data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL,...)
    } else {
        boxplot(data, range=0, ylab=ylab, xlab=xlab, at=at, main=main, col=NULL, border="white", ...)
    }
        
    if(is.null(at)){
        xx=jitter(rep.int(1:length(data), times=sapply(data, length)))
    } else{
        xx=jitter(at[rep.int(1:length(data), times=sapply(data, length))])
    }    
    points(xx, unlist(data), cex=cex,pch=pch,col=col)
    
    # inference
    if (length(test)>0 & length(data)==2) {
        y.1=data[[1]]
        y.2=data[[2]]
        sub="p-value: "
        if ("t" %in% test) sub=sub%+%signif(t.test(y.1, y.2)$p.value,2) 
        if (length(test)==2) sub=sub%+%"|"
        if ("w" %in% test) sub=sub%+%signif(wilcox.test(y.1, y.2)$p.value,2)
        title(sub=sub)
    }    
    
}


# an old function I wrote
#myboxplot=function(formula,data,showNames=TRUE, subset=NULL,... ){
#    if (!is.null(subset)) 
#        data = subset (data, subset)
#
#    names = boxplot(formula, data=data)$names
#
#    resp = as.character ( formula )[[2]]
#    splittedforms = splitFormula( formula , "+")
#    groupvars = list (length(splittedforms ))
#    for (i in 1:length(splittedforms ) ) {
#        name = as.character ( splittedforms [[i]] )[2]
#        groupvars[[i]] = data [,name]
#    }
#    groupeddata = split (data, groupvars)
#    summarydata = sapply ( groupeddata , function (y) { sd1=sd(y[resp]); m1=mean(y[resp]); c ( m1-2*sd1, m1-1*sd1, m1, m1+1*sd1, m1+2*sd1 ) } )
#    
#    ifelse(showNames,
#      bxp( list(stats= summarydata, names= names ),  ...),
#      bxp( list(stats= summarydata, names=rep("", length (groupeddata ) ) ), ...)
#    )   
#
#    names
#}
#


# called butterfly.plot, because it is meant to plot two treatment arms at two time points, the two arms are plotted in a mirror fashion, see "by analyte.pdf" for an example
# if dat2 is null: dat is matrix with four columns. each row is one subject, the columns will be plotted side by side, with lines connecting values from one ptid
# if dat2 is not null, dat has two columns, which are plotted side by side with lines connecting them, same for dat2
# if add is true, no plot function will be called
butterfly.plot=function (dat, dat2=NULL, add=FALSE, xaxislabels=rep("",4), x.ori=0, xlab="", ylab="", cex.axis=1, ...){
    if (!add) plot(0,0,type="n",xlim=c(1,4),ylim=range(dat), xaxt="n", xlab=xlab, ylab=ylab, ...)
    for (i in 1:nrow(dat)) {
        lines (1:2+x.ori, dat[i,1:2], lwd=.25, col=ifelse(dat[i,1]<=dat[i,2],"red","black"))
        if (is.null(dat2)) {
            lines (2:3+x.ori, dat[i,2:3], lwd=.25, col="lightgreen")
            lines (3:4+x.ori, dat[i,3:4], lwd=.25, col=ifelse(dat[i,3]<=dat[i,4],"black","red"))
        }
    }
    if (!is.null(dat2)) {
        for (i in 1:nrow(dat2)) {
            lines (3:4+x.ori, dat2[i,1:2], lwd=.25, col=ifelse(dat2[i,1]<=dat2[i,2],"black","red"))
        }
    }
    axis(side=1, at=1:4+x.ori, labels=xaxislabels, cex.axis=cex.axis)
}

# dat is a matrix of two columns
merge.overlap=function(dat) {
    out=NULL
    cursor=1
    j=0
    for (i in 1:nrow(dat)) {
        if (dat[i,1]<=cursor+1) {
            # overlapping out
            out[j,2]=dat[i,2]
            cursor=dat[i,2]
        } else {
            # new epitope
            out=rbind (out, c(dat[i,1], dat[i,2]))
            cursor=dat[i,2]
            j=j+1
        }
    }
    out
}

##  From the R-help e-mail by Ted Harding: http://tolstoy.newcastle.edu.au/R/e2/help/07/03/12853.html
##  See also http://tolstoy.newcastle.edu.au/R/help/05/05/4254.html
pava <- function(x, wt = rep(1, length(x)))
{
    n <- length(x)
    if (n <= 1) return(x)
    lvlsets <- 1:n
    repeat 
    {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) break
        i <- min((1:(n-1))[viol])
    
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i+1]
        ilvl <- ( (lvlsets == lvl1) | (lvlsets == lvl2) )
    
        x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])     
        lvlsets[ilvl] <- lvl1
    }
    x
} 

# the integratd density of a normal random variable, whose mean and precision follow a normal gamma distribution. It is a three-parameter t distribution. 
# keywords: normgamma
# when x has length greater than 1 and same.distr is true, the data are considered to be from the same mean and variance
dnorm.norm.gamma = function(x, p, same.distr=FALSE, log=FALSE) {
    mu.0=p[1]; lambda=p[2]; a=p[3]; beta=p[4]
    n=length(x)
    if (!same.distr) {
        if (n>1) {
            mu.0=rep(mu.0, n)
            lambda=rep(lambda, n)
            a=rep(a, n)
            beta=rep(beta, n)
        }
        ll= 1/2*log(lambda/(2*pi*(1+lambda))) + log(gamma(a+1/2)/gamma(a)) + a*log(beta) - (a+1/2)*log(beta+lambda*(x-mu.0)**2/(2*(1+lambda)))
    } else{
        s2=var(x)
        ll= 1/2*log(lambda/((2*pi)**n*(n+lambda))) + log(gamma(a+n/2)/gamma(a)) + a*log(beta) - (a+n/2)*log(beta + n*s2/2 +n*lambda*(mean(x)-mu.0)**2/(2*(1+lambda)))
    }
    ll=unname(ll)
    ifelse(log, ll, exp(ll))
}


# simulation samples from a normal random variable, whose mean and precision follow a normal gamma distribution
rnorm.norm.gamma = function(n, mu.0, lambda, alpha, beta) {
    tao=rgamma(n, shape=alpha, rate=beta)
    mu=rnorm(n, mean=mu.0, sd=sqrt(1/(lambda*tao)))
    rnorm(n, mean=mu, sd=tao**-.5)
}

# simulate correlated normal random variables, correlation is sqrt(alpha)
rnorm.cor = function (n, mu, sd, alpha) {
    out=numeric(n)
    out[1]=rnorm(1, mu, sd)
    if (n>1) 
        for (i in 2:n) {
            out[i]=sqrt(alpha)*out[i-1] + sqrt(1-alpha)*rnorm(1, mu, sd)
        }
    out
}

# copied from pairs help page
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, method="pearson", ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, method=method, use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.nothing=function(x, ...) {}
mypairs=function(dat, method="pearson"){
    pairs(dat, lower.panel=panel.smooth, upper.panel=panel.cor, method=method, diag.panel=panel.hist)
}


# if lty is specified, a line will be drawn
mylegend=function(legend, x, lty=NULL,bty="n", ...) {
    x=switch(x, "topleft", "top", "topright", "left", "center" , "right", "bottomleft", "bottom", "bottomright")
    legend(bty=bty,x=x, legend=legend, lty=lty, ...)
}

# mixture normal density funtion 
# mix.p: proportion of first component
dmixnorm=function (x, mix.p, sd1, sd2, log=FALSE){
    out=mix.p*dnorm (x,0,sd1) + (1-mix.p)*dnorm (x,0,sd2)
    if(log) out=log(out)
    out
}

getFileStem=function(file.name){
    strsplit(file.name,"\\.")[[1]][1]
}
getStem=getFileStem
fileStem=getFileStem

# get first prinipal component
pr.1=function(x){
    x.s=scale(x)
    pr.s = prcomp(x.s)
    out=c(x.s %*% pr.s$rotation[,1])
    attr(out, "rotation")=pr.s$rotation[,1]
    out
}

# category.var is 
myreshapewide=function(dat, category.var, outcome.var, idvar=NULL){
    if (is.null(idvar)) {
        idvar=setdiff(names(dat), c(category.var,outcome.var))
        # if any of idvar has NA then it is a problem if not treated, here we opt to remove such columns
        idvar=idvar[sapply(idvar, function (idvar.) all(!is.na(dat[,idvar.])) )]
    } else {
        # remove those columns with changing values within an id
        tmp=apply(aggregate (x=dat[,!names(dat) %in% idvar,drop=FALSE], by=dat[,names(dat) %in% idvar,drop=FALSE], function(y) length(rle(y)$values)==1), 2, all)
        varying.var=names(tmp)[!tmp]
        dat=dat[,!names(dat) %in% setdiff(varying.var, c(category.var,outcome.var)),drop=FALSE]
    }
    reshape(dat, direction="wide", v.names=outcome.var, timevar=category.var, idvar=idvar )
}

plot.cor <- function(object, ...) UseMethod("plot.cor") 

plot.cor.default=function(x,y,...){
    dat=data.frame(x,y)
    names(dat)=c(names(x), names(y))
    plot.cor(as.formula(names(y)%+%"~"%+%names(x)), dat, ...)
}

# col can be used to highlight some points
plot.cor.formula=function(formula,data,main="",method=c("pearson","spearman"),col=1,cex=.5,plot.abline=TRUE,...){
    vars=dimnames(attr(terms(formula),"factors"))[[1]]
    cor.=sapply (method, function (method) {
        cor(data[,vars[1]],data[,vars[2]],method=method,use="p")
    })
    main=main%+%ifelse(main=="","",", ")
    main=main%+%"cor: "%+%concatList(round(cor.,2),"|")
#    if (length(cor.)==1) {
#        main=main%+%"cor: "%+%round(cor.,2) 
#    } else {
#        main=main%+%"cor: "%+%round(cor.,2) 
#    }
    plot(formula,data,main=main,col=col,cex=cex,...)
    if(plot.abline) abline(0,1)
    cor.
}

# compute covariability as defined in Clarke 1995 page 2271, second formula
# x is a vector of characters, so is y
covariability=function(x,y){
    tmp=table(x)
    p.x=data.frame(tmp/(sum(tmp)))
    tmp=table(y)
    p.y=data.frame(tmp/(sum(tmp)))
    
    # tabulate all pairwise combination of aa
    pair = paste(x,y,sep="")
    tmp=table(pair)
    dat=data.frame(tmp/(sum(tmp)))
    
    ps=numeric(nrow(dat))
    for (i in 1:nrow(dat)) {
        row.=dat[i,]
        a2=as.character(row.[1,1])
        p.i=p.x$Freq[p.x$x==substr(a2,1,1)]
        p.j=p.y$Freq[p.y$y==substr(a2,2,2)]
        p.ij=row.[1,2]
        ps[i]=p.ij**2*log(p.ij/p.i/p.j)        
    }
    sum(ps)
}
## test
#dat=readFastaFile ("D:/1CHMM/domains/SET/SETpfamseed/seq/SETpfamseed_aligned.fasta")
#seqs.mat.a=mysapply(dat, s2c)
#covar=matrix(0,ncol(seqs.mat.a),ncol(seqs.mat.a))
#for (i in 1:(ncol(seqs.mat.a)-1)){
#    for (j in (i+1):ncol(seqs.mat.a)){
#        covar[i,j]<-covar[j,i]<-covariability(seqs.mat.a[,i], seqs.mat.a[,j])
#    }
#}
#sort(round(covar,3), decre=TRUE)[1:100]

# output format:
#+1 1:0.708333 2:1 3:1 4:-0.320755 5:-0.105023 6:-1 7:1 8:-0.419847 9:-1 10:-0.225806 12:1 13:-1 
# y: vector of 1 and non-1
# z: matrix or data.frame of covariates
# file.name: e.g. test. There should be no file extension
# ws: vector of weights
write.svm.input=function(y, z, file.name="svm_input", ws=NULL){
    rows=sapply(1:length(y), function(i){
        paste(c(ifelse(y[i]==1,"+1","-1"),(1:ncol(z))%+%":"%+%z[i,]), collapse=" ")
    })
    write(paste(rows,collapse="\n"), file = file.name%+%".dat", append = FALSE, sep = " ")        
    # write weight file
    if (!is.null(ws)) write(paste(ws,collapse="\n"), file = file.name%+%".wgt", append = FALSE, sep = " ")    
}

# matrix trace function
tr=function(m) sum(diag(m))

# keep column names
# return data frame if input is data frame
myscale=function(x){
    flag=is.data.frame(x)
    nam=dimnames(x)[[2]]
    aux=scale(x)
    dimnames(aux)[[2]]=nam
    if(flag) {
        tmp=as.data.frame(aux)
        mostattributes(tmp)=attributes(x)
        attr(tmp,"scaled:center")=attr(aux,"scaled:center")
        attr(tmp,"scaled:scale")=attr(aux,"scaled:scale")
#        
#        # add old attr
#        tmp1=attributes(aux)
#        for(a in tmp1){
#            #print(a)
#            print(names(a))
#            print("a")
#            #attr(tmp,"scaled:scale")=attr(x,"scaled:scale")
#        }
        aux=tmp
        names(aux)=names(x)
    }
    aux
}

# m is data frame or matrix, each row is an observation
# w is a named vector, the names may not match m, may not have same number even
# return a vector
weighted.ave=function(m, w){
    if(is.data.frame(m)){
        nam=names(m)
        m=as.matrix(m)
    } else {
        nam=colnames(m)
    }
    ret=w[match(nam, names(w))]
    ret[is.na(ret)]=0
    names(ret)=nam
    c(m %*% ret )
}




abline.pts=function(pt1, pt2=NULL){
    if (is.null(pt2)) {
        if (nrow(pt1)>=2) {
            pt2=pt1[2,]
            pt1=pt1[1,]
        } else {
            stop("wrong input")
        }
    }
    slope=(pt2-pt1)[2]/(pt2-pt1)[1]
    intercept = pt1[2]-slope*pt1[1]
    abline(intercept, slope)
}
#abline.pts(c(1,1), c(2,2))

abline.pt.slope=function(pt1, slope,...){
    intercept = pt1[2]-slope*pt1[1]
    abline(intercept, slope,...)
}
#abline.pt.slope(c(1,1), 1)


mean.med=function(x, na.rm=FALSE){
    c(mean=mean(x, na.rm=na.rm), sd=sd(x, na.rm=na.rm), median=median(x, na.rm=na.rm), iqr=IQR(x, na.rm=na.rm), var=var(x, na.rm=na.rm))
}


# st.fit is an object returned by sn::st.mle
get.st.mean=function(st.fit) {
    p=st.fit$dp
    
    xi=p["location"];
    alpha=p["shape"];
    df=p["df"]
    #omega=sd(X);
    omega=p["scale"]
    delta=alpha/sqrt(1+alpha^2);
    mu=delta*sqrt(df/pi)*gamma(0.5*(df-1))/gamma(0.5*df);
    
    # the 4 first moments
    moment1=xi+omega*mu;
    moment2=xi^2 + 2*omega*mu*xi + omega^2 * df/(df-2);
    moment3=xi^3 + 3*omega*mu*xi^2 + 3*omega^2*df/(df-2)*xi +
    omega^3*mu*(3-delta^2)*df/(df-3);
    moment4=xi^4 + 4*omega*mu*xi^3 + 6*omega^2*df/(df-2)*xi^2 +
    4*omega^3*mu*(3-delta^2)*df/(df-3)*xi+omega^4*3*df^2/((df-2)*(df-4));
    
    # the 4 useful measures
    mean=moment1;
    var=moment2;
    skew=moment3/var^(3/2);
    kurt=moment4/var^2 - 3;
    
    c(mean)
}

# default row.names to FALSE
# file name needs no file extension
mywrite.csv = function(x, file="tmp", row.names=FALSE, digits=NULL, ...) {  
    if (!is.null(digits)) {
        if(length(digits)==1) {
            x=round(x,digits)
        } else {
            for (i in 1:ncol(x)) {
                x[,i]=round(x[,i], digits[i])
            }                
        }
    }
    print("Writing csv file to "%+%getwd()%+%"/"%+%file%+%".csv", quote=FALSE)
    write.csv(x, file=file%+%".csv", row.names=row.names, ...)
}



read.tsv=function (file, header = TRUE, sep = "\t", ...) {
    read.csv(file, header = header, sep = sep, ...)
}

read.sv=function (file, header = TRUE, ...) {
    sep=","
    if (getExt(file)=="tsv") sep="\t"
    read.csv(file, header = header, sep = sep, ...)
}

keepWarnings <- function(expr) { 
    localWarnings <- list() 
    value <- withCallingHandlers(expr, 
    warning = function(w) { 
    localWarnings[[length(localWarnings)+1]] <<- w 
    invokeRestart("muffleWarning") 
    }) 
    list(value=value, warnings=localWarnings) 
} 

mymatplot=function(x, make.legend=TRUE, legend=NULL, legend.x=9, lty=1:5, pch=NULL, col=1:6, legend.title=NULL, xlab=NULL, ylab="", draw.x.axis=TRUE, bg=NULL, ...) {
    if (is.null(xlab)) xlab=names(dimnames(x))[1]
    if (is.null(legend.title)) legend.title=names(dimnames(x))[2]
    matplot(x=x, lty=lty, pch=pch, col=col, xlab=xlab, xaxt="n", ylab=ylab, bg=bg, ...)
    if(draw.x.axis) axis(side=1, at=1:nrow(x), labels=rownames(x))
    if (make.legend) {
        if (is.null(legend)) legend=colnames(x)
        mylegend(legend, x=legend.x, lty=lty, title=legend.title, pch=pch, col=col, pt.bg=bg)
    }
}


#Hosmer-Lemeshow goodness of fit test
hosmerlem = function(y, yhat, g=10) {
  cutyhat = cut(yhat,
     breaks = quantile(yhat, probs=seq(0,
       1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq,p.value=P))
}

# s="prefix_name" is changed to "name"
remove.prefix=function(s,sep="_"){
    tmp=strsplit(s,sep)    
    sapply(tmp, function (x) {
        if (length(x)==1) return (x)
        concatList(x[-1],sep)
    })
}

# make table that shows both counts/frequency and proportions
# style 1: count only; 2: count + percentage; 3: percentage only
table.prop=function (x,y,digit=1,style=1) {
    tbl=table(x,y)
    prop = apply (tbl, 2, prop.table)
    if (style==2) {
        res = tbl %+% " (" %+% round(prop*100,1) %+% ")"
    } else if (style==3) {
        res = round(prop*100,digit)
    } else res=tbl
    res=matrix(res, nrow=2)    
    dimnames(res)=dimnames(tbl)
    names(dimnames(res))=NULL
    res
}


# from Thomas, found on R user group
methods4<-function(classes, super=FALSE, ANY=FALSE){ 
    if (super) classes<-unlist(sapply(classes, function(cl) getAllSuperClasses(getClass(cl)))) 
    if (ANY) classes<-c(classes,"ANY") 
    gens<-getGenerics()@.Data 
    sigs<-lapply(gens, function(g) linearizeMlist(getMethods(g))@classes) 
    names(sigs)<-gens@.Data 
    sigs<-lapply(sigs, function(gen){ gen[unlist(sapply(gen, function(sig) any(sig %in% classes)))]}) 
    sigs[sapply(sigs,length)>0] 
} 
 

# this function contains "\""), which makes all the following line be miss-interpreted as being in quotes
myprint.default = function (..., newline=TRUE, digits=3) {     
    digits.save=getOption("digits")
    options(digits=digits)
    object <- as.list(substitute(list(...)))[-1]
    x=list(...)
    for (i in 1:length(x)) {
        tmpname <- deparse(object[[i]])[1]
        #str(tmpname)
        #str(gsub("\\\\","\\",gsub("\"", "", tmpname)))
        #str(x[[i]])
        #if (gsub("\\\\","\\",gsub("\"", "", tmpname))!=x[[i]]) {
        if (contain(tmpname, "\"") | contain(tmpname, "\\")) {
            for (a in x[[i]]) cat(a," ")
        } else {
            cat (tmpname %+% " = ")
            for (a in x[[i]]) cat(a," ")
            if (i!=length(x)) cat ("; ")
        }
    }
    if (newline)  cat("\n")
    options(digits=digits.save)
}
#a="hello"; b="you"; myprint (a, b); myprint ("test"); myprint.default ("\t")
