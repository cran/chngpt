performance.unit.test=function(formula.1, formula.2, family, data, B, I){
    
    y=model.frame(formula.1, data)[,1]
    Z=model.matrix(formula.1, data)
    tmp=model.matrix(formula.2, data)[,-1,drop=F]
    chngpt.var.name=setdiff(colnames(tmp), colnames(Z))[1]
    if(is.na(chngpt.var.name)) stop("Something is wrong. Check the formulal. ")
    z.1.name=intersect(colnames(tmp), colnames(Z))
    chngpt.var = tmp[,chngpt.var.name]
    chngpt.var.sorted=sort(chngpt.var)
    #print(chngpt.var.sorted); chngpt.var.sorted=chngpt.var[order(chngpt.var)]; print(chngpt.var.sorted)
    Z.sorted=Z[order(chngpt.var),,drop=FALSE]
    y.sorted=y[order(chngpt.var)]
    data.sorted=data[order(chngpt.var),]
    z.1 = tmp[,z.1.name] # if the intersection is a null set, z.1 is a matrix of n x 0 dimension
    has.itxn = length(z.1.name)>0
    
    type="segmented"
    n=nrow(Z)
    p.z=ncol(Z)
    p.2=switch(type, step=1, hinge=1, segmented=2, segmented2=2, stegmented=3)
    p.2.itxn=p.2*ifelse(has.itxn,2,1)
    p=p.z+p.2.itxn+1 #total number of paramters, including threshold
    
    # make formula that includes all parameters but threshold
    formula.new = if (type %in% c("segmented","segmented2","stegmented")) update(formula.1, as.formula("~.+"%.%chngpt.var.name)) else formula.1
    f.alt=get.f.alt(type, has.itxn, z.1.name, chngpt.var.name)
    formula.new=update(formula.new, as.formula(f.alt))
    
    if (FALSE) {
        myprint(type, est.method, has.itxn)     
        print(formula.new)
        myprint(p.z, p.2.itxn, p)
    }
        
    -.Call("performance_unit_test", cbind(Z.sorted,chngpt.var.sorted, if(type=="segmented") chngpt.var.sorted), as.double(y.sorted), B, I)
}
