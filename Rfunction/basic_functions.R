####################
##### ML Method
####################
# Calculating log likelihood in under-report case
logit.model.under = function(param, data, eta = c(0,0)){
  # data
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  beta = param[1:m]
  gamma = param[(m+1):(2*m)]
  gamma_t = param[(2*m+1)]
  
  # sensitivity parameter
  eta0 = eta[1]
  eta1 = eta[2]
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  # outcome models + propensity score model
  m0 = exp((x %*% gamma))/(1 + exp((x %*% gamma)))
  m1 = exp((x %*% gamma) + t.vec*gamma_t)/(1 + exp((x %*% gamma) + t.vec*gamma_t))
  ps = exp((x %*% beta))/(1 + exp((x %*% beta)))
  
  # p_{ab} = pr(Y = a, Tstar = b | X)
  p11 = (1-eta1)*m1*ps
  p10 = eta1*m1*ps + m0*(1-ps)
  p01 = (1-eta0)*(1-m1)*ps
  p00 = eta0*(1-m1)*ps + (1-m0)*(1-ps)
  
  like.val = log(p11)*y*t + log(p01)*(1-y)*t + log(p10)*y*(1-t) + log(p00)*(1-y)*(1-t)
  sum.like.val = sum(like.val)
  
  return(-sum.like.val)
}

# Calculating log likelihood in over-report case
logit.model.over = function(param, data, zeta = c(0,0)){
  # data
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  beta = param[1:m]
  gamma = param[(m+1):(2*m)]
  gamma_t = param[(2*m+1)]
  
  # sensitivity parameter
  zeta0 = zeta[1]
  zeta1 = zeta[2]
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  # outcome models + propensity score model
  m0 = exp((x %*% gamma))/(1 + exp((x %*% gamma)))
  m1 = exp((x %*% gamma) + t.vec*gamma_t)/(1 + exp((x %*% gamma) + t.vec*gamma_t))
  ps = exp((x %*% beta))/(1 + exp((x %*% beta)))
  
  # p_{ab} = pr(Y = a, Tstar = b | X)
  p11 = m1*ps + zeta1*m0*(1-ps)
  p01 = (1-m1)*ps + zeta0*(1-m0)*(1-ps)
  p10 = (1-zeta1)*m0*(1-ps)
  p00 = (1-zeta0)*(1-m0)*(1-ps)
  
  like.val = log(p11)*y*t + log(p01)*(1-y)*t + log(p10)*y*(1-t) + log(p00)*(1-y)*(1-t)
  sum.like.val = sum(like.val)
  
  return(-sum.like.val)
}

#################################
## Maximum Likelihood Method
#################################
ML.logistic.under = function(data, eta = c(0,0)){
  # data
  require(optimx)
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  init.ps = glm(t ~ x, family="binomial")
  init.m = glm(y ~ x + t, family = "binomial")
  
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  
  init.param = c(init.ps$coefficients, init.m$coefficients)
  
  # opt.res = optim(init.param, logit.model, method = "BFGS", data = data, eta = eta)
  opt.res =suppressWarnings(optimx(init.param, logit.model.under, method = "BFGS", data = data, eta = eta))
  
  est.beta = t(opt.res[1, 1:m])
  est.gamma = t(opt.res[1, (m+1):(2*m)])
  est.gamma_t = t(opt.res[1, (2*m+1)])
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  est.e =  exp((x %*% est.beta))/(1 + exp((x %*% est.beta)))
  est.m0 =  exp((x %*% est.gamma))/(1 + exp((x %*% est.gamma)))
  est.m1 =  exp((x %*% est.gamma) + t.vec %*% est.gamma_t)/(1 + exp((x %*% est.gamma) + t.vec %*% est.gamma_t))
  
  est.p0 = mean(est.m0)
  est.p1 = mean(est.m1)
  
  risk.diff = est.p1 - est.p0
  odds.ratio = (est.p1*(1-est.p0))/(est.p0*(1-est.p1))
  
  return(list(potential.outcome = c(est.p1, est.p0), risk.diff = risk.diff, est.OR = odds.ratio, beta = est.beta, gamma = est.gamma, gamma.t = est.gamma_t, prop = est.e, m0 = est.m0, m1 = est.m1, likeli.val = opt.res$value))
}

ML.logistic.over = function(data, zeta = c(0,0)){
  # data
  require(optimx)
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  init.ps = glm(t ~ x, family="binomial")
  init.m = glm(y ~ x + t, family = "binomial")
  
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  
  init.param = c(init.ps$coefficients, init.m$coefficients)
  
  # opt.res = optim(init.param, logit.model, method = "BFGS", data = data, eta = eta)
  opt.res =suppressWarnings(optimx(init.param, logit.model.over, method = "BFGS", data = data, zeta = zeta))
  
  est.beta = t(opt.res[1, 1:m])
  est.gamma = t(opt.res[1, (m+1):(2*m)])
  est.gamma_t = t(opt.res[1, (2*m+1)])
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  est.e =  exp((x %*% est.beta))/(1 + exp((x %*% est.beta)))
  est.m0 =  exp((x %*% est.gamma))/(1 + exp((x %*% est.gamma)))
  est.m1 =  exp((x %*% est.gamma) + t.vec %*% est.gamma_t)/(1 + exp((x %*% est.gamma) + t.vec %*% est.gamma_t))
  
  est.p0 = mean(est.m0)
  est.p1 = mean(est.m1)
  
  risk.diff = est.p1 - est.p0
  odds.ratio = (est.p1*(1-est.p0))/(est.p0*(1-est.p1))
  
  return(list(potential.outcome = c(est.p1, est.p0), risk.diff = risk.diff, est.OR = odds.ratio, beta = est.beta, gamma = est.gamma, gamma.t = est.gamma_t, prop = est.e, m0 = est.m0, m1 = est.m1, likeli.val = opt.res$value))
}

#################################
## Stratification 
#################################
# Producing count matrix
SP.set = function(newdata, length, score){
  y = colnames(newdata)[1]
  trt_s = colnames(newdata)[2]
  ps.seq = quantile(score, prob = seq(0, 1, length.out = length))
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
  for(i in 1:(length(ps.seq)-1)){
    lower = ps.seq[i]
    upper = ps.seq[i+1]
    sbg = newdata[(score > lower & score <= upper),]
    
    if(i == 1){
      sbg = newdata[(score >= lower & score <= upper),]
    }
    
    astar.vec[i] = sum(sbg[trt_s] == 1 & sbg[y] == 1)
    bstar.vec[i] = sum(sbg[trt_s] == 1 & sbg[y] == 0)
    cstar.vec[i] = sum(sbg[trt_s] == 0 & sbg[y] == 1)
    dstar.vec[i] = sum(sbg[trt_s] == 0 & sbg[y] == 0)
  }
  n.vec = astar.vec + bstar.vec + cstar.vec + dstar.vec
  
  count.mat = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  return (count.mat)
}

# Stratification for under-reported case
SP.inference.under = function(count.mat, eta){
  eta0 = eta[1]
  eta1 = eta[2]
  
  ## The following corresponds to the cells in Table 1. 
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  
  c = c.star - (eta1/(1-eta1))*a.star
  d = d.star - (eta0/(1-eta0))*b.star
  a = c.star + a.star - c
  b = b.star + d.star - d
  
  # modify the matrix
  mat = Mod.sp(cbind(a,b,c,d))
  if(is.vector(mat)){
    a = mat[1]
    b = mat[2]
    c = mat[3]
    d = mat[4]
  }else{
    a = mat[,1]
    b = mat[,2]
    c = mat[,3]
    d = mat[,4]
  }
  
  n.vec = a+b + c+d
  
  p1 = a/(a+b)
  p0 = c/(c+d)
  
  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)
  
  risk.diff = weighted.p1 - weighted.p0
  
  return(list(risk.diff = risk.diff))
}

# Stratification for over-reported case
SP.inference.over = function(count.mat, zeta){
  zeta0 = zeta[1]
  zeta1 = zeta[2]

  ## The following corresponds to the cells in Table 1.
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control

  a = a.star - (zeta1/(1-zeta1))*c.star
  b = b.star - (zeta0/(1-zeta0))*d.star
  c = a.star + c.star - a
  d = b.star + d.star - b

  # modify the matrix
  mat = Mod.sp(cbind(a,b,c,d))
  if(is.vector(mat)){
    a = mat[1]
    b = mat[2]
    c = mat[3]
    d = mat[4]
  }else{
    a = mat[,1]
    b = mat[,2]
    c = mat[,3]
    d = mat[,4]
  }
  
  n.vec = a+b+c+d
  n1 = a+b
  n0 = c+d

  p1 = a/(a+b)
  p0 = c/(c+d)

  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)

  risk.diff = weighted.p1 - weighted.p0

  return(list(risk.diff = risk.diff))
}

#############################
## Blocking
#############################
# "smahal" and "makeblock" functions from the "blockingChallenge package"
# Calculating distance matrix
smahal <- function (X){
  units <- rownames(X)
  X <- as.matrix(X)
  s<-apply(X,2,sd)
  stopifnot(!all(s==0))
  X<-X[,s>0,drop=FALSE]
  m <- dim(X)[1]
  rownames(X) <- 1:m
  k <- dim(X)[2]
  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- cov(X)
  vuntied <- var(1:m)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  out <- matrix(NA, m, m)
  rownames(out) <- units#rownames(X)
  colnames(out) <- units
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please load package 'MASS'")
  }
  icov<-MASS::ginv(cv)
  for (i in 1:m) {
    out[i, ] <- sqrt(mahalanobis(X, X[i, ], icov, inverted = T))
  }
  return (out)
}

# Make blocks
makeblocks = function(distmat, bsize, Itr=30, .data=NULL, vars, maxItr = 200, verbose=0, ...){
  #if (!requireNamespace("optmatch", quietly = TRUE)) {
  #stop("Error: package optmatch (>= 0.9-1) not loaded. You must install optmatch first and agree to the terms of its license.")
  #}
  require(optmatch)
  stopifnot(length(bsize)==1)
  stopifnot(bsize==round(bsize))
  stopifnot(bsize>=3)
  
  if(missing(distmat)){
    stopifnot(!is.null(.data))
    stopifnot(is.data.frame(.data)|is.matrix(.data))
    stopifnot((dim(.data)[1])>(2*bsize))
    if(is.null(rownames(.data)))
      rownames(.data) = 1:nrow(.data)
    if(!missing(vars)){
      stopifnot(all(vars %in% colnames(.data)))
      
      distmat = smahal(.data[,vars])
      rownames(distmat) = rownames(.data)
      colnames(distmat) = rownames(.data)
    } 
    if(missing(vars)){
      distmat = smahal(.data)
      rownames(distmat) = rownames(.data)
      colnames(distmat) = rownames(.data)
    }
  }
  
  if(!missing(distmat)){
    stopifnot(is.data.frame(distmat)|is.matrix(distmat))
    if(dim(distmat)[1]!=dim(distmat)[2]){
      if(!is.null(.data)){
        .data = .data
      } else { .data = distmat }
      
      if(is.null(rownames(.data)))
        rownames(.data) = 1:nrow(.data)
      if(!missing(vars)){
        stopifnot(all(vars %in% colnames(.data)))
        
        distmat = smahal(.data[,vars])
        rownames(distmat) = rownames(.data)
        colnames(distmat) = rownames(.data)
      }
      if(missing(vars)){
        distmat = smahal(.data)
        rownames(distmat) = rownames(.data)
        colnames(distmat) = rownames(.data)
      }
    }
    if(is.null(rownames(distmat))) rownames(distmat) = 1:nrow(distmat)
    stopifnot((dim(distmat)[1])>(2*bsize))
  }
  
  nobs<-dim(distmat)[1] #number of observations
  nblocks<-floor(nobs/bsize) #number of blocks
  units <- rownames(distmat)
  
  if(!missing(maxItr) & !missing(Itr)) {
    if(maxItr==Itr){
      Itr4est = -1
    } else if(Itr > 30){
      Itr4est = 30
    } else Itr4est = Itr
    
    runallItr = TRUE
  } else { 
    if(missing(maxItr)){
      runallItr = FALSE
      
      if(missing(Itr)){
        Itr = 30
        Itr4est <- Itr
      } else if(Itr > 30){
        Itr4est = 30
      } else Itr4est = Itr
    } else {
      runallItr = TRUE
      
      Itr = min(maxItr, 30)
      Itr4est <- Itr
    }
  }
  
  prev.strat = NA
  prev.strat.dist = 10*bsize*nblocks*max(distmat)
  
  strat.dist <- c()
  
  
  itr <- 1
  Itr_est <- -1
  while(1){
    #for(itr in 1:Itr){
    cat("Iteration no.", itr,":: ")
    # Randomly select the first unit in each block
    if(verbose>0) cat("\n\tSelecting random starts: ")
    
    first <- sample(units, nblocks, replace=FALSE)
    
    current.strat <- matrix(NA, nblocks, bsize)
    rownames(current.strat) <- first
    
    distmat.temp1 <- distmat[first, setdiff(units, first)]
    
    if(eval(parse(text="requireNamespace('optmatch', quietly = TRUE)"))){
      om <- eval(parse(text="suppressWarnings(pairmatch(distmat.temp1, controls=bsize-1))"))
    } else {
      if(itr==1) { 
        cat("\n")
        cat("Error: package optmatch (>= 0.9-1) not loaded. You must install optmatch first and agree to the terms of its license.")
        cat("\nRunning code check mode.")
        cat("\n")
      }
      
      om <- paste0("1.", c(sample(1:length(first)), sample(rep(1:length(first), bsize-1))))
      om <- c(om, rep(NA, (length(units)%%bsize)))
      om = factor(om)
      names(om) = c(first, colnames(distmat.temp1))
    }
    
    invisible(gc())
    stopifnot(all(units %in% names(om)))
    stopifnot(length(levels(om)) == nblocks)
    
    rm(distmat.temp1)
    for(b in rownames(current.strat))
      current.strat[b,] = na.omit(names(om)[om==om[b]])
    
    sinks <- names(om)[is.na(om)]
    
    current.strat.dist = sum(sapply(first, function(b) sum(distmat[b, current.strat[b,]], na.rm=TRUE)))
    
    if(verbose>0) cat("\n\t First match created.  Average within strata distance = ", current.strat.dist/nblocks,"\n")
    if(verbose>0) cat("\tTrying local optimization.")
    
    ntry = 0
    nrogue = 1
    
    while(1){
      ##identify the furthest block unit in each block
      first.nxt.itr = first
      rogue.unit = c()
      rogue.unitid = list()
      rogue.dist = 0
      
      ntry = ntry+1
      if(verbose>1) cat("\n\t\t ",ntry, ": ( nrogue=",nrogue,")")
      for(b in rownames(current.strat)){
        unitid = order(distmat[b, current.strat[b,]], decreasing=TRUE)[1:nrogue]#which.max(distmat[b, current.strat[b,]])[1]
        names(unitid) = current.strat[b,unitid]
        if(distmat[b, tail(names(unitid),1)] == 0){
          first.nxt.itr = setdiff(first.nxt.itr, b)
        } else {
          rogue.unit = c(rogue.unit, names(unitid))
          rogue.dist = rogue.dist + sum(distmat[b, names(unitid)])
          rogue.unitid[[b]] = unitid
        }
      }
      
      if(length(first.nxt.itr)==0) break;#########
      
      rogue.unit = c(rogue.unit, sinks)# add the unmatched units in the last try to the rogue ones 
      
      #distmat.temp2 <- distmat[first.nxt.itr, rogue.unit]
      distmat.temp2 <- distmat[first.nxt.itr, rogue.unit]
      for(b in first.nxt.itr)
        distmat.temp2[b,rogue.unit] = colMeans(distmat[current.strat[b,-rogue.unitid[[b]]], rogue.unit])
      
      
      if(eval(parse(text="requireNamespace('optmatch', quietly = TRUE)"))){
        om2 <- eval(parse(text="suppressWarnings(pairmatch(distmat.temp2, controls=nrogue))"))
        
      } else { 
        om2 <- paste0("1.", c(sample(1:length(first.nxt.itr)), sample(rep(1:length(first.nxt.itr), nrogue))))
        om2 <- c(om2, rep(NA, (sum(dim(distmat.temp2))-length(first.nxt.itr)*(1+nrogue)) ))
        om2 = factor(om2)
        names(om2) = c(first.nxt.itr, colnames(distmat.temp2))
      }
      
      invisible(gc())
      stopifnot(names(om2)[1:length(first.nxt.itr)] == first.nxt.itr)
      stopifnot(length(levels(om2)) == length(first.nxt.itr))
      rm(distmat.temp2)
      
      new.dist = 0
      matched.unit = list()
      for(b in first.nxt.itr){
        matched.unit[[b]] = names(na.omit(om2[om2==om2[b]]))[-1]
        new.dist = new.dist + sum(distmat[b, matched.unit[[b]]])
      }
      
      if(verbose>0) cat(" ", new.dist/nblocks, rogue.dist/nblocks)
      if(new.dist < rogue.dist){
        cat(" . ")
        for(b in first.nxt.itr)
          current.strat[b, rogue.unitid[[b]]] = matched.unit[[b]]
        sinks = names(om2[is.na(om2)])
        current.strat.dist = current.strat.dist - rogue.dist + new.dist
      }
      
      if(new.dist >= rogue.dist){
        if(nrogue >= 2){ 
          nrogue = 1
          break;
        } else{ nrogue = nrogue+1 }
      }
    }
    
    if(verbose>0) cat("\n\tFinal  average within block distance of iteration", itr,"::",current.strat.dist/nblocks,".\n")
    strat.dist <- c(strat.dist,current.strat.dist)
    
    #print(prev.strat.dist)
    #print(current.strat.dist)
    if(current.strat.dist < prev.strat.dist){
      prev.strat.dist = current.strat.dist
      #print(prev.strat.dist)
      prev.strat = current.strat
    }
    cat("best average within block distance::", prev.strat.dist/nblocks,"\n") 
    
    
    if(itr == Itr4est){
      C_K <- 1 + (max(strat.dist)/min(strat.dist)-1)/(2*(max(strat.dist)/min(strat.dist)))^2
      pC <- mean(strat.dist <= C_K*min(strat.dist))
      pC <- pC/2
      eps <- 0.1
      Itr_est <- ceiling( (eps^(-1)-1)*(pC^(-1)-1) )
      #print(Itr_est )
    }
    itr <- itr + 1
    
    if(itr > Itr & Itr_est > Itr){
      Itr <- Itr_est
      cat("\nRunning", Itr-(itr-1), "more iteration(s).\n\n")
      if(Itr4est<30){ 
        Itr4est = ifelse(Itr>30, 30, Itr)
      } else if(Itr4est>=30) Itr_est = Itr-1
    } else if(itr > Itr) break;
  }
  
  cat("\n")
  ## summary F statistics
  bk = rep(NA, nobs)
  names(bk) = units
  for(b in 1:nrow(prev.strat))
    bk[prev.strat[b,]] = b
  
  unblock = sum(is.na(bk))
  if (unblock>0) bk[is.na(bk)]<-(-(1:unblock))
  
  who<-outer(bk,bk,"==")
  inside<-sum(distmat[who])
  outside<-sum(distmat[!who])
  
  ndist<-bsize*(bsize-1)*nblocks
  
  #print(c(mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist))))
  #return(strat.dist)
  
  if(!is.null(.data)){
    .data = cbind(.data, strata=NA)
    .data[names(bk),'strata'] = bk
    
    cat("A new column, strata, is added to the data with strata ids.\n")
    return(list(strata = prev.strat, .data = .data)) 
    ##mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist)) ))
    
  }
  
  if(is.null(.data))
    return(list(strata = prev.strat))##, mean.inside=inside/ndist, mean.outside = outside/(nobs*nobs-(nobs+ndist))))
  
}

dist_ftn <- function(a, b) {sqrt(sum((a - b)^2))}

# Producing block matrix
Block.set <- function(newdata, bsize){
  y = colnames(newdata)[1]
  trt_s = colnames(newdata)[2]
  X = newdata[,-c(1,2)]
  X = data.frame(scale(X)[1:length(X[,1]),1:length(X[1,])])
  newdata[,-c(1,2)] = X
  distmat1 = smahal(X)
  res = makeblocks(distmat1, bsize=bsize, data=X)
  strnum = length(newdata[,1])/bsize
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, strnum)
  X1 = matrix(NA, ncol = length(X[1,]), nrow = strnum)
  for (i in 1:strnum){
    sbg = newdata[as.integer(res$strata[i,]),]
    astar.vec[i] = sum(sbg[trt_s] == 1 & sbg[y] == 1)
    bstar.vec[i] = sum(sbg[trt_s] == 1 & sbg[y] == 0)
    cstar.vec[i] = sum(sbg[trt_s] == 0 & sbg[y] == 1)
    dstar.vec[i] = sum(sbg[trt_s] == 0 & sbg[y] == 0)
    X1[i,] = apply(sbg[,-c(1,2)], 2, mean)
  }
  n.vec = astar.vec + bstar.vec + cstar.vec + dstar.vec
  block.mat = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec, X1, n.vec)
  return (block.mat)
}

# Blocking for under-reported case
Block.inference.under = function(block.mat, eta){
  eta0 = eta[1]
  eta1 = eta[2]
  c = block.mat[,3] - (eta1/(1-eta1))*block.mat[,1]
  d = block.mat[,4] - (eta0/(1-eta0))*block.mat[,2]
  a = block.mat[,3] + block.mat[,1] - c
  b = block.mat[,2] + block.mat[,4] - d
  block.mat[,1:4] = cbind(a,b,c,d)
  
  block.mat = Mod.block(block.mat)
  
  if(is.vector(block.mat)){
    a = block.mat[1]
    b = block.mat[2]
    c = block.mat[3]
    d = block.mat[4]
  }else{
    a = block.mat[,1]
    b = block.mat[,2]
    c = block.mat[,3]
    d = block.mat[,4]
  }
  
  n.vec = a+b+c+d
  
  p1 = a/(a+b)
  p0 = c/(c+d)
  
  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)
  
  risk.diff = weighted.p1 - weighted.p0
  
  return(list(risk.diff = risk.diff))
}

# Blocking for over-reported case
Block.inference.over = function(block.mat, zeta){
  zeta0 = zeta[1]
  zeta1 = zeta[2]
  
  a = block.mat[,1] - (zeta1/(1-zeta1))*block.mat[,3]
  b = block.mat[,2] - (zeta0/(1-zeta0))*block.mat[,4]
  c = block.mat[,1] + block.mat[,3] - a
  d = block.mat[,2] + block.mat[,4] - b
  block.mat[,1:4] = cbind(a,b,c,d)
  
  block.mat = Mod.block(block.mat)
  
  if(is.vector(block.mat)){
    a = block.mat[1]
    b = block.mat[2]
    c = block.mat[3]
    d = block.mat[4]
  }else{
    a = block.mat[,1]
    b = block.mat[,2]
    c = block.mat[,3]
    d = block.mat[,4]
  }
  
  n.vec = a+b+c+d
  
  p1 = a/(a+b)
  p0 = c/(c+d)
  
  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)
  
  risk.diff = weighted.p1 - weighted.p0
  
  return(list(risk.diff = risk.diff))
}

#################################
## Nearest Neighbor Combination
################################# 
# Check if there is any problem
is_prob = function(vec){
  val = FALSE
  if (vec[1]+vec[2]<=0){val = TRUE}
  if (vec[3]+vec[4]<=0){val = TRUE}
  if (vec[1]<0||vec[2]<0||vec[3]<0||vec[4]<0){val = TRUE}
  return (val)
}

# Modify count matrix in propensity and prognostic score stratifications
Mod.sp = function(mat){
  while (TRUE){
    if(is.vector(mat)) {break;}
    l = length(mat[,1])
    ind = TRUE
    for (i in 1:(l-1)){
      if(is_prob(mat[i,])){
        ind = FALSE
      }
    }
    if(ind == TRUE){
      if(is_prob(mat[l,])){
        mat[l-1,] = mat[l-1,]+mat[l,]
        mat = mat[-l,]
        next;
      }else{
        break;
      }
    }
    for (i in 1:(l-1)){
      if(is_prob(mat[i,])){
        mat[i,] = mat[i,]+mat[i+1,]
        mat = mat[-(i+1),]
        break;
      }
    }
  }
  return(mat)
}

# Modify block matrix in blocking method
Mod.block = function(block.mat){
  while(TRUE){
    if(is.vector(block.mat)) {break;}
    l = length(block.mat[,1])
    ind = TRUE
    loc = 0
    for (i in 1:l){
      if(is_prob(block.mat[i,1:4])){
        loc = i
        ind = FALSE
        break;
      }
    }
    if (ind){break;}
    num = ifelse(loc==l,1,loc+1)
    dis = dist_ftn(block.mat[loc,1:4],block.mat[num,1:4])
    for (i in 1:l){
      if(l==loc){next;}
      if(dist_ftn(block.mat[loc,1:4],block.mat[l,1:4])<dis){
        num = l
        dis = dist_ftn(block.mat[loc,1:4],block.mat[l,1:4])
      }
    }
    block.mat[num,1:4] = block.mat[num,1:4]+block.mat[loc,1:4]
    l = length(block.mat[1,])
    block.size = block.mat[,l]
    for(i in 5:(l-1)){
      block.mat[num,i] = (block.mat[num,i]*block.mat[num,l]+block.mat[loc,i]*block.mat[loc,l])/(block.mat[num,l]+block.mat[loc,l])
    }
    block.mat[num,l] = block.mat[num,l]+block.mat[loc,l]
    block.mat = block.mat[-loc,]
  }
  
  return (block.mat)
}
