
## time-varying coef DistCox (ODACT)

# Distributed time-varying Cox regression
# C. Jason Liang
# 8/14/2019


# Local constant partial likelihood
# Cai and Sun (2003) propose maximizing a local constant partial likelihood to estimate the time-dependent log-hazard ratio β(t) from a Cox regression model with time-dependent β(t):
# 
#   λ(t|x)=λ0(t)exp(β(t)x)
# Note that the “local” in “local constant partial likelihood” refers to “local in time” rather than a “local site” in the context of a distributed algorithm.
# 
# The local constant log partial likelihood is defined as
# 
# L(β,t)=(nhn)−1∑i=1nK(yi−thn)δi[βTxi−log{∑j∈R(ti)exp(βTxj)}],
# where solving for β provides an estimate of β(t). The K(⋅) function is a kernel, such as the Epanechnikov kernel, centered around t with data-adaptive bandwidth hn.
# 
# Note that L(β,t) is similar to the usual partial log likelihood L(β), with the main difference being that each likelihood contribution is now weighted by a kernel centered around t. It then follows that the gradient and hessian of L(β,t) can be written as locally weighted versions of ∇L(β) and ∇2L(β):
# 
#   ∇L(β,t)∇2L(β,t)=(nhn)−1∑i=1nK(yi−thn)δi{xi−∑j∈R(ti)exp(βTxj)xj∑j∈R(ti)exp(βTxj)}=(nhn)−1∑i=1nK(yi−thn)δi⎡⎣⎢⎢{∑j∈R(ti)exp(βTxj)xj}⨂2−∑j∈R(ti)exp(βTxj)∑j∈R(ti)exp(βTxj)x⨂2j{∑j∈R(ti)exp(βTxj)}2⎤⎦⎥⎥
# ## Functions for local in time estimation:


require(RcppArmadillo)
require(Rcpp)
sourceCpp('utils.cpp')

## DO NOT use order.time=T, always order time before calling 
## llpl llplg llplh sllpl sllplg D12

# log-lik
llpl <- function(beta, times, status, covars, tt=0.5, h=500,
                 order.time=F) {  
  # locally constant log partial likelihood
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  if (is.null(ncol(covars))) 
    lp <- covars * beta
  else lp <- covars %*% beta
  
  # -sum((K * status * (lp - log(cumsum(exp(lp)[N:1])[N:1])))[nn])
  ## cum_sum is rcpp function, ~ 20 times faster
  -sum((K * status * (lp - log(cum_sum(exp(lp), reversely = T))))[nn])
}

# gradient of ll
llplg <- function(beta, times, status, covars, tt=0.5, h=500,
                  order.time=F) {
  # locally constant log partial likelihood gradient
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  if (is.null(ncol(covars))) 
    lp <- covars * beta
  else lp <- covars %*% beta
  elp <- c(exp(lp))
  # num <- apply(elp*covars, 2, function(x) rev(cumsum(rev(x))))
  # den <- cumsum(elp[N:1])[N:1]   # %o% rep(1, ncol(covars))
  ## cum_sum_col is rcpp function, ~ 20 times faster
  num <- cum_sum_cols(elp*covars, reversely = T)
  den <- c(cum_sum(elp, reversely = T))
  
  -colSums(((covars - num/den) * (K * status))[nn, , drop = FALSE])
}

# hessian of ll
llplh <- function(beta, times, status, covars, tt=0.5, h=500,
                  order.time=F) {
  # locally constant log partial likelihood gradient
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  nn <- 1 - ((times - tt)/h)^2 >= 0
  K <- (1 - ((times - tt)/h)^2) * 0.75/h
  N <- length(times)
  px = ncol(covars)
  if (is.null(ncol(covars))){
    lp <- covars * beta
  }else{
    lp <- covars %*% beta
  }
  
  elp <- c(exp(lp))
  
  # den <- ((cumsum(elp[N:1])[N:1])^2)
  # num1 <- apply(elp*covars, 2, function(x) rev(cumsum(rev(x))))
  # num1 = t(apply(num1, 1, function(x) c(x %o% x ) ))   # matrix(, nrow=N)
  # num2 = elp * t(apply(covars, 1, function(x) c(x %o% x) ) )
  # num2 <- sqrt(den) * apply(num2, 2, function(x) rev(cumsum(rev(x))) )
  
  ## cum_sum_col Xotimes2 are rcpp function, ~ 10 times faster
  den <- c(cum_sum(elp, reversely = T)^2)
  num1 <- cum_sum_cols(elp*covars, reversely = T)
  num1 <- Xotimes2(num1)
  num2 <- elp * Xotimes2(covars)
  num2 <- sqrt(den) * cum_sum_cols(num2, reversely = T)
  
  hh = -colSums((((num1-num2)/den) * (K * status))[nn, , drop = FALSE])
  hh = matrix(hh, px, px)
  
  return(hh)
}

# surrogate ll 
sllpl <- function(beta, times, status, covars, betabar, ind1, order=1, tt=0.5, h=.25,
                  order.time=F) {
  # surrogate locally constant log partial likelihood
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  if(order==1){
    # first order version (set 'order=1')
    llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
      sum(beta * (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) -
                    llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])))
  } else{
    # second order version (set 'order=2')
    llpl(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
      sum(beta * (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
                    llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]))) +
      0.5 * t(beta - betabar) %*% (llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
                                     llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
    
  }
}

# gradient of sll
sllplg <- function(beta, times, status, covars, betabar, ind1, order=1, tt=0.5, h=.25,
                   order.time=F) {
  # surrogate locally constant log partial likelihood
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  if(order==1){
    # first order version (set 'order=1')
    llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
      (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) -
         llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) 
  } else{
    # second order version (set 'order=2')
    llplg(beta, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) +
      (llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
         llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) +
      0.5 * (llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times) - 
               llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1])) %*% (beta - betabar)
    
  }
}

# 1st and 2nd order derivatives of llpl, using local or pooled data
D12 <- function(times, status, covars, betabar, ind1,  tt=0.5, h=.25,
                order.time=F) {
  # output derivatives at betabar
  if(order.time==T){
    times = times[order(times)]
    status = status[order(times)]
    covars = covars[order(times), ,drop=FALSE]
  }
  
  # second order version (set 'order=2')
  D1.N = llplg(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times)  
  D1.1 = llplg(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) 
  D2.N = llplh(betabar, times=times, status=status, covars=covars, tt=tt, h=h)/length(times)  
  D2.1 = llplh(betabar, times=times[ind1], status=status[ind1], covars=covars[ind1, , drop=FALSE], tt=tt, h=h)/length(times[ind1]) 
  
  return(list(D1.N=D1.N, D1.1=D1.1, D2.N=D2.N, D2.1=D2.1))
}
 
seq_fit <- function(da, fn=llpl, times=seq(0,1,0.1), h=0.1, betabar, ...){
  # wrapper function to make it easier to call optim() multiple times
  apply(as.matrix(times), 1, function(x)
    optim(par=betabar, fn=fn, method="Nelder-Mead", times=da[,1], status=da[,2], covars=da[,-c(1:2),drop=FALSE], 
          tt=x, h=h, ...)$par)
}


seq_fit_list <- function(da, fn=llpl, times=seq(0,1,0.1), h=0.1, betabar, ...){
  # wrapper function to make it easier to call optim() multiple times
  lapply(as.list(times), function(x)
    optim(par=betabar, fn=fn, method="Nelder-Mead", times=da[,1], status=da[,2], covars=da[,-c(1:2),drop=FALSE], 
          tt=x, h=h, ...) )
}


# generate time-varying cox data
# lambda(t|m) = exp(tm+m) w/ baseline haz = 1
# M ~ Unif(0,1)
Stm <- function(t,m,c0=-1){
  #S(t|m)
  #exp((1/m) - (exp(t*m)/m))
  # exp(exp(m)/m - exp(t*m + m)/m)
  exp(exp(m*c0)/m*(1- exp(t*m)) )
} 

# generate random times given m
Stm_inv <- function(zo, m,c0=-1){
  # helper function for generating random times given m
  #zo for random number from Unif(0, 1)
  #  log(1 - m * log(zo))/m
  (log(exp(m*c0) - m*log(zo)) - m*c0)/m
} 

ftm <- function(t,m){
  # spot check to make sure random numbers are valid
  # f(t, m)
  #  if (m < 1 | m > 3) return(0)
  #  (1/2) * exp((1/m) - (exp(t*m)/m) + t*m)
  #(1/2) * exp(exp(m)/m - exp(t*m + m)/m + t*m + m)
  exp(exp(m)/m - exp(t*m + m)/m + t*m + m)
}

ft <- function(t){
  # f(t), numerically integrated
  # again, just spot checking
  ms <- seq(0,1,.001)
  sum(ftm(rep(t, length(ms)), ms))*.001
}

coxlcgen <- function(n=10000, cens=4.5, c0=-1){
  # generate data from Cox model with time-varying coefficient
  # higher 'cens' values means more censoring
  # default of 4.5 translates to roughly 80% censoring, or 20% events
  mo <- runif(n, 0, 1) # covariate values
  m2 <- runif(n, 0, 1) # another covariate with beta(t)==0
  to <- Stm_inv(runif(n), mo, c0=c0) # event times
  mc <- runif(n, 0, cens) # covariate values for generating censoring times
  # co <- Stm_inv(runif(n), mc, c0=c0) # censoring times
  # co = mc
  # plot(density(rweibull(100, shape=2, scale=15)))
  # co = rweibull(n, shape=cens, scale=3)
  co = rbeta(n, 2, cens)*5
  do <- 1*(to < co) # event indicator (1: event; 0: censor)
  yo <- co; 
  yo[to<co] <- to[to<co] # min(event time, censor time)
  
  mo <- mo[order(yo)]
  do <- do[order(yo)]
  # to = to[order(yo)]
  # co = co[order(yo)]  
  yo <- yo[order(yo)]
  da <- cbind(yo, do, mo, m2) # , to=to, co=co
  return(da)
}



do.one.lc <- function(n=10000, K=10, n.site=rep(n/K, K), cens=10, c0=-1,
                      evalt = seq(0, 1, 0.2), h=0.2, init='local'){  # cens smaller, more events
  # estimate the 
  # 1) pooled estimator, 
  # 2) local estimators, 
  # 3) mean of local estimators,
  # 4) 1st order surrogate estimator, and 
  # 5) 2nd order surrogate estimator
  
  px <- 2
  da <- coxlcgen(n=n, cens=cens, c0=c0) 
  id.site = sample(rep(1:K, n.site))
  inddist <- lapply(1:K, function(a)  which(id.site==a))
  
  # gold <- seq_fit(da=da, fn=llpl, times=evalt, h=h, betabar = c(0,0))
  gold.list <- seq_fit_list(da=da, fn=llpl, times=evalt, h=h, betabar = rep(0,px), hessian=T)
  gold <- matrix(unlist(lapply(gold.list, function(a) a$par) ), nrow=px)
  gold.sd <- matrix(unlist(lapply(gold.list, function(a) sqrt(diag(solve(a$hessian))*0.6/h)) ), nrow=px)
  
  locals.list <- lapply(inddist, function(x){
    fit <- seq_fit_list(da=da[x,], fn=llpl, times=evalt, h=h, betabar = rep(0,px), hessian=T)
  }) 
  locals <- sapply(locals.list, function(a){sapply(a, function(a) a$par) }, simplify="array") 
  meta = apply(locals, c(1,2), function(a) sum(a*n.site)/sum(n.site) )
  
  b.init = apply(as.matrix(evalt), 1, function(x){
    if(init=='true'){
      b.init <- c(1+x,0)
    }else if(init=='local'){
      b.init <- locals[,findInterval(x, evalt),1]
    } else if(init=='meta'){
      b.init <- meta[,findInterval(x, evalt)]
    } else if(init=='zero'){
      b.init <- rep(0, px)
    }
  })
 
  sl2 = array(NA,c(px, length(evalt), K))
  for(k in 1:K){ 
    sl2.list <- lapply(as.list(evalt), function(x)
      optim(par= b.init[,findInterval(x, evalt)], 
            fn=sllpl, method="Nelder-Mead", # hessian=T,
            times=da[,1], status=da[,2], covars=da[,-c(1,2),drop=FALSE],
            tt=x, h=h, betabar=b.init[,findInterval(x, evalt)], order=2, ind1=inddist[[k]]))
    sl2[,,k] <- matrix(unlist(lapply(sl2.list, function(a) a$par) ), nrow=px)
  }
  sl2avg = apply(sl2,1:2,mean)
  
  cat(date(), "\n")
  return(list(evalt=evalt, gold=gold, locals=locals, meta=meta, sl2=sl2, sl2avg=sl2avg))
}
 