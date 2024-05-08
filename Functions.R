##############################################################
#### Load Packages ##################################

library(DEoptim)
library(VineCopula)
library(rvinecopulib)
library(vinereg)
library(parallel)
library(foreach)
library(CopReg)

##############################################################


##############################################################
#### Conditional densities
##############################################################
## This function get the conditional density f(xvec|L = U(1/p))
## Inputs: xvec: vector of matrix of quantiles
##         RVM: An RVineMatrix() object of (X,L) including the 
##              structure and the pair-copula families and parameters
##         Up: U(1/p) of L
##         p: probability level of upper quantile
##         skt_result: fitted parameters of skew-t marginals of X
##         set: vector of components of X where marginals are fitted
##              using semi-parametric method
##         set_margin: fitted parameters of the marginals for components
##                     in set
##         PnL_margin: fitted parameters of the marginal for L
################################################################
con_den <- function(xvec, RVM, Up=NULL, p=NULL, skt_result, 
                    set, set_margin, PnL_margin){
  d = length(RVM$names)
  
  if(is.null(p)){
    if(is.numeric(PnL_margin))
      p = 1-pstjf(Up, PnL_margin)
    if(is.list(PnL_margin))
      p = 1-pmargin(x=Up, method = "edfun", ul = PnL_margin$ul,
                    sigmal = PnL_margin$sigmal, 
                    xil = PnL_margin$xil, 
                    ur = PnL_margin$ur, 
                    sigmar = PnL_margin$sigmar, 
                    xir = PnL_margin$xir, 
                    pfun = PnL_margin$data.dist$pfun)
  }
  
  xvec = matrix(xvec, ncol = (d-1))
  uvec = xvec
  dvec = xvec
  for(i in 1:(d-1)){
    if(!rownames(skt_result)[i] %in% set){
      uvec[,i] = pstjf(xvec[,i], skt_result[i,])
      dvec[,i] = dstjf((xvec[,i] - skt_result[i,1])/skt_result[i,2],
                       a = skt_result[i,3],
                       b = skt_result[i,4])/skt_result[i,2]
    }
    else{
      j = which(set == rownames(skt_result)[i])
      uvec[,i] = pmargin(x = xvec[,i], method = "edfun", 
                         ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                         xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                         sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                         pfun = set_margin[[j]]$data.dist$pfun)
      dvec[,i] = dmargin(xvec[,i], method = "edfun", 
                         ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                         xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                         sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                         dfun = set_margin[[j]]$data.dist$dfun, 
                         pfun = set_margin[[j]]$data.dist$pfun)
    }
  }
  
  prod_margin = apply(log(dvec), MARGIN = 1, FUN = sum)
  uvec = cbind(rep(1-p, nrow(uvec)), uvec)
  re = log(RVinePDF(uvec, RVM)) + prod_margin ## at lof scale
  return(re)
}



##############################################################
## This function get the conditional density f(xvec|L >= U(1/p))
## Inputs: xvec: vector of matrix of quantiles
##         RVM: An RVineMatrix() object of (X,L) including the 
##              structure and the pair-copula families and parameters
##         Up: U(1/p) of L
##         p: probability level of upper quantile
##         skt_result: fitted parameters of skew-t marginals of X
##         set: vector of components where marginals are fitted using
##              semi-parametric method
##         set: vector of components of X where marginals are fitted
##              using semi-parametric method
##         PnL_margin: fitted parameters of the marginal for L
################################################################
con_den_gt <- function(xvec, RVM_X, RVM_all, p=NULL, Up=NULL, 
                       skt_result, set, set_margin, PnL_margin){
  d = length(RVM_all$names)
  xvec = matrix(xvec, ncol = (d-1))
  
  if(is.null(Up)) {
    if(is.numeric(PnL_margin))
      Up = qstjf(1-p, PnL_margin)
    if(is.list(PnL_margin))
      Up = qmargin(u = 1-p, method = "edfun", ul = PnL_margin$ul,
                   sigmal = PnL_margin$sigmal,
                   xil = PnL_margin$xil,
                   ur = PnL_margin$ur, 
                   sigmar = PnL_margin$sigmar, 
                   xir = PnL_margin$xir,
                   qfun = PnL_margin$data.dist$qfun,  
                   pfun = PnL_margin$data.dist$pfun)
  }
  
  if(is.null(p)){
    if(is.numeric(PnL_margin))
      p = 1-pstjf(Up, PnL_margin)
    if(is.list(PnL_margin))
      p = 1-pmargin(x=Up, method = "edfun", ul = PnL_margin$ul,
                    sigmal = PnL_margin$sigmal, 
                    xil = PnL_margin$xil, 
                    ur = PnL_margin$ur, 
                    sigmar = PnL_margin$sigmar, 
                    xir = PnL_margin$xir, 
                    pfun = PnL_margin$data.dist$pfun)
  }
  
  uvec = xvec
  dvec = xvec
  for(i in 1:(d-1)){
    if(!rownames(skt_result)[i] %in% set){
      uvec[,i] = pstjf(xvec[,i], skt_result[i,])
      dvec[,i] = dstjf((xvec[,i] - skt_result[i,1])/skt_result[i,2],
                       a = skt_result[i,3],
                       b = skt_result[i,4])/skt_result[i,2]
    }
    else{
      j = which(set == rownames(skt_result)[i])
      uvec[,i] = pmargin(x = xvec[,i], method = "edfun", 
                         ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                         xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                         sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                         pfun = set_margin[[j]]$data.dist$pfun)
      dvec[,i] = dmargin(xvec[,i], method = "edfun", 
                         ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                         xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                         sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                         dfun = set_margin[[j]]$data.dist$dfun, 
                         pfun = set_margin[[j]]$data.dist$pfun)
    }
  }
  xvec = cbind(rep(Up, nrow(xvec)), xvec)
  uvec = cbind(rep(1-p, nrow(uvec)), uvec)
  
  joint_x = log(RVinePDF(uvec[,-1], RVM_X)) + rowSums(log(dvec))
  
  uvec = uvec[,c(2:d,1)]  ## put L to the last column
  M = RVM_all[["Matrix"]][d:1,d:1]
  fam = RVM_all[["family"]][d:1,d:1]
  par = RVM_all[["par"]][d:1,d:1]
  par2 = RVM_all[["par2"]][d:1,d:1]
  perm = diag(M)
  for(j in 1:d){
    M[which(RVM_all[["Matrix"]][d:1,d:1] == perm[j])] = j
  }
  uvec = matrix(uvec, ncol = d)
  con_prob = RVinePCond(udat = matrix(uvec[,perm], ncol = d), A=M, ntrunc = (d-1), fam, par, par2)[,d]
  return(joint_x + log(1 - con_prob)) ## at lof scale
}



##############################################################
## This function get the joint density f(x)
## Inputs: x: vector of matrix of quantiles
##         RVM: An RVineMatrix() object of X including the structure 
##              and the pair-copula families and parameters
##         skt.par: fitted parameters of skew-t marginals of X
##         set: vector of components of X where marginals are fitted
##              using semi-parametric method
##         set_margin: fitted parameters of the marginals for components
##                     in set
################################################################
joint_den_x = function(x, RVM, skt.par, set, set_margin){
  d = length(RVM$names)
  x = matrix(x, ncol = d)
  uvec = dmargin = x
  for(i in 1:d){
    if(!rownames(skt.par)[i] %in% set){
      uvec[,i] = pstjf(x[,i], skt.par[i,])
      dmargin[,i] = dstjf((x[,i] - skt.par[i,1])/skt.par[i,2],
                          a = skt.par[i,3],
                          b = skt.par[i,4])/skt.par[i,2]
    }
    else{
      j = which(set == rownames(skt.par)[i])
      uvec[,i] = pmargin(x = x[,i], method = "edfun", 
                         ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                         xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                         sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                         pfun = set_margin[[j]]$data.dist$pfun)
      dmargin[,i] = dmargin(x[,i], method = "edfun", 
                            ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                            xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                            sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                            dfun = set_margin[[j]]$data.dist$dfun, 
                            pfun = set_margin[[j]]$data.dist$pfun)
    }
  }
  prod_margin = apply(dmargin, MARGIN = 1, FUN = prod)
  re = RVinePDF(uvec, RVM)*prod_margin
  return(re)
}


##############################################################
## This function get the conditional density f(x)*if{gest(x)>=thred}
## Inputs: x: vector of matrix of quantiles
##         thred: threshold of L
##         RVM: An RVineMatrix() object of (X,L) including the 
##              structure and the pair-copula families and parameters
##         RVM_X: An RVineMatrix() object of X including the structure 
##              and the pair-copula families and parameters
##         skt.par: fitted parameters of skew-t marginals of X
##         set: vector of components where marginals are fitted using
##              semi-parametric method
##         set_margin: fitted parameters of the marginals for components
##                     in set
##         PnL_margin: fitted parameters of the marginal for L
################################################################
con_den_g_est <- function(x, thred, RVM, RVM_X, skt.par, 
                          set, set_margin, PnL_margin){
  if(is.numeric(PnL_margin))
    L_qfunc = function(u){
      qstjf(p = u, param = PnL_margin)
    }
  if(is.list(PnL_margin))
    L_qfunc = function(u){
      qmargin(u, method = "edfun", ul = PnL_margin$ul, 
                 sigmal = PnL_margin$sigmal,xil = PnL_margin$xil,
                 ur = PnL_margin$ur, sigmar = PnL_margin$sigmar, 
                 xir = PnL_margin$xir,
                 qfun = PnL_margin$data.dist$qfun, 
                 pfun = PnL_margin$data.dist$pfun)
    }
  u = x
  for(i in 1:length(x)){
    if(!rownames(skt.par)[i] %in% set){
      u[i] = pstjf(x[i], param = skt.par[i,])
    }
    else{
      j = which(set == rownames(skt.par)[i])
      u[i] = pmargin(x = x[i], method = "edfun", 
                     ul = set_margin[[j]]$ul, sigmal = set_margin[[j]]$sigmal,
                     xil = set_margin[[j]]$xil, ur = set_margin[[j]]$ur, 
                     sigmar = set_margin[[j]]$sigmar, xir = set_margin[[j]]$xir,
                     pfun = set_margin[[j]]$data.dist$pfun)
    }
  }
  u = matrix(u, nrow = 1)
  colnames(u) = rownames(skt.par)
  pred_L = predict(RVM, newdata = u, y_marginal_qfunc = L_qfunc, mean = F)
  if(pred_L < thred) f = 0
  if(pred_L >= thred) {
    f = joint_den_x(x, RVM_X, skt.par, set, set_margin)
  }
  return(f)
}

##############################################################
## This function get the conditional density f(x)*if{g(x)>=thred}
## Inputs: x: vector of matrix of quantiles
##         thred: threshold of L
##         RVM: An RVineMatrix() object of X including the structure 
##              and the pair-copula families and parameters
##         skt.par: fitted parameters of skew-t marginals of X
##         set: vector of components where marginals are fitted using
##              semi-parametric method
##         set_margin: fitted parameters of the marginals for components
##                     in set
##         weight: weight vector in function g
################################################################
con_den_known <- function(x, thred, RVM, skt.par = NULL, set, 
                          set_margin, weight){
  if(sum(weight*x) < thred) f = 0
  if(sum(weight*x) >= thred) {
    f = joint_den_x(x, RVM_X, skt.par, set, set_margin)
  }
  return(f)
}




##############################################################
#### Marginal fitting (skew-t distribution)
##############################################################




##############################################################
#### Marginal fitting (semi-parametric model)
#### Only used in Applications
##############################################################


##############################################################
## This function estimate parameter xi in GPD dist
## Inputs: x: vector of data
##         thred: threshold used to estimate GPD dist
##         sigma: value of parameter sigma in GPD dist
##############################################################
fitxi <- function(x, thred, sigma){
  xsub = x[x > thred]
  z = (xsub - thred)/sigma
  llk = function(par){
    f = (1+par*z)^(-1/par - 1)/sigma
    return(sum(log(f)))
  }
  xi = optimize(llk, interval = c(0,10), maximum = TRUE)$maximum
  return(xi)
}


##############################################################
## This function estimate KDE for central data
## Inputs: data: vector of data
##############################################################
fit_kde <- function(data){
  h = bw.nrd0(data)
  dfun = function(x){
    f = x
    for(i in 1:length(x)){
      f[i] = mean(dnorm(x[i], mean = data, sd = h))
    }
    return(f)
  }
  pfun = function(x){
    P = x
    for(i in 1:length(x)){
      P[i] = mean(pnorm(x[i], mean = data, sd = h))
    }
    return(P)
  }
  qfun = function(x){
    q = x
    for(i in 1:length(x)){
      sfun = function(y){(pfun(y) - x[i])^2}
      q[i] = optimize(sfun, interval = 2*range(data),
                      tol = .Machine$double.eps)$minimum
    }
    return(q)
  }
  return(list(dfun = dfun, pfun = pfun, qfun = qfun))
}


##############################################################
## This function fit semi-parametric marginal to data
## Inputs: data: vector of data
##         method: method used to fit the central data:
##                 "edfun": KDE method
##                 "skt": skew-t method
##         ul: upper threshold to estimate GPD
##         ur: lower threshold to estimate GPD
##         con: add constraint to remove jump point or not
##############################################################
fitmix <- function(data, method, ul = NULL, ur = NULL, con = TRUE){
  if(is.null(ul)) ul = quantile(data,0.1)
  if(is.null(ur)) ur = quantile(data,0.9)
  
  if(!con){
    ## right tail
    par.mle = ismev::gpd.fit(data[data>0], threshold = ur, show = FALSE)$mle
    sigmar = par.mle[1]
    xir = par.mle[2]
    
    ## left tail
    par.mle = ismev::gpd.fit(abs(data[data<=0]), threshold = abs(ul), 
                             show = FALSE)$mle
    sigmal = par.mle[1]
    xil = par.mle[2]
  }
  
  if(method == "edfun"){
    data.dist = fit_kde(data)
    skt.par = NULL
    if(con){
      sigmar = (1-data.dist$pfun(ur))/data.dist$dfun(ur)
      xir = fitxi(data[data>0], ur, sigmar)
      sigmal = data.dist$pfun(ul)/data.dist$dfun(ul)
      xil = fitxi(abs(data[data<=0]), abs(ul), sigmal)
    }
  }
  if(method == "skt"){
    set.seed(100)
    skt.par = DEoptim(nllkstjf, lower = c(10*min(data),0,0,0), 
                      upper = c(10*max(data),10*sd(data),10,10), 
                      ydat = data, 
                      control = DEoptim.control(trace = F))$optim$bestmem
    data.dist = NULL
    if(con){
      sigmar = (1-pstjf(ur, skt.par))/dstjf((ur - skt.par[1])/skt.par[2], 
                                            a = skt.par[3], 
                                            b =skt.par[4])*skt.par[2]
      xir = fitxi(data[data>0], ur, sigmar)
      sigmal = pstjf(ul, skt.par)/dstjf((ul - skt.par[1])/skt.par[2], 
                                        a = skt.par[3], 
                                        b =skt.par[4])*skt.par[2]
      xil = fitxi(abs(data[data<=0]), abs(ul), sigmal)
    }
  }
  
  return(list(ul = ul, sigmal = sigmal, xil = xil, ur = ur, sigmar = sigmar,
              xir = xir, data.dist = data.dist, skt.par = skt.par))
}


##############################################################
## Thes functions are density, probability and quantile functions
## for a fitted object from fitmix()
## Inputs: x: vector of quantile
##         u: vector of probabilities
##         method: method used to fit the central data:
##                 "edfun": KDE method
##                 "skt": skew-t method
##         ul: upper threshold to estimate GPD
##         sigmal: fitted sigma for left tail
##         sigmar: fitted sigma for right tail
##         xil: fitted xi for left tail
##         xir: fitted xi for right tail
##         dfun: density function for central part data
##         pfun: probability function for central part data
##         qfun: quantile function for central part data
##         skt.par: fitted skew-t parameters for central part data
##############################################################
dmargin <- function(x, method, ul, sigmal, xil, ur, sigmar, xir, 
                    dfun = NULL, pfun = NULL, skt.par = NULL){
  if(method == "edfun") {
    f = dfun(x)
    phi = pfun(c(ul,ur))
  }
  if(method == "skt") {
    f = dstjf((x - skt.par[1])/skt.par[2], a = skt.par[3],
              b = skt.par[4])/skt.par[2]
    phi = pstjf(x = c(ul,ur), param = skt.par)
  }
  
  ## right tail
  ind <- which(x > ur)  
  y <- 1 + xir*(x[ind] - ur)/sigmar
  f[ind] <- (1-phi[2])/sigmar*y^(-1/xir - 1)
  
  ## left tail
  ind <- which(x < ul)
  y <- (1 + xil*(-x[ind] + ul)/sigmal)
  f[ind] <- phi[1]/sigmal*y^(-1/xil - 1)
  
  return(f)
}


pmargin <- function(x, method, ul, sigmal, xil, ur, sigmar, xir, 
                    pfun = NULL, skt.par = NULL){
  if(method == "edfun") {
    Fx = pfun(x)
    phi = pfun(c(ul,ur))
  }
  if(method == "skt") {
    Fx = pstjf(x = x, param = skt.par)
    phi = pstjf(x = c(ul,ur), param = skt.par)
  }
  
  ind <- which(x > ur)
  y <- (1 + xir*(x[ind] - ur)/sigmar)
  y = pmax(y,0)
  Fx[ind] <- 1 - (1-phi[2])*y^(-1/xir)
  
  ind <- which(x < ul)
  y <- (1 + xil*(-x[ind] + ul)/sigmal)
  y = pmax(y,0)
  Fx[ind] <- phi[1]*y^(-1/xil)
  
  return(Fx)
}


qmargin <- function(u, method, ul, sigmal, xil, ur, sigmar, xir, 
                    qfun = NULL, pfun = NULL, skt.par = NULL){
  if(method == "edfun"){
    X <- qfun(u)
    phi = pfun(c(ul,ur))
  }
  if(method == "skt"){
    X = qstjf(p = u, param = skt.par)
    phi = pstjf(x = c(ul,ur), param = skt.par)
  }
  
  ind <- which(u > phi[2])  
  y <- (1-u[ind])/(1-phi[2])
  X[ind] <- (y^(-xir) - 1)/xir*sigmar + ur
  
  ## left tail
  ind <- which(u < phi[1])  
  y <- u[ind]/phi[1]
  X[ind] <- -((y^(-xil) - 1)/xil*sigmal - ul)
  
  return(X)
}