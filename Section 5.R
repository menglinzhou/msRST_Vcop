####################################################################
#### Section 5 #####################################################

path = ".../"
source(paste(path, "Functions.R", sep = ""))

#### Take Meta-vine distribution for example #######################

## import fitted models for Portfolio A
## including RVineMatrix() object of risk factor changes X
## skew-t marginals for risk factor changes and portfolio losses
load(paste(path, "model.RData")) 

## fitted coefficient from linear regression for Portfolio A
weight = c(-1.944, -0.018, 0.011, 0.011)
thred = 0.03 ## threshold at which to estimate the stress scenario
####################################################################


#### Data generation ###############################################

##############################################################
## This function generated random samples from meta-vine model
## Inputs: RVM: An RVineMatrix() object of X including the 
##              structure and the pair-copula families and parameters
##         skt_result: fitted parameters of skew-t marginals of X
##         w: weight vector to determine L
##         seed: seed of the random generation
################################################################
generate_data <- function(RVM, skt_result = NULL, w, seed = 2022){
  set.seed(seed)
  simudata = RVineSim(3000, RVM)
  colnames(simudata) = RVM$names
  simxdata = simudata
  for(i in 1:ncol(simudata)){
    name = colnames(simudata)[i]
    simxdata[,i] = qstjf(simudata[,i], param = skt_result[name,])
  }
  L = simxdata%*%matrix(w, ncol = 1)
  data = data.frame(L, cbind(simxdata))
  colnames(data) = c("L", colnames(simudata))
  return(data)
}


#### Take the first iteration as example #############
#### change seed to generate different samples########

sim = generate_data(RVM = RVM_A, skt_result = A_skt, 
                    w = weight, seed = 1)
####################################################################



#### Fit margins ###################################################
sim_skt = matrix(NA, nrow = ncol(sim), ncol = 4)
for(i in 1:ncol(sim)){
  set.seed(100)
  sim_skt[i,] = DEoptim(nllkstjf, lower = c(10*min(sim[,i]),0,0,0), 
                        upper = c(10*max(sim[,i]),10*sd(sim[,i]),10,10), 
                        ydat = sim[,i], 
                        control = DEoptim.control(trace = F))$optim$bestmem}
rownames(sim_skt) = colnames(sim)
####################################################################



#### CM1 estimator #################################################
udata = sim
for(i in 1:ncol(udata)){
  udata[,i] = (rank(sim[,i]) - 0.5)/length(sim[,i])
} ## transform data to uscale

RVM_fit = RVineStructureSelect(data = udata, cores = 10,
                               familyset = c(0,1,2,3,4,5,6,7,8,10))
## fit RVM to udata

optim_fun = function(x){
  den = con_den(x, RVM = RVM_fit[[i]], Up = thred, 
                skt_result = sim_skt[-1,], set = c(), 
                set_margin = NULL, PnL_margin = sim_skt[1,])
  return(-den)}

CM1 = DEoptim(optim_fun, lower = 10*apply(sim[,-1], 2, min),
              upper = 10*apply(sim[,-1], 2, max),
              control = DEoptim.control(trace = F, 
                                        itermax = 400))$optim$bestmem
## optimize conditional density function with DEoptim
####################################################################




#### CM2 estimator #################################################
udata = sim
for(i in 1:ncol(udata)){
  udata[,i] = (rank(sim[,i]) - 0.5)/length(sim[,i])
} ## transform data to uscale

udata = data.frame(udata[,-1], L = udata$L)
## put L as the last column

VA = getVinearray(vWcor(as.matrix(udata)), n = nrow(udata), 
                  method = "MSTleaf", cfi_bd = 1.01)
RVM_reg = CopReg(data = udata, vine_array = VA$vine_array, 
                 family_set = c(1,2,3), cores = core)
RVM_reg$family = RVM_fit$family[d:1,d:1]
RVM_reg$par = RVM_fit$par[d:1,d:1]
RVM_reg$par2 = RVM_fit$par2[d:1,d:1]
RVM_reg$logLik = RVM_fit$logLik
RVM_reg$pair.logLik = RVM_fit$pair.logLik[d:1,d:1]
## fit vine regression to udata

optim_fun = function(x){
  den = con_den_g_est(x, RVM = RVM_reg, RVM_X = RVM_X,
                      thred=thred, skt.par = sim_skt[-1,], 
                      set = c(), set_margin = NULL, 
                      PnL_margin = sim_skt[1,])
  return(-log(den))}

CM2 = DEoptim(optim_fun, lower = 10*apply(sim[,-1], 2, min),
              upper = 10*apply(sim[,-1], 2, max),
              control = DEoptim.control(trace = F, 
                                        itermax = 400))$optim$bestmem
## optimize conditional density function with DEoptim
####################################################################



#### CM3 estimator #################################################
udata = sim
for(i in 1:ncol(udata)){
  udata[,i] = (rank(sim[,i]) - 0.5)/length(sim[,i])
} ## transform data to uscale

udata = data.frame(udata[,-1], L = udata$L)
## put L as the last column

VA = getVinearray(vWcor(as.matrix(udata)), n = nrow(udata), 
                  method = "MSTleaf", cfi_bd = 1.01)
RVM_fit = RVineCopSelect(udata, familyset = c(0,1,2,3,4,5,6,7,8,10), 
                         Matrix = VA$vine_array, cores = cores)
RVM_X = RVineMatrix(Matrix = RVM_fit$Matrix[-1,-1], 
                    family = RVM_fit$family[-1,-1],
                    par = RVM_fit$par[-1,-1], par2 = RVM_fit$par2[-1,-1],
                    names = RVM_fit$names[-length(RVM_fit$names)])
## fit RVM to udata

optim_fun = function(x){
  den = con_den_gt(x, RVM_X = RVM_X,
                   RVM_all = RVM_fit, Up = thred, 
                   skt_result = sim_skt[-1,], 
                   set = c(), set_margin = NULL,
                   PnL_margin = sim_skt[1,])
  return(-den)}

CM3 = DEoptim(optim_fun, lower = 10*apply(sim[,-1], 2, min),
              upper = 10*apply(sim[,-1], 2, max),
              control = DEoptim.control(trace = F, 
                                        itermax = 400))$optim$bestmem
## optimize conditional density function with DEoptim
####################################################################



#### CM_star estimator #############################################
optim_fun = function(x){
  den = con_den_known(x, RVM = RVM_X, thred =thred,
                      skt_result = sim_skt[-1,], 
                      set = c(), set_margin = NULL,
                      PnL_margin = sim_skt[1,])
  return(-log(den))}

CM_star = DEoptim(optim_fun, lower = 10*apply(sim[,-1], 2, min),
              upper = 10*apply(sim[,-1], 2, max),
              control = DEoptim.control(trace = F, 
                                        itermax = 400))$optim$bestmem
## optimize conditional density function with DEoptim
####################################################################