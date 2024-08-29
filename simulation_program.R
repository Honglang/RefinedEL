# Filename: simulation_program.R
# Authors: Xiang Wang, Honglang Wang
# Article: "Empirical Likelihood Based Inference for Functional Mean Models Accounting for Within-Subject Correlation"
# Sections in the article: main (Section 4. Simulation Studies) + supplementary material (Section 14. Simulation Studies Using Estimated Covariance)
# Required Packages: fdapace, PLRModels, locpol, emplik, doParallel, foreach


# Read functions ####

require(fdapace)
require(PLRModels)
require(locpol)
require(emplik)
require(doParallel)
require(foreach)

# load instrumental functions
source("../functions/Functions.R")


# Settings ####

detectCores()
cl <- makeCluster(20) # 20 is the number of cores (the number should be less than detectCores() for parallel computing)
registerDoParallel(cl)

args<-commandArgs(TRUE)
parameters = as.numeric(args) # e.g. parameters = c(100,5,201,1.5,0.75,3,2,1,1)
SampleSize = parameters[1]
mi = parameters[2]; Mi = rep(mi,SampleSize)
Ng = parameters[3]
Lambda = parameters[4]
VarianceCoef = parameters[5]
CorrelationCoef = parameters[6]
MeanFunction = parameters[7]
WorkingCov = parameters[8] # 0 is for true covariance, and 1 is for estimated covariance
RandomSeed = parameters[9]

set.seed(RandomSeed)
TimeRange = 1
TimePoints = seq(0, TimeRange, length.out = Ng)


# Generate data ####

print(Sys.time())

DataInfo = try(
  GenerateData(SampleSize, Mi, TimeRange, Lambda, VarianceCoef, CorrelationCoef, MeanFunction)
  , silent=FALSE)
# length(DataInfo)

while (length(DataInfo) == 1) {
  print('Generate Wrong Data')
  DataInfo = try(
    GenerateData(SampleSize, Mi, TimeRange, Lambda, VarianceCoef, CorrelationCoef, MeanFunction)
    , silent=FALSE)
}

D = DataInfo$D

UniqueObservedTime = sort(unique(D$Time))

if( WorkingCov == 1) {
  
  print(Sys.time())
  
  # FPCA
  fpca.res = FPCA(DataInfo$Y, DataInfo$Time, list(methodBwMu = 'GCV', methodBwCov = 'GCV', dataType='Sparse', error=TRUE, kernel='epan', verbose=TRUE, usergrid = TRUE))
  # CreateCovPlot(fpca.res, covPlotType = "Fitted", isInteractive = TRUE)
  
  # using estimated covariance as working covariance
  DataInfo$V_cov = lapply(DataInfo$Time, function(T_i){
    Cov_i = fpca.res$fittedCov[match(T_i, UniqueObservedTime),match(T_i, UniqueObservedTime)] + diag(fpca.res$sigma2, length(T_i))
    V_i = try(solve(Cov_i), silent = FALSE)
    if (length(V_i) == 1) {
      Var_i = diag(Cov_i)
      V_i = diag(x=1/Var_i, nrow=length(Var_i))
    }
    if( !isSymmetric.matrix(V_i) ) {
      V_i = (t(V_i) + V_i) / 2
    }
    return(V_i)
  })
  
  # worrking independence
  DataInfo$V_ind = lapply(DataInfo$Time, function(T_i){
    Cov_i = fpca.res$fittedCov[match(T_i, UniqueObservedTime),match(T_i, UniqueObservedTime)] + diag(fpca.res$sigma2, length(T_i))
    Var_i = diag(Cov_i)
    Var_i_inverse = diag(x=1/Var_i, nrow=length(Var_i))
    return(Var_i_inverse)
  })
  
  # worrking I.I.D.
  DataInfo$V_iid = lapply(DataInfo$Time, function(T_i){
    return( diag(1, nrow = length(T_i), ncol = length(T_i)) )
  })
  
  print(Sys.time())
}


# Bandwidth selection ####

print(Sys.time())

# criterions
criterions = c('ACV','ROT','LOSCV')

# Bandwidth selection for initial estimators

h_initial_range <- c(0.015, 1)

K = SampleSize # K = SampleSize for LOSCV, or K = 5 for 5-fold CV

CV_NW = function(h){
  ACVfunction(h, method = 'NW', D=D)
}
(h_NW = optimize(CV_NW, interval=h_initial_range, tol = 0.03)$minimum)

CV_ll = function(h){
  ACVfunction(h, method = 'll', D=D)
}
(h_ll = optimize(CV_ll, interval=h_initial_range, tol = 0.03)$minimum)

# Fan and Gijbels(1996)'s Rule of thumb is better for small bandwidth of local cubic estimator with small bias
(h_cubic = thumbBw(D$Time, D$Y, deg=3, kernel=EpaK,weig=rep(1,length(D$Y))))

print(paste('h_NW, h_ll, h_cubic = ',h_NW,h_ll,h_cubic))

print(Sys.time())

# Bandwidth selection for estimators and inference

print(Sys.time())

h_RAEL_range <- c(0.015, 0.25*SampleSize^(-1/5))

CV_RAEL = function(h){
  PE(h=h, D, DataInfo, K=K, method='RAEL', h_NW, h_ll, h_cubic)
}
(h_RAEL = optimize(CV_RAEL, interval=h_RAEL_range, tol = 0.03)$minimum)

print(paste('h_RAEL =',h_RAEL))

print(Sys.time())

h_RBEL_range <- c(0.015, 1)

CV_h_RBEL = function(h){
  PE(h=h, D, DataInfo, K=K, method='RBEL', h_NW, h_ll, h_cubic)
}
(h_RBEL = optimize(CV_h_RBEL, interval=h_RBEL_range, tol = 0.03)$minimum)

print(paste('h_RBEL =',h_RBEL))

print(Sys.time())

CV_h_rl_cov = function(h){
  PE(h=h, D, DataInfo, K=K, method='rl_cov', h_NW, h_ll, h_cubic)
}
(h_rl_cov = optimize(CV_h_rl_cov, interval=h_RBEL_range, tol = 0.03)$minimum)

print(Sys.time())

CV_h_gl_cov = function(h){
  PE(h=h, D, DataInfo, K=K, method='gl_cov', h_NW, h_ll, h_cubic)
}
(h_gl_cov = optimize(CV_h_gl_cov, interval=h_RBEL_range, tol = 0.03)$minimum)

print(paste('h_rl_cov, h_gl_cov =', h_rl_cov, h_gl_cov))

print(Sys.time())


# Estimators ####

print(Sys.time())

# NW estimator ####
timestart<-Sys.time()
mu_NW = apply(matrix(TimePoints, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_NW = data.frame(Method='NW',h_ini=h_NW,h_NW=h_NW,Bandwidth=h_NW,
                           SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                           runningtime=runningtime,criterion=criterions[1])

# Local linear estimator ####
timestart<-Sys.time()
mu_ll = apply(matrix(TimePoints, ncol = 1), 1, local_linear, D$Time, D$Y, h = h_ll)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_ll = data.frame(Method='ll',h_ini=h_ll,h_NW=h_NW,Bandwidth=h_ll,
                           SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                           runningtime=runningtime,criterion=criterions[1])

# Local cubic estimator ####
timestart<-Sys.time()
mu_cubic = locpol(Y~Time, data=D, weig=rep(1,nrow(D)), bw=h_cubic,kernel=EpaK,deg=3, xeval=TimePoints)$lpFit$Y
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_cubic = data.frame(Method='cubic',h_ini=h_cubic,h_NW=h_NW,Bandwidth=h_cubic,
                              SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                              runningtime=runningtime,criterion=criterions[3])

Estimators_LP = t(matrix(c(mu_NW,Parameters_NW,mu_ll,Parameters_ll,mu_cubic,Parameters_cubic), nrow = Ng+length(Parameters_NW)))
colnames(Estimators_LP) = c(1:Ng,names(Parameters_NW))

# Residual-Adjusted Empirical Likelihood (RAEL) estimator (Bias-corrected local constant estimator) ####
timestart<-Sys.time()
mu_initial_vec = apply(matrix(TimePoints, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
mu_initial_list = apply(matrix(D$Time, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
mu_RAEL = apply(matrix(TimePoints, ncol = 1), 1, BiasCorrected_local_cons, D$Time, D$Y, h = h_RAEL, mu_initial_vec, mu_initial_list, TimePoints)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_RAEL = data.frame(Method='mu_RAEL',h_ini=h_NW,h_NW=h_NW,Bandwidth=h_RAEL,
                             SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                             runningtime=runningtime,criterion=criterions[3])

# Refined Bias-corrected Empirical Likelihood (RBEL) estimator ####
timestart<-Sys.time()
mu_local_cubic_initial_list = lapply(DataInfo$Time, function(x) locpol(Y~Time, data=D, weig=rep(1,nrow(D)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
mu_initial_vec = apply(matrix(TimePoints, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW) )
mu_RBEL = apply(matrix(TimePoints, ncol = 1), 1, BiasCorrected_refined_cons_cov, DataInfo$D, DataInfo$Time, DataInfo$Y, DataInfo$V_cov, h = h_RBEL, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_RBEL = data.frame(Method='mu_RBEL',h_ini=h_cubic,h_NW=h_NW,Bandwidth=h_RBEL,
                             SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                             runningtime=runningtime,criterion=criterions[3])

# Refined linear Estimators ####
timestart<-Sys.time()
mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x, ncol = 1), 1, local_linear, D$Time, D$Y, h = h_ll) )
mu_rl_cov = apply(matrix(TimePoints, ncol = 1), 1, refined_linear_cov, DataInfo$D, DataInfo$Time, DataInfo$Y, DataInfo$V_cov, h = h_rl_cov, mu_initial_list)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_rl_cov = data.frame(Method='mu_rl_cov',h_ini=h_ll,h_NW=h_NW,Bandwidth=h_rl_cov,
                               SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                               runningtime=runningtime,criterion=criterions[3])

# Global linear Estimators ####
timestart<-Sys.time()
mu_gl_cov = global_linear_cov(TimePoints, DataInfo$D, DataInfo$Time, DataInfo$Y, DataInfo$V_cov, h_ll, h_gl_cov, Threshold = 0.001)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_gl_cov = data.frame(Method='mu_gl_cov',h_ini=h_ll,h_NW=h_NW,Bandwidth=h_gl_cov,
                               SampleSize=SampleSize,mi=mi,Ng=Ng,Lambda=Lambda,VarianceCoef=VarianceCoef,CorrelationCoef=CorrelationCoef,MeanFunction=MeanFunction,WorkingCov=WorkingCov,RandomSeed=RandomSeed,
                               runningtime=runningtime,criterion=criterions[3])

Estimators_RAEL_RBEL_Refined_Global = t(matrix(c(mu_RAEL, Parameters_RAEL, 
                                                 mu_RBEL,Parameters_RBEL,
                                                 mu_rl_cov,Parameters_rl_cov,
                                                 mu_gl_cov,Parameters_gl_cov),
                                               nrow = Ng+length(Parameters_RBEL)))
colnames(Estimators_RAEL_RBEL_Refined_Global) = c(1:Ng,names(Parameters_RBEL))

# save Estimators ### 

Estimators = rbind(Estimators_LP, Estimators_RAEL_RBEL_Refined_Global)

Estimators

saveRDS(Estimators, file=paste(getwd(),'/intermediate_results/Estimator_SampleSize',parameters[1],'mi',parameters[2],'Ng',parameters[3],'Lambda',parameters[4],'VarianceCoef',parameters[5],'CorrelationCoef',parameters[6],'MeanFunction',parameters[7],'WorkingCov',parameters[8],'RandomSeed',parameters[9],'.RDS',sep=''))


# RAEL CI ####

print(Sys.time())

(TimePoints.CI = seq(0.1,0.9,by=0.2))
(Ng = parameters[3] = length(TimePoints.CI))

RAEL_CI = foreach(t_scalar = TimePoints.CI, .combine = 'rbind', .packages = c('PLRModels','locpol','emplik')) %dopar% {
  
  mu_initial_vec = apply(matrix(t_scalar, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar)
  mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar) )
  
  (MELE_t = optimise(f = function(d) el.test( local_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, mu_initial_vec, mu_initial_list, h=h_RAEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR`, 
                     interval = range(D$Y), tol = 0.001)$minimum)
  
  (CILB_t = optimise(f = function(d) abs(el.test( local_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, mu_initial_vec, mu_initial_list, h=h_RAEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(range(D$Y)[1],MELE_t), tol = 0.001)$minimum)
  
  (CIUB_t = optimise(f = function(d) abs(el.test( local_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, mu_initial_vec, mu_initial_list, h=h_RAEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(MELE_t,range(D$Y)[2]), tol = 0.001)$minimum)
  
  is.rael.ci.cover = as.numeric( (CILB_t < mu(t_scalar, MeanFunction)) & (mu(t_scalar, MeanFunction) < CIUB_t) )
  
  rael.ci.length = CIUB_t - CILB_t
  
  CI_return_matrix = cbind(MELE_t,CILB_t,CIUB_t,h_NW,h_RAEL,is.rael.ci.cover, rael.ci.length)
  rownames(CI_return_matrix) = paste('t=',t_scalar,sep='')
  
  return(CI_return_matrix)
}

RAEL_CI

saveRDS(RAEL_CI, file=paste(getwd(),'/intermediate_results/RAEL_CI_SampleSize',parameters[1],'mi',parameters[2],'Ng',parameters[3],'Lambda',parameters[4],'VarianceCoef',parameters[5],'CorrelationCoef',parameters[6],'MeanFunction',parameters[7],'WorkingCov',parameters[8],'RandomSeed',parameters[9],'.RDS',sep=''))


# RBEL CI ####

print(Sys.time())

(TimePoints.CI = seq(0.1,0.9,by=0.2))
(Ng = parameters[3] = length(TimePoints.CI))

RBEL_CI = foreach(t_scalar = TimePoints.CI, .combine = 'rbind', .packages = c('PLRModels','locpol','emplik')) %dopar% {
  
  mu_initial_vec = apply(matrix(t_scalar, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar)
  mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar) )
  
  (MELE_t = optimise(f = function(d) el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR`, 
                     interval = range(D$Y), tol = 0.001)$minimum)
  
  (CILB_t = optimise(f = function(d) abs(el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(range(D$Y)[1],MELE_t), tol = 0.001)$minimum)
  
  (CIUB_t = optimise(f = function(d) abs(el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(MELE_t,range(D$Y)[2]), tol = 0.001)$minimum)
  
  is.rbel.ci.cover = as.numeric( (CILB_t < mu(t_scalar, MeanFunction)) & (mu(t_scalar, MeanFunction) < CIUB_t) )
  
  rbel.ci.length = CIUB_t - CILB_t
  
  CI_return_matrix = cbind(MELE_t,CILB_t,CIUB_t,h_NW,h_RBEL, is.rbel.ci.cover, rbel.ci.length)
  rownames(CI_return_matrix) = paste('t=',t_scalar,sep='')
  
  return(CI_return_matrix)
}

RBEL_CI

saveRDS(RBEL_CI, file=paste(getwd(),'/intermediate_results/RBEL_CI_SampleSize',parameters[1],'mi',parameters[2],'Ng',parameters[3],'Lambda',parameters[4],'VarianceCoef',parameters[5],'CorrelationCoef',parameters[6],'MeanFunction',parameters[7],'WorkingCov',parameters[8],'RandomSeed',parameters[9],'.RDS',sep=''))


# RAEL and RBEL Test ####

print(Sys.time())

(deviate.seq = seq(0,0.4,by=0.1))
(TimePoints.Test = seq(0.1, 0.9, by=0.2))
(Ng = parameters[3] = length(TimePoints.Test))

RAEL_RBEL_Test = foreach(d = deviate.seq, .packages = c('PLRModels','locpol','emplik','doParallel','foreach')) %:% 
  
  foreach(t_scalar = TimePoints.Test, .combine = 'rbind', .packages = c('PLRModels','locpol','emplik','doParallel','foreach')) %dopar% {
    
    for(i in 1:SampleSize){
      DataInfo$Y[[i]] = (DataInfo$Y)[[i]] + d*(1+(DataInfo$Time)[[i]]) # y+d*(1+t)
      
      DataInfo$GenerateTrueMean[[i]] = (DataInfo$GenerateTrueMean)[[i]] + d*(1+(DataInfo$Time)[[i]]) # 2*(t-0.5)^2+d*(1+t)
    }
    
    DataInfo$D$Y = (DataInfo$D)$Y + d*(1+(DataInfo$D)$Time) # y+d*(1+t)
    D = DataInfo$D
    
    # RAEL test ####
    mu_initial_vec = apply(matrix(t_scalar, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar)
    mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar) )
    
    ( el_RAEL = el.test( local_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, mu_initial_vec, mu_initial_list, h=h_RAEL, TimePoints=t_scalar, h_cubic, d = mu(t_scalar, MeanFunction)), c(0)) )
    
    (n2llr = el_RAEL$`-2LLR`)
    (pval = el_RAEL$Pval)
    (is.reject.qchisq = (n2llr >= qchisq(1-0.05,df=1)))
    (is.reject.pval = (pval <= 0.05))
    
    rael_test_return_matrix = cbind(n2llr,pval,h_NW,h_RAEL,is.reject.qchisq,is.reject.pval)
    rownames(rael_test_return_matrix) = paste('RAEL, t=',t_scalar,sep='')
    
    # RBEL test ####
    
    mu_initial_vec = apply(matrix(t_scalar, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar)
    mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar) )
    
    ( el_RBEL = el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d = mu(t_scalar, MeanFunction)), c(0)) )
    
    (n2llr = el_RBEL$`-2LLR`)
    (pval = el_RBEL$Pval)
    (is.reject.qchisq = (n2llr >= qchisq(1-0.05,df=1)))
    (is.reject.pval = (pval <= 0.05))
    
    rbel_test_return_matrix = cbind(n2llr,pval,h_cubic,h_RBEL,is.reject.qchisq,is.reject.pval)
    rownames(rbel_test_return_matrix) = paste('RBEL, t=',t_scalar,sep='')
    
    test_return_matrix = rbind(rael_test_return_matrix, rbel_test_return_matrix)
    
    return(test_return_matrix)
  }

(Ng = parameters[3] = length(TimePoints.Test))
names(RAEL_RBEL_Test) = paste('d=',deviate.seq,sep='')

RAEL_RBEL_Test

saveRDS(RAEL_RBEL_Test, file=paste(getwd(),'/intermediate_results/RAEL_RBEL_Test_SampleSize',parameters[1],'mi',parameters[2],'Ng',parameters[3],'Lambda',parameters[4],'VarianceCoef',parameters[5],'CorrelationCoef',parameters[6],'MeanFunction',parameters[7],'WorkingCov',parameters[8],'RandomSeed',parameters[9],'.RDS',sep=''))

print(Sys.time())