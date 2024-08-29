# Filename: ADNI_AD_data_application.R
# Authors: Xiang Wang, Honglang Wang
# Article: "Empirical Likelihood Based Inference for Functional Mean Models Accounting for Within-Subject Correlation"
# Sections in the article: main article (Section 5. Real Data Applications)
# Required Packages: face, PLRModels, locpol, emplik, doParallel, foreach


# Read functions ####

require(face)
require(PLRModels)
require(locpol)
require(emplik)
require(doParallel)
require(foreach)

# load instrumental functions
source("../functions/Functions.R")


# Settings ####

set.seed(1)

detectCores()
cl <- makeCluster(20) # 20 is the number of cores (the number should be less than detectCores() for parallel computing)
registerDoParallel(cl)


# ADNI (AD) data ####

print(Sys.time())

ADNIMERGE <- read.csv("ADNIMERGE.csv")

summary(ADNIMERGE[ADNIMERGE$AGE > 0,c('RID','DX_bl','DX','AGE','Years_bl','MMSE','MMSE_bl')])

table(ADNIMERGE$DX)
table(ADNIMERGE$DX_bl)

DF = na.omit(ADNIMERGE[ADNIMERGE$AGE >= 55,c('RID','DX_bl','AGE','Years_bl','MMSE')])

table(DF$DX_bl)

summary(DF)

SingleObsIndex = names(table(DF$RID))[table(DF$RID)==1]

DF_AD = DF[(DF$DX_bl=='AD') & (!(DF$RID %in% SingleObsIndex)),]

table(DF_AD$DX_bl)

I_vec = DF_AD$RID

Age_bl = DF_AD$AGE + DF_AD$Years_bl
# boxplot(Age_bl)
summary(Age_bl)

(TimeRange = range(Age_bl))

t_vec = (Age_bl - TimeRange[1]) / (TimeRange[2]-TimeRange[1])

Y_vec = DF_AD$MMSE


D = data.frame(Time=as.numeric( t_vec ), Y= Y_vec, Index=as.integer( I_vec ))
summary(D);dim(D)

DataInfo = list()

DataInfo$Y = lapply(unique(D$Index), function(i){
  return(D$Y[D$Index==i])
})
DataInfo$Time = lapply(unique(D$Index), function(i){
  return(D$Time[D$Index==i])
})

# Remove subjects that have repeated measurements at same tij. Or we can take average of the repeated measurements.
D = D[!(D$Index %in% which( unlist( lapply(DataInfo$Time, function(x) length(x) != length(unique(x)))))),]

D$Index = as.integer(factor(D$Index))

DataInfo = list()

DataInfo$D = D

DataInfo$Y = lapply(unique(D$Index), function(i){
  return(D$Y[D$Index==i])
})
DataInfo$Time = lapply(unique(D$Index), function(i){
  return(D$Time[D$Index==i])
})


(SampleSize = length(DataInfo$Time))
Mi = sapply(DataInfo$Time, length)
table(Mi)
summary(Mi)
(mi = round(mean(Mi), digits = 2))


# FACE ####

# require(face)

face.res = face.sparse(data = data.frame(argvals=D[,c(1)], subj=D[,c(3)], y=D[,c(2)]), newdata = data.frame(argvals=D[,c(1)], subj=D[,c(3)], y=D[,c(2)]),
                       center = TRUE, argvals.new = D[,1],
                       knots = 7,
                       p = 3, m = 2, lambda = NULL, lambda_mean = NULL,
                       search.length = 14,
                       lower = -3, upper = 10, lower2 = -3, upper2 = 5,
                       calculate.scores = FALSE,pve=0.99,two_step=FALSE)

UniqueGrids = sort(unique(D$Time))

DataInfo$V_cov = lapply(DataInfo$Time, function(T_i){
  COV_i = face.res$Chat.new[match(T_i, UniqueGrids),match(T_i, UniqueGrids)] + diag(face.res$sigma2, nrow = length(T_i))
  V_i = solve(COV_i)
  V_i = (t(V_i) + V_i) / 2
  return(V_i)
})


# Bandwidth selection ####

print(Sys.time())

# criterions
criterions = c('ACV','ROT','LOSCV')

# Bandwidth selection for initial estimators

h_initial_range <- c(0.1, 1)

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

(h.seq = seq(0.1,1,by=0.05))

Sys.time()

pe.seq = foreach(h = h.seq, .combine = 'rbind', .packages = c('PLRModels','locpol','doParallel','foreach')) %dopar% {
  pe.rc.bc.cov = PE(h=h, D, DataInfo, K=SampleSize, method='RBEL', h_NW, h_ll, h_cubic)
  return( c( pe.rc.bc.cov = pe.rc.bc.cov ) )
}
(h_RBEL = h.seq[which.min(pe.seq)])

Sys.time()


# Estimators ####

print(Sys.time())

(TimePoints <- (60:90- TimeRange[1]) / (TimeRange[2]-TimeRange[1]) )
(Ng = length(TimePoints))

# NW estimator ####
timestart<-Sys.time()
mu_NW = apply(matrix(TimePoints, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_NW = data.frame(Method='NW',h_ini=h_NW,h_NW=h_NW,Bandwidth=h_NW,
                           runningtime=runningtime,criterion=criterions[1])

# Local linear estimator ####
timestart<-Sys.time()
mu_ll = apply(matrix(TimePoints, ncol = 1), 1, local_linear, D$Time, D$Y, h = h_ll)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_ll = data.frame(Method='ll',h_ini=h_ll,h_NW=h_NW,Bandwidth=h_ll,
                           runningtime=runningtime,criterion=criterions[1])

# Local cubic estimator ####
timestart<-Sys.time()
mu_cubic = locpol(Y~Time, data=D, weig=rep(1,nrow(D)), bw=h_cubic,kernel=EpaK,deg=3, xeval=TimePoints)$lpFit$Y
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_cubic = data.frame(Method='cubic',h_ini=h_cubic,h_NW=h_NW,Bandwidth=h_cubic,
                              runningtime=runningtime,criterion=criterions[2])

Estimators_LP = t(matrix(c(mu_NW,Parameters_NW,mu_ll,Parameters_ll,mu_cubic,Parameters_cubic), nrow = Ng+length(Parameters_NW)))
colnames(Estimators_LP) = c(1:Ng,names(Parameters_NW))

# Refined Bias-corrected Empirical Likelihood (RBEL) estimator ####
timestart<-Sys.time()
mu_local_cubic_initial_list = lapply(DataInfo$Time, function(x) locpol(Y~Time, data=D, weig=rep(1,nrow(D)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
mu_initial_vec = apply(matrix(TimePoints, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW)
mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW) )
mu_RBEL = apply(matrix(TimePoints, ncol = 1), 1, BiasCorrected_refined_cons_cov, DataInfo$D, DataInfo$Time, DataInfo$Y, DataInfo$V_cov, h = h_RBEL, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list)
timeend<-Sys.time(); runningtime<-timeend-timestart
Parameters_RBEL = data.frame(Method='mu_RBEL',h_ini=h_cubic,h_NW=h_NW,Bandwidth=h_RBEL,
                             runningtime=runningtime,criterion=criterions[3])

Estimators_RBEL = t(matrix(c(mu_RBEL,Parameters_RBEL),
                           nrow = Ng+length(Parameters_RBEL)))
colnames(Estimators_RBEL) = c(1:Ng,names(Parameters_RBEL))

# Save estimators ### 

Estimators = rbind(Estimators_LP, Estimators_RBEL)

Estimators

saveRDS(Estimators, file=paste(getwd(),'/intermediate_results/ADNI_AD_Estimator.RDS',sep=''))


# RBEL working unstructured covariance CI ####

print(Sys.time())

(TimePoints.CI <- (60:90- TimeRange[1]) / (TimeRange[2]-TimeRange[1]) )
(Ng = length(TimePoints.CI))

RBEL_CI = foreach(t_scalar = TimePoints.CI, .combine = 'rbind', .packages = c('PLRModels','locpol','emplik')) %dopar% {
  
  mu_initial_vec = apply(matrix(t_scalar, ncol = 1), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar)
  mu_initial_list = lapply(DataInfo$Time, function(x) apply(matrix(x), 1, local_cons, D$Time, D$Y, h = h_NW, TimePoints=t_scalar) )
  
  (MELE_t = optimise(f = function(d) el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR`, 
                     interval = range(D$Y), tol = 0.01)$minimum)
  
  (CILB_t = optimise(f = function(d) abs(el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(range(D$Y)[1],MELE_t), tol = 0.01)$minimum)
  
  (CIUB_t = optimise(f = function(d) abs(el.test( refined_cons_bias_corrected_estimating_equation(t_scalar, Time=DataInfo$Time, Y=DataInfo$Y, V_cov=DataInfo$V_cov, mu_initial_vec, mu_initial_list, h=h_RBEL, TimePoints=t_scalar, h_cubic, d), c(0))$`-2LLR` - qchisq(1-0.05,df=1)),
                     interval = c(MELE_t,range(D$Y)[2]), tol = 0.01)$minimum)
  
  rbel.ci.length = CIUB_t - CILB_t
  
  CI_return_matrix = cbind(MELE_t, CILB_t, CIUB_t, h_NW, h_cubic, h_RBEL, rbel.ci.length)
  rownames(CI_return_matrix) = paste('t=',t_scalar,sep='')
  
  return(CI_return_matrix)
}

print(Sys.time())

RBEL_CI

saveRDS(RBEL_CI, file=paste(getwd(),'/intermediate_results/ADNI_AD_RBEL_CI.RDS',sep=''))


# Draw figures ####

# getwd()
# RBEL_CI <- readRDS(paste(getwd(),'/intermediate_results/ADNI_AD_RBEL_CI.RDS',sep=''))

ADNI_AD_grids = 60:90

colnames(RBEL_CI)

pdf("ADNI_AD_RBEL.pdf",width=5, height=5)

par(mfrow=c(1,1), mar = c(5, 5, 5/5, 1) * 0.5, mgp = c(1.5, 0.5, 0), ps = 12, cex = 0.9)

plot(RBEL_CI[,'MELE_t'] ~ ADNI_AD_grids, type='l', xlab='Age', ylab='MMSE', ylim=c(20,30))
lines(RBEL_CI[,'CILB_t'] ~ ADNI_AD_grids, lty=2)
lines(RBEL_CI[,'CIUB_t'] ~ ADNI_AD_grids, lty=2)
legend('top', c(paste('AD (h = ', round(h_RBEL, digits = 1), ')',sep=''),'MELE','95% CI'), lty=c(NA,1:2))

dev.off()


# The selected bandwidth

print('The selected ADNI_AD bandwidth h_RBEL')
h_RBEL
h_RBEL*(TimeRange[2]-TimeRange[1])
TimeRange