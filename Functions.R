# Data setting ####

# true f (density) function ####
f_T = dunif
rf_T = runif

# true mean function ####
mu = function(t, MeanFunction){
  if(MeanFunction==1){
    return(20*(sin(t*2*pi*3)))
  }
  if(MeanFunction==2){
    return(2*(t-0.5)^2)
  }
}

# H0: mean function mu0(t) ####
mu0 = function(d){
  return(d)
}

# H1: mean function mu1(t) ####
mu1 = function(t, MeanFunction, d=0){
  if(MeanFunction==1){
    return(20*(sin(t*2*pi*3))+d*20*(t-0.5))
  }
  if(MeanFunction==2){
    return(2*(t-0.5)^2+d*(1+t))
  }
}

# true sigma function ####
sigma = function(t, VarianceCoef){
  return(VarianceCoef * log(5+t))
}

# true rho function ####
rho = function(tij, tik, TimeRange, CorrelationCoef){
  return( 1 - (abs(tij - tik)/TimeRange)^(1/CorrelationCoef) )
}

# function to compute the inverse square root of a matrix
fnMatSqrt = function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}

# generate data function ####
GenerateData = function(SampleSize, Mi, TimeRange, Lambda, VarianceCoef, CorrelationCoef, MeanFunction){
  
  Time = lapply(Mi, function(m){
    sort(rf_T(m))
  })
  GenerateTrueMean = lapply(Time, mu, MeanFunction=MeanFunction)
  
  GenerateTrueSigma = lapply(Time, sigma, VarianceCoef)
  GenerateTrueRho = lapply(Time, function(x) {
    rhomatrix = matrix(nrow = length(x), ncol = length(x))
    for (j in 1:length(x)) {
      for (k in 1:length(x)) {
        rhomatrix[j, k] = rho(x[j], x[k], TimeRange, CorrelationCoef)
      }
    }
    return(rhomatrix)
  } )
  GenerateTrueCov = lapply(c(1:SampleSize), function(x) diag(GenerateTrueSigma[[x]], nrow = Mi[x]) %*% GenerateTrueRho[[x]] %*% diag(GenerateTrueSigma[[x]], nrow = Mi[x]) + diag(Lambda^2, nrow = Mi[x]))
  
  V_cov = lapply(GenerateTrueCov, solve)
  
  GenerateTrueVar = lapply(GenerateTrueCov, diag)
  V_ind = lapply(GenerateTrueVar, function(x) diag(1/x))
  
  V_iid = lapply(Mi, function(x) diag(1, nrow=x, ncol=x))
  
  Y = list()
  for (i in 1:SampleSize) {
    N_0_I = rnorm(Mi[i], mean = 0, sd = 1)
    A = fnMatSqrt(GenerateTrueCov[[i]])# t(A) %*% N_0_I is subject-specific random trajectory
    Y[[i]] = GenerateTrueMean[[i]] + t(A) %*% N_0_I + rnorm(n=Mi[i], mean=0, sd=Lambda) # mvrnorm(1, GenerateTrueMean[[i]], GenerateTrueCov[[i]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  }
  
  D = cbind(unlist(Time), unlist(Y), rep(1:SampleSize, Mi))
  
  colnames(D) = c("Time", "Y", "Index")
  (D = as.data.frame(D))
  
  return(list(D = D, Time = Time, Y = Y, GenerateTrueMean = GenerateTrueMean, V_cov = V_cov, V_ind = V_ind, V_iid = V_iid))
}


# Kernal based estimators ####

# kernal function ####
kerf = function(x){
  kij = .75*(1-x^2)*(abs(x) <=  1) + 0*(abs(x)>1)
  return(kij)
}

# kernal function ( . / h) /h ####
Kerf = function(x, h){
  Kij = kerf(x/h)/h
  return(Kij)
}

# kernal matrix function ( . / h) /h ####
Ki = function(x, h){
  Ki = diag(kerf(x/h)/h, nrow = length(x)) # better to claim the nrow / ncol of diag
  return(Ki)
}

# ones matrix function ####
one = function(Mi){
  return(as.matrix(rep(1, Mi)))
}

# w_ij vector function ####
wij = function(t0, tij){
  return(matrix(c(1, tij-t0), nrow = 2))
}

# W_i matrix function ####
Wi = function(t0, tvec){
  return(cbind(rep(1, length(tvec)), tvec-t0))
}

# Local constant estimator ####
# t0: grid time point t to estimate
# tvec: all T_ij as a vector
# Yvec: all Y_ij as a vector
# h: bandwidth
local_cons = function(t0, tvec, Yvec, h,TimePoints){
  Kvec = Kerf(tvec-t0, h)
  return(sum(Yvec*Kvec)/sum(Kvec))
}

# Bias-corrected local constant estimator ####
# t0: grid time point t to estimate
# tvec: all T_ij as a vector
# Yvec: all Y_ij as a vector
# h: bandwidth
# mu_initial_vec: initial estimator as a vector for bias correction
# mu_initial_list: initial estimator as a list for bias correction
BiasCorrected_local_cons = function(t0, tvec, Yvec, h, mu_initial_vec, mu_initial_list,TimePoints){
  Kvec = Kerf(tvec-t0, h)
  est = sum((Yvec - (mu_initial_list - mu_initial_vec[which.min(abs(TimePoints-t0))]) )*Kvec)/sum(Kvec)
  return(est)
}

# Local linear estimator ####
# t0: grid time point t to estimate
# tvec: all T_ij as a vector
# Yvec: all Y_ij as a vector
# h: bandwidth
S = function(t0, tvec, h, r){
  Kvec = Kerf(tvec-t0, h)
  Tvec = ((tvec-t0)/h)^r
  return(mean(Kvec*Tvec))
}
R = function(t0, tvec, Yvec, h, r){
  Kvec = Kerf(tvec-t0, h)
  Tvec = ((tvec-t0)/h)^r
  return(mean(Kvec*Tvec*Yvec))
}

local_linear = function(t0, tvec, Yvec, h){
  return( (R(t0, tvec, Yvec, h, 0)*S(t0, tvec, h, 2) - R(t0, tvec, Yvec, h, 1)*S(t0, tvec, h, 1)) / (S(t0, tvec, h, 0)*S(t0, tvec, h, 2)-(S(t0, tvec, h, 1))^2) )
}


# ACV bandwidth selection criterion for initial local constant / linear estimator ####

ACVfunction = function(h, method, D){
  
  if(method=='NW'){
    H = t(sapply(D$Time, function(t){
      K = Kerf(D$Time - t, h)
      Hij = K / sum(K)
      return(Hij)
    }))
  }
  
  if(method=='ll'){
    H = t(sapply(D$Time, function(t){
      One = rep(1,length(D$Time))
      W = cbind(One, D$Time - t)
      K = Ki(W[,2], h)
      Hij = (solve(t(W) %*% K %*% W) %*% t(W) %*% K)[1,]
      return(Hij)
    }))
  }
  
  if(method=='cubic'){
    H = t(sapply(D$Time, function(t){
      One = rep(1,length(D$Time))
      W = cbind(One, D$Time - t, (D$Time - t)^2/2, (D$Time - t)^3/6)
      K = Ki(W[,2], h)
      Hij = (solve(t(W) %*% K %*% W) %*% t(W) %*% K)[1,]
      return(Hij)
    }))
  }
  
  ACV_L2norm = sapply(unique(D$Index), function(i){
    Ii = diag(1, nrow = sum(D$Index==i))
    Hii = H[D$Index==i, D$Index==i]
    ACV_i = solve(Ii-Hii) %*% (D$Y[D$Index==i] - (H %*% D$Y)[D$Index==i])
    return(sqrt(sum(ACV_i^2)))
  })
  
  ACV = mean( ACV_L2norm )
  return(ACV)
}


# Refined constant estimator (using working unstructured covariance) ####
# t0: grid time point t to estimate
# D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
# Time list (each element is a m_i-dim vector for time record)
# Y list (each element is a m_i-dim vector for response record)
# V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
# h: bandwidth
# mu_initial_list: initial estimator as a list
refined_cons_cov = function(t0, D, Time, Y, V_cov, h, mu_initial_list){
  
  SampleSize = length(Time)
  Mi = sapply(Time, length)
  
  Kmat =  lapply(Time, function(x) Ki(x-t0, h))
  
  Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_initial_list[[x]]) ) )
  Numerator = unlist(Numerator)
  
  Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
  Denominator = unlist(Denominator)
  
  return(sum(Numerator)/sum(Denominator))
}

# Bias-corrected Refined constant estimator (using working unstructured covariance) ####
# t0: grid time point t to estimate
# D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
# Time list (each element is a m_i-dim vector for time record)
# Y list (each element is a m_i-dim vector for response record)
# V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
# h: bandwidth
# mu_initial_vec: initial estimator as a vector for bias correction
# mu_initial_list: initial estimator as a list for bias correction
# mu_local_cubic_initial_list: initial local cubic estimator as a list
BiasCorrected_refined_cons_cov = function(t0, D, Time, Y, V_cov, h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list){
  
  SampleSize = length(Time)
  Mi = sapply(Time, length)
  
  Kmat =  lapply(Time, function(x) Ki(x-t0, h))
  
  Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + (h * Kmat[[x]]) %*% matrix(mu_initial_vec[which.min(abs(TimePoints-t0))] - mu_initial_list[[x]]) ) )
  Numerator = unlist(Numerator)
  
  Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
  Denominator = unlist(Denominator)
  
  return(sum(Numerator)/sum(Denominator))
}


# Refined linear estimator (using working unstructured covariance) ####
# t0: grid time point t to estimate
# D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
# Time list (each element is a m_i-dim vector for time record)
# Y list (each element is a m_i-dim vector for response record)
# V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
# h: bandwidth
# mu_initial_list: initial estimator as a list
refined_linear_cov = function(t0, D, Time, Y, V_cov, h, mu_initial_list){
  
  SampleSize = length(Time)
  Mi = sapply(Time, length)
  
  Kmat =  lapply(Time, function(x) Ki(x-t0, h))
  
  Numerator = lapply(c(1:SampleSize), function(x) t(Wi(t0, Time[[x]])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_initial_list[[x]]) ) )
  Numerator <- Reduce("+", Numerator)
  
  Denominator = lapply(c(1:SampleSize), function(x) h * t(Wi(t0, Time[[x]])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% Wi(t0, Time[[x]]) )
  Denominator <- Reduce("+", Denominator)
  
  return((solve(Denominator) %*% Numerator)[1, ])
}


# Global linear estimator (using working unstructured covariance) with iteration ####
# TimePoints: grid time points t's to estimate
# D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
# Time list (each element is a m_i-dim vector for time record)
# Y list (each element is a m_i-dim vector for response record)
# V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
# OptimalBandwidth_local_cons: the bandwidth used for initial local_cons estimator
# Bandwidth: the bandwidth used for global linear estimator
# Threshold to break iteration
global_linear_cov = function(TimePoints, D, Time, Y, V_cov, OptimalBandwidth_local_linear, Bandwidth, Threshold){
  
  ut = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, local_linear, D$Time, D$Y, OptimalBandwidth_local_linear) )
  UT = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = ut) )
  
  max_diff = max(abs(unlist(ut)-unlist(UT)))
  
  # iteration = 0
  while (max_diff>Threshold) {
    ut = UT
    UT = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = ut) )
    max_diff = max(abs(unlist(ut)-unlist(UT)))
    # iteration = iteration+1
    # print(paste("iteration & max_diff: ", c(iteration, max_diff)))
  }
  
  UT = apply(matrix(TimePoints, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = UT)
  return(UT)
}


# The optimal bandwidth h that minimizes the prediction error PE(h) ####

PE = function(h, D, DataInfo, K, method, h_NW, h_ll, h_cubic){
  
  print(paste('h =',h))
  
  n_full = length(unique(D$Index))
  folds_i <- sample(rep(1:K, length.out = n_full))
  
  SSE = foreach(k = 1:K, .combine = '+', .packages = c('PLRModels','locpol')) %dopar% {
    
    # kernal function ####
    kerf = function(x){
      kij = .75*(1-x^2)*(abs(x) <=  1) + 0*(abs(x)>1)
      return(kij)
    }
    
    # kernal function ( . / h) /h ####
    Kerf = function(x, h){
      Kij = kerf(x/h)/h
      return(Kij)
    }
    
    # kernal matrix function ( . / h) /h ####
    Ki = function(x, h){
      Ki = diag(kerf(x/h)/h, nrow = length(x)) # better to claim the nrow / ncol of diag
      return(Ki)
    }
    
    # ones matrix function ####
    one = function(Mi){
      return(as.matrix(rep(1, Mi)))
    }
    
    # w_ij vector function ####
    wij = function(t0, tij){
      return(matrix(c(1, tij-t0), nrow = 2))
    }
    
    # W_i matrix function ####
    Wi = function(t0, tvec){
      return(cbind(rep(1, length(tvec)), tvec-t0))
    }
    
    # Local constant estimator ####
    # t0: grid time point t to estimate
    # tvec: all T_ij as a vector
    # Yvec: all Y_ij as a vector
    # h: bandwidth
    local_cons = function(t0, tvec, Yvec, h,TimePoints){
      Kvec = Kerf(tvec-t0, h)
      return(sum(Yvec*Kvec)/sum(Kvec))
    }
    
    # Bias-corrected local constant estimator ####
    # t0: grid time point t to estimate
    # tvec: all T_ij as a vector
    # Yvec: all Y_ij as a vector
    # h: bandwidth
    # mu_initial_vec: initial estimator as a vector for bias correction
    # mu_initial_list: initial estimator as a list for bias correction
    BiasCorrected_local_cons = function(t0, tvec, Yvec, h, mu_initial_vec, mu_initial_list,TimePoints){
      Kvec = Kerf(tvec-t0, h)
      est = sum((Yvec - (mu_initial_list - mu_initial_vec[which.min(abs(TimePoints-t0))]) )*Kvec)/sum(Kvec)
      return(est)
    }
    
    # Local linear estimator ####
    # t0: grid time point t to estimate
    # tvec: all T_ij as a vector
    # Yvec: all Y_ij as a vector
    # h: bandwidth
    S = function(t0, tvec, h, r){
      Kvec = Kerf(tvec-t0, h)
      Tvec = ((tvec-t0)/h)^r
      return(mean(Kvec*Tvec))
    }
    R = function(t0, tvec, Yvec, h, r){
      Kvec = Kerf(tvec-t0, h)
      Tvec = ((tvec-t0)/h)^r
      return(mean(Kvec*Tvec*Yvec))
    }
    
    local_linear = function(t0, tvec, Yvec, h){
      return( (R(t0, tvec, Yvec, h, 0)*S(t0, tvec, h, 2) - R(t0, tvec, Yvec, h, 1)*S(t0, tvec, h, 1)) / (S(t0, tvec, h, 0)*S(t0, tvec, h, 2)-(S(t0, tvec, h, 1))^2) )
    }
    
    
    # Refined constant estimator (using working unstructured covariance) ####
    # t0: grid time point t to estimate
    # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
    # Time list (each element is a m_i-dim vector for time record)
    # Y list (each element is a m_i-dim vector for response record)
    # V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
    # h: bandwidth
    # mu_initial_list: initial estimator as a list
    refined_cons_cov = function(t0, D, Time, Y, V_cov, h, mu_initial_list){
      
      SampleSize = length(Time)
      Mi = sapply(Time, length)
      
      Kmat =  lapply(Time, function(x) Ki(x-t0, h))
      
      Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_initial_list[[x]]) ) )
      Numerator = unlist(Numerator)
      
      Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
      Denominator = unlist(Denominator)
      
      return(sum(Numerator)/sum(Denominator))
    }
    
    # Bias-corrected Refined constant estimator (using working unstructured covariance) ####
    # t0: grid time point t to estimate
    # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
    # Time list (each element is a m_i-dim vector for time record)
    # Y list (each element is a m_i-dim vector for response record)
    # V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
    # h: bandwidth
    # mu_initial_vec: initial estimator as a vector for bias correction
    # mu_initial_list: initial estimator as a list for bias correction
    # mu_local_cubic_initial_list: initial local cubic estimator as a list
    BiasCorrected_refined_cons_cov = function(t0, D, Time, Y, V_cov, h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list){
      
      SampleSize = length(Time)
      Mi = sapply(Time, length)
      
      Kmat =  lapply(Time, function(x) Ki(x-t0, h))
      
      Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + (h * Kmat[[x]]) %*% matrix(mu_initial_vec[which.min(abs(TimePoints-t0))] - mu_initial_list[[x]]) ) )
      Numerator = unlist(Numerator)
      
      Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
      Denominator = unlist(Denominator)
      
      return(sum(Numerator)/sum(Denominator))
    }
    
    # Refined linear estimator (using working unstructured covariance) ####
    # t0: grid time point t to estimate
    # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
    # Time list (each element is a m_i-dim vector for time record)
    # Y list (each element is a m_i-dim vector for response record)
    # V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
    # h: bandwidth
    # mu_initial_list: initial estimator as a list
    refined_linear_cov = function(t0, D, Time, Y, V_cov, h, mu_initial_list){
      
      SampleSize = length(Time)
      Mi = sapply(Time, length)
      
      Kmat =  lapply(Time, function(x) Ki(x-t0, h))
      
      Numerator = lapply(c(1:SampleSize), function(x) t(Wi(t0, Time[[x]])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_initial_list[[x]]) ) )
      Numerator <- Reduce("+", Numerator)
      
      Denominator = lapply(c(1:SampleSize), function(x) h * t(Wi(t0, Time[[x]])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% Wi(t0, Time[[x]]) )
      Denominator <- Reduce("+", Denominator)
      
      return((solve(Denominator) %*% Numerator)[1, ])
    }
    
    
    # Global linear estimator (using working unstructured covariance) with iteration ####
    # TimePoints: grid time points t's to estimate
    # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
    # Time list (each element is a m_i-dim vector for time record)
    # Y list (each element is a m_i-dim vector for response record)
    # V_cov list (each element is a m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
    # OptimalBandwidth_local_cons: the bandwidth used for initial local_cons estimator
    # Bandwidth: the bandwidth used for global linear estimator
    # mu_initial_list: initial estimator as a list
    # Threshold to break iteration
    global_linear_cov = function(TimePoints, D, Time, Y, V_cov, OptimalBandwidth_local_linear, Bandwidth, Threshold){
      
      ut = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, local_linear, D$Time, D$Y, OptimalBandwidth_local_linear) )
      UT = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = ut) )
      
      max_diff = max(abs(unlist(ut)-unlist(UT)))
      
      # iteration = 0
      while (max_diff>Threshold) {
        ut = UT
        UT = lapply(Time, function(x) apply(matrix(x, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = ut) )
        max_diff = max(abs(unlist(ut)-unlist(UT)))
        # iteration = iteration+1
        # print(paste("iteration & max_diff: ", c(iteration, max_diff)))
      }
      
      UT = apply(matrix(TimePoints, ncol = 1), 1, refined_linear_cov, D, Time, Y, V_cov, Bandwidth, mu_initial_list = UT)
      return(UT)
    }
    
    SE = 0
    
    # print(paste('h =',h,', k =',k))
    
    test_i <- which(folds_i == k)
    trainD <- D[!(D$Index %in% test_i), ]
    testD <- D[D$Index %in% test_i, ]
    
    if(method == 'NW') {
      # NW estimator ####
      mu_NW = apply(matrix(testD$Time, ncol = 1), 1, local_cons, trainD$Time, trainD$Y, h = h)
      SE = SE + sum( (testD$Y - mu_NW)^2 )
    }
    
    if(method == c('RAEL')) {
      
      mu_initial_vec = apply(matrix(testD$Time, ncol = 1), 1, local_cons, trainD$Time, trainD$Y, h = h_NW)
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x), 1, local_cons, trainD$Time, trainD$Y, h = h_NW) )
      # Bias-corrected local constant estimator ####
      # t0: grid time point t to estimate
      # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
      # Time list (each element is a m_i-dim vector for time record)
      # Y list (each element is a m_i-dim vector for response record)
      # h: bandwidth
      # mu_initial_vec: initial estimator as a vector for bias correction
      # mu_initial_list: initial estimator as a list for bias correction
      BiasCorrected_local_cons_pe = function(t0, Time, Y, h, mu_initial_vec, mu_initial_list){
        SampleSize = length(Time)
        Mi = sapply(Time, length)
        Kvec = lapply(Time, function(x) diag(Ki(x-t0, h)))
        Numerator = lapply(c(1:SampleSize), function(x) ( Y[[x]] - c(mu_initial_list[[x]] - mu_initial_vec[which.min(abs(testD$Time-t0))]) ) * Kvec[[x]] )
        Numerator = unlist(Numerator)
        Denominator = unlist(Kvec)
        return(sum(Numerator)/sum(Denominator))
      }
      mu_RAEL = apply(matrix(testD$Time, ncol = 1), 1, BiasCorrected_local_cons_pe, DataInfo$Time[-test_i], DataInfo$Y[-test_i], h = h, mu_initial_vec, mu_initial_list)
      SE = SE + sum( (testD$Y - mu_RAEL)^2 )
    }
    
    if(method == 'll') {
      mu_ll = apply(matrix(testD$Time, ncol = 1), 1, local_linear, trainD$Time, trainD$Y, h = h)
      SE = SE + sum( (testD$Y - mu_ll)^2 )
    }
    
    if(method == 'cubic') {
      mu_cubic = locpol(Y~Time, data=trainD, weig=rep(1,nrow(trainD)), bw=h,kernel=EpaK,deg=3, xeval=testD$Time)$lpFit$Y
      SE = SE + sum( (testD$Y - mu_cubic)^2 )
    }
    
    if(method == c('rc_cov')) {
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x, ncol = 1), 1, local_cons, tvec = trainD$Time, Yvec = trainD$Y, h = h_NW))
      mu_rc_cov = apply(matrix(testD$Time, ncol = 1), 1, refined_cons_cov, trainD, DataInfo$Time[-c(test_i)], DataInfo$Y[-c(test_i)], DataInfo$V_cov[-c(test_i)], h = h, mu_initial_list)
      SE = SE + sum( (testD$Y - mu_rc_cov)^2 )
    }
    
    if(method == c('RBEL')) {
      
      mu_local_cubic_initial_list = lapply(DataInfo$Time[-test_i], function(x) locpol(Y~Time, data=trainD, weig=rep(1,nrow(trainD)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
      
      mu_initial_vec = apply(matrix(testD$Time, ncol = 1), 1, local_cons, trainD$Time, trainD$Y, h = h_NW)
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x), 1, local_cons, trainD$Time, trainD$Y, h = h_NW) )
      
      # Bias-corrected Refined constant estimator (using working unstructured covariance) ####
      # t0: grid time point t to estimate
      # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
      # Time list (each element is a m_i-dim vector for time record)
      # Y list (each element is a m_i-dim vector for response record)
      # V_cov list (each element is an m_i by m_i inverse matrix of cov_i+1(j ==  k)*lambda^2)
      # h: bandwidth
      # mu_initial_vec: initial estimator as a vector for bias correction
      # mu_initial_list: initial estimator as a list for bias correction
      # mu_local_cubic_initial_list: initial local cubic estimator as a list
      BiasCorrected_refined_cons_cov_pe = function(t0, D, Time, Y, V_cov, h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list){
        
        SampleSize = length(Time)
        Mi = sapply(Time, length)
        
        Kmat =  lapply(Time, function(x) Ki(x-t0, h))
        
        Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + (h * Kmat[[x]]) %*% matrix(mu_initial_vec[which.min(abs(testD$Time-t0))] - mu_initial_list[[x]]) ) )
        Numerator = unlist(Numerator)
        
        Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
        Denominator = unlist(Denominator)
        
        return(sum(Numerator)/sum(Denominator))
      }
      
      mu_RBEL = apply(matrix(testD$Time, ncol = 1), 1, BiasCorrected_refined_cons_cov_pe, trainD, DataInfo$Time[-test_i], DataInfo$Y[-test_i], DataInfo$V_cov[-test_i], h = h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list)
      
      SE = SE + sum( (testD$Y - mu_RBEL)^2 )
    }
    
    if(method == c('RBEL_ind')) {
      
      mu_local_cubic_initial_list = lapply(DataInfo$Time[-test_i], function(x) locpol(Y~Time, data=trainD, weig=rep(1,nrow(trainD)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
      
      mu_initial_vec = apply(matrix(testD$Time, ncol = 1), 1, local_cons, trainD$Time, trainD$Y, h = h_NW)
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x), 1, local_cons, trainD$Time, trainD$Y, h = h_NW) )
      
      # Bias-corrected Refined constant estimator (using working independence) ####
      # t0: grid time point t to estimate
      # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
      # Time list (each element is a m_i-dim vector for time record)
      # Y list (each element is a m_i-dim vector for response record)
      # V_ind list (each element is an m_i by m_i inverse matrix of var_i+1(j ==  k)*lambda^2)
      # h: bandwidth
      # mu_initial_vec: initial estimator as a vector for bias correction
      # mu_initial_list: initial estimator as a list for bias correction
      # mu_local_cubic_initial_list: initial local cubic estimator as a list
      BiasCorrected_refined_cons_ind_pe = function(t0, D, Time, Y, V_ind, h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list){
        
        SampleSize = length(Time)
        Mi = sapply(Time, length)
        
        Kmat =  lapply(Time, function(x) Ki(x-t0, h))
        
        Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_ind[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + (h * Kmat[[x]]) %*% matrix(mu_initial_vec[which.min(abs(testD$Time-t0))] - mu_initial_list[[x]]) ) )
        Numerator = unlist(Numerator)
        
        Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_ind[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
        Denominator = unlist(Denominator)
        
        return(sum(Numerator)/sum(Denominator))
      }
      
      mu_RBEL = apply(matrix(testD$Time, ncol = 1), 1, BiasCorrected_refined_cons_ind_pe, trainD, DataInfo$Time[-test_i], DataInfo$Y[-test_i], DataInfo$V_ind[-test_i], h = h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list)
      
      SE = SE + sum( (testD$Y - mu_RBEL)^2 )
    }
    
    if(method == c('RBEL_iid')) {
      
      mu_local_cubic_initial_list = lapply(DataInfo$Time[-test_i], function(x) locpol(Y~Time, data=trainD, weig=rep(1,nrow(trainD)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
      
      mu_initial_vec = apply(matrix(testD$Time, ncol = 1), 1, local_cons, trainD$Time, trainD$Y, h = h_NW)
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x), 1, local_cons, trainD$Time, trainD$Y, h = h_NW) )
      
      # Bias-corrected Refined constant estimator (using working independence) ####
      # t0: grid time point t to estimate
      # D: Data set with 3 columns: [ Time | Y | Index ], each column is a n-dim vector
      # Time list (each element is a m_i-dim vector for time record)
      # Y list (each element is a m_i-dim vector for response record)
      # V_iid list (each element is an m_i by m_i identity matrix)
      # h: bandwidth
      # mu_initial_vec: initial estimator as a vector for bias correction
      # mu_initial_list: initial estimator as a list for bias correction
      # mu_local_cubic_initial_list: initial local cubic estimator as a list
      BiasCorrected_refined_cons_iid_pe = function(t0, D, Time, Y, V_iid, h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list){
        
        SampleSize = length(Time)
        Mi = sapply(Time, length)
        
        Kmat =  lapply(Time, function(x) Ki(x-t0, h))
        
        Numerator = lapply(c(1:SampleSize), function(x) t(one(Mi[x])) %*% Kmat[[x]] %*% V_iid[[x]] %*% (matrix(Y[[x]]) - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + (h * Kmat[[x]]) %*% matrix(mu_initial_vec[which.min(abs(testD$Time-t0))] - mu_initial_list[[x]]) ) )
        Numerator = unlist(Numerator)
        
        Denominator = lapply(c(1:SampleSize), function(x) h * t(one(Mi[x])) %*% Kmat[[x]] %*% V_iid[[x]] %*% Kmat[[x]] %*% one(Mi[x]) )
        Denominator = unlist(Denominator)
        
        return(sum(Numerator)/sum(Denominator))
      }
      
      mu_RBEL = apply(matrix(testD$Time, ncol = 1), 1, BiasCorrected_refined_cons_iid_pe, trainD, DataInfo$Time[-test_i], DataInfo$Y[-test_i], DataInfo$V_iid[-test_i], h = h, mu_local_cubic_initial_list, mu_initial_vec, mu_initial_list)
      
      SE = SE + sum( (testD$Y - mu_RBEL)^2 )
    }
    
    if(method == c('rl_cov')) {
      mu_initial_list = lapply(DataInfo$Time[-test_i], function(x) apply(matrix(x, ncol = 1), 1, local_linear, tvec = trainD$Time, Yvec = trainD$Y, h = h_ll))
      mu_rl_cov = apply(matrix(testD$Time, ncol = 1), 1, refined_linear_cov, trainD, DataInfo$Time[-c(test_i)], DataInfo$Y[-c(test_i)], DataInfo$V_cov[-c(test_i)], h = h, mu_initial_list)
      SE = SE + sum( (testD$Y - mu_rl_cov)^2 )
    }
    
    if(method == c('gl_cov')) {
      mu_gl_cov = global_linear_cov(testD$Time, trainD, DataInfo$Time[-c(test_i)], DataInfo$Y[-c(test_i)], DataInfo$V_cov[-c(test_i)], h_ll, h, Threshold = 0.005)
      SE = SE + sum( (testD$Y - mu_gl_cov)^2 )
    }
    
    if (is.na(SE)){
      SE = 99999
    }
    return(SE)
  }
  
  CV = SSE/sum(sapply(DataInfo$Time, length))
  if (is.na(CV)){
    CV = 99999
  }
  
  print(paste('CV =',CV))
  
  return(CV)
}


# local constant bias-corrected estimating equations with NW initial for bias correction and h_RAEL ####

local_cons_bias_corrected_estimating_equation = function(t0, Time, Y, mu_initial_vec, mu_initial_list, h, TimePoints, h_cubic, d){
  
  mu_0_t0 = mu0(d)
  
  SampleSize = length(Time)
  Mi = sapply(Time, length)
  
  Kmat =  lapply(Time, function(x) Ki(x-t0, h))
  
  Z_i_s = sapply(c(1:SampleSize), function(x){
    sum( ( matrix(Y[[x]]) - mu_0_t0 - ( mu_initial_list[[x]] - mu_initial_vec[which.min(abs(t0-TimePoints))] ) ) * diag(x = Kmat[[x]]) )
  } )
  
  names(Z_i_s) = paste('i=',1:length(Time),sep='')
  
  return(Z_i_s)
}


# refined constant bias-corrected estimating equations with NW initial for bias correction, cubic initial for refining procedure, and h_RBEL ####

refined_cons_bias_corrected_estimating_equation = function(t0, Time, Y, V_cov, mu_initial_vec, mu_initial_list, h, TimePoints, h_cubic, d){
  
  mu_0_t0 = mu0(d)
  
  SampleSize = length(Time)
  Mi = sapply(Time, length)
  
  Kmat =  lapply(Time, function(x) Ki(x-t0, h))
  
  mu_local_cubic_initial_list = lapply(Time, function(x) locpol(Y~Time, data=D, weig=rep(1,nrow(D)), bw=h_cubic,kernel=EpaK,deg=3, xeval=x)$lpFit$Y)
  
  Z_i_s = sapply(c(1:SampleSize), function(x){
    t(one(Mi[x])) %*% Kmat[[x]] %*% V_cov[[x]] %*% (matrix(Y[[x]]) - h*Kmat[[x]]%*%one(Mi[x])*mu_0_t0 - ( diag(1, nrow=Mi[x]) - (h * Kmat[[x]]) ) %*% matrix(mu_local_cubic_initial_list[[x]]) + h*Kmat[[x]]%*%one(Mi[x])*(mu_initial_vec[which.min(abs(t0-TimePoints))]-mu_initial_list[[x]]))
  } )
  
  names(Z_i_s) = paste('i=',1:length(Time),sep='')
  
  return(Z_i_s)
}