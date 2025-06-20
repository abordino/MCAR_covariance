source("MCARtest/computeRcycle.R")
source("MCARtest/find_SigmaS.R")
source("MCARtest/indexConsistency.R")

library(missMethods)
library(MASS)
library(norm)
library(corpcor)
library(pracma)
library(future.apply)
library(furrr)

mySqrtm = function(A){
  E = eigen(A)
  V = E$vectors; U = solve(V)
  if(length(E$values) == 1){
    D = as.matrix(E$values)
  }
  else{
    D = diag(E$values)
  }
  
  return(V %*% abs(D)^(1/2) %*% U)
}

corr.compTest.cycle = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### compute R^0
  ####--------------------------------------------------------------------------
  
  tmp = computeRcycle(SigmaS0)
  
  QS_hat = list()
  for (i in 1:d0){
    QS_hat[[i]] = matrix(c(1, cos(tmp$x[i]), cos(tmp$x[i]), 1), ncol = 2)
  }
  
  R_hat_0 = tmp$R
  
  ####--------------------------------------------------------------------------
  # Early return if R_hat_0 exceeds threshold
  ####--------------------------------------------------------------------------
  
  if (R_hat_0 >= 3 / 4) {
    return(as.numeric(0))
  }
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  rot_data_pattern = list()
  for (i in 1:d0){
    rot_data_pattern[[i]] = scale(data_pattern0[[i]]) %*%
      solve(mySqrtm(as.matrix(SigmaS0[[i]]))) %*%
      mySqrtm(as.matrix(QS_hat[[i]]))
  }
  
  compute_R_hat_b = function(b){
    
    SigmaS_b = list()
    for (i in 1:d0){
      n_S = 200
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      SigmaS_b[[i]] = cor(tmp_data)
    }
    
    return(computeRcycle(SigmaS_b)$R)
  }
  
  plan(multicore)  
  R_hats = future_sapply(seq_len(B), compute_R_hat_b, future.seed = TRUE)
  
  plan(sequential)   
  sum_indicator = sum(R_hats >= R_hat_0)
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
  
}

mean.consTest.cycle = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### compute M^0
  ####--------------------------------------------------------------------------
  
  M_hat_0 = M(mu_S0, patterns0) 
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  
  avmu = av(mu_S0, patterns0)
  
  rot_data_pattern = list()
  for (i in 1:d0){
    rot_data_pattern[[i]] = data_pattern0[[i]] - 
      do.call(rbind, replicate(200, mu_S0[[i]], simplify = FALSE)) +
      do.call(rbind, replicate(200, avmu[patterns0[[i]]], simplify = FALSE))
  }
  
  compute_M_hat_b = function(b){
    mu_S_b = list()
    
    for (i in 1:d0){
      n_S = 200
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      mu_S_b[[i]] = colMeans(tmp_data)
    }
    
    return(M(mu_S_b, patterns0))
  }
  
  plan(multicore)  
  M_hats = future_sapply(seq_len(B), compute_M_hat_b, future.seed = TRUE)
  
  plan(sequential)   
  sum_indicator = sum(M_hats >= M_hat_0)
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
}

var.consTest.cycle = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### rescale the data
  ####--------------------------------------------------------------------------
  
  avsigma = av(sigma_squared_S0, patterns0)
  
  X[,1] = X[,1]/sqrt(avsigma[1])
  data_pattern0[[1]][,1] = data_pattern0[[1]][,1]/sqrt(avsigma[1])
  data_pattern0[[d0]][,1] = data_pattern0[[d0]][,1]/sqrt(avsigma[1])
  
  for (j in 2:(d0-1)){
    X[,j] = X[,j]/sqrt(avsigma[j])
    data_pattern0[[j]][,1] = data_pattern0[[j]][,1]/sqrt(avsigma[j])
    data_pattern0[[j-1]][,2] = data_pattern0[[j-1]][,2]/sqrt(avsigma[j])
  }
  
  X[,d0] = X[,d0]/sqrt(avsigma[1])
  data_pattern0[[d0]][,2] = data_pattern0[[d0]][,2]/sqrt(avsigma[d0])
  data_pattern0[[d0-1]][,2] = data_pattern0[[d0-1]][,2]/sqrt(avsigma[d0])
  
  for (j in 1:d){
    sigma_squared_S0[[j]] = diag(cov(data_pattern0[[j]]))
  }
  
  V_hat_0 = V(sigma_squared_S0, patterns0) 
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  
  avsigma = av(sigma_squared_S0, patterns0)
  
  rot_data_pattern = list()
  for (i in 1:d){
    rot_data_pattern[[i]] = data_pattern0[[i]]%*%
      diag(1/sqrt(sigma_squared_S0[[i]]), ncol = 2)%*% 
      diag(sqrt(avsigma[patterns0[[i]]]), ncol = 2)
  }
  
  compute_V_hat_b = function(b){
    
    sigma_squared_S_b = list()
    data_b = list()
    for (i in 1:d0){
      n_S = 200
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      sigma_squared_S_b[[i]] = diag(cov(tmp_data))
      data_b[[i]] = tmp_data
    }
    
    avsigma = av(sigma_squared_S_b, patterns0)
    
    data_b[[1]][,1] = data_b[[1]][,1]/sqrt(avsigma[1])
    data_b[[d0]][,1] = data_b[[d0]][,1]/sqrt(avsigma[1])
    
    for (j in 2:(d0-1)){
      data_b[[j]][,1] = data_b[[j]][,1]/sqrt(avsigma[j])
      data_b[[j-1]][,2] = data_b[[j-1]][,2]/sqrt(avsigma[j])
    }
    
    data_b[[d0]][,2] = data_b[[d0]][,2]/sqrt(avsigma[d0])
    data_b[[d0-1]][,2] = data_b[[d0-1]][,2]/sqrt(avsigma[d0])
    
    for (j in 1:d0){
      sigma_squared_S_b[[j]] = diag(cov(data_b[[j]]))
    }
    
    return(V(sigma_squared_S_b, patterns0))
  }
  
  plan(multicore)  
  V_hats = future_sapply(seq_len(B), compute_V_hat_b, future.seed = TRUE)
  
  plan(sequential)   
  sum_indicator = sum(V_hats >= V_hat_0)
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
}