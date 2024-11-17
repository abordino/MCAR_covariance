source("computeR.R")
source("find_SigmaS.R")
source("indexConsistency.R")
library(missMethods)
library(MASS)
library(norm)
library(corpcor)
library(pracma)

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

corr.compTest = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### compute R^0
  ####--------------------------------------------------------------------------
  
  result = get_SigmaS(X, min_diff = 1)
  d = result$ambient_dimension; N_S = result$n_S; n = sum(unlist(N_S))
  n_pattern = result$n_pattern; patterns = result$pattern
  
  SigmaS = result$SigmaS 
  
  data_pattern = result$data_pattern
  
  ####--------------------------------------------------------------------------
  #### compute c
  ####--------------------------------------------------------------------------
  c = 1
  for(i in 1:n_pattern){
    cand = min(eigen(SigmaS[[i]])$values)
    if (cand < c){
      c = cand
    }
  }
  
  ####--------------------------------------------------------------------------
  
  tmp = computeR(patterns, SigmaS)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  R_hat_0 = tmp$R
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = scale(data_pattern[[i]]) %*%
      solve(mySqrtm(as.matrix(SigmaS[[i]]))) %*%
      mySqrtm(as.matrix(QS_hat[[i]]))
  }
  
  sum_indicator = 0
  for (b in 1:B){
    r_ind = 0
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:n_pattern){
      n_S = dim(rot_data_pattern[[i]])[1]
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      while(dim(unique(tmp_data))[1] <= dim(tmp_data)[2]){
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      }
      X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
      r_ind = r_ind + n_S
    }
    X = as.matrix(X[1:r_ind,])
    
    ####--------------------------------------------------------------------------
    #### compute R^b
    ####--------------------------------------------------------------------------
    result = get_SigmaS(X, min_diff = 1)
    SigmaS_b = result$SigmaS
    
    R_hat_b = computeR.reg(patterns, SigmaS_b, 4/c)$R
    
    if (R_hat_b >= R_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
  
}

mean.consTest = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### compute M^0
  ####--------------------------------------------------------------------------
  
  result = get_SigmaS(X, min_diff = 10)
  d = result$ambient_dimension; N_S = result$n_S; n = sum(unlist(N_S))
  n_pattern = result$n_pattern; patterns = result$pattern
  
  mu_S = result$mu_S
  
  data_pattern = result$data_pattern
  
  ####--------------------------------------------------------------------------
  
  M_hat_0 = M(mu_S, patterns) 
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  
  avmu = av(mu_S, patterns)
  
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = data_pattern[[i]] - 
      do.call(rbind, replicate(dim(data_pattern[[i]])[1], mu_S[[i]], simplify = FALSE)) +
      do.call(rbind, replicate(dim(data_pattern[[i]])[1], avmu[patterns[[i]]], simplify = FALSE))
  }
  
  sum_indicator = 0
  for (b in 1:B){
    
    r_ind = 0
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:n_pattern){
      n_S = dim(rot_data_pattern[[i]])[1]
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      while(dim(unique(tmp_data))[1] <= dim(tmp_data)[2]){
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      }
      X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
      r_ind = r_ind + n_S
    }
    
    X = as.matrix(X[1:r_ind,])
    
    ####--------------------------------------------------------------------------
    #### compute M^b
    ####--------------------------------------------------------------------------
    
    result = get_SigmaS(X, min_diff = 10)
    mu_S_b = result$mu_S
    
    M_hat_b = M(mu_S_b, patterns)
    
    if (M_hat_b >= M_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
}

var.consTest = function(X, B){
  
  ####--------------------------------------------------------------------------
  #### rescale the data
  ####--------------------------------------------------------------------------
  
  result = get_SigmaS(X, min_diff = 10)
  d = result$ambient_dimension; N_S = result$n_S; n = sum(unlist(N_S))
  n_pattern = result$n_pattern; patterns = result$pattern
  avsigma = av(result$sigma_squared_S, patterns)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(avsigma[j])
  }
  
  ####--------------------------------------------------------------------------
  #### compute V^0
  ####--------------------------------------------------------------------------
  result = get_SigmaS(X, min_diff = 10)
  sigma_squared_S = result$sigma_squared_S
  
  data_pattern = result$data_pattern
  
  V_hat_0 = V(sigma_squared_S, patterns) 
  
  ####--------------------------------------------------------------------------
  #### rotate X, to make it look like it's from H0
  ####--------------------------------------------------------------------------
  
  avsigma = av(sigma_squared_S, patterns)
  
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = data_pattern[[i]]%*%
      diag(1/sqrt(sigma_squared_S[[i]]), ncol = length(sigma_squared_S[[i]]))%*% 
      diag(sqrt(avsigma[patterns[[i]]]), ncol = length(avsigma[patterns[[i]]]))
  }
  
  sum_indicator = 0
  for (b in 1:B){
    
    r_ind = 0
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:n_pattern){
      n_S = dim(rot_data_pattern[[i]])[1]
      tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      while(dim(unique(tmp_data))[1] <= dim(tmp_data)[2]){
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
      }
      X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
      r_ind = r_ind + n_S
    }
    
    X = as.matrix(X[1:r_ind,])
    
    ####--------------------------------------------------------------------------
    #### rescale data and compute V^b
    ####--------------------------------------------------------------------------
    
    result = get_SigmaS(X, min_diff = 10)
    avsigma = av(result$sigma_squared_S, result$patterns)
    
    for (j in 1:d){
      X[,j] = X[,j]/sqrt(avsigma[j])
    }
    
    result = get_SigmaS(X, min_diff = 10)
    sigma_squared_S_b = result$sigma_squared_S
    
    V_hat_b = V(sigma_squared_S_b, patterns)
    
    if (V_hat_b >= V_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  return(p_hat)
}
