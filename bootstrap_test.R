source("computeR.R")
source("little_test.R")
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

MCAR_corrTest = function(X, alpha, B){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension
  n_pattern = result$n_pattern; patterns = result$pattern
  C_S = result$C_S; SigmaS = result$SigmaS
  data_pattern = result$data_pattern
  
  tmp = computeR(patterns, SigmaS)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  R_hat_0 = tmp$R
  
  #### rotate X, to make it look like it's from H0
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = t((mySqrtm(as.matrix(QS_hat[[i]])))%*%(solve(mySqrtm(as.matrix(SigmaS[[i]]))))%*%
                                t(scale(data_pattern[[i]])))
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
    
    result = get_SigmaS(X)
    SigmaS_b = result$SigmaS
    
    # R_hat_b = computeR(patterns, SigmaS_b)$R
    R_hat_b = computeR(patterns, SigmaS_b)$R
    
    if (R_hat_b >= R_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
  
}

MCAR_covTest = function(X, alpha, B){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension; N_S = result$n_S
  n_pattern = result$n_pattern; patterns = result$pattern
  av_sigma = av_sigma(result$sigma_squared_S, patterns)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  result = get_SigmaS(X)
  C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  data_pattern = result$data_pattern
  
  tmp = computeR(patterns, SigmaS)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  T_hat_0 = tmp$R + V(sigma_squared_S, patterns)
  
  #### rotate X, to make it look like it's from H0
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = t((mySqrtm(as.matrix(QS_hat[[i]])))%*%
                                (solve(mySqrtm(as.matrix(C_S[[i]]))))%*%
                                t(data_pattern[[i]]))
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
    
    result = get_SigmaS(X)
    av_sigma = av_sigma(result$sigma_squared_S, result$patterns)
    
    for (j in 1:d){
      X[,j] = X[,j]/sqrt(av_sigma[j])
    }
    
    result = get_SigmaS(X)
    SigmaS_b = result$SigmaS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S
    patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
    
    tmp = computeR(patterns, SigmaS_b)
    
    T_hat_b = tmp$R + V(sigma_squared_S_b, patterns)
    
    if (T_hat_b >= T_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}

MCAR_meanTest = function(X, alpha, B){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension; N_S = result$n_S
  n_pattern = result$n_pattern; patterns = result$pattern
  mu_S = result$mu_S; av_mu = av_mu(mu_S, patterns, N_S)
  C_S = result$C_S
  data_pattern = result$data_pattern
  
  T_hat_0 = M(mu_S, patterns, C_S, N_S) 
  
  #### rotate X, to make it look like it's from H0
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = data_pattern[[i]] - mu_S[[i]] +  av_mu[patterns[[i]]]
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
    
    result = get_SigmaS(X)
    C_S_b = result$C_S
    mu_S_b = result$mu_S
    patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
    
    T_hat_b = M(mu_S_b, patterns, C_S_b, N_S)
    
    if (T_hat_b >= T_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}

MCAR_meancovTest = function(X, alpha, B){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension;  N_S = result$n_S
  n_pattern = result$n_pattern; patterns = result$pattern
  av_sigma = av_sigma(result$sigma_squared_S, patterns)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  result = get_SigmaS(X)
  mu_S = result$mu_S; av_mu = av_mu(mu_S, patterns, N_S)
  C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  data_pattern = result$data_pattern
  
  tmp = computeR(patterns, SigmaS)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  T_hat_0 = tmp$R + V(sigma_squared_S, patterns, C_S, N_S) + M(mu_S, patterns, C_S, N_S) 
  
  #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((mySqrtm(as.matrix(QS_hat[[i]])))%*%
                                  (solve(mySqrtm(as.matrix(C_S[[i]]))))%*%
                                  t(data_pattern[[i]] - mu_S[[i]] +  av_mu[patterns[[i]]]))
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
        
      result = get_SigmaS(X)
      av_sigma = av_sigma(result$sigma_squared_S, result$patterns)
        
      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }
      
      result = get_SigmaS(X)
      mu_S_b = result$mu_S
      SigmaS_b = result$SigmaS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      
      # tmp = computeR(patterns, SigmaS_b)
      tmp = computeR(patterns, SigmaS_b)
      
      T_hat_b = tmp$R + V(sigma_squared_S_b, patterns, C_S_b, N_S) +
        M(mu_S_b, patterns, C_S_b, N_S)
      
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
    }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}
