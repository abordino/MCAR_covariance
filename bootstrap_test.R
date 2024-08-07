source("computeR.R")
source("little_test.R")
source("find_SigmaS.R")
source("indexConsistency.R")
library(missMethods)
library(MASS)
library(norm)
library(corpcor)

MCAR_corr_test = function(X, alpha, B, type="np"){
  
  result = get_SigmaS(X)
  SigmaS = result$SigmaS; patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
  
  for (i in 1:n_pattern){
    if (det(as.matrix(SigmaS[[i]])) < 1e-10){
      print("SigmaS is singular!")
      return(NA)
    }
  }
  
  tmp = computeR(patterns, SigmaS)
  R_hat_0 = tmp$R
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
    if (det(as.matrix(QS_hat[[i]])) < 1e-10){
      return(NA)
    }
  }
  
  if (type=="np"){
    
    #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((sqrtm(as.matrix(QS_hat[[i]]))$B)%*%(solve(sqrtm(as.matrix(SigmaS[[i]]))$B))%*%
                                  t(scale(data_pattern[[i]])))
    }

    sum_indicator = 0
    for (b in 1:B){
      SigmaS_b = list()
      for (i in 1:n_pattern){
        n_S = dim(rot_data_pattern[[i]])[1]
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        if(dim(unique(tmp_data))[1] < 2){
          tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        }
        SigmaS_b[[i]] = cov2cor(var(tmp_data))
      }
      
      R_hat_b = computeR(patterns, SigmaS_b)$R
      
      if (R_hat_b >= R_hat_0){
        sum_indicator = sum_indicator + 1
      }
    }
  }
else if (type=="p"){
  sum_indicator = 0
  
  for (b in 1:B){
    SigmaS_b = list()
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(QS_hat[[i]])[1]
      tmp_data = as.matrix(mvrnorm(n_S, rep(0, card_S), QS_hat[[i]]))
      SigmaS_b[[i]] = cov2cor(var(tmp_data))
    }
    
    R_hat_b = computeR(patterns, SigmaS_b)$R
    
    if (R_hat_b >= R_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
}
else {
  print("An error occured")
}

  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}

MCAR_meancovTest = function(X, alpha, B, type="np"){
  
  result = get_SigmaS(X); d = result$ambient_dimension; av_sigma = compute_av("var", X)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  result = get_SigmaS(X)
  av_mu = compute_av("mean", X)
  muS = result$muS; C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
  
  for (i in 1:n_pattern){
    if (det(as.matrix(C_S[[i]])) < 1e-10){
      print("SigmaS is singular!")
      return(NA)
    }
  }
  
  tmp = computeR(patterns, SigmaS)
  T_hat_0 = tmp$R + M(muS, patterns) + V(sigma_squared_S, patterns)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
    if (det(as.matrix(QS_hat[[i]])) < 1e-10){
      print("QS is singular!")
      return(NA)
    }
  }
  
  if (type=="np"){
    
    #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((sqrtm(as.matrix(QS_hat[[i]]))$B)%*%(solve(sqrtm(as.matrix(C_S[[i]]))$B))%*%
                                  t(data_pattern[[i]] - muS[[i]] +  av_mu[patterns[[i]]]))
    }
    
    sum_indicator = 0
    for (b in 1:B){
      
      r_ind = 0
      X = data.frame(matrix(nrow = d*n, ncol = d))
      for (i in 1:n_pattern){
        n_S = dim(rot_data_pattern[[i]])[1]
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        while(dim(unique(tmp_data))[1] < 2){
          tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        }
        X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
        r_ind = r_ind + n_S
      }
      X = as.matrix(X[1:r_ind,])
        
      result = get_SigmaS(X); av_sigma = compute_av("var", X)
        
      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }
      
      result = get_SigmaS(X)
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      muS_b = result$muS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S; SigmaS_b = result$SigmaS
      
      T_hat_b = computeR(patterns, SigmaS_b)$R + M(muS_b, patterns) + V(sigma_squared_S_b, patterns) 
      
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
    }
  }
  else if (type=="p"){

    sum_indicator = 0
    for (b in 1:B){
      
      r_ind = 0
      X = data.frame(matrix(nrow = d*n, ncol = d))
      for (i in 1:n_pattern){
        n_S = dim(data_pattern[[i]])[1]
        tmp_data = as.matrix(mvrnorm(n_S, av_mu[patterns[[i]]], as.matrix(QS_hat[[i]])))
        X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
        r_ind = r_ind + n_S
      }
      X = as.matrix(X[1:r_ind,])
      
      result = get_SigmaS(X); av_sigma = compute_av("var", X)
      
      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }
      
      result = get_SigmaS(X)
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      muS_b = result$muS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S; SigmaS_b = result$SigmaS
      
      T_hat_b = computeR(patterns, SigmaS_b)$R + M(muS_b, patterns) + V(sigma_squared_S_b, patterns) 
      
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
    }
  }
  else {
    print("An error occured")
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}

if (sys.nframe() == 0){
  set.seed(7)
  n = 200
  MC = 30
  d = 5
  
  # Select the copula
  cp = claytonCopula(param = c(1), dim = d)
  # Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
  P = mvdc(copula = cp, margins = c(rep("lnorm",d)),
           paramMargins = rep(list(1),d) )
  data = rMvdc(n, P)
  
  
  #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
  X = delete_MCAR(data, .12, c(1,2))
  result = get_SigmaS(X); d = result$ambient_dimension; av_sigma = compute_av("var", X)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  result = get_SigmaS(X)
  av_mu = compute_av("mean", X)
  muS = result$muS; C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
  
  for (i in length(SigmaS)){
    if (min(eigen(SigmaS[[i]])$values) < 1e-10){
      print("SigmaS is singular!")
      return(NA)
    }
  }
  
  tmp = computeR(patterns, SigmaS)
  T_hat_0 = tmp$R + M(muS, patterns) + V(sigma_squared_S, patterns)
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
    #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((sqrtm(as.matrix(QS_hat[[i]]))$B)%*%(sqrtm(as.matrix(C_S[[i]]))$Binv)%*%
                                  t(data_pattern[[i]] - muS[[i]] +  av_mu[patterns[[i]]]))
    }
  
    sum_indicator = 0
      r_ind = 0
      X = data.frame(matrix(nrow = d*n, ncol = d))
      for (i in 1:n_pattern){
        n_S = dim(rot_data_pattern[[i]])[1]
        tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        if(dim(unique(tmp_data))[1] < 2){
          tmp_data = as.matrix(rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),])
        }
        X[(1+r_ind):(n_S+r_ind), patterns[[i]]] = tmp_data
        r_ind = r_ind + n_S
      }
      X = as.matrix(X[1:r_ind,])
  
      result = get_SigmaS(X); av_sigma = compute_av("var", X)
  
      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }
  
      result = get_SigmaS(X)
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      muS_b = result$muS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S; SigmaS_b = result$SigmaS
  
      T_hat_b = computeR(patterns, SigmaS_b)$R + M(muS_b, patterns) + V(sigma_squared_S_b, patterns)
  
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
}
