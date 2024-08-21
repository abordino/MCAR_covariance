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
  D = diag(E$values)
  return(V %*% abs(D)^(1/2) %*% U)
}

MCAR_covTest = function(X, alpha, B){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension
  n_pattern = result$n_pattern; patterns = result$pattern
  av_sigma = compute_av(result$sigma_squared_S, patterns)
  
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
    av_sigma = compute_av(result$sigma_squared_S, result$patterns)
    
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
  d = result$ambient_dimension
  n_pattern = result$n_pattern; patterns = result$pattern
  mu_S = result$mu_S; av_mu = compute_av(mu_S, patterns)
  C_S = result$C_S
  data_pattern = result$data_pattern
  
  T_hat_0 = M(mu_S, patterns, C_S) 
  
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
    
    T_hat_b = M(mu_S_b, patterns, C_S_b)
    
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
  d = result$ambient_dimension
  n_pattern = result$n_pattern; patterns = result$pattern
  av_sigma = compute_av(result$sigma_squared_S, patterns)
  
  for (j in 1:d){
    X[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  result = get_SigmaS(X)
  mu_S = result$mu_S; av_mu = compute_av(mu_S, patterns)
  C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  data_pattern = result$data_pattern
  
  tmp = computeR(patterns, SigmaS)
  Q_hat = tmp$Sigma/(1-tmp$R)
  
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  T_hat_0 = tmp$R + V(sigma_squared_S, patterns) + M(mu_S, patterns, C_S) 
  
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
      av_sigma = compute_av(result$sigma_squared_S, result$patterns)
        
      for (j in 1:d){
        X[,j] = X[,j]/sqrt(av_sigma[j])
      }
      
      result = get_SigmaS(X)
      mu_S_b = result$mu_S
      SigmaS_b = result$SigmaS; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S
      patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
      
      tmp = computeR(patterns, SigmaS_b)
      
      T_hat_b = tmp$R + V(sigma_squared_S_b, patterns) +
        M(mu_S_b, patterns, C_S_b)
      
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
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
  mu_S = result$mu_S; C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
  patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
  
  for (i in length(SigmaS)){
    if (min(eigen(SigmaS[[i]])$values) < 1e-10){
      print("SigmaS is singular!")
      return(NA)
    }
  }
  
  tmp = computeR(patterns, SigmaS)
  T_hat_0 = tmp$R + M(mu_S, patterns) + V(sigma_squared_S, patterns)
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
    #### rotate X, to make it look like it's from H0
    rot_data_pattern = list()
    for (i in 1:n_pattern){
      rot_data_pattern[[i]] = t((mySqrtm(as.matrix(QS_hat[[i]]))$B)%*%(mySqrtm(as.matrix(C_S[[i]]))$Binv)%*%
                                  t(data_pattern[[i]] - mu_S[[i]] +  av_mu[patterns[[i]]]))
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
      mu_S_b = result$mu_S; C_S_b = result$C_S; sigma_squared_S_b = result$sigma_squared_S; SigmaS_b = result$SigmaS
  
      T_hat_b = computeR(patterns, SigmaS_b)$R + M(mu_S_b, patterns, N_S) + V(sigma_squared_S_b, patterns)
  
      if (T_hat_b >= T_hat_0){
        sum_indicator = sum_indicator + 1
      }
      
      
      #### another example
      X = delete_MAR_1_to_x(data, 0.05, c(1,2), cols_ctrl = c(3, 3), x = 9)
        result = get_SigmaS(X)
        d = result$ambient_dimension
        
        av_sigma = compute_av("var", X)

        for (j in 1:d){
          X[,j] = X[,j]/sqrt(av_sigma[j])
        }

        result = get_SigmaS(X)
        av_mu = compute_av("mean", X)
        mu_S = result$mu_S; C_S = result$C_S; sigma_squared_S = result$sigma_squared_S; SigmaS = result$SigmaS
        patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
        n_S = result$n_S

        tmp = computeR(patterns, SigmaS)
        Q_hat = tmp$Sigma/(1-tmp$R)

        QS_hat = list()
        for (i in 1:n_pattern){
          QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
        }

        T_hat_0 = tmp$R + V(sigma_squared_S, patterns) + M(mu_S, patterns, C_S, n_S)
        
      
        #another example
        alpha = 0.05
        n = 300
        MC = 300
        d = 5
        
        # Select the copula
        cp = claytonCopula(param = c(1), dim = d)
        # Generate the multivariate distribution (in this case it is just bivariate) with chisqal and t marginals
        P = mvdc(copula = cp, margins = c(rep(yyy,d)),
                 paramMargins = rep(list(c(mean = 0, sd = 1)),d))
        
        data = rMvdc(n, P)
        X = delete_MAR_rank(data, p, c(1,2), cols_ctrl = c(3, 4))
        
        
        result = get_SigmaS(X)
        SigmaS = result$SigmaS; patterns = result$pattern; n_pattern = result$n_pattern; data_pattern = result$data_pattern
        
        tmp = computeR(patterns, SigmaS)
        R_hat_0 = tmp$R
        Q_hat = tmp$Sigma/(1-tmp$R)
        QS_hat = list()
        for (i in 1:n_pattern){
          QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
        }
        
        rot_data_pattern = list()
        for (i in 1:n_pattern){
          rot_data_pattern[[i]] = t((mySqrtm(as.matrix(QS_hat[[i]])))%*%(solve(mySqrtm(as.matrix(SigmaS[[i]]))))%*%
                                      t(scale(data_pattern[[i]])))
        }

}
