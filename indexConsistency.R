source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)
library(foreach)
library(doSNOW)
library(doParallel)
library(future)
library(parallel)
library(foreach)
library(doRNG) 
library(doFuture)
library(future.apply)

#### compute av_\mu and av_\sigma
compute_av = function(type, X){
  
  result = get_SigmaS(X)
  d = result$ambient_dimension
  
  if (type=="mean"){
    value_S = result$muS
  }
  else if (type=="var"){
    value_S = result$sigma_squared_S
  }
  else {
    print("An error occured")
  }
  
  patterns = result$pattern
  card_patterns = length(patterns)
  n_pattern = result$n_pattern
  data_pattern = result$data_pattern
  
  res = numeric(length = d)
  for (j in 1:d){
    tmp = 0
    card_Sj = 0
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]){
        tmp = tmp + value_S[[i]][which(patterns[[i]] == j)[[1]]]
        card_Sj = card_Sj+1
      }
    }
    
    res[j] = tmp/card_Sj
  }
  
  return(res)
}

M = function(mu_S, patterns){
  n_pattern = length(mu_S)
  
  max = 0
  couples = t(combn(1:n_pattern, 2))
  n_couples = dim(t(combn(1:n_pattern, 2)))[1]
  couples = as.list(data.frame(t(combn(1:n_pattern, 2))))
  
  for (j in 1:d){
    for (i in 1:n_couples){
      if ((j %in% patterns[[ couples$X1[i] ]])&(j %in% patterns[[ couples$X2[i] ]])) {
        candidate = abs(mu_S[[ couples$X1[i] ]][which(patterns[[couples$X1[i]]] == j)[[1]]] 
                        - mu_S[[ couples$X2[i] ]][which(patterns[[couples$X2[i]]] == j)[[1]]])
        if (candidate > max){
          max = candidate
        }
      }
    }
  }
  
  return(max)
}

V = function(sigma_squared_S, patterns){
  n_pattern = length(sigma_squared_S)
  
  min = 1
  
  for (j in 1:d){
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]) {
        candidate = sigma_squared_S[[i]][which(patterns[[i]] == j)[[1]]]
        if (candidate < min){
          min = candidate
        }
      }
    }
  }
  
  return(1-min)
}



#### test this function
if (sys.nframe() == 0){
  n = 2
  t1 = pi/2
  t3 = pi/4
  t2 = pi/3
  d = 3
  
  SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }
  
  SigmaS[[1]][1,2] = cos(t1); SigmaS[[1]][2,1] = cos(t1)
  SigmaS[[2]][1,2] = cos(t2); SigmaS[[2]][2,1] = cos(t2)
  SigmaS[[3]][1,2] = cos(t3); SigmaS[[3]][2,1] = cos(t3)
  
  X = data.frame(matrix(nrow = 3*n, ncol = 3))
  X[1:n, c(1,2)] = mvrnorm(n, c(0,0), SigmaS[[1]])
  X[(n+1):(2*n), c(2, 3)] = mvrnorm(n, c(0,0), SigmaS[[2]])
  X[(2*n+1):(3*n), c(1, 3)] = mvrnorm(n, c(0,0), SigmaS[[3]])
  X = as.matrix(X)
  
  X
  xxx = get_SigmaS(X)$patterns
  
  get_SigmaS(X)$muS
  M(get_SigmaS(X)$muS, xxx)
  av_mu = compute_av("mean", X)
  
  av_sigma = compute_av("var", X); av_sigma
  
  X_new = X
  for (j in 1:3){
    X_new[,j] = X[,j]/sqrt(av_sigma[j])
  }
  
  X_new
  get_SigmaS(X_new)$sigma_squared_S
  compute_av("var", X_new)
  
  V(get_SigmaS(X_new)$sigma_squared_S, xxx)
}


# 
# alpha = 0.05
# n = 200
# MC = 10
# d = 5
# 
# # Select the copula
# cp = claytonCopula(param = c(1), dim = d)
# # Generate the multivariate distribution (in this case it is just bivariate) with normal and t marginals
# P = mvdc(copula = cp, margins = c(rep("exp",d)),
#          paramMargins = rep(list(1),d) )
# data = rMvdc(n, P)
# 
# X = delete_MCAR(data, 0.03, c(1,3,5))
# 
# X
# patterns = get_SigmaS(X)$patterns
# 
# mu_S = get_SigmaS(X)$muS
# 
# n_pattern = length(mu_S)
# 
# max = 0
# couples = t(combn(1:n_pattern, 2))
# n_couples = dim(t(combn(1:n_pattern, 2)))[1]
# couples = as.list(data.frame(t(combn(1:n_pattern, 2))))
# 
# for (j in 1:d){
#   for (i in 1:n_couples){
#     print(paste("j = ", j))
#     print(paste("couple number ", i))
#     print(couples$X1[i])
#     print(couples$X2[i])
#     print(patterns[[ couples$X1[i] ]])
#     print(patterns[[ couples$X2[i] ]])
#     if ((j %in% patterns[[ couples$X1[i] ]])&(j %in% patterns[[ couples$X2[i] ]])) {
#       candidate = abs(mu_S[[ couples$X1[i] ]][which(patterns[[couples$X1[i]]] == j)[[1]]] 
#                       - mu_S[[ couples$X2[i] ]][which(patterns[[couples$X2[i]]] == j)[[1]]])
#       if (candidate > max){
#         max = candidate
#       }
#     }
#   }
# }
# 
# M(get_SigmaS(X)$muS, xxx)
# av_mu = compute_av("mean", X)
# 
# av_sigma = compute_av("var", X); av_sigma
# 
# X_new = X
# for (j in 1:3){
#   X_new[,j] = X[,j]/sqrt(av_sigma[j])
# }
# 
# X_new
# get_SigmaS(X_new)$sigma_squared_S
# compute_av("var", X_new)
# 
# V(get_SigmaS(X_new)$sigma_squared_S, xxx)
