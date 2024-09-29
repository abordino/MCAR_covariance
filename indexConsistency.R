source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)

#### compute av_\mu and av_\sigma
av_sigma = function(value_S, patterns){
  
  n_pattern = length(patterns)
  d = max(unlist(patterns))
  
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

av_mu = function(mu_S, patterns, n_S){
  n_pattern = length(patterns)
  d = max(unlist(patterns))
  
  res = numeric(length = d)
  for (j in 1:d){
    num = 0
    den = 0
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]){
        num = num + n_S[[i]]*mu_S[[i]][which(patterns[[i]] == j)[[1]]]
        den = den + n_S[[i]]
      }
    }
    
    res[j] = num/den
  }
  
  return(res)
}

M = function(mu_S, patterns, MahS, n_S){

  av_mu = av_mu(mu_S, patterns, n_S)
  n_pattern = length(patterns)
  d = max(unlist(patterns))

  sum = 0
  for(i in 1:n_pattern){
    dim_S = length(patterns[[i]])
    MahS[[i]] = diag(dim_S)
    sum = sum + n_S[[i]]*sum(abs(mu_S[[i]] - av_mu[patterns[[i]]]))
  }

  return(sum/sum(unlist(n_S)))
}

V = function(sigma_squared_S, patterns){

  n_pattern = length(patterns)
  d = max(unlist(patterns))

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
