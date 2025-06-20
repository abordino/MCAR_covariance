source("MCARtest/find_SigmaS.R")

library(missMethods)
library(MASS)
library(norm)
library(pracma)

av = function(x_S, patterns){
  n_pattern = length(patterns)
  d = max(unlist(patterns))
  
  res = numeric(length = d)
  for (j in 1:d){
    num = 0
    den = 0
    for (i in 1:n_pattern){
      if (j %in% patterns[[i]]){
        num = num + x_S[[i]][which(patterns[[i]] == j)[[1]]]
        den = den + 1
      }
    }
    
    res[j] = num/den
  }
  
  return(res)
}

M = function(mu_S, patterns){
  
  av_x = av(mu_S, patterns)
  n_pattern = length(patterns)
  d = max(unlist(patterns))
  
  num = 0
  for(i in 1:n_pattern){
    cand = Norm(mu_S[[i]] - av_x[patterns[[i]]], p = 1)
    if (cand > num){
      num = cand
    }
  } 
  
  den = 0
  for(i in 1:n_pattern){
    cand = length(patterns[[i]])
    if (cand > den){
      den = cand
    }
  }
  
  return(num/den)
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
