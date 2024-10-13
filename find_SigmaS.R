library(copula)
library(missMethods)
library(misty)

get_SigmaS = function(X, min_diff){
  
  #----------------------------------------------------------------------------------------
  # create vector with indicators of NA's
  #----------------------------------------------------------------------------------------
  d = dim(X)[2]
  v = 1:d
  pattern_indicator = unique(na.indicator(X))
  n_pattern = dim(pattern_indicator)[1]
  
  #----------------------------------------------------------------------------------------
  ### create sequence of patterns e.g. c(2,3) if d=3 and just X1 is missing
  #----------------------------------------------------------------------------------------
  patterns = list()
  for (row in 1:n_pattern){
    patterns[[row]] = v[as.logical(pattern_indicator[row,])]
  }
  
  #----------------------------------------------------------------------------------------
  ### add a column to X, corresponding to the pattern of missingness
  #----------------------------------------------------------------------------------------
  tmp = na.indicator(X)
  extra_col = numeric(length = dim(tmp)[1])
  for (row in 1:dim(tmp)[1]){
    extra_col[row] = paste(v[as.logical(tmp[row,])], collapse = "")
  }
  
  X = cbind(X, extra_col)
  
  #----------------------------------------------------------------------------------------
  ### split the data X into multiple subset, each of which corresponds to a certain pattern
  #----------------------------------------------------------------------------------------
  data_pattern = list()
  for (i in 1:n_pattern){
    tmp = matrix(X[X[,dim(X)[2]] == paste(patterns[[i]], collapse = ""),], ncol = dim(X)[2])
    data_pattern[[i]] = apply(matrix(tmp[, patterns[[i]]], 
                                     ncol = length(patterns[[i]])), c(1,2), as.numeric)
  }
  
  #----------------------------------------------------------------------------------------
  ####### remove pattern with sample size too small
  #----------------------------------------------------------------------------------------
  deletion = c()
  for (i in 1:n_pattern){
    if (dim(data_pattern[[i]])[1] <= dim(data_pattern[[i]])[2] + min_diff){
      deletion = c(deletion, i)
    }
  }

  if (length(deletion) > 0){
    data_pattern = data_pattern[-deletion]
    patterns = patterns[-deletion]
  }

  n_pattern = n_pattern - length(deletion)
  
  #----------------------------------------------------------------------------------------
  # return output as a list
  #----------------------------------------------------------------------------------------
  mu_S = list()
  C_S = list()
  sigma_squared_S = list()
  SigmaS = list()
  n_S = list()
  for (i in 1:n_pattern){
    n_S[[i]] = dim(data_pattern[[i]])[1]
    mu_S[[i]] = colMeans(data_pattern[[i]])
    C_S[[i]] = var(data_pattern[[i]])
    sigma_squared_S[[i]] = diag(C_S[[i]])
    SigmaS[[i]] = cov2cor(C_S[[i]])
  }
  
  my_list = list("n_S" = n_S, "patterns" = patterns, "n_pattern" = n_pattern,
                 "data_pattern" = data_pattern, "mu_S" = mu_S,
                 "C_S" = C_S, "sigma_squared_S" = sigma_squared_S, 
                 "SigmaS" = SigmaS, "ambient_dimension" = d)
  return(my_list)
}

if (sys.nframe() == 0){
  n = 1000
  cp = claytonCopula(param = c(1), dim = 5)
  P = mvdc(copula = cp, margins = c("exp", "exp", "exp", "exp", "exp"),
           paramMargins = list(list(1), list(1), list(1), list(1), list(1)))
  X = rMvdc(n, P)
  X = delete_MCAR(X, 0.1, c(1,4,5))
  
  get_SigmaS(X)
}
