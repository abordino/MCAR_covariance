library(copula)
library(missMethods)

get_SigmaS = function(X){
  
  #### create vector with indicators of NA's
  v = 1:dim(X)[2]
  pattern_indicator = unique(na.indicator(X))
  n_pattern = dim(pattern_indicator)[1]
  
  ### create sequence of patterns e.g. c(2,3) if d=3 and just X1 is missing
  patterns = list()
  for (row in 1:n_pattern){
    patterns[[row]] = v[as.logical(pattern_indicator[row,])]
  }
  
  ### add a column to X, corresponding to the pattern of missingness
  extra_col = c()
  tmp = na.indicator(X)
  for (row in 1:dim(tmp)[1]){
    extra_col = c(extra_col, paste(v[as.logical(tmp[row,])], collapse = ""))
  }
  
  X = cbind(X, extra_col)
  
  ### split the data X into multiple subset, each of which corresponds to a certain pattern
  data_pattern = list()
  for (i in 1:n_pattern){
    tmp = matrix(X[X[,dim(X)[2]] == paste(patterns[[i]], collapse = ""),], ncol = dim(X)[2])
    data_pattern[[i]] = apply(matrix(tmp[, patterns[[i]]], ncol = length(patterns[[i]])), c(1,2), as.numeric)
  }
  
  deletion = c()
  for (i in 1:n_pattern){
    if (dim(data_pattern[[i]])[1] < 3){ # put 3 in order for Little's test to work
      deletion = c(deletion, i)
    }
  }
  
  if (length(deletion) > 0){
    data_pattern = data_pattern[-deletion]
    patterns = patterns[-deletion]
  }
  
  n_pattern = n_pattern - length(deletion)
  
  SigmaS = list()
  for (i in 1:n_pattern){
    SigmaS[[i]] = cor(data_pattern[[i]])
  }
  
  my_list = list("pattern" = patterns, "n_pattern" = n_pattern,
                 "data_pattern" = data_pattern, "SigmaS" = SigmaS)
  return(my_list)
}

if (sys.nframe() == 0){
  n = 1000
  cp = claytonCopula(param = c(1), dim = 5)
  P = mvdc(copula = cp, margins = c("exp", "exp", "exp", "exp", "exp"),
           paramMargins = list(list(1), list(1), list(1), list(1), list(1)))
  X = rMvdc(n, P)
  X = delete_MCAR(X, 0.1, c(1,4,5))
  
  result = get_SigmaS(X)
  print(result$data_pattern)
  print(result$SigmaS)
  print(result$pattern)
}
