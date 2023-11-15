MCAR_corr_test = function(X, alpha, B){
  source("computeR.R")
  source("find_SigmaS.R")
  
  result = get_SigmaS(X)
  SigmaS = result$SigmaS
  patterns = result$pattern
  n_pattern = result$n_pattern
  data_pattern = result$data_pattern
  
  for (i in length(SigmaS)){
    if (min(eigen(SigmaS[[i]])$values) < 10^-7){
      print("SigmaS is singular!")
      return(SigmaS)
    }
  }
  
  tmp = computeR(patterns, SigmaS)
  
  R_hat_0 = tmp$R
  
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  for (i in 1:n_pattern){
    QS_hat[[i]] = Q_hat[patterns[[i]], patterns[[i]]]
  }
  
  #### rotate X, to make it look like it's from H0
  rot_data_pattern = list()
  for (i in 1:n_pattern){
    rot_data_pattern[[i]] = t((sqrtm(QS_hat[[i]])$B)%*%(solve(sqrtm(SigmaS[[i]])$B))%*%t(data_pattern[[i]]))
  }
  
  #### create bootstrap samples, and compute R^(b)
  sum_indicator = 0
  
  for (b in 1:B){
    SigmaS_b = list()
    for (i in 1:n_pattern){
      n_S = dim(rot_data_pattern[[i]])[1]
      tmp_data = rot_data_pattern[[i]][sample(1:n_S, n_S, replace = T),]
      SigmaS_b[[i]] = cor(tmp_data)
    }
    
    R_hat_b = computeR(patterns, SigmaS_b)$R
    
    if (R_hat_b >= R_hat_0){
      sum_indicator = sum_indicator + 1
    }
  }
  
  p_hat = (1+sum_indicator)/(B+1)
  decision = p_hat < alpha
  return(decision)
}