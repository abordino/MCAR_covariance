little_test = function(X, alpha){
  s = prelim.norm(as.matrix(X))
  thetahat = em.norm(s)
  mu_true = getparam.norm(s,thetahat,corr=TRUE)$mu
  Sigma_true = getparam.norm(s,thetahat,corr=TRUE)$r
  
  result = get_SigmaS(X)
  SigmaS = result$SigmaS
  patterns = result$pattern
  n_pattern = result$n_pattern
  data_pattern = result$data_pattern
  
  d_aug = 0
  for (i in 1:n_pattern){
    n_S = dim(data_pattern[[i]])[1]
    x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
    L_S = Sigma_true[patterns[[i]],patterns[[i]]]
    Sigma_S = cor(data_pattern[[i]])
    d_aug = d_aug + n_S*t(x_S)%*%solve(n_S*L_S/(n_S-1))%*%t(t(x_S)) + 
      n_S*(sum(diag(Sigma_S%*%solve(L_S))) - 2 - log(det(Sigma_S)) + log(det(L_S)))
  }

  little_d = (d_aug > qchisq(1-alpha, 6))
  return(little_d)
}
