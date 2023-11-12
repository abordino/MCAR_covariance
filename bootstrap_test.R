MCAR_corr_test = function(X, alpha, B){
  source("simul_general.R")
  
  tmp = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS)
  
  Q_hat = tmp$Sigma/(1-tmp$R)
  QS_hat = list()
  QS_hat[[1]] = Q_hat[c(1,2),c(1,2)]
  QS_hat[[2]] = Q_hat[c(2,3),c(2,3)]
  QS_hat[[3]] = Q_hat[c(1,3),c(1,3)]
  
  #### rotate X
  X_tilde = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
  X_tilde[1:n, c("X1", "X2")] = t((sqrtm(QS_hat[[1]])$B)%*%(solve(sqrtm(SigmaS[[1]])$B))%*%t(X[1:n, c("X1", "X2")]))
  X_tilde[(n+1):(2*n), c("X2", "X3")] = t(sqrtm(QS_hat[[2]])$B%*%solve(sqrtm(SigmaS[[2]])$B)%*%t(X[(n+1):(2*n), c("X2", "X3")]))
  X_tilde[(2*n+1):(3*n), c("X1", "X3")] = t(sqrtm(QS_hat[[3]])$B%*%solve(sqrtm(SigmaS[[3]])$B)%*%t(X[(2*n+1):(3*n), c("X1", "X3")]))
  
  
  SigmaS_0 = list()
  SigmaS_0[[1]] = cor(X_tilde[1:n, c("X1", "X2")])
  SigmaS_0[[2]] = cor(X_tilde[(n+1):(2*n), c("X2", "X3")])
  SigmaS_0[[3]] = cor(X_tilde[(2*n+1):(3*n), c("X1", "X3")])
  
  R_hat_0 = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS_0)$R
  
  #### create bootstrap samples, and compute R^(b)
  sum_indicator = 0
  
  for (b in 1:B){
    SigmaS_b = list()
    X1_b = X_tilde[1:n, c("X1", "X2")][sample(1:n, n, replace = T),]
    SigmaS_b[[1]] = cor(X1_b)
    
    X2_b = X_tilde[(n+1):(2*n), c("X2", "X3")][sample(1:n, n, replace = T),]
    SigmaS_b[[2]] = cor(X2_b)
    
    X3_b = X_tilde[(2*n+1):(3*n), c("X1", "X3")][sample(1:n, n, replace = T),]
    SigmaS_b[[3]] = cor(X3_b)
    
    R_hat_b = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS_b)$R
    sum_indicator = sum_indicator + as.numeric(R_hat_b <= R_hat_0)
  }
  
  our_d = (1 + sum_indicator <= alpha*(B+1))
  return(our_d)
}

