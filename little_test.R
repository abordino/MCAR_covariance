little_test = function(X, alpha){
  s = prelim.norm(as.matrix(X))
  thetahat = em.norm(s)
  mu_true = getparam.norm(s,thetahat,corr=TRUE)$mu
  Sigma_true = getparam.norm(s,thetahat,corr=TRUE)$r     # possibly also right-multiply by diag(getparam.norm(s,thetahat,corr=TRUE)$sdv)%*%
  
  d_squared = n*t((colMeans(X1) - mu_true[c(1,2)]))%*%solve(n*Sigma_true[c(1,2),c(1,2)]/(n-1))%*%t(t((colMeans(X1) - mu_true[c(1,2)]))) +
    n*t((colMeans(X2) - mu_true[c(2,3)]))%*%solve(n*Sigma_true[c(2,3),c(2,3)]/(n-1))%*%t(t((colMeans(X2) - mu_true[c(2,3)]))) +
    n*t((colMeans(X3) - mu_true[c(1,3)]))%*%solve(n*Sigma_true[c(1,3),c(1,3)]/(n-1))%*%t(t((colMeans(X3) - mu_true[c(1,3)])))
  
  d_cov = n*(sum(diag(cor(X1)%*%solve(Sigma_true[c(1,2),c(1,2)]))) - 2 - log(det(cor(X1))) + log(det(Sigma_true[c(1,2),c(1,2)]))) +
    n*(sum(diag(cor(X2)%*%solve(Sigma_true[c(2,3),c(2,3)]))) - 2 - log(det(cor(X2))) + log(det(Sigma_true[c(2,3),c(2,3)]))) +
    n*(sum(diag(cor(X3)%*%solve(Sigma_true[c(1,3),c(1,3)]))) - 2 - log(det(cor(X3))) + log(det(Sigma_true[c(1,3),c(1,3)])))
  
  d_aug = d_squared + d_cov
  little_d = (d_aug > qchisq(1-alpha, 6))
  return(little_d)
}