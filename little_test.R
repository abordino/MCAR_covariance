little_test = function(X, alpha, type="mean&cov"){
  s = prelim.norm(as.matrix(X))
  thetahat = em.norm(s)
  mu_true = getparam.norm(s,thetahat,corr=TRUE)$mu
  Sigma_true = getparam.norm(s,thetahat,corr=TRUE)$r
  
  result = get_SigmaS(X)
  SigmaS = result$SigmaS
  patterns = result$pattern
  n_pattern = result$n_pattern
  data_pattern = result$data_pattern
  d = result$ambient_dimension
  
  if(type == "mean"){
    d_squared = 0
    df = -d
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])
      
      d_squared = d_squared + n_S*t(x_S)%*%solve(n_S*L_S/(n_S-1))%*%t(t(x_S))
      df = df + card_S
      print(df)
    }
    little_d = (d_squared > qchisq(1-alpha, df))
  }
  
  else if(type == "cov"){
    d_cov = 0
    df = -d*(d+1)/2
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])
      
      d_cov = d_cov + n_S*(sum(diag(Sigma_S%*%solve(L_S))) - card_S - log(det(Sigma_S)) + log(det(L_S)))
      df = df + card_S*(card_S+1)/2
      print(df)
    }
    
    little_d = (d_cov > qchisq(1-alpha, df))
  }
    
  else{
    d_aug = 0
    df = -d*(d+3)/2
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = Sigma_true[patterns[[i]],patterns[[i]]]
      Sigma_S = cor(data_pattern[[i]])
      
      d_aug = d_aug + n_S*t(x_S)%*%solve(n_S*L_S/(n_S-1))%*%t(t(x_S)) + 
        n_S*(sum(diag(Sigma_S%*%solve(L_S))) - card_S - log(det(Sigma_S)) + log(det(L_S)))
      df = df + card_S*(card_S+3)/2
      print(df)
    }
    
    little_d = (d_aug > qchisq(1-alpha, df))
  }
  
  return(little_d)

}
