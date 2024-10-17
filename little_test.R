little_test = function(X, type="mean&cov"){
  s = prelim.norm(as.matrix(X))
  thetahat = em.norm(s)
  mu_true = getparam.norm(s,thetahat,corr=FALSE)$mu
  Sigma_true = getparam.norm(s,thetahat,corr=FALSE)$sigma
  
  result = get_SigmaS(X, min_diff = 0)
  Omega_S = result$C_S
  patterns = result$pattern; data_pattern = result$data_pattern
  n_pattern = result$n_pattern; d = result$ambient_dimension
  
  if(type == "mean"){
    d_squared = 0
    df = -d
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = as.matrix(Sigma_true[patterns[[i]], patterns[[i]]])
      
      d_squared = d_squared + n_S*t(x_S)%*%solve(L_S)%*%x_S
      df = df + card_S
    }
    p_L = pchisq(d_squared, df, lower.tail = FALSE)
  }
  
  else if(type == "cov"){
    
    d_cov = 0
    df = -d*(d+1)/2
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = as.matrix(Sigma_true[patterns[[i]], patterns[[i]]])
      
      d_cov = d_cov + n_S*(sum(diag(Omega_S[[i]]%*%solve(L_S))) - card_S -
                                  log(det(Omega_S[[i]])) + log(det(L_S)))
      df = df + card_S*(card_S+1)/2
      print(df)
    }
    
    p_L = pchisq(d_cov, df, lower.tail = FALSE)
  }
  
  else{
    
    d_aug = 0
    df = -d*(d+3)/2
    
    for (i in 1:n_pattern){
      n_S = dim(data_pattern[[i]])[1]
      card_S = dim(data_pattern[[i]])[2]
      
      x_S = colMeans(data_pattern[[i]]) - mu_true[patterns[[i]]]
      L_S = as.matrix(Sigma_true[patterns[[i]], patterns[[i]]])
      
      d_aug = d_aug + n_S*t(x_S)%*%solve(L_S)%*%x_S + 
        n_S*(sum(diag(Omega_S[[i]]%*%solve(L_S))) - card_S -
               log(det(Omega_S[[i]])) + log(det(L_S)))
      df = df + card_S*(card_S+3)/2
    }
    
    p_L = pchisq(d_aug, df, lower.tail = FALSE)
  }
  
  return(p_L)
  
}
