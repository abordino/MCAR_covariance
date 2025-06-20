library(Rcsdp)
library(Matrix)
library(npmr)
library(matrixcalc)
library(pracma)


############ FUNCTION ###############
computeR = function(patterns=list(), SigmaS=list()) {
  
  #----------------------------------------------------------------------------------------
  ##### DEFINE PATTERN IF NOT SPECIFIED
  #----------------------------------------------------------------------------------------
  add_pattern = TRUE
  if (length(patterns) > 0){
    add_pattern = FALSE
  }
  
  i = 1
  while(add_pattern == TRUE){
    prompt1 = "Enter variable numbers (space-separated list): \n"
    tmp = as.vector(as.integer(strsplit(readline(prompt1), " ")[[1]]))
    patterns[[i]] = tmp
    
    prompt2 = "Do you want to add another pattern? (Write TRUE or FALSE): \n"
    add_pattern = as.logical((readline(prompt2)))
    i = i+1
  }
  
  #----------------------------------------------------------------------------------------
  # remove duplicates
  patterns = unique(patterns)
  card_patterns = length(patterns)
  
  d = 0 # computing d
  for (S in patterns){
    if (max(S) > d){
      d = max(S)
    }
  }
  
  
  
  #----------------------------------------------------------------------------------------
  ########## GENERATE SIGMA IF NOT SPECIFIED
  #----------------------------------------------------------------------------------------
  if (length(SigmaS) == 0){
    prompt2 = "Do you want a random SigmaS or a compatible one? (Write 1 or 2 respectively): \n"
    comp = as.numeric((readline(prompt2)))
    
    if (comp == 1){
      for(j in 1:card_patterns){
        card_S = length(patterns[[j]])
        A = matrix(runif((card_S)^2)*2-1, ncol=card_S) 
        tmp_0 = t(A) %*% A
        SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
      }
    } else if(comp == 2){
      SigmaS=list()
      A = matrix(runif(d^2)*2-1, ncol=d)
      Sigma = cov2cor(t(A) %*% A)
      for(j in 1:card_patterns){
        SigmaS[[j]]=as.matrix(Sigma[patterns[[j]],patterns[[j]]])
      }
    } else {
      return("An error occured.")
    }
  }
  
  #----------------------------------------------------------------------------------------
  ########## WRITING THE SDP PROBLEM ##########################
  #----------------------------------------------------------------------------------------
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)))
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    C[[1+j]] = matrix(0, card_S, card_S)
  }
  
  # number of constraints
  m = d-1
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    m = m + card_S*(card_S+1)/2
  }
  
  b=rep(0,m)
  
  A=rep(list(0),m)
  Ablank=c(list(matrix(0,d,d)))
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    Ablank[[1+j]] = matrix(0, card_S, card_S)
  }
  
  
  ind = 0
  for(k in 1:card_patterns){
    card_S = length(patterns[[k]])
    for (i in 1:card_S){
      for (j in i:card_S){
        b[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i] = SigmaS[[k]][i,j]
        
        e_i = rep(0, card_S); e_i[i]=1
        e_j = rep(0, card_S); e_j[j]=1
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]] = Ablank
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]][[k+1]] = e_i%*%t(e_j)/2 + e_j%*%t(e_i)/2
        
        tmp_1 = e_i%*%t(e_j)/2 +  e_j%*%t(e_i)/2
        tmp_2 = matrix(0,d,d)
        tmp_2[patterns[[k]],patterns[[k]]] = tmp_1
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]][[1]] = tmp_2
      }
    }
    ind = ind + card_S*(card_S+1)/2
  }
  
  for (i in (d-2):0){
    A[[m-i]] = Ablank
    A[[m-i]][[1]][1,1]=1
    A[[m-i]][[1]][2+i,2+i]=-1
  }
  
  sizes = numeric(length = 1+card_patterns)
  sizes[1] = d; ind = 2
  for (S in patterns){
    sizes[ind] = length(S)
    ind = ind+1
  }
  
  K = list(type=rep("s",1+card_patterns),size=sizes) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  #----------------------------------------------------------------------------------------
  echo=FALSE
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj

  
  #----------------------------------------------------------------------------------------
  # Optimal XS 
  #----------------------------------------------------------------------------------------
  tmp0=c(list())
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    tmp0[[j]] = matrix(0, card_S, card_S)
  }
  
  YS=tmp0
  ind = 0
  for (k in 1:card_patterns){
    card_S = length(patterns[[k]])
    for (i in 1:card_S){
      for (j in i:card_S){
        tmp = as.numeric(abs(SDP$y[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]/2) > 0.0001)*
          SDP$y[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]
        YS[[k]][i,j] = tmp/2
        YS[[k]][j,i] = YS[[k]][j,i] + tmp/2
        
      }
    }
    ind = ind + card_S*(card_S+1)/2
  }
  
  XS0 = tmp0
  count = rep(0,d)
  for (k in 1:card_patterns){
    for(i in 1:d){
      if(i %in% patterns[[k]]){
        count[i] = count[i] + 1
      }
    }
  }
  
  for (k in 1:card_patterns){
    pattern = patterns[[k]]
    if (length(pattern) == 1){
      XS0[[k]] = as.matrix(1/count[pattern])
    } else{
      XS0[[k]] = diag(1/count[pattern])
    }
    
  }
  
  XS = c(list())
  for (k in 1:card_patterns){
    XS[[k]] = d*YS[[k]] - XS0[[k]]
  }
  
  #----------------------------------------------------------------------------------------
  # Optimal Sigma
  #----------------------------------------------------------------------------------------
  Sigma = SDP$X[[1]]
  
  SigmaSprime = list()
  for (i in 1:card_patterns){
    SigmaSprime[[i]] = SDP$X[[i+1]]
  }
  
  
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  my_list = list("R" = R, "XS" = XS, "XS0" = XS0, "YS" = YS,
                 "Sigma" = Sigma, "SigmaS" = SigmaS, 
                 "SigmaSprime" = SigmaSprime)
  return(my_list)
  
}

# for the regularised version alpha = 1/c corresponds to the d/c in the paper
computeR.reg = function(patterns=list(), SigmaS=list(), alpha) {
  
  #----------------------------------------------------------------------------------------
  ##### DEFINE PATTERN IF NOT SPECIFIED
  #----------------------------------------------------------------------------------------
  add_pattern = TRUE
  if (length(patterns) > 0){
    add_pattern = FALSE
  }
  
  i = 1
  while(add_pattern == TRUE){
    prompt1 = "Enter variable numbers (space-separated list): \n"
    tmp = as.vector(as.integer(strsplit(readline(prompt1), " ")[[1]]))
    patterns[[i]] = tmp
    
    prompt2 = "Do you want to add another pattern? (Write TRUE or FALSE): \n"
    add_pattern = as.logical((readline(prompt2)))
    i = i+1
  }
  
  #----------------------------------------------------------------------------------------
  # remove duplicates
  patterns = unique(patterns)
  card_patterns = length(patterns)
  
  d = 0 # computing d
  for (S in patterns){
    if (max(S) > d){
      d = max(S)
    }
  }
  
  
  
  #----------------------------------------------------------------------------------------
  ########## GENERATE SIGMA IF NOT SPECIFIED
  #----------------------------------------------------------------------------------------
  if (length(SigmaS) == 0){
    prompt2 = "Do you want a random SigmaS or a compatible one? (Write 1 or 2 respectively): \n"
    comp = as.numeric((readline(prompt2)))
    
    if (comp == 1){
      for(j in 1:card_patterns){
        card_S = length(patterns[[j]])
        A = matrix(runif((card_S)^2)*2-1, ncol=card_S) 
        tmp_0 = t(A) %*% A
        SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
      }
    } else if(comp == 2){
      SigmaS=list()
      A = matrix(runif(d^2)*2-1, ncol=d)
      Sigma = cov2cor(t(A) %*% A)
      for(j in 1:card_patterns){
        SigmaS[[j]]=as.matrix(Sigma[patterns[[j]],patterns[[j]]])
      }
    } else {
      return("An error occured.")
    }
  }
  
  #----------------------------------------------------------------------------------------
  ########## WRITING THE SDP PROBLEM ##########################
  #----------------------------------------------------------------------------------------
  
  # Matrix giving objective function (1/d)*tr(Sigma)
  C = c(list((1/d)*diag(d)))
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    C[[1+j]] = matrix(0, card_S, card_S)
  }
  C[[2+card_patterns]] = as.matrix(-alpha)
  
  # number of constraints
  m = d-1
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    m = m + card_S*(card_S+1)/2
  }
  
  b=rep(0,m)
  
  A=rep(list(0),m)
  Ablank=c(list(matrix(0,d,d)))
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    Ablank[[1+j]] = matrix(0, card_S, card_S)
  }
  Ablank[[2+card_patterns]] = as.matrix(0)
  
  
  ind = 0
  for(k in 1:card_patterns){
    card_S = length(patterns[[k]])
    for (i in 1:card_S){
      for (j in i:card_S){
        b[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i] = SigmaS[[k]][i,j]
        
        e_i = rep(0, card_S); e_i[i]=1
        e_j = rep(0, card_S); e_j[j]=1
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]] = Ablank
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]][[k+1]] = e_i%*%t(e_j)/2 + e_j%*%t(e_i)/2
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]][[2+card_patterns]] = -1*as.matrix(as.numeric(i == j))
        
        tmp_1 = e_i%*%t(e_j)/2 +  e_j%*%t(e_i)/2
        tmp_2 = matrix(0,d,d)
        tmp_2[patterns[[k]],patterns[[k]]] = tmp_1
        A[[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]][[1]] = tmp_2
      }
    }
    ind = ind + card_S*(card_S+1)/2
  }
  
  for (i in (d-2):0){
    A[[m-i]] = Ablank
    A[[m-i]][[1]][1,1]=1
    A[[m-i]][[1]][2+i,2+i]=-1
  }
  
  sizes = numeric(length = 2+card_patterns)
  sizes[1] = d; ind = 2
  for (S in patterns){
    sizes[ind] = length(S)
    ind = ind+1
  }
  sizes[2+card_patterns] = 1
  
  K = list(type=rep("s",2+card_patterns),size=sizes) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  #----------------------------------------------------------------------------------------
  echo=FALSE
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  
  
  #----------------------------------------------------------------------------------------
  # Optimal XS 
  #----------------------------------------------------------------------------------------
  tmp0=c(list())
  for(j in 1:card_patterns){
    card_S = length(patterns[[j]])
    tmp0[[j]] = matrix(0, card_S, card_S)
  }
  
  YS=tmp0
  ind = 0
  for (k in 1:card_patterns){
    card_S = length(patterns[[k]])
    for (i in 1:card_S){
      for (j in i:card_S){
        tmp = as.numeric(abs(SDP$y[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]/2) > 0.0001)*
          SDP$y[ind + (i-1)*(card_S+1)-(i-1)*i/2 + j+1-i]
        YS[[k]][i,j] = tmp/2
        YS[[k]][j,i] = YS[[k]][j,i] + tmp/2
        
      }
    }
    ind = ind + card_S*(card_S+1)/2
  }
  
  XS0 = tmp0
  count = rep(0,d)
  for (k in 1:card_patterns){
    for(i in 1:d){
      if(i %in% patterns[[k]]){
        count[i] = count[i] + 1
      }
    }
  }
  
  for (k in 1:card_patterns){
    pattern = patterns[[k]]
    if (length(pattern) == 1){
      XS0[[k]] = as.matrix(1/count[pattern])
    } else{
      XS0[[k]] = diag(1/count[pattern])
    }
    
  }
  
  XS = c(list())
  for (k in 1:card_patterns){
    XS[[k]] = d*YS[[k]] - XS0[[k]]
  }
  
  #----------------------------------------------------------------------------------------
  # Optimal Sigma
  #----------------------------------------------------------------------------------------
  Sigma = SDP$X[[1]]
  
  SigmaSprime = list()
  for (i in 1:card_patterns){
    SigmaSprime[[i]] = SDP$X[[i+1]]
  }
  
  z = SDP$X[[2+card_patterns]]
  
  
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  my_list = list("R" = R, "XS" = XS, "XS0" = XS0, "z" = z, "YS" = YS,
                 "Sigma" = Sigma, "SigmaS" = SigmaS, 
                 "SigmaSprime" = SigmaSprime)
  return(my_list)
  
}
