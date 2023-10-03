library(Rcsdp)
library(Matrix)
library(npmr)
library(matrixcalc)
library(pracma)


############ FUNCTION ###############
computeR = function(patterns=list(), SigmaS=list()) {
  library(Rcsdp)
  library(Matrix)
  
  #----------------------------------------------------------------------------------------
  ##### DEFINE PATTERN IF NOT SPECIFIED
  #----------------------------------------------------------------------------------------
  add_pattern = T
  if (length(patterns) > 0){
    add_pattern = F
  }
  
  i = 1
  while(add_pattern == T){
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
    } else if(comp ==2){
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
  
  sizes = c(d)
  for (S in patterns){
    card_S = length(S)
    sizes = c(sizes, card_S)
  }
  K = list(type=rep("s",card_patterns+1),size=sizes) # Variables are the PSD matrices Sigma,Z_{1,2},...,Z_{d,1}
  
  #----------------------------------------------------------------------------------------
  echo=F
  SDP=csdp(C, A, b, K) #Running the function
  R = 1-SDP$pobj
  
  #----------------------------------------------------------------------------------------
  ####### THE FINAL OUTPUT
  #----------------------------------------------------------------------------------------
  print("------------------------------------------------")
  print("COMPUTE R FOR A SEQUENCE OF CORRELATION MATRICES")
  print("------------------------------------------------")
  print(paste("The dimensionality of the problem is",d))
  print(paste("There are", card_patterns, "patterns:"))
  for (S in patterns){
    print(S)
  }
  print("------------------------------------------------")
  print("The sequence of correlation matrices is:")
  print(SigmaS)
  print("------------------------------------------------")
  print(paste("R is equal to ", R))
  is_compatibile = (R < 0.0001)
  if (is_compatibile == T){
    print("Hence, SigmaS is compatible")
  } else{
    print("Hence, SigmaS is not compatible")
  }
  
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
    card_S = length(patterns[[k]])
    pattern = patterns[[k]]
    XS0[[k]] = diag(1/count[pattern])
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
    SigmaSprime[[1]] = SDP$X[[i+1]]
  }
  
  
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  my_list = list("R" = R, "XS" = XS, "XS0" = XS0, "Sigma" = Sigma, "SigmaS" = SigmaS, 
                 "SigmaSprime" = SigmaSprime)
  return(my_list)
  
}


##### TESTING ON A SPECIFIC EXAMPLE (ALL d-1)
d = 5
SigmaS=list() 
for(j in 1:d){
  A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
  tmp_0 = t(A) %*% A
  SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
}
result = computeR(patterns = list(c(2,3,4,5),c(1,3,4,5),c(1,2,4,5),c(1,2,3,5),c(1,2,3,4)), SigmaS = SigmaS)

AAA = result$SigmaS_prime
computeR(patterns = list(c(2,3,4,5),c(1,3,4,5),c(1,2,4,5),c(1,2,3,5),c(1,2,3,4)), SigmaS = AAA)$R

############# test on 3-cycle ####################
R_rho = c()
d = 3
SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(0)
SigmaS[[1]][2,1] = cos(0)
SigmaS[[2]][1,2] = cos(0)
SigmaS[[2]][2,1] = cos(0)

for (rho in seq(0, pi, length.out = 1000)){
  SigmaS[[3]][1,2] = cos(rho)
  SigmaS[[3]][2,1] = cos(rho)
  
  result = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS = SigmaS)
  R_rho = c(R_rho, result$R)
}

plot(seq(0, pi, length.out = 1000), R_rho, col="red", pch=19, cex=0.1)
curve(sin(x/2)^2, add=T)


R_rho = c()
theta_2 = pi/4
d = 3
SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
for(j in 1:d){
  x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
}

SigmaS[[1]][1,2] = cos(0)
SigmaS[[1]][2,1] = cos(0)
SigmaS[[2]][1,2] = cos(theta_2)
SigmaS[[2]][2,1] = cos(theta_2)

for (rho in seq(0, pi, length.out = 1000)){
  SigmaS[[3]][1,2] = cos(rho)
  SigmaS[[3]][2,1] = cos(rho)
  
  result = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS = SigmaS)
  R_rho = c(R_rho, result$R)
}

plot(seq(0, pi, length.out = 1000), R_rho, col="red", pch=19, cex=0.1)
curve(abs((cos(x)-cos(theta_2)))/2, add=T)

################## 1-Schatten + maximum modulus LOWER and UPPER BOUND. ########

m_k_plus_s_k = c()
R_k = c()
d=3

for (kkk in 1:3000){
  SigmaS=list() 

  a = matrix(runif((d)^2)*2-1, ncol=d)
  tmp_0 = t(a) %*% a
  SigmaS[[2]]=round(as.matrix(cov2cor(tmp_0)),10)
  
  a = matrix(runif((d-1)^2)*2-1, ncol=d-1)
  tmp_0 = t(a) %*% a
  SigmaS[[1]]=round(as.matrix(cov2cor(tmp_0)),10)
  
  if(is.positive.definite(SigmaS[[1]])==F || is.positive.definite(SigmaS[[2]])==F){
    print("NOT PSD")
    next
  }
  result = computeR(list(c(1,2),c(1,2,3)), SigmaS = SigmaS)
  
  R_k = c(R_k, result$R)
  
  tmp = abs(SigmaS[[1]][1,2]-SigmaS[[2]][1,2])/2
  m_k_plus_s_k = c(m_k_plus_s_k, tmp)
  
}

plot(m_k_plus_s_k, R_k, col="red", pch=19, cex=0.1, main="d = 3, |x-y|/2")
curve(x^1, add=T)

##### EVEN FOR d=3, WHERE R=0 IFF |x-y|=0, THE LOWER BOUND SEEMS TO BE LOOSE. 
####  WHAT ARE WE MISSING HERE? ANY Patterns? For varying y?
  R_y = c()
  for (y in seq(-1, 1, length.out=1000)){
    SigmaS[[1]][1,2] = y
    SigmaS[[1]][2,1] = y
    result = computeR(list(c(1,2),c(1,2,3)), SigmaS = SigmaS)
    R_y = c(R_y, result$R)
    
  }
  plot(seq(-1, 1, length.out=1000), R_y, col="red", pch=19, cex=0.1, 
       main=paste("d = 3, for varying y, and x = ", SigmaS[[2]][1,2]))
  curve(abs(x-SigmaS[[2]][1,2])/2, add=T)
  curve(1-(1-x)/(1-SigmaS[[2]][1,2]), add=T)
  curve(1-(1+x)/(1+SigmaS[[2]][1,2]), add=T)
  
 #### \mathbb{S} = {[d-2]U{d}, [d-1]} ###############
  R_n = c()
  theta_n = c()
  ub_n = c()
  d = 7
  
  for (i in 1:1000){

  SigmaS=list() 
  for(j in 1:2){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1)
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=round(as.matrix(cov2cor(tmp_0)),10)
  }
  
  result = computeR(list(c(1,2,3,4,5,7),c(1,2,3,4,5,6)), SigmaS = SigmaS)
  R_n = c(R_n, result$R)
  
  O = abs(SigmaS[[1]][1:(d-2),1:(d-2)] - SigmaS[[2]][1:(d-2),1:(d-2)])
  i = which(O==max(O), arr.ind=T)[1,1]
  j = which(O==max(O), arr.ind=T)[1,2]
  
  tmp = abs(SigmaS[[1]][i,j] - SigmaS[[2]][i,j])/2
  theta_n = c(theta_n, tmp)
  
  tmp = 1 - min((1-min(SigmaS[[1]][i,j],SigmaS[[2]][i,j]))/(1-max(SigmaS[[1]][i,j],SigmaS[[2]][i,j])), 
                (1+min(SigmaS[[1]][i,j],SigmaS[[2]][i,j]))/(1+max(SigmaS[[1]][i,j],SigmaS[[2]][i,j])))
  ub_n = c(ub_n, tmp)
  
}
  
  plot(1:1000, sort(R_n), col="red", pch=19, cex=0.1, ylim = c(0,1), type="l",
       main=paste("d = 7, for varying y, and x = ", SigmaS[[1]][1,2]))
  lines(theta_n[sort(R_n, index.return=TRUE)$ix], col="blue")
  lines(ub_n[sort(R_n, index.return=TRUE)$ix], col="green")

  
  ############## block 3-cycle ############
  
  ######## SV lower bound
  for (m in 2:7){
    R_P = c()
    R_SV = c()
    normP = c()
    R_cycle = c()
    
    rho = 0.8
    
    for (kkk in 1:1000){
      
      P = 0.5*matrix(runif((m)^2)*2-1, ncol=m)
      
      while(norm(P, type="2") > 1){
        P = 0.5*matrix(runif((m)^2)*2-1, ncol=m)
      }
      
      P = runif(1,0,1)*P
      sqSingValues=eigen(t(P)%*%P)$values
      guess=(1/m)*sum(pmax(sqSingValues-(1-rho)/2,0))
      normP = c(normP, guess)
      
      SigmaS=list()
      
      SigmaS[[1]]=diag(2*m)
      SigmaS[[1]][1:m, (m+1):(2*m)] = P
      SigmaS[[1]][(m+1):(2*m), 1:m] = t(P)
      
      SigmaS[[2]]=diag(2*m)
      SigmaS[[2]][1:m, (m+1):(2*m)] = -P
      SigmaS[[2]][(m+1):(2*m), 1:m] = -t(P)
      
      SigmaS[[3]]=diag(2*m)
      SigmaS[[3]][1:m, (m+1):(2*m)] = rho*diag(m)
      SigmaS[[3]][(m+1):(2*m), 1:m] = rho*diag(m)
      
      result = computeR(list(c(1:(2*m)),c(1:m,(2*m+1):(3*m)),c((m+1):(3*m))),
                        SigmaS = SigmaS)
      R_P = c(R_P, result$R)
      
      # Q = diag(sqrt(sqSingValues))
      # SigmaS=list()
      # 
      # SigmaS[[1]]=diag(2*m)
      # SigmaS[[1]][1:m, (m+1):(2*m)] = Q
      # SigmaS[[1]][(m+1):(2*m), 1:m] = t(Q)
      # 
      # SigmaS[[2]]=diag(2*m)
      # SigmaS[[2]][1:m, (m+1):(2*m)] = -Q
      # SigmaS[[2]][(m+1):(2*m), 1:m] = -t(Q)
      # 
      # SigmaS[[3]]=diag(2*m)
      # SigmaS[[3]][1:m, (m+1):(2*m)] = rho*diag(m)
      # SigmaS[[3]][(m+1):(2*m), 1:m] = rho*diag(m)
      # 
      # result = computeR(list(c(1:(2*m)),c(1:m,(2*m+1):(3*m)),c((m+1):(3*m))),
      #                   SigmaS = SigmaS)
      # 
      # R_SV = c(R_SV, result$R)
      
      tmp = 0
      for (j in 1:m){
        cycle = list()
        cycle[[1]] = diag(2)
        cycle[[1]][1,2] = sqrt(sqSingValues[j]); cycle[[1]][2,1] = sqrt(sqSingValues[j])
  
        cycle[[2]] = diag(2)
        cycle[[2]][1,2] = rho; cycle[[2]][2,1] = rho
  
        cycle[[3]] = diag(2)
        cycle[[3]][1,2] = -sqrt(sqSingValues[j]); cycle[[3]][2,1] = -sqrt(sqSingValues[j])
        
        resultC = computeR(list(c(1,2), c(2,3), c(1,3)), SigmaS = cycle)
        tmp = tmp + resultC$R/m
        print("AAAAAAAA")
        print(resultC$XS)
      }
      R_cycle = c(R_cycle, tmp)

    # if (R_P[kkk]-R_P[kkk]>0.001){
    #   break
    # }
      
    # if ((R_P[kkk]-R_cycle[kkk]<0.001)&R_cycle[kkk]>0.001){
    #     break
    #   }
    }
    
    
    plot(normP, R_P, col="red", pch=19, cex=0.2,
         main=paste("m = ", m, "; rho = ", rho), ylim=c(0,1))
    points(normP, R_SV, col="blue", pch=19, cex=0.2)
    points(normP, R_cycle, col="cyan", pch=19, cex=0.2)
    curve(x^1, add=T, col="green")
    
    # if (R_P[kkk]-R_P[kkk]>0.001){
    #   break
    # }
    

  } 
  
  result$XS
  resultC$XS
  R_P[kkk]
  R_cycle[kkk]
  result$SigmaS
  resultC$SigmaS
  
  
  #### 3-cycle lower bound 
  R_l = c()
  rho = 0.55
  for (l in seq(0,1,length.out=1000)){
    SigmaS = list()
    
    SigmaS[[1]] = diag(2)
    SigmaS[[1]][1,2] = l; SigmaS[[1]][2,1] = l
    
    SigmaS[[2]] = diag(2)
    SigmaS[[2]][1,2] = -l; SigmaS[[2]][2,1] = -l
    
    SigmaS[[3]] = diag(2)
    SigmaS[[3]][1,2] = rho; SigmaS[[3]][2,1] = rho
    
    result = computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS = SigmaS)
    
    R_l = c(R_l, result$R)
    
  }
  
  plot(seq(0,1,length.out=1000), R_l, col="red", pch=19, cex=0.1,
       main=paste("rho = ", rho), ylim=c(0,1))
  curve(x^2-(1-rho)/2, add=T, col="green")
  d = function(x){
    return(acos(((1-rho)-sqrt(8*x*rho + rho^2-8*x-10*rho+9))/(4*(1-x))))
  }
  curve(1-(1-x)/(1+cos(d(x))), col="blue", add = T)
