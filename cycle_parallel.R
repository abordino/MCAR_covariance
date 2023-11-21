library(missMethods)
library(MASS)
library(norm)
library(parallel)
library(copula)
library(missMethods)
library(misty)
library(Rcsdp)
library(Matrix)
library(npmr)
library(matrixcalc)
library(pracma)
library(foreach)
library(doSNOW)
library(doParallel)

computeR = function(patterns=list(), SigmaS=list()) {
  
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
  # ####### THE FINAL OUTPUT
  # #----------------------------------------------------------------------------------------
  # print("------------------------------------------------")
  # print("COMPUTE R FOR A SEQUENCE OF CORRELATION MATRICES")
  # print("------------------------------------------------")
  # print(paste("The dimensionality of the problem is",d))
  # print(paste("There are", card_patterns, "patterns:"))
  # for (S in patterns){
  #   print(S)
  # }
  # print("------------------------------------------------")
  # print("The sequence of correlation matrices is:")
  # print(SigmaS)
  # print("------------------------------------------------")
  # print(paste("R is equal to ", R))
  # is_compatibile = (R < 0.0001)
  # if (is_compatibile == T){
  #   print("Hence, SigmaS is compatible")
  # } else{
  #   print("Hence, SigmaS is not compatible")
  # }
  
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

get_SigmaS = function(X){
  
  #### create vector with indicators of NA's
  d = dim(X)[2]
  v = 1:d
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
    if (dim(data_pattern[[i]])[1] < 10){ # put 3 in order for Little's test to work
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
                 "data_pattern" = data_pattern, "SigmaS" = SigmaS, 
                 "ambient_dimension" = d)
  return(my_list)
}

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

MCAR_corr_test = function(X, alpha, B){
  
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

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 200
M = 50
t3 = pi/6
t2 = pi/4
n_cores = detectCores()-1

little_power = c()
little_power_cov = c()
R = c()

bootstrap_test_t1 = function(t1){
    SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
    for(j in 1:3){
      x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
    }
    
    SigmaS[[1]][1,2] = cos(t1)
    SigmaS[[1]][2,1] = cos(t1)
    SigmaS[[2]][1,2] = cos(t2)
    SigmaS[[2]][2,1] = cos(t2)
    SigmaS[[3]][1,2] = cos(t3)
    SigmaS[[3]][2,1] = cos(t3)
    
    ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
    our_decisions = c()
    for (i in 1:M){
      
      #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
      X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
      X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
      X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
      
      columns = c("X1","X2","X3") 
      X = data.frame(matrix(nrow = 3*n, ncol = 3))
      X[1:n, c("X1", "X2")] = X1
      X[(n+1):(2*n), c("X2", "X3")] = X2
      X[(2*n+1):(3*n), c("X1", "X3")] = X3
      X = as.matrix(X)
      
      ### run our tests
      our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
    }
    
    res = mean(our_decisions)
    return(res)
}

our_power = mclapply(seq(t2-t3, pi - t3, length.out = n_cores), bootstrap_test_t1, mc.cores = n_cores)

##### easy part
for(t1 in seq(t2-t3, pi - t3, length.out = n_cores)){
    
    #### POPULATION LEVEL ######
    SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
    for(j in 1:3){
      x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
    }
    
    SigmaS[[1]][1,2] = cos(t1)
    SigmaS[[1]][2,1] = cos(t1)
    SigmaS[[2]][1,2] = cos(t2)
    SigmaS[[2]][2,1] = cos(t2)
    SigmaS[[3]][1,2] = cos(t3)
    SigmaS[[3]][2,1] = cos(t3)
    
    R = c(R, computeR(list(c(1,2),c(2,3), c(1,3)), SigmaS)$R)
    
    
    ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
    little_decisions = c()
    little_decisions_cov = c()
    for (i in 1:M){
      
      #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
      X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
      X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
      X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
      
      columns = c("X1","X2","X3") 
      X = data.frame(matrix(nrow = 3*n, ncol = 3))
      X[1:n, c("X1", "X2")] = X1
      X[(n+1):(2*n), c("X2", "X3")] = X2
      X[(2*n+1):(3*n), c("X1", "X3")] = X3
      X = as.matrix(X)
      
      ### run little's test
      little_decisions = c(little_decisions, little_test(X, alpha))
      little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))
      
    }
    
    little_power = c(little_power, mean(little_decisions))
    little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  }

plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(R, little_power_cov, col="orange", pch=18)
points(R, our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))

