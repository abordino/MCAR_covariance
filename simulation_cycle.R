setwd("~/Documents/phd/MCAR_covariance/MCAR/")
source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)
library(parallel)
library(doParallel)

MC_rep = function(i){
  
  setwd("~/Documents/phd/MCAR_covariance/MCAR/")
  source("computeR.R")
  source("little_test.R")
  source("bootstrap_test.R")
  source("find_SigmaS.R")
  library(missMethods)
  library(MASS)
  library(norm)
  library(parallel)
  library(doParallel)
  
  #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
  X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
  X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
  X3 = mvrnorm(n, c(0,0), SigmaS[[3]])
  
  columns = c("X1","X2","X3") 
  X = data.frame(matrix(nrow = 3*n, ncol = length(columns)))
  X[1:n, c("X1", "X2")] = X1
  X[(n+1):(2*n), c("X2", "X3")] = X2
  X[(2*n+1):(3*n), c("X1", "X3")] = X3
  X = as.matrix(X)
  
  result_vect = c(little_test(X, alpha),
                     little_test(X, alpha, "cov"), 
                     MCAR_corr_test(X, alpha, B = 1))
  return(result_vect)
}

#### set up the cluster
n_cores = detectCores()
cl = makeCluster(n_cores-1, outfile = "Log.txt")
registerDoParallel(cl)

# load all the libraries
clusterCall(cl, function() library(MASS))
clusterCall(cl, function() library(norm))
clusterCall(cl, function() library(missMethods))
clusterCall(cl, function() library(misty))
clusterCall(cl, function() library(Rcsdp))
clusterCall(cl, function() library(Matrix))
clusterCall(cl, function() library(npmr))
clusterCall(cl, function() library(matrixcalc))
clusterCall(cl, function() library(pracma))

# load all the variables
clusterExport(cl, varlist = ls())

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 200
M = 5000
t3 = pi/6
t2 = pi/4

little_power = c()
little_power_cov = c()
our_power = c()
R = c()

for(t1 in seq(t2-t3, pi - t3, length.out = 1)){

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

  
  decisions = foreach(i = 1:M, .combine = "rbind",.inorder = F) %dopar% {
    MC_rep(i)
  }
  
  
  # little_decisions = c()
  # little_decisions_cov = c()
  # our_decisions = c()
  # 
  # for (i in 1:M){
  #   little_decisions = c(little_decisions, decisions[[i]]$little_dec)
  #   little_decisions_cov = c(little_decisions_cov, decisions[[i]]$little_dec_cov)
  #   our_decisions = c(our_decisions, decisions[[i]]$our_dec)
  # }
  # 
  # little_power = c(little_power, mean(little_decisions))
  # little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  # our_power = c(our_power, mean(our_decisions))
}

plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(R, little_power_cov, col="orange", pch=18)
points(R, our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))

stopCluster(cl)


