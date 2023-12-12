source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 500
M = 50
d = 50

computeR_kkk = function(beta){
  #### POPULATION LEVEL ######
  SigmaS=list()
  
  SigmaS[[1]]=diag(2*d)
  SigmaS[[1]][1:d, (d+1):(2*d)] = P
  SigmaS[[1]][(d+1):(2*d), 1:d] = t(P)
  
  SigmaS[[2]]=diag(2*d)
  SigmaS[[2]][1:d, (d+1):(2*d)] = -P
  SigmaS[[2]][(d+1):(2*d), 1:d] = -t(P)
  
  SigmaS[[3]]=diag(2*d)
  SigmaS[[3]][1:d, (d+1):(2*d)] = beta*diag(d)
  SigmaS[[3]][(d+1):(2*d), 1:d] = beta*diag(d)
  
  return(computeR(list(c(1:(2*d)),c(1:d,(2*d+1):(3*d)),c((d+1):(3*d))),
                  SigmaS = SigmaS)$R)
}

little_power = function(beta){
  #### POPULATION LEVEL ######
  SigmaS=list()
  
  SigmaS[[1]]=diag(2*d)
  SigmaS[[1]][1:d, (d+1):(2*d)] = P
  SigmaS[[1]][(d+1):(2*d), 1:d] = t(P)
  
  SigmaS[[2]]=diag(2*d)
  SigmaS[[2]][1:d, (d+1):(2*d)] = -P
  SigmaS[[2]][(d+1):(2*d), 1:d] = -t(P)
  
  SigmaS[[3]]=diag(2*d)
  SigmaS[[3]][1:d, (d+1):(2*d)] = beta*diag(d)
  SigmaS[[3]][(d+1):(2*d), 1:d] = beta*diag(d)
  
  
  little_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X1 = mvrnorm(n, rep(0, 2*d), SigmaS[[1]])
    X2 = mvrnorm(n, rep(0, 2*d), SigmaS[[2]])
    X3 = mvrnorm(n, rep(0, 2*d), SigmaS[[3]])
    
    X = data.frame(matrix(nrow = 3*n, ncol = 3*d))
    X[1:n, c(1:(2*d))] = X1
    X[(n+1):(2*n), c(1:d,(2*d+1):(3*d))] = X2
    X[(2*n+1):(3*n), c((d+1):(3*d))] = X3
    X = as.matrix(X)
    
    ### run little's test
    little_decisions[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decisions))
}

little_power_cov = function(beta){
  #### POPULATION LEVEL ######
  SigmaS=list()
  
  SigmaS[[1]]=diag(2*d)
  SigmaS[[1]][1:d, (d+1):(2*d)] = P
  SigmaS[[1]][(d+1):(2*d), 1:d] = t(P)
  
  SigmaS[[2]]=diag(2*d)
  SigmaS[[2]][1:d, (d+1):(2*d)] = -P
  SigmaS[[2]][(d+1):(2*d), 1:d] = -t(P)
  
  SigmaS[[3]]=diag(2*d)
  SigmaS[[3]][1:d, (d+1):(2*d)] = beta*diag(d)
  SigmaS[[3]][(d+1):(2*d), 1:d] = beta*diag(d)
  
  
  little_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X1 = mvrnorm(n, rep(0, 2*d), SigmaS[[1]])
    X2 = mvrnorm(n, rep(0, 2*d), SigmaS[[2]])
    X3 = mvrnorm(n, rep(0, 2*d), SigmaS[[3]])
    
    X = data.frame(matrix(nrow = 3*n, ncol = 3*d))
    X[1:n, c(1:(2*d))] = X1
    X[(n+1):(2*n), c(1:d,(2*d+1):(3*d))] = X2
    X[(2*n+1):(3*n), c((d+1):(3*d))] = X3
    X = as.matrix(X)
    
    ### run little's test
    little_decisions[i] = little_test(X, alpha = 0.05, type = "cov")
  }
  return(mean(little_decisions))
}

bootstrap_power = function(beta){
  #### POPULATION LEVEL ######
  SigmaS=list()
  
  SigmaS[[1]]=diag(2*d)
  SigmaS[[1]][1:d, (d+1):(2*d)] = P
  SigmaS[[1]][(d+1):(2*d), 1:d] = t(P)
  
  SigmaS[[2]]=diag(2*d)
  SigmaS[[2]][1:d, (d+1):(2*d)] = -P
  SigmaS[[2]][(d+1):(2*d), 1:d] = -t(P)
  
  SigmaS[[3]]=diag(2*d)
  SigmaS[[3]][1:d, (d+1):(2*d)] = beta*diag(d)
  SigmaS[[3]][(d+1):(2*d), 1:d] = beta*diag(d)
  
  our_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X1 = mvrnorm(n, rep(0, 2*d), SigmaS[[1]])
    X2 = mvrnorm(n, rep(0, 2*d), SigmaS[[2]])
    X3 = mvrnorm(n, rep(0, 2*d), SigmaS[[3]])
    
    X = data.frame(matrix(nrow = 3*n, ncol = 3*d))
    X[1:n, c(1:(2*d))] = X1
    X[(n+1):(2*n), c(1:d,(2*d+1):(3*d))] = X2
    X[(2*n+1):(3*n), c((d+1):(3*d))] = X3
    X = as.matrix(X)
    
    ### run our tests
    our_decisions[i] = MCAR_corr_test(X, alpha = 0.05, B = 99, type = "p")
  }
  
  return(mean(our_decisions))
}


######## USING FOREACH AND DORNG

start.time = Sys.time()

# definition of P
P = 2*matrix(runif((d)^2)*2-1, ncol=d)/d

while(norm(P, type="2") > 0.75){
  P = matrix(runif((d)^2)*2-1, ncol=d)/d
}

normP = norm(P, type="2")
xxx = seq(1-2*normP^2, 0.9, length.out = 8)


### work in parallel
registerDoFuture()
plan(multicore)
RNGkind("L'Ecuyer-CMRG")
set.seed(232)

start.time = Sys.time()

R = foreach(kkk = xxx, .combine = 'c') %dorng% computeR_kkk(kkk)
little_power = foreach(kkk = xxx, .combine = 'c') %dorng% little_power(kkk)
little_power_cov = foreach(kkk = xxx, .combine = 'c') %dorng% little_power_cov(kkk)
our_power = foreach(kkk = xxx, .combine = 'c') %dorng% bootstrap_power(kkk)

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

### plot the simulation

png(paste("pictures/block_3_cycle-",d, "P.png"))
plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type = "b")
lines(R, little_power_cov, col="orange", pch=18, type = "b")
lines(R, our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("bottomright",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))
dev.off()