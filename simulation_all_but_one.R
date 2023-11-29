source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
library(missMethods)
library(MASS)
library(norm)

######### 3-cycle: setting 1 ############
alpha = 0.05
n = 200
M = 100
d = 10

##### easy part
computeR_kkk = function(t1){
  #### POPULATION LEVEL ######
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  v = 1:d
  pattern = list()
  for (i in 1:d){
    pattern[[i]] = v[-i]
  }
  
  return(computeR(patterns = pattern, SigmaS = SigmaS)$R)
}

little_power = function(t1){
  #### POPULATION LEVEL ######
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  v = 1:d
  pattern = list()
  for (i in 1:d){
    pattern[[i]] = v[-i]
  }
  
  little_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), pattern[[i]]] = mvrnorm(n, rep(0, d-1), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    ### run little's test
    little_decisions[i] = little_test(X, alpha = 0.05)
  }
  
  return(mean(little_decisions))
}

little_power_cov = function(t1){
  #### POPULATION LEVEL ######
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  v = 1:d
  pattern = list()
  for (i in 1:d){
    pattern[[i]] = v[-i]
  }
  
  little_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), pattern[[i]]] = mvrnorm(n, rep(0, d-1), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    ### run little's test
    little_decisions[i] = little_test(X, alpha = 0.05, type = "cov")
  }
  return(mean(little_decisions))
}

#### parallel using parallel
bootstrap_power = function(t1){
  #### POPULATION LEVEL ######
  SigmaS=list() 
  for(j in 1:d){
    A = matrix(runif((d-1)^2)*2-1, ncol=d-1) 
    tmp_0 = t(A) %*% A
    SigmaS[[j]]=as.matrix(cov2cor(tmp_0))
  }
  
  v = 1:d
  pattern = list()
  for (i in 1:d){
    pattern[[i]] = v[-i]
  }
  
  little_decisions = logical(length = M)
  for (i in 1:M){
    
    #### generate dataset
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), pattern[[i]]] = mvrnorm(n, rep(0, d-1), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    ### run our tests
    our_decisions[i] = MCAR_corr_test(X, alpha = 0.05, B = 99)
  }
  
  return(mean(our_decisions))
}


xxx = 1:6

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

# plot the simulation
png(paste("all_but_one-", d, ".png"))
plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "", type = "b")
lines(R, little_power_cov, col="orange", pch=18, type = "b")
lines(R, our_power, col="blue", pch=18, type = "b")
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))
dev.off()