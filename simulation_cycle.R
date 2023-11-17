setwd("~/Documents/phd/MCAR_covariance/MCAR/")
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
M = 1000 # with 4000 I expect 4 hours
t3 = pi/6
t2 = pi/4

little_power = c()
little_power_cov = c()
our_power = c()
R = c()

start.time = Sys.time()

for(t1 in seq(0.1, pi - 0.1, length.out = 20)){
  
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
  our_decisions = c()
  for (i in 1:M){
    
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
    
    ### run little's test
    little_decisions = c(little_decisions, little_test(X, alpha))
    little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))
    
    ### run our tests
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  our_power = c(our_power, mean(our_decisions))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

plot(R, little_power, col="green", ylim = c(0,1), pch=18, xlab = "", ylab = "")
points(R, little_power_cov, col="orange", pch=18)
points(R, our_power, col="blue", pch=18)
abline(h = alpha, col="red")
legend("center",
       legend = c("Little's power", "Little's power cov", "Our power"),
       col = c("green", "orange", "blue"),
       pch = c(18, 18, 18))

######### d-cycle: high-dimensional ############
d = 20

alpha = 0.05
n = 200
M = 50# with 4000 I expect 4 hours
angle = pi/(2*(d-1))

our_power = c()
R = c()

start.time = Sys.time()

for(t1 in seq(pi/2, 3*pi/4, length.out = 10)){
  
  #### POPULATION LEVEL ######
  SigmaS=list() #Random 2x2 correlation matrices (necessarily consistent)
  for(j in 1:d){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
    SigmaS[[j]][1,2] = cos(angle)
    SigmaS[[j]][2,1] = cos(angle)
  }
  
  SigmaS[[d]][1,2] = cos(t1)
  SigmaS[[d]][2,1] = cos(t1)
  
  patterns = list()
  v = 1:d
  for (i in 1:(d-1)){
    patterns[[i]] = c(v[i], v[i+1])
  }
  patterns[[d]] = c(1,d)
  R = c(R, computeR(patterns, SigmaS)$R)
  
  
  ###### SAMPLE LEVEL, REPEATING THE TEST M TIMES #######
  our_decisions = c()
  for (i in 1:M){
    
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), patterns[[i]]] = mvrnorm(n, rep(0, 2), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    ### run our tests
    our_decisions = c(our_decisions, MCAR_corr_test(X, alpha, B = 99))
  }
  
  our_power = c(our_power, mean(our_decisions))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

plot(R, our_power, col="blue", ylim = c(0,1), pch=18, xlab = "", ylab = "")
abline(h = alpha, col="red")
legend("topleft", legend = "Our power", col = "blue", pch = 18)


