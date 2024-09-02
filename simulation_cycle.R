source("computeR.R")
source("little_test.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
source("indexConsistency.R")
library(missMethods)
library(norm)
library(naniar)
library(MASS)
library(naniar)
library(latex2exp)
library(compositions)

#----------------------------------------------------------------------------------------
######### 3-cycle: lognormal ############
#----------------------------------------------------------------------------------------
alpha = 0.05
n = 1000
MC = 800
t3 = pi/6
t2 = pi/3
d = 3

little_power_mean = c()
little_power = c()
little_power_cov = c()
our_power = c()
our_power_corr = c()
our_power_mean = c()

start.time = Sys.time()

R = c()
for(t1 in seq(t2+t3, pi-0.01, length.out = 8)){

  #----------------------------------------------------------------------------------------
  #### ESTIMATE R FROM DATA, WITH INDEPENDENT SAMPLE ######
  #----------------------------------------------------------------------------------------
  SigmaS=list()
  for(j in 1:3){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }

  SigmaS[[1]][1,2] = cos(t1)
  SigmaS[[1]][2,1] = cos(t1)
  SigmaS[[2]][1,2] = cos(t2)
  SigmaS[[2]][2,1] = cos(t2)
  SigmaS[[3]][1,2] = cos(t3)
  SigmaS[[3]][2,1] = cos(t3)

  tmp = 0
  for (i in 1:MC){

    #----------------------------------------------------------------------------------------
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    #----------------------------------------------------------------------------------------
    X1 = rlnorm.rplus(n, c(0,0), SigmaS[[1]])
    X2 = rlnorm.rplus(n, c(0,0), SigmaS[[2]])
    X3 = rlnorm.rplus(n, c(0,0), SigmaS[[3]])

    S = list()
    S[[1]] = cor(X1)
    S[[2]] = cor(X2)
    S[[3]] = cor(X3)

    tmp = tmp + computeR(patterns = list(c(1,2), c(2,3), c(1,3)), SigmaS = S)$R/MC
  }
  R = c(R, tmp)
}

for(t1 in seq(t2+t3, pi-0.01, length.out = 8)){

  SigmaS=list()
  for(j in 1:3){
    x=runif(2,min=-1,max=1); y=runif(2,min=-1,max=1); SigmaS[[j]]=cov2cor(x%*%t(x) + y%*%t(y))
  }

  SigmaS[[1]][1,2] = cos(t1)
  SigmaS[[1]][2,1] = cos(t1)
  SigmaS[[2]][1,2] = cos(t2)
  SigmaS[[2]][2,1] = cos(t2)
  SigmaS[[3]][1,2] = cos(t3)
  SigmaS[[3]][2,1] = cos(t3)

  #----------------------------------------------------------------------------------------
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  #----------------------------------------------------------------------------------------
  little_decisions_mean = c()
  little_decisions = c()
  little_decisions_cov = c()
  our_decisions = c()
  our_decisions_corr = c()
  our_decisions_mean = c()

  for (i in 1:MC){

    #----------------------------------------------------------------------------------------
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    #----------------------------------------------------------------------------------------
    X1 = rlnorm.rplus(n, c(0,0), SigmaS[[1]])
    X2 = rlnorm.rplus(n, c(0,0), SigmaS[[2]])
    X3 = rlnorm.rplus(n, c(0,0), SigmaS[[3]])

    columns = c("X1","X2","X3")
    X = data.frame(matrix(nrow = 3*n, ncol = 3))
    X[1:n, c("X1", "X2")] = X1
    X[(n+1):(2*n), c("X2", "X3")] = X2
    X[(2*n+1):(3*n), c("X1", "X3")] = X3
    X = as.matrix(X)

    #----------------------------------------------------------------------------------------
    ### run little's test
    #----------------------------------------------------------------------------------------
    little_decisions_mean = c(little_decisions_mean, mcar_test(data.frame(X))$p.value < alpha)
    little_decisions = c(little_decisions, little_test(X, alpha))
    little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))

    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    our_decisions = c(our_decisions, MCAR_meancovTest(X, alpha, B = 99))
    our_decisions_corr = c(our_decisions_corr, MCAR_corrTest(X, alpha, B = 99))
    our_decisions_mean = c(our_decisions_mean, MCAR_meanTest(X, alpha, B = 99))
  }

  little_power_mean = c(little_power_mean, mean(little_decisions_mean))
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  our_power = c(our_power, mean(our_decisions))
  our_power_corr = c(our_power_corr, mean(our_decisions_corr))
  our_power_mean = c(our_power_mean, mean(our_decisions_mean))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken


png("3_cycle_lognormal_1.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(sort(R), little_power_mean, col="green", ylim = c(0,1), pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(sort(R), little_power, col="orange", pch=19, type = "b")
lines(sort(R), little_power_cov, col="brown", pch=20, type = "b")
lines(sort(R), our_power, col="blue", pch=21, type = "b")
lines(R, our_power_corr, col="darkviolet", pch=22, type = "b")
lines(R, our_power_mean, col="aquamarine2", pch=23, type = "b")
lines(sort(R), rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{cov}$)'),
                  TeX(r'($d^2_{aug}$)'), "Omnibus"),
       col = c("green", "brown", "orange", "blue"),
       pch = c(18, 19, 20, 21))
dev.off()
# 
#----------------------------------------------------------------------------------------
# ######### 3-cycle: setting 1 ############
#----------------------------------------------------------------------------------------
alpha = 0.05
n = 200
MC = 200
t3 = pi/6
t2 = pi/3

little_power_mean = c()
little_power = c()
little_power_cov = c()
our_power = c()
our_power_corr = c()
our_power_mean = c()

R = c()

start.time = Sys.time()

for(t1 in seq(t2+t3, (pi + t2 + t3)/2, length.out = 8)){

  #----------------------------------------------------------------------------------------
  #### POPULATION LEVEL ######
  #----------------------------------------------------------------------------------------
  SigmaS=list()
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

  #----------------------------------------------------------------------------------------
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  #----------------------------------------------------------------------------------------
  little_decisions_mean = c()
  little_decisions = c()
  little_decisions_cov = c()
  our_decisions = c()
  our_decisions_corr = c()
  our_decisions_mean = c()

  for (i in 1:MC){
    #----------------------------------------------------------------------------------------
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    #----------------------------------------------------------------------------------------
    X1 = mvrnorm(n, c(0,0), SigmaS[[1]])
    X2 = mvrnorm(n, c(0,0), SigmaS[[2]])
    X3 = mvrnorm(n, c(0,0), SigmaS[[3]])

    columns = c("X1","X2","X3")
    X = data.frame(matrix(nrow = 3*n, ncol = 3))
    X[1:n, c("X1", "X2")] = X1
    X[(n+1):(2*n), c("X2", "X3")] = X2
    X[(2*n+1):(3*n), c("X1", "X3")] = X3
    X = as.matrix(X)

    #----------------------------------------------------------------------------------------
    ### run little's tests
    #----------------------------------------------------------------------------------------
    little_decisions_mean = c(little_decisions_mean, mcar_test(data.frame(X))$p.value < alpha)
    little_decisions = c(little_decisions, little_test(X, alpha))
    little_decisions_cov = c(little_decisions_cov, little_test(X, alpha, "cov"))

    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    our_decisions = c(our_decisions, MCAR_meancovTest(X, alpha, B = 99))
    our_decisions_corr = c(our_decisions_corr, MCAR_corrTest(X, alpha, B = 99))
    our_decisions_mean = c(our_decisions_mean, MCAR_meanTest(X, alpha, B = 99))
  }

  little_power_mean = c(little_power_mean, mean(little_decisions_mean))
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  our_power = c(our_power, mean(our_decisions))
  our_power_corr = c(our_power_corr, mean(our_decisions_corr))
  our_power_mean = c(our_power_mean, mean(our_decisions_mean))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("3_cycle_1.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, little_power_mean, col="green", ylim = c(0,1), pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, little_power, col="orange", pch=19, type = "b")
lines(R, little_power_cov, col="brown", pch=20, type = "b")
lines(R, our_power, col="blue", pch=21, type = "b")
lines(R, our_power_corr, col="darkviolet", pch=22, type = "b")
lines(R, our_power_mean, col="aquamarine2", pch=23, type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{cov}$)'),
                  TeX(r'($d^2_{aug}$)'), "Omnibus"),
       col = c("green", "brown", "orange", "blue"),
       pch = c(18, 19, 20, 21))
dev.off()

# #----------------------------------------------------------------------------------------
# ######### d-cycle: high-dimensional ############
# #----------------------------------------------------------------------------------------
d = 200

alpha = 0.05
n = 200
MC = 1 # with 4000 I expect 4 hours
angle = pi/(2*(d-1))

our_power_corr = c()
R = c()

start.time = Sys.time()

for(t1 in seq(pi/2, 5*pi/8, length.out = 8)){
  
  #----------------------------------------------------------------------------------------
  #### POPULATION LEVEL ######
  #----------------------------------------------------------------------------------------
  SigmaS=list()
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
  
  
  #----------------------------------------------------------------------------------------
  ###### SAMPLE LEVEL, REPEATING THE TEST MC TIMES #######
  #----------------------------------------------------------------------------------------
  our_decision_corr = c()
  for (i in 1:MC){
    
    #----------------------------------------------------------------------------------------
    #### generate dataset from patter S = {{1,2},{2,3},{1,3}}
    #----------------------------------------------------------------------------------------
    X = data.frame(matrix(nrow = d*n, ncol = d))
    for (i in 1:d){
      X[(1+(i-1)*n):(i*n), patterns[[i]]] = mvrnorm(n, rep(0, 2), SigmaS[[i]])
    }
    X = as.matrix(X)
    
    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    our_decision_corr = c(our_decision_corr, MCAR_covTest(X, alpha, B = 99))
  }
  
  our_power_corr = c(our_power_corr, mean(our_decision_corr))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("100_cycle.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, our_power_corr, col="darkviolet", ylim = c(0,1), pch=22,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.3,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("Corr"),
       col = c("darkviolet"),
       pch = c(22))
dev.off()
