source("computeR.R")
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
n = 60
MC = 30
t3 = pi/6
t2 = pi/6
d = 3

little_power = c()
combined_power = c()
our_power = c()

R = c()
for(t1 in seq(t2+t3-pi/6, pi-pi/6, length.out = 8)){

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

for(t1 in seq(t2+t3-pi/6, pi-pi/6, length.out = 8)){

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
  little_decisions = c()
  combined_decisions = c()
  our_decisions = c()

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
    p_L = mcar_test(data.frame(X))$p.value
    little_decisions = c(little_decisions,  p_L < alpha)

    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    p_R = corr.compTest(X, B= 99)
    p_M = mean.consTest(X, B= 99)
    p_V = bootstrapTestV(X, B = 99)
    combined_decisions = c(combined_decisions, 
                           -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6))
    our_decisions = c(our_decisions, 
                      -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6))
  }

  little_power = c(little_power, mean(little_decisions))
  combined_power = c(combined_power, mean(combined_decisions))
  our_power = c(our_power, mean(our_decisions))
}


png("pictures/3_cycle_lognormal_1.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, little_power, col="green", ylim = c(0,1), pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, our_power, col="blue", pch=21, type = "b")
lines(R, combined_power, col="orange", pch=20, type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
       col = c("green", "orange", "blue"),
       pch = c(18, 20, 21))
dev.off()

#----------------------------------------------------------------------------------------
# ######### 3-cycle: setting 1 ############
#----------------------------------------------------------------------------------------
alpha = 0.05
n = 60
MC = 30
t3 = pi/6
t2 = pi/3

little_power = c()
our_power = c()
combined_power = c()

R = c()
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
  our_decisions = c()
  combined_decisions = c()

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
    ### run little's test
    #----------------------------------------------------------------------------------------
    p_L = mcar_test(data.frame(X))$p.value
    little_decisions = c(little_decisions,  p_L < alpha)
    
    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    p_R = corr.compTest(X, B= 99)
    p_M = mean.consTest(X, B= 99)
    p_V = bootstrapTestV(X, B = 99)
    combined_decisions = c(combined_decisions, 
                           -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6))
    our_decisions = c(our_decisions, 
                      -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6))
  }

  little_power = c(little_power, mean(little_decisions))
  combined_power = c(combined_power, mean(combined_decisions))
  our_power = c(our_power, mean(our_decisions))
}

png("pictures/3_cycle_1.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, little_power, col="green", ylim = c(0,1), pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, our_power, col="blue", pch=21, type = "b")
lines(R, combined_power, col="orange", pch=20, type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
       col = c("green", "orange", "blue"),
       pch = c(18, 20, 21))
dev.off()
