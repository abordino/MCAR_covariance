# Load external functions and libraries
source("computeR.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
source("indexConsistency.R")
source("little_test.R")

library(missMethods)
library(MASS)
library(norm)
library(latex2exp)
library(naniar)
library(compositions)

#----------------------------------------------------------------------------------------
######### 3-cycle: lognormal ############
#----------------------------------------------------------------------------------------
alpha = 0.05
n = 200
MC = 1000
t3 = pi/4
t2 = 3*pi/4
d = 3

little_power = c()
little_power_cov = c()
our_power = c()
combined_power = c()
our_power_corr = c()

R = c()
for(t1 in seq(pi/2 + pi/12, pi/12, length.out = 12)){

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
R

for(t1 in seq(pi/2 + pi/12, pi/12, length.out = 12)){

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
  little_decisions_cov = c()
  our_decisions = c()
  combined_decisions = c()
  our_decisions_corr = c()

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
    p_aug = little_test(X, alpha)
    little_decisions = c(little_decisions,  p_L < alpha)
    little_decisions_cov = c(little_decisions_cov, p_aug < alpha)
    
    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    p_R = corr.compTest(X, B = 99)
    p_M = mean.consTest(X, B = 99)
    p_V = var.consTest(X, B = 99)
    
    combined_decisions = c(combined_decisions, max(c(p_R, p_M, p_V) < alpha/3))
    our_decisions = c(our_decisions, max(c(p_R, p_L, p_V) < alpha/3))
    our_decisions_corr = c(our_decisions_corr, p_R < alpha)
  }

  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  combined_power = c(combined_power, mean(combined_decisions))
  our_power = c(our_power, mean(our_decisions))
  our_power_corr = c(our_power_corr, mean(our_decisions_corr))
}

save(R, little_power, little_power_cov, our_power, combined_power, our_power_corr, alpha, 
     file = "pictures/3_cycle_Lnorm.RData")
png("pictures/3_cycle_Lnorm.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, little_power, col="green", pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylim = c(0,1), ylab = "Power", type = "b")
lines(R, little_power_cov, col="black", pch=19, type = "b")
lines(R, our_power, col="blue", pch=21, type = "b")
lines(R, combined_power, col="orange", pch=20, type = "b")
lines(R, our_power_corr, col="darkviolet", pch=25, type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{aug}$)'), "Combined", "Omnibus", TeX(r'($p_R$)')),
       col = c("green", "black", "orange", "blue", "darkviolet"),
       pch = c(18, 19, 20, 21, 25))
dev.off()
