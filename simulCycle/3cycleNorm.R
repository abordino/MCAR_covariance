rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Documents/phd/MCAR_covariance/simulation_final/")

# Load external functions and libraries
source("MCARtest/computeR.R")
source("MCARtest/find_SigmaS.R")
source("MCARtest/indexConsistency.R")
source("bootstrap_test.R")
source("bootstrap_test_cycle.R")
source("little_test.R")

library(missMethods)
library(MASS)
library(norm)
library(latex2exp)
library(naniar)

#----------------------------------------------------------------------------------------
# ######### 3-cycle: setting 1 ############
#----------------------------------------------------------------------------------------
set.seed(200801)

alpha = 0.05
n = 200
MC = 500
t3 = pi/6
t2 = pi/3

little_power = c()
little_power_cov = c()
our_power = c()
combined_power = c()
our_power_corr = c()


R = c()
for(t1 in seq(t2+t3, (pi + t2 + t3)/2, length.out = 8)){
  
  print(t1)
  
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
  little_decisions = c()
  little_decisions_cov = c()
  our_decisions = c()
  combined_decisions = c()
  our_decisions_corr = c()
  
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
    p_aug = little_test(X, alpha)
    little_decisions = c(little_decisions,  p_L < alpha)
    little_decisions_cov = c(little_decisions_cov, p_aug < alpha)
    
    #----------------------------------------------------------------------------------------
    ### run our tests
    #----------------------------------------------------------------------------------------
    p_R = corr.compTest(X, B = 99)
    p_M = mean.consTest(X, B = 99)
    p_V = var.consTest(X, B = 99)
    
    combined_decisions = c(combined_decisions, max(c(p_R, p_L, p_V) < alpha/3))
    our_decisions = c(our_decisions, max(c(p_R, p_M, p_V) < alpha/3))
    our_decisions_corr = c(our_decisions_corr, p_R < alpha)
  }
  
  little_power = c(little_power, mean(little_decisions))
  little_power_cov = c(little_power_cov, mean(little_decisions_cov))
  combined_power = c(combined_power, mean(combined_decisions))
  our_power = c(our_power, mean(our_decisions))
  our_power_corr = c(our_power_corr, mean(our_decisions_corr))
}

#----------------------------------------------------------------------------------------
# ######### save data ############
#----------------------------------------------------------------------------------------
data_to_save = data.frame(R, little_power, little_power_cov, our_power, 
                          combined_power, our_power_corr, alpha)
write.csv(data_to_save, "simulCycle/results/3_cycle_Norm.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------
# ######### plot results ############
#----------------------------------------------------------------------------------------
data = read.csv("simulCycle/results/3_cycle_Norm.csv")

png("simulCycle/pictures/3_cycle_Norm.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(data$R, data$little_power, col="green", pch=18,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylim = c(0,1), ylab = "Power", type = "b")
lines(data$R, data$little_power_cov, col="black", pch=19, type = "b")
lines(data$R, data$our_power, col="blue", pch=21, type = "b")
lines(data$R, data$combined_power, col="orange", pch=20, type = "b")
lines(data$R, data$our_power_corr, col="darkviolet", pch=25, type = "b")
lines(data$R, rep(alpha, length(data$R)), lty = 2, col = "red")
legend("right", inset = c(-0.4,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{aug}$)'), "Combined", "Omnibus", TeX(r'($p_R$)')),
       col = c("green", "black", "orange", "blue", "darkviolet"),
       pch = c(18, 19, 20, 21, 25))
dev.off()