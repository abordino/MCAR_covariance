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

# #----------------------------------------------------------------------------------------
# ######### d-cycle: high-dimensional ############
# #----------------------------------------------------------------------------------------
d = 100

alpha = 0.05
n = 200
MC = 50 # with 4000 I expect 4 hours
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
    our_decision_corr = c(our_decision_corr, MCAR_meancovTest(X, alpha, B = 99))
  }
  
  our_power_corr = c(our_power_corr, mean(our_decision_corr))
}

end.time = Sys.time()
time.taken = round(end.time - start.time,2)
time.taken

png("100_cycle.png")
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(R, our_power_corr, col="blue", ylim = c(0,1), pch=21,
     xlab = TeX(r'($R(\Sigma_\$)$)'), ylab = "Power", type = "b")
lines(R, rep(alpha, length(R)), lty = 2, col = "red")
legend("right", inset = c(-0.3,0), xpd = TRUE,
       horiz = FALSE, lty = 1, bty = "n",
       legend = c("Omnibus"),
       col = c("blue"),
       pch = c(21))
dev.off()
