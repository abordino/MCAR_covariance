source("computeR.R")
source("bootstrap_test.R")
source("find_SigmaS.R")
source("indexConsistency.R")
library(missMethods)
library(MASS)
library(norm)
library(latex2exp)
library(future.apply)
library(naniar)

alpha = 0.05
n = 100
MC = 40
xxx =  seq(0.05, 0.4, length.out=7)

# #----------------------------------------------------------------------------------------
# #### miss methods
# #----------------------------------------------------------------------------------------

d = 3
yyy = "exp"

  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------

  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(rate = 1),d)) #list(c(mean = 0, sd = 1))

  data = rMvdc(n, P)

  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------

  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)

  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)

    for (i in 1:MC){
      X = delete_MAR_rank(data, p, c(1,2), cols_ctrl = c(3, 3))

      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }

    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)

    ind = ind+1
  }

  png(paste("pictures/", yyy, "_", "np", "_MAR1.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()

  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(rate = 1),d))

  data = rMvdc(n, P)

  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------

  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3, 3), x = 9)
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }

  png(paste("pictures/", yyy, "_", "np", "_MAR2.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()

  #----------------------------------------------------------------------------------------
  ######### MCAR ############
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(rate = 1),d))

  data = rMvdc(n, P)

  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MCAR(data, p, c(1,2))
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }

  png(paste("pictures/", yyy,"_","np","_MCAR.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Size", type="b", main = "")
  lines(xxx, rep(alpha, length(xxx)), lty = 2, col = "red")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
  
  
  
d = 3
yyy = "norm"

  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MAR_rank(data, p, c(1,2), cols_ctrl = c(3, 3))
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy, "_", "np", "_MAR1.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
  
  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3, 3), x = 9)
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy, "_", "np", "_MAR2.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
  
  #----------------------------------------------------------------------------------------
  ######### MCAR ############
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MCAR(data, p, c(1,2))
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy,"_","np","_MCAR.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Size", type="b", main = "")
  lines(xxx, rep(alpha, length(xxx)), lty = 2, col = "red")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()

  
  d = 5
  yyy = "lnorm"
  
  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MAR_rank(data, p, c(1,2), cols_ctrl = c(3, 4))
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy, "_", "np", "_MAR1.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
  
  #----------------------------------------------------------------------------------------
  ######### MAR ##############################
  #----------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MAR_1_to_x(data, p, c(1,2), cols_ctrl = c(3, 4), x = 9)
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy, "_", "np", "_MAR2.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Power", type="b", main = "")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
  
  #----------------------------------------------------------------------------------------
  ######### MCAR ############
  #----------------------------------------------------------------------------------------
  
  #----------------------------------------------------------------------------------------
  # Select the copula
  #----------------------------------------------------------------------------------------
  cp = claytonCopula(param = c(1), dim = d)
  #----------------------------------------------------------------------------------------
  # Generate the multivariate distribution
  #----------------------------------------------------------------------------------------
  P = mvdc(copula = cp, margins = c(rep(yyy,d)),
           paramMargins = rep(list(c(mean = 0, sd = 1)),d))
  
  data = rMvdc(n, P)
  
  #----------------------------------------------------------------------------------------
  # run little's and our tests
  #----------------------------------------------------------------------------------------
  
  little_power = numeric(length = 7)
  combined_power = numeric(length = 7)
  our_power = numeric(length = 7)
  
  ind = 1
  for (p in xxx){
    little_decision = logical(length = MC)
    combined_decision = logical(length = MC)
    our_decision = logical(length = MC)
    
    for (i in 1:MC){
      X = delete_MCAR(data, p, c(1,2))
      
      p_L = mcar_test(data.frame(X))$p.value
      p_R = corr.compTest(X, B= 99)
      p_M = mean.consTest(X, B= 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] =  p_L < alpha
      
      # our_decision[i] = p_V < alpha
     our_decision[i] = -2*(log(p_R)+log(p_M)+log(p_V)) > qchisq(1-2*alpha/3, 6)
      combined_decision[i] = -2*(log(p_R)+log(p_L)+log(p_V)) > qchisq(1-2*alpha/3, 6)
    }
    
    little_power[ind] = mean(little_decision)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    
    ind = ind+1
  }
  
  png(paste("pictures/", yyy,"_","np","_MCAR.png"))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(xxx, little_power, col="green", ylim = c(0,1), pch=18,
       xlab = "Missingness probability p", ylab = "Size", type="b", main = "")
  lines(xxx, rep(alpha, length(xxx)), lty = 2, col = "red")
  lines(xxx, our_power, col="blue", pch=21, type = "b")
  lines(xxx, combined_power, col="orange", pch=20, type = "b")
  legend("right", inset = c(-0.4,0), xpd = TRUE,
         horiz = FALSE, lty = 1, bty = "n",
         legend = c(TeX(r'($d^2_\mu$)'), "Combined", "Omnibus"),
         col = c("green", "orange", "blue"),
         pch = c(18, 20, 21))
  dev.off()
