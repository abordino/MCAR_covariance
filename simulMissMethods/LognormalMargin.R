rm(list = ls())  # Clear environment
gc()             # Free memory

setwd("~/Desktop/simulation_regularised/")

# Load external functions and libraries
source("MCARtest/computeR.R")
source("MCARtest/find_SigmaS.R")
source("MCARtest/indexConsistency.R")
source("bootstrap_test.R")
source("little_test.R")

library(missMethods)
library(MASS)
library(norm)
library(latex2exp)
library(naniar)

set.seed(221198)

# Constants and parameters
alpha = 0.05
n = 200
MC = 500
xxx = seq(0.05, 0.4, length.out = 7)

#----------------------------------------------------------------------------------------
# Helper function to run the tests and calculate power
#----------------------------------------------------------------------------------------
run_tests <- function(data_gen_func, file_suffix, method, x = NULL, under_null = FALSE) {
  little_power = numeric(length(xxx))
  little_power_cov = numeric(length(xxx))
  combined_power = numeric(length(xxx))
  our_power = numeric(length(xxx))
  our_power_corr = numeric(length(xxx))
  
  
  for (ind in seq_along(xxx)) {
    p = xxx[ind]
    little_decision = logical(MC)
    little_decision_cov = logical(MC)
    combined_decision = logical(MC)
    our_decision = logical(MC)
    our_decision_corr = logical(MC)
    
    for (i in 1:MC) {
      X = data_gen_func(p, x)
      
      p_L = mcar_test(data.frame(X))$p.value
      p_aug = little_test(X, alpha)
      
      p_R = corr.compTest(X, B = 99)
      p_M = mean.consTest(X, B = 99)
      p_V = var.consTest(X, B = 99)
      
      little_decision[i] = p_L < alpha
      little_decision_cov[i] = p_aug < alpha
      
      our_decision[i] = max(c(p_R, p_M, p_V) < alpha/3)
      combined_decision[i] = max(c(p_R, p_L, p_V) < alpha/3)
      our_decision_corr[i] = p_R < alpha
    }
    
    little_power[ind] = mean(little_decision)
    little_power_cov[ind] = mean(little_decision_cov)
    our_power[ind] = mean(our_decision)
    combined_power[ind] = mean(combined_decision)
    our_power_corr[ind] = mean(our_decision_corr)
  }
  
  if (under_null == TRUE){
    # Save the data
    data.to.save = data.frame(xxx, little_power, little_power_cov, our_power, 
                              combined_power, our_power_corr)
    write.csv(data.to.save, 
              file = paste0("simulMissMethods/results/", yyy, "_", file_suffix, ".csv"),
              row.names = FALSE)
    
    # Plot the results
    png(paste0("simulMissMethods/pictures/", yyy, "_", file_suffix, ".png"))
    par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    plot(xxx, little_power, col = "green", ylim = c(0, 1), pch = 18,
         xlab = "Missingness probability p", ylab = "Size", type = "b", main = "")
    lines(xxx, little_power_cov, col = "black", pch = 19, type = "b")
    lines(xxx, rep(alpha, length(xxx)), lty = 2, col = "red")
    lines(xxx, our_power, col = "blue", pch = 21, type = "b")
    lines(xxx, combined_power, col = "orange", pch = 20, type = "b")
    lines(xxx, our_power_corr, col="darkviolet", pch=25, type = "b")
    legend("right", inset = c(-0.4,0), xpd = TRUE,
           horiz = FALSE, lty = 1, bty = "n",
           legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{aug}$)'), "Combined", "Omnibus", TeX(r'($p_R$)')),
           col = c("green", "black", "orange", "blue", "darkviolet"),
           pch = c(18, 19, 20, 21, 25))
    dev.off()
  } else{
    # Save the data
    data.to.save = data.frame(xxx, little_power, little_power_cov, our_power, 
                              combined_power, our_power_corr)
    write.csv(data.to.save, 
              file = paste0("simulMissMethods/results/", yyy, "_", file_suffix, ".csv"),
              row.names = FALSE)
    
    # Plot the results
    png(paste0("simulMissMethods/pictures/", yyy, "_", file_suffix, ".png"))
    par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    plot(xxx, little_power, col = "green", ylim = c(0, 1), pch = 18,
         xlab = "Missingness probability p", ylab = "Power", type = "b", main = "")
    lines(xxx, little_power_cov, col = "black", pch = 19, type = "b")
    lines(xxx, our_power, col = "blue", pch = 21, type = "b")
    lines(xxx, combined_power, col = "orange", pch = 20, type = "b")
    lines(xxx, our_power_corr, col="darkviolet", pch=25, type = "b")
    legend("right", inset = c(-0.4,0), xpd = TRUE,
           horiz = FALSE, lty = 1, bty = "n",
           legend = c(TeX(r'($d^2_\mu$)'), TeX(r'($d^2_{aug}$)'), "Combined", "Omnibus", TeX(r'($p_R$)')),
           col = c("green", "black", "orange", "blue", "darkviolet"),
           pch = c(18, 19, 20, 21, 25))
    dev.off()
  }
  
}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Run simulations for d = 5 and lognormal data
#----------------------------------------------------------------------------------------
d = 5
yyy = "lnorm"


#----------------------------------------------------------------------------------------
# Data generation functions
#----------------------------------------------------------------------------------------
generate_data_MAR1 <- function(p, x = NULL) {
  delete_MAR_rank(data, p, c(1, 2), cols_ctrl = c(3, 4))
}

generate_data_MAR2 <- function(p, x) {
  delete_MAR_1_to_x(data, p, c(1, 2), cols_ctrl = c(3, 4), x = x)
}

generate_data_MCAR <- function(p, x = NULL) {
  delete_MCAR(data, p, c(1, 2))
}

#----------------------------------------------------------------------------------------
# Generate base data (Clayton Copula)
#----------------------------------------------------------------------------------------
cp = claytonCopula(param = c(1), dim = d)
P = mvdc(copula = cp, margins = rep(yyy, d), paramMargins = rep(list(c(mean = 0, sd = 1)),d))
data = rMvdc(n, P)

#----------------------------------------------------------------------------------------
# Run tests for MAR (Method 1)
#----------------------------------------------------------------------------------------
run_tests(generate_data_MAR1, "MAR1", "MAR")

#----------------------------------------------------------------------------------------
# Run tests for MAR (Method 2)
#----------------------------------------------------------------------------------------
run_tests(generate_data_MAR2, "MAR2", "MAR", x = 9)

#----------------------------------------------------------------------------------------
# Run tests for MCAR
#----------------------------------------------------------------------------------------
run_tests(generate_data_MCAR, "MCAR", "MCAR", under_null = TRUE)