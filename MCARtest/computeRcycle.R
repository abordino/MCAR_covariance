library(Rcsdp)
library(Matrix)
library(npmr)
library(matrixcalc)
library(pracma)


computeRcycle = function(SigmaS, tol = 1e-12, return_x = TRUE) {
  ## ------------------------------------------------------------
  ##  Input checks and extraction of cos(theta_i)
  ## ------------------------------------------------------------
  if (!is.list(SigmaS) || length(SigmaS) < 2)
    stop("SigmaS must be a list of at least two 2x2 matrices.")
  
  cos_theta =vapply(SigmaS, function(M) {
    if (!is.matrix(M) || any(dim(M) != 2))
      stop("Every element of SigmaS must be a 2x2 matrix.")
    M[1, 2]
  }, numeric(1))
  
  d =length(cos_theta)
  cos_theta1 =cos_theta[1]
  
  ## ------------------------------------------------------------
  ##  Helper: given s = cos(x1), compute f(s) = x1 - sum_{j>=2} xj
  ## ------------------------------------------------------------
  f =function(s) {
    if (s <= -1 || s >= 1) return(NaN)    
    x1 =acos(s)
    cos_xj =1 - ((1 - cos_theta[-1]) / (1 + cos_theta1)) * (1 + s)
    cos_xj =pmin(1, pmax(-1, cos_xj))
    xj =acos(cos_xj)
    x1 - sum(xj)
  }
  
  ## ------------------------------------------------------------
  ##  Root-finding: locate s* such that f(s*) = 0
  ## ------------------------------------------------------------
  lower =-1 + 1e-9
  upper = 1 - 1e-9
  fl =f(lower)
  fu =f(upper)
  
  if (is.nan(fl) || is.nan(fu) || fl * fu > 0)
    stop("Could not bracket a root.  Check that the input angles satisfy the theory assumptions.")
  
  root =uniroot(f, interval = c(lower, upper), tol = tol)$root
  x1   =acos(root)
  
  ## complete x-vector
  if (return_x) {
    cos_xj =1 - ((1 - cos_theta[-1]) / (1 + cos_theta1)) * (1 + root)
    cos_xj =pmin(1, pmax(-1, cos_xj))
    x =c(x1, acos(cos_xj))
  }
  
  ## ------------------------------------------------------------
  ##  Final value  R_cycle = 1 + (1 - cos(theta1)) / (1 + cos(x1))
  ## ------------------------------------------------------------
  Rcycle = 1 - (1 + cos_theta1) / (1 + root)
  
  if (return_x) {
    return(list(R = Rcycle, x = x))
  } else {
    return(Rcycle)
  }
}