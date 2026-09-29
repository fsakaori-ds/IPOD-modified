# Modified version of IPOD() from the R package "leapp".
#
# Modification by Fumitake Sakaori, 2026-09-29.
#
# The original IPOD() function has a length.out argument, but this
# argument is not used when constructing the lambda sequence.
# This modified version replaces the fixed step size (by = -0.1)
# with length.out = length.out.
#
# This file is distributed under the GNU General Public License
# version 3 or later.
#
# Reference:
# She, Y. and Owen, A. B. (2011).
# Outlier Detection Using Nonconvex Penalized Regression.
# Journal of the American Statistical Association, 106, 626-639.

IPOD2 <- function (X, Y, H, method = "hard", TOL = 1e-04, length.out = 50) 
{
  if (is.null(X)) {
    r = 0
  }
  else {
    r = ncol(X)
  }
  N = length(Y)
  ress = NULL
  gammas = NULL
  if (is.null(X)) {
    betaInit = NULL
    lambdas = seq(
      round(norm(matrix(Y, N, 1), "I") / 1 + 1),
      0,
      length.out = length.out
    )
  }
  else {
    betaInit = rlm(Y ~ X - 1)$coefficients
    tmp = t(diag(ncol(H)) - H) %*% Y / sqrt(1 - diag(H))
    lambdas = seq(
      round(norm(tmp, "I") / 1 + 1),
      0,
      length.out = length.out
    )
  }
  for (sigma in lambdas) {
    sigma = sigma/sqrt(2 * log(N))
    result = IPODFUN(X, Y, H, sigma, betaInit, method = method, 
                     TOL = TOL)
    gammas = cbind(gammas, result$gamma)
    ress = cbind(ress, result$ress)
  }
  DF = colSums(abs(gammas) > 1e-05)
  if (!is.null(X)) {
    Q = qr.Q(qr(X), complete = TRUE)
    X.unscaled = t(Q[, (r + 1):N])
    Y.eff = X.unscaled %*% Y
    sigmaSqEsts = colSums((Y.eff %*% matrix(rep(1, ncol(gammas)), 
                                            1, ncol(gammas)) - X.unscaled %*% gammas)^2)/(length(Y.eff) - 
                                                                                            DF)
  }
  else {
    sigmaSqEsts = colSums((Y %*% matrix(rep(1, ncol(gammas)), 
                                        1, ncol(gammas)) - gammas)^2)/(length(Y) - DF)
  }
  sigmaSqEsts[sigmaSqEsts < 0] = 0
  if (!all(DF <= N - r)) {
    gammas = gammas[, DF <= N - r]
    sigmaSqEsts = sigmaSqEsts[DF <= N - r]
    lambdas = lambdas[DF <= N - r]
    ress = ress[, DF <= N - r]
    DF = DF[DF <= N - r]
  }
  redundantInds = which(sum(abs(gammas[, -1] - gammas[, 1:(ncol(gammas) - 
                                                             1)])) < 0.001) + 1
  if (length(redundantInds) > 0) {
    gammas = gammas[, -redundantInds]
    lambdas = lambdas[-redundantInds]
    sigmaSqEsts = sigmaSqEsts[-redundantInds]
    ress = ress[, -redundantInds]
    DF = DF[-redundantInds]
  }
  redundantInds = which(abs(sigmaSqEsts[-1] - sigmaSqEsts[1:(length(sigmaSqEsts) - 
                                                               1)]) < 0.001) + 1
  if (length(redundantInds) > 0) {
    gammas = gammas[, -redundantInds]
    lambdas = lambdas[-redundantInds]
    sigmaSqEsts = sigmaSqEsts[-redundantInds]
    ress = ress[, -redundantInds]
    DF = DF[-redundantInds]
  }
  if (!all(DF <= N/2)) {
    gammas = gammas[, DF <= N/2]
    sigmaSqEsts = sigmaSqEsts[DF <= N/2]
    lambdas = lambdas[DF <= N/2]
    ress = ress[, DF <= N/2]
    DF = DF[DF <= N/2]
  }
  riskEst = ((N - r) * log(sigmaSqEsts * (N - r - DF)/(N - 
                                                         r)) + (log(N - r) + 1) * (DF + 1))/(N - r)
  optSet = which(riskEst == min(riskEst))
  gammaOptInd = optSet[which(DF[optSet] == min(DF[optSet]))[1]]
  gammaOpt = gammas[, gammaOptInd]
  resOpt = ress[, gammaOptInd]
  tau = mad(ress[gammas[, gammaOptInd] == 0, gammaOptInd])
  resOpt.scale = resOpt/tau
  p = 2 * pnorm(-abs(resOpt.scale))
  return(list(p = p, resOpt.scale = resOpt.scale, gamma = gammaOpt))
}
