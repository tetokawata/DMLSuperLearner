#' @title Sensitivity Analysis for AIPW estimators
#'
#' @description \code{SenstivityAIPW} is used to implement Robinson estimation with Super Learner
#' Y, M, G, Alpha, RY2, RD2
#' @param Y Outcome vector
#' @param M Plugin score vector
#' @param G Outcome predictors
#' @param Alpha Reiz representer
#' @param RY2 R2 between prediction error of Y and unobservables
#' @param RD2 R2 between estimated and ideal Reiz representer
#'
#' @return A list containing follows:
#' \item{Nuisance}{ estimated nuisance}
#' \item{Score}{ Score}
#' \item{Main}{ Estimation results}
#'
#' @export
#'
#' @examples
#'
#' # Example data
#' X <- matrix(rnorm(200,0,1),100)
#' D <- 2*X[,1] + runif(100,-10,10)
#' Y <- 3*X[,1] + rnorm(100,0,1)
#'
#' # Estimation
#' fit <- Robinson_SLearner(X = X, D = D, Y = Y)
#' fit$parameter_estimation
#'
#' @references
#'
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending omitted variable bias. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(1), 39-67.
#'
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2021). Omitted Variable Bias in Machine Learned Causal Models. arXiv preprint arXiv:2112.13398.
#'


EstPartialLinear <- function(Y, M, G, Alpha, RY2, RD2) {
  result <- list()
  TempN <- length(Y)
  TempRY2 <- RY2
  TempRD2 <- RD2

  CY2 <- TempRY2
  CD2 <- (1 - TempRD2) / TempRD2

  TempY <- Y
  TempX <- X
  TempD <- D
  TempM <- M
  TempAlpha <- Alpha
  TempG <- G

  result$Main$Theta <- mean(TempM + (TempY - TempG) * TempAlpha)
  result$Nuisance$Sigma2 <- mean((TempY - TempG)^2)
  result$Nuisance$V2 <- mean(TempAlpha^2)
  result$Nuisance$S2 <- result$Nuisance$Sigma2 * result$Nuisance$V2
  result$Nuisance$CY2 <- TempRY2
  result$Nuisance$CD2 <- (1 - TempRD2) / TempRD2

  result$Main$Lower <- result$Main$Theta -
    sqrt(result$Nuisance$S2) * sqrt(result$Nuisance$CY2) * sqrt(result$Nuisance$CD2)
  result$Main$Upper <- result$Main$Theta +
    sqrt(result$Nuisance$S2) * sqrt(result$Nuisance$CY2) * sqrt(result$Nuisance$CD2)

  result$Score$Theta <- result$Main$Theta - TempM - (TempY - TempG) * TempAlpha
  result$Score$Sigma2 <- result$Nuisance$Sigma2 - (TempY - TempG)^2
  result$Score$V2 <- result$Nuisance$V2 - TempAlpha^2
  result$Score$Lower <- result$Score$Theta -
    0.5*
    ((sqrt(result$Nuisance$CY2)*sqrt(result$Nuisance$CD2))/sqrt(result$Nuisance$S2))*
    (result$Nuisance$Sigma2*result$Score$V2 + result$Nuisance$V2*result$Score$Sigma2)
  result$Score$Upper <- result$Score$Theta +
    0.5*
    ((sqrt(result$Nuisance$CY2)*sqrt(result$Nuisance$CD2))/sqrt(result$Nuisance$S2))*
    (result$Nuisance$Sigma2*result$Score$V2 + result$Nuisance$V2*result$Score$Sigma2)
  result$Main$ThetaSD <- sd(result$Score$Theta)/sqrt(TempN)
  result$Main$UpperSD <- sd(result$Score$Upper)/sqrt(TempN)
  result$Main$LowerSD <- sd(result$Score$Lower)/sqrt(TempN)
  return(result)
}
