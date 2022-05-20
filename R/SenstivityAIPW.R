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
#' N <- 1000
#' X <- runif(N,0,1)
#' D <- sample(0:1,N,replace = TRUE)
#' Y <- D + X + rnorm(N,0,1)
#'
#' Fit <- AIPW$new(Y = Y,
#' A = D,
#' W = X,
#' Q.SL.library = c("SL.lm"),
#' g.SL.library = c("SL.lm")
#' )
#' Fit$fir()$Fit$summary()
#'
#' EstG <- Fit$obs_est$mu
#' EstM <- Fit$obs_est$mu1 - Fit$obs_est$mu0
#' HatD <- Fit$obs_est$p_score
#' EstAlpha <- D/HatD - (1-D)/(1-HatD)
#' FitTest <- EstPartialLinear(Y = Y, M = EstM, G = EstG, Alpha = EstAlpha, RY2 = 0.05, RD2 = 0.95)
#' tibble(PointATE = FitTest$Main$Theta,
#' PointUpper95th = FitTest$Main$Theta + 1.96*FitTest$Main$ThetaSD,
#' PointLower95th = FitTest$Main$Theta - 1.96*FitTest$Main$ThetaSD,
#' BoundUpper = FitTest$Main$Upper,
#' BoundLower = FitTest$Main$Lower,
#' BoundUpper95th = FitTest$Main$Upper + 1.96*FitTest$Main$UpperSD,
#' BoundLower95th = FitTest$Main$Lower - 1.96*FitTest$Main$LowerSD)
#'
#'
#' @references
#'
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending omitted variable bias. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(1), 39-67.
#'
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2021). Omitted Variable Bias in Machine Learned Causal Models. arXiv preprint arXiv:2112.13398.
#'


SenstivityAIPW <- function(Y, M, G, Alpha, RY2, RD2) {
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
