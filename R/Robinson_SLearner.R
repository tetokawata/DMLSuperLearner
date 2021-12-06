#' @title Robinson estimation with Super Learner
#'
#' @description \code{Robinson_SLearner} is used to implement Robinson estimation with Super Learner
#'
#' @param Y Outcome vector
#' @param D Treatment vector
#' @param X Control as data frame
#' @param SL.DML.library SL.library to estimate nuisance functions. Use first method if all SL coefficient are zero.
#' @param DML.V Number of folds
#' @param SL.V Number of folds in each SuperLearner
#' @param weights Sample weights
#'
#' @return A list containing follows:
#' \item{data}{ Data frame including estimated nuisance}
#' \item{parameter_estimation}{ Estimation result}
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
#' Robinson, P. M. (1988). Root-N-consistent semiparametric regression. Econometrica, 931-954.
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters.
#'
#' van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007) Super Learner. Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#'


Robinson_SLearner <- function(X,
                      D,
                      Y,
                      SL.DML.library = c("SL.ranger","SL.mean"),
                      DML.V = 2,
                      SL.V = 2,
                      core  = 1,
                      process = TRUE,
                      weights = NULL){
  require(magrittr)
  require(furrr)
  require(SuperLearner)
  # Weight
  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  # Data split
  DML.group <- sample(1:DML.V, length(Y), replace = TRUE)
  # Estimate nuisance
  eachDML <- function(g){
    temp = SuperLearner(X = X[DML.group != g,],
                                      Y = Y[DML.group != g],
                                      newX = X[DML.group == g,],
                                      SL.library = SL.DML.library,
                                      cvControl = SuperLearner.CV.control(V = SL.V))
    Y.hat = temp$library.predict[,1]
    if (is_empty(temp$SL.predict) == FALSE) {Y.hat = temp$SL.predict}
    Y.all <- temp$library.predict |> as.data.frame()
    colnames(Y.all) <- stringr::str_c("Y",SL.DML.library,sep = ".")
    temp = SuperLearner(X = X[DML.group != g,],
                                      Y = D[DML.group != g],
                                      newX = X[DML.group == g,],
                                      SL.library = SL.DML.library,
                                      cvControl = SuperLearner.CV.control(V = SL.V))
    D.hat = temp$library.predict[,1]
    if (is_empty(temp$SL.predict) == FALSE) {D.hat = temp$SL.predict}
    D.all <- temp$library.predict|> as.data.frame()
    colnames(D.all) <- stringr::str_c("D",SL.DML.library,sep = ".")
    tibble::tibble(Y = Y[DML.group == g],
                   D = D[DML.group == g],
                   Y.hat = Y.hat,
                   D.hat = D.hat,
                   X[DML.group == g,],
                   D.all,
                   Y.all,
                   weights = weights[DML.group == g])
  }
  # Iteration
  plan(multisession, workers = core)
  result <- list(data = future_map_dfr(1:DML.V,
                 eachDML,
                 .options = furrr_options(seed = 1L),
                 .progress = process)
                 )
  # Estimate R-learner
  result$parameter_estimation <-
    result$data |>
    dplyr::mutate(Y.oht = Y - Y.hat,
                  D.oht = D - D.hat) %>%
    estimatr::lm_robust(Y.oht ~ 0 + D.oht,
              data = .,
              weights = weights) |>
    generics::tidy() |>
    dplyr::mutate(outcome = "Y",
           term = "Rlearner")
  result
}

