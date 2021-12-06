#' @title AIPW with Super Learner
#'
#' @description \code{AIPW_SLearner} is used to implement AIPW with Super Learner. Current version trim sample with propensity <= 0.05 or >= 0.95
#'
#' @param Y Outcome vector
#' @param D Treatment vector (Binary)
#' @param X Control as data frame
#' @param SL.DML.library SL.library to estimate nuisance functions. Use first method if all SL coefficient are zero.
#' @param DML.V Number of folds
#' @param SL.V Number of folds in each SuperLearner
#' @param weights Sample weights
#' @param D.family gaussian() or binomial() to describe the error distribution.
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
#' fit <- AIPW_SLearner(X = X, D = D, Y = Y)
#' fit$parameter_estimation
#'
#' @references
#' Robins, J. M., Rotnitzky, A., & Zhao, L. P. (1994). Estimation of regression coefficients when some regressors are not always observed. Journal of the American statistical Association, 89(427), 846-866.
#'
#' van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007) Super Learner. Statistical Applications of Genetics and Molecular Biology, 6, article 25.
#'
#' Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., & Newey, W. (2017). Double/debiased/neyman machine learning of treatment effects. American Economic Review, 107(5), 261-65.
#'

AIPW_SLearner <- function(X,
                          D,
                          Y,
                          SL.DML.library = c("SL.ranger","SL.mean"),
                          DML.V = 2,
                          SL.V = 2,
                          core  = 1,
                          process = TRUE,
                          weights = NULL,
                          D.family = binomial()){
  require(magrittr)
  require(furrr)
  require(SuperLearner)
  # weight
  if(is.null(weights)) {
    weights <- rep(1, length(Y))
  }
  # Grouping
  DML.group <- sample(1:DML.V, length(Y), replace = TRUE)
  # Each nuisance
  eachDML <- function(g){
    # Y1
    temp = SuperLearner(X = X[DML.group != g & D == 1,],
                                      Y = Y[DML.group != g & D == 1],
                                      newX = X[DML.group == g,],
                                      SL.library = SL.DML.library,
                                      cvControl = SuperLearner.CV.control(V = SL.V)
                                      )
    Y1.hat = temp$library.predict[,1]
    if (is_empty(temp$SL.predict) == FALSE) {Y1.hat = temp$SL.predict}
    Y1.all <- temp$library.predict |> as.data.frame()
    colnames(Y1.all) <- stringr::str_c("Y1",SL.DML.library,sep = ".")
    # Y0
    temp = SuperLearner(X = X[DML.group != g & D == 0,],
                                      Y = Y[DML.group != g & D == 0],
                                      newX = X[DML.group == g,],
                                      SL.library = SL.DML.library,
                                      cvControl = SuperLearner.CV.control(V = SL.V)
                                      )
    Y0.hat = temp$library.predict[,1]
    if (is_empty(temp$SL.predict) == FALSE) {Y0.hat = temp$SL.predict}
    Y0.all <- temp$library.predict |> as.data.frame()
    colnames(Y0.all) <- stringr::str_c("Y0",SL.DML.library,sep = ".")
    # D
    temp = SuperLearner(X = X[DML.group != g,],
                                      Y = D[DML.group != g],
                                      newX = X[DML.group == g,],
                                      SL.library = SL.DML.library,
                                      cvControl = SuperLearner.CV.control(V = SL.V),
                                      family = D.family
    )
    D.hat = temp$library.predict[,1]
    if (is_empty(temp$SL.predict) == FALSE) {D.hat = temp$SL.predict}
    D.all <- temp$library.predict |> as.data.frame()
    colnames(D.all) <- stringr::str_c("D",SL.DML.library,sep = ".")
    tibble::tibble(Y = Y[DML.group == g],
                   D = D[DML.group == g],
                   Y1.hat = Y1.hat[,1],
                   Y0.hat = Y0.hat[,1],
                   D.hat = D.hat[,1],
                   X[DML.group == g,],
                   D.all,
                   Y1.all,
                   Y0.all,
                   weights = weights[DML.group == g])
  }
  # Iteration
  plan(multisession, workers = core)
  result <- list(data = future_map_dfr(1:DML.V,
                                       eachDML,
                                       .options = furrr_options(seed = 1L),
                                       .progress = process)
  )
  # Store results
  result$parameter_estimation <-
    result$data |>
    dplyr::filter(D.hat >= 0.05 & D.hat <= 0.95) |>
    dplyr::mutate(Y1.potential = Y1.hat + (D/D.hat)*(Y - Y1.hat),
           Y0.potential = Y0.hat + ((1-D)/(1-D.hat))*(Y - Y0.hat),
           aipw = Y1.potential - Y0.potential) %>%
    estimatr::lm_robust(aipw ~ 1,
              data = .,
              weights = weights) |>
    generics::tidy() |>
    dplyr::mutate(outcome = "ATE",
           term = "AIPW")
  result
}

