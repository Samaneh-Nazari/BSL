#' Prepare Data for DAG-Probit Model
#'
#' This function takes a data matrix \code{Y} and preprocesses it for the DAG-probit model.
#' It separates the first column as the binary response variable \code{y} and the remaining
#' columns as the covariate matrix \code{X}. The response variable is then converted
#' \code{to binary (0/1)} based on a cutoff value.
#'
#' @param Y An \eqn{n,q} data matrix where the first column is a continuous variable
#'   that will be converted to a binary response.
#' @param cutoff A numeric value. Values greater than the cutoff are set to 1,
#'   and values less than or equal to the cutoff are set to 0. Defaults to 0.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{y}: An \eqn{n,1} numeric vector of the binary response variable.
#'   \item \code{X}: An \eqn{n,q-1} numeric matrix of covariates.
#' }
#'
#' @export
#'
prepare_data <- function(Y, cutoff = 0) {
  y <- Y[,1]
  X <- Y[,-1]

  # Convert the continuous y to binary
  y[y > cutoff] <- 1
  y[y < cutoff] <- 0

  return(list(y = y, X = X))
}
