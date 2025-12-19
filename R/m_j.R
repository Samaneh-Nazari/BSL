#' Computes node-marginal likelihood of a DAG model (j in {2,...,q})
#'
#' This function computes the logarithm marginal likelihood for a specific node j (where j is in {2,...,q}) in a DAG,
#' given the data X. This calculation is based on a DAG-Probit model with a Normal-Gamma
#' prior distribution and utilizes a Laplace approximation for calculating marginal likelihood.
#'
#'
#' @param j An integer representing the label of the \code{j} in \code{dag} (must be > 1).
#' @param dag A \eqn{q,q} adjacency matrix representing the \code{dag} structure.
#' @param X An \eqn{n,q} data matrix of observations.
#' @param n number of observations (rows) in the data matrix \eqn{X}.
#' @param a A hyperparameter of the prior distribution.
#' @param g A hyperparameter of the prior distribution, related to the prior variance.
#'
#' @return The logarithm of the marginal likelihood for the specified \code{j}.
#'
#' @details
#' The function handles two cases: when the node has no parents and when it has one or more parents.
#' The calculation involves terms from the Normal-Gamma prior and approximations based on the data.
#' The formula for \code{a} hyperparameter has been adjusted to be \code{nrow(X)-1}.
#'
#' @export
m_j = function(j, dag, X, n, a, g){

  ## This function computes the marginal likelihood of a dag given the data X relative to a generic node j in {2,...,q}

  # n is the sample size
  # a,g are hyperparameters of the prior distribution

  pa_j = BSL::pa(j, DAG=dag)

  y  = X[,j]
  XX = as.matrix(X[,pa_j])
  q=ncol(X)

  p_j   = length(pa_j)
  a=n-1
  j_pos = q - p_j
  a_j_star = (a + q - 2*j_pos + 3)/2 - p_j/2 - 1

  if(length(pa_j) == 0){
    a_bessel=a_j_star-n/2
    beta=0.5*g
    alpha=0.5*sum(y^2)
    m = - 0.5*n*log(2*pi)-lgamma(a_j_star) + a_j_star*log(0.5*g)+log(sqrt(pi))+0.25*log((alpha^(2*a_bessel-1))/(beta^(2*a_bessel+1)))-2*sqrt(alpha*beta)
  }
  else{

    T_j = g*diag(1, p_j)

    T_j_hat = T_j + t(XX)%*%XX
    L_j_hat = solve(T_j_hat)%*%t(XX)%*%y
    a_bessel=a_j_star-n/2
    beta=0.5*g
    alpha=0.5*(sum(y^2)-t(L_j_hat)%*%T_j_hat%*%L_j_hat)
    m = - 0.5*n*log(2*pi) +0.5*log(det(T_j)) - 0.5*log(det(T_j_hat)) +a_j_star*log(0.5*g) - lgamma(a_j_star)+log(sqrt(pi))+0.25*log((alpha^(2*a_bessel-1))/(beta^(2*a_bessel+1)))-2*sqrt(alpha*beta)
  }

  return(m)

}

