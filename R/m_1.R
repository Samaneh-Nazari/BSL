#' Computes marginal likelihood for node 1
#'
#' Computes the log marginal likelihood for node 1 (the 0-1 response variable) in a DAG,
#' given the data X. This function uses a different formula tailored for the Probit
#' model's response variable.
#'
#' @param dag A \eqn{q,q} adjacency matrix representing the DAG structure.
#' @param X An \eqn{n,q} data matrix of observations.
#' @param n number of observations (rows) in the data matrix \code{X}.
#' @param a A hyperparameter of the prior distribution.
#' @param g A hyperparameter of the prior distribution.
#'
#' @return The log of the marginal likelihood for node 1.
#'
#' @export
m_1 = function(dag, X, n, a, g){

  ## This function computes the marginal likelihood of a dag given the data X relative to node 1 (the 0-1 response variable)

  # n is the sample size
  # a,g are hyperparameters of the prior distribution

  pa_1 = BSL::pa(1, DAG=dag)

  y  = X[,1]
  XX = as.matrix(X[,pa_1])

  if(length(pa_1) == 0){

    m = -0.5*n*log(2*pi) -0.5*sum(y^2)

  } else{

    p_1 = length(pa_1)

    DD = t(XX)%*%XX

    V = DD + diag(g, p_1)

    b_hat = solve(V)%*%t(XX)%*%y
    e_hat = y - XX%*%b_hat

    m = - 0.5*n*log(2*pi) - 0.5*sum(y^2) +
      0.5*p_1*log(g) - 0.5*log(det(V)) + 0.5*t(b_hat)%*%V%*%b_hat


  }

  return(m)

}
