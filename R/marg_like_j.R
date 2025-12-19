#' Compute log marginal likelihood for a Node in a DAG
#'
#' A wrapper function that calculates the log marginal likelihood for a specified node
#' in a DAG. It automatically determines whether to use the formula for the general
#' node or the special formula for the response node (node 1).
#'
#' @param j An integer representing the label of the \code{j} for which the likelihood is calculated.
#' @param dag A \eqn{q,q} adjacency matrix representing the DAG structure.
#' @param X An \eqn{n,q} data matrix of observations.
#' @param a A hyperparameter of the prior distribution.
#' @param g A hyperparameter of the prior distribution.
#' @param n number of observations (rows) in the data matrix \eqn{X}.
#'
#' @return A single numeric value: the log of the marginal likelihood for the specified node \code{j}.
#'
#' @export
marg_like_j = function(j, dag, X, a, g, n){

  ## General function to compute the marginal likelihood relative to a node j in the dag

  # n is the sample size
  # a,g are hyperparameters of the prior distribution

  if(j != 1){

    out_m = m_j(j = j, dag = dag, X = X, n = n, a = a, g = g)

  }else{

    out_m = m_1(dag = dag, X = X, n = n, a = a, g = g)

  }

  return(out_m)

}
