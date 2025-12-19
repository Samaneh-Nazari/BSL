#' Sample from the Posterior Distribution of the Covariance Matrix
#'
#' This function generates a single sample from the posterior distribution of the
#' covariance matrix Sigma, given a data matrix \code{Y} and a \code{DAG} structure A. The sampling
#' is performed by first sampling from the posterior of the Cholesky decomposition
#' parameters (\code{L} and \code{D}).
#'
#' @param Y An \eqn{n,q} data matrix of observations.
#' @param A A \eqn{q,q} adjacency matrix representing the \code{DAG} structure.
#' @param g A hyperparameter of the prior distribution.
#' @param a A hyperparameter of the prior distribution.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{Sigma_post}: A \eqn{q,q} matrix, which is a single sample from the posterior distribution of the covariance matrix Sigma.
#'   \item \code{L_post}: A \eqn{q,q} lower-triangular matrix, the Cholesky factor \code{L} from the reparameterization.
#'   \item \code{D_post}: A \eqn{q,q} diagonal matrix, the diagonal matrix \code{D} from the reparameterization.
#' }
#'
#' @export
#'
#'
posterior_sigma = function(Y, A, g, a){

  # This is a function to sample from the posterior of (D,L), a reparameterization of the covariance matrix Sigma

  ###########
  ## INPUT ##
  ###########

  # Y :   the (n,q) data matrix
  # A :   a DAG, represented by its (q,q) adjacency matrix
  # g,a : hyperparameters of the DAG prior distribution

  ############
  ## OUTPUT ##
  ############

  # Sigma_post : (q,q) matrix representing a sample from the posterior of Sigma
  # L_post, D_post : (q,q) matrices, reparameterization of Sigma_post

  n = nrow(Y)
  q = ncol(Y)
  a = n-1
  g = 2


  # Next function samples from the posterior of node-parameters L_j, D_jj

  post_DL_j = function(A = A, j, Y = Y, n = n, q = q, g = g, a = a){

    # A :   (q,q) adjacency matrix of the DAG
    # j :   a node (j = 1, ..., q)
    # Y :   (n,q) data matrix
    # n :   number of observations (rows of Y)
    # q :   number of variables (columns of Y)
    # g,a : hyperparameters of the DAG prior distribution

    pa_j  = BSL::pa(j, A)

    XX = Y[,pa_j]
    y  = Y[,j]

    p_j  = length(pa_j)
    j_pos = q - length(pa_j)
    a=n-1
    a_j  = (a + q - 2*j_pos + 3)/2 - p_j/2 - 1

    L_j = matrix(0, 1, ncol = p_j)

    if(length(pa_j) == 0){
      D_jj =(rgamma(1,n/2-a_j,sum(y^2)/2))^(-1)*rgamma(1,1,0.5*g)

    } else{

      b_hat = solve(t(XX)%*%XX + diag(g,p_j))%*%t(XX)%*%y

      d = t(b_hat)%*%(diag(g, p_j) + t(XX)%*%XX)%*%b_hat

      D_jj =(rgamma(1,n/2-a_j,sum(y^2)/2 - d/2))^(-1)*rgamma(1,1,0.5*g)

      L_j = rmvnorm(1, -b_hat, D_jj*solve(t(XX)%*%XX + diag(g, p_j)))

    }

    return(list(L_j_post = L_j, D_jj_post = D_jj))

  }

  # For node 1 (fix conditional variance to 1 and sample L_1)

  pa_1 = BSL::pa(1, A)
  p_1 = length(pa_1)

  XX = Y[,pa_1]
  y  = Y[,1]

  L_1 = matrix(0, 1, ncol = p_1)

  if(length(pa_1) == 0){

    D_11 = 1

  } else{

    D_11 = 1

    b_hat = solve(t(XX)%*%XX + diag(g,p_1))%*%t(XX)%*%y

    L_1 = rmvnorm(1, -b_hat, solve(t(XX)%*%XX + diag(g, p_1)))

  }

  # For the other nodes (2,...,q) sample from the Cholesky parameters

  DL_post_js = lapply(X = 1:q, FUN = post_DL_j, A = A, Y = Y, n = n, q = q, g = g, a = a)

  DL_post_js[[1]] = list(L_j_post = L_1, D_jj_post = D_11)


  j_parents = sapply(X = 1:q, FUN = BSL::pa, DAG = A)

  Omega  = list()
  Sigma  = list()
  Beta   = list()
  D_post = list()

  L = matrix(0, q, q)
  D = matrix(0, q, q)

  for(j in 1:q){

    L[unlist(j_parents[j]),j] = DL_post_js[[j]][[1]]
    D[j,j] = DL_post_js[[j]][[2]]

  }

  diag(L) = 1

  Sigma  = t(solve(L))%*%D%*%solve(L)


  return(list(Sigma_post = Sigma, L_post = L, D_post = D))

}
