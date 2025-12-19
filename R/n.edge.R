#' Count the Number of Edges in a DAG
#'
#' This function counts the number of directed edges in a DAG represented by an adjacency matrix.
#'
#' @param A A \eqn{q,q} adjacency matrix of the DAG.
#'
#' @return An integer representing the total number of edges.
#' @keywords internal
#'
n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}
