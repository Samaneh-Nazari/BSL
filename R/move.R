#' Perform local moves given a DAG
#'
#' This function locally proposes a DAG by inserting (\code{id}), deleting (\code{dd}) or reversing (\code{rd}) an edge between two \code{nodes} on a given DAG's adjacency matrix.
#' It randomly selects a possible move from a predefined set of actions, applies it, and checks if the resulting graph is still a DAG.
#' The process is repeated until a valid move that results in a DAG is found.
#'
#' @param A A \eqn{q,q} adjacency matrix of the current \code{DAG}.
#' @param q The number of \code{nodes} in the \code{DAG}.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{A_new}: The adjacency matrix after the valid move has been applied.
#'   \item \code{type.operator}: An integer \eqn{{1,2,3}} indicating the type of move performed (1: insert, 2: delete, 3: reverse).
#'   \item \code{nodes}: A numeric vector of \eqn{(2,1)}, indicating the pair of nodes involved in the move.
#' }
#'
#' @details
#' The function defines three types of moves internally:
#' \describe{
#'   \item{1 (Insertion)}{Adds a directed edge from node x to node y.}
#'   \item{2 (Deletion)}{Deletes the existing directed edge from node x to node y.}
#'   \item{3 (Reversal)}{Deletes the edge from node x to node y and adds a new edge from node y to node x.}
#' }
#' This function relies on helper functions \code{id()}, \code{dd()}, and \code{rd()}, as well as a function called \code{is.DAG()} to validate the graph structure.
#'
#'
#' @export
move = function(A, q = q){

  # Define the action names locally to avoid global variables
  actions = c("id", "dd", "rd")

  # A: adjacency matrix of the DAG
  # q: number of vertices

  # Output:
  # A direct successor DAG
  # the (estimated) number of direct successors

  A_na = A_na_na = A

  A_na[1,] = NA

  A_na_na[1,] = A_na_na[,1] = NA

  diag(A_na) = diag(A_na_na) = NA

  id_set = c()
  dd_set = c()
  rd_set = c()

  # set of nodes for id

  set_id = which(A_na == 0, TRUE)

  if(length(set_id) != 0){
    id_set = cbind(1, set_id)
  }

  # set of nodes for dd

  set_dd = which(A == 1, TRUE)

  if(length(set_dd != 0)){
    dd_set = cbind(2, set_dd)
  }

  # set of nodes for rd

  set_rd = which(A_na_na == 1, TRUE)

  if(length(set_rd != 0)){
    rd_set = cbind(3, set_rd)
  }

  O = rbind(id_set, dd_set, rd_set)

  repeat {

    i = sample(dim(O)[1],1)

    act_to_exe  = paste0(actions[O[i,1]],"(A=A,c(",as.vector(O[i,2]),",",as.vector(O[i,3]),"))")
    A_succ      = eval(parse(text = act_to_exe))
    act_to_eval = paste0("is.DAG(A_succ)")
    val = eval(parse(text = act_to_eval))

    if (val != 0){
      break
    }
  }

  A_new = A_succ

  return(list(A_new = A_new, type.operator = O[i,1], nodes = O[i,2:3]))

}

#' Insert a Directed Edge
#'
#' Inserts a directed edge from node x to node y in the adjacency matrix.
#'
#' @param A A \eqn{q,q} adjacency matrix.
#' @param nodes A numeric vector of \eqn{(2,1)}, \code{c(x, y)}, representing the nodes.
#'
#' @return The modified adjacency matrix.
#' @keywords internal
id = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 1
  return(A)
}

#' Delete a Directed Edge
#'
#' Deletes the directed edge from node x to node y in the adjacency matrix.
#'
#' @param A A \eqn{q,q} adjacency matrix.
#' @param nodes A numeric vector of \eqn{(2,1)}, \code{c(x, y)}, representing the nodes.
#'
#' @return The modified adjacency matrix.
#' @keywords internal
dd = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  return(A)
}

#' Reverse a Directed Edge
#'
#' Reverses the direction of an edge from x -> y to y -> x in the adjacency matrix.
#'
#' @param A A \eqn{q,q} adjacency matrix.
#' @param nodes A numeric vector of \eqn{(2,1)}, \code{c(x, y)}, representing the nodes.
#'
#' @return The modified adjacency matrix.
#' @keywords internal
rd = function(A, nodes){ # reverse D x -> y
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  A[y,x] = 1
  return(A)
}
