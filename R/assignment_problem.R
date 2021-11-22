

#' Compute the cost of a given labeling
#'
#' @param cost   Square cost matrix
#' @param label  Integer labeling. Should be some permutation of
#'               \code{[1, ..., nrow(cost)]}.
#'
#' @return Cost associated with the labeling (scalar)
#'
#' @examples
#' C <- simulate.difficult.cost(10)  # Cost matrix
#' f <- ctc(C)
#' assignment.cost(C, f) == attr(f, "cost")  # [1] TRUE
#'
#' @seealso \code{\link{ctc}}, \code{\link{simulate.difficult.cost}}
#'
assignment.cost <- function(cost, label) {
  if ( length(label) != nrow(cost) )
    stop("Length of label must match dimension of cost matrix")
  if ( nrow(cost) != ncol(cost) )
    stop("Cost matrix should be square")
  label <- as.integer(label)
  r <- range(label)
  if ( r[1L] < 1L || r[2L] > nrow(cost) )
    stop("Label range does not match dimension of cost matrix")
  x <- 0
  for ( j in 1:ncol(cost) ) {
    x <- x + cost[label[j], j]
  }
  x
}


#' Compute cost matrix for label pairs
#'
#' Input \code{a} taken as the reference lables.
#' Implementation uses \code{base::table}.
#'
#' @param a  Reference labels
#' @param b  (Partially) Paired labels. Inputs \code{a} and \code{b}
#'           must be of equal length
#' @param .unpaired  Cost associated with assigning a
#'           non-overlapping pair of labels
#' @param .names  Include table row/column names
#'
#' @return Square cost matrix where the columns are associated with
#'   the reference labels, \code{a}, and the rows are associated with
#'   the compared labels, \code{b}
#'
#' @examples
#' x <- sample.int(10, 100, replace = TRUE)  # Reference labels
#' map <- sample(1:10)
#' y <- map[x]                               # Renamed labels
#' abcost(x, y)
#'
#' @seealso \code{\link[base]{table}}
#'
abcost <- function(
  a, b,
  .unpaired = 2 * length(a),
  .names = TRUE
) {
  if ( length(a) != length(b) ) {
    stop("length(", deparse(substitute(a)), ") != length(",
         deparse(substitute(b)), ")")
  }
  .unpaired <- as.integer(.unpaired)
  ua <- sort(unique(a))
  ub <- sort(unique(b))
  n <- max(length(ua), length(ub))
  cost <- matrix(.unpaired, n, n)
  ct <- -as.matrix(table(b, a))
  ct[ct == 0] <- .unpaired
  cost[1:nrow(ct), 1:ncol(ct)] <- ct
  if ( .names ) {
    colnames(cost)[1:length(ua)] <- as.character(ua)
    rownames(cost)[1:length(ub)] <- as.character(ub)
  }
  cost
}



#' Solve the Assignment Problem
#'
#' Given a cost function \eqn{C(a, b)}, \eqn{a \in A} and
#' \eqn{b \in B}, the "Assignment Problem" is to find a bijective map
#' \eqn{f : A \to B} that minimizes,
#' \deqn{\sum_{a \in A} C(a, f(a)).}{\sum_a C(a, f(a)).}
#' In the solution, we assume that A is the index set of the *columns*
#' of the cost matrix \eqn{C}.
#'
#' @param cost  Square cost matrix
#' @param .tol  Tolerance for equivalence of non-integer costs.
#'              If not provided, costs are coerced to integer values
#'
#' @return The map \eqn{f(A)}
#'
#' @examples
#' C <- simulate.difficult.cost(10)  # Cost matrix
#' ctc(C)
#'
#' @seealso \code{\link{abcost}}, \code{\link{simulate.difficult.cost}}
#'
#' @references
#' Carpaneto, Giorgio, and Paolo Toth.
#' "Primal-dual algrorithms for the assignment problem."
#' Discrete Applied Mathematics 18.2 (1987): 137-153.
#' (\href{https://www.sciencedirect.com/science/article/pii/0166218X87900163}{sciencedirect.com})
#'
ctc <- function(cost, .tol = NA) {
  if ( nrow(cost) != ncol(cost) )
    stop("Cost matrix should be square")
  if ( is.na(.tol) || is.null(.tol) ) {
    cost <- matrix(as.integer(cost), nrow(cost))
    x <- .Call("ctc_integer_", cost, PACKAGE = "assignment") + 1L
  } else {
    .tol <- abs(.tol[1L])
    x <- .Call("ctc_double_", cost, .tol, PACKAGE = "assignment") + 1L
  }
  structure(x, cost = assignment.cost(cost, x))
}




#' Simulate difficult to solve cost matrix
#'
#' @usage
#' simulate.difficult.cost(n)
#'
#' @param n  Dimension of simulated cost matrix
#' @return \code{n x n} cost matrix
#'
simulate.difficult.cost <- function(n) {
  n <- as.integer(n)
  apply(outer(1:n, 1:n), 2, function(x) sapply(x, sample.int, 1))
}


