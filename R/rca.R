#' Relevant Component Analysis
#'
#' Performs relevant component analysis on the given data.
#'
#' The RCA function takes a data set and a set of positive constraints
#' as arguments and returns a linear transformation of the data space
#' into better representation, alternatively, a Mahalanobis metric
#' over the data space.
#'
#' Relevant component analysis consists of three steps:
#' \enumerate{\item locate the test point
#' \item compute the distances between the test points
#' \item find \eqn{k} shortest distances and the bla}
#' The new representation is known to be optimal in an information
#' theoretic sense under a constraint of keeping equivalent data
#' points close to each other.
#'
#' @param x matrix or data frame of original data.
#'          Each row is a feature vector of a data instance.
#' @param chunks list of \code{k} numerical vectors.
#'               Each vector represents a chunklet, the elements
#'               in the vectors indicate where the samples locate
#'               in \code{x}. See examples for more information.
#' @param useD optional. When not given, RCA is done in the
#' original dimension and \code{B} is full rank. When \code{useD} is given,
#' RCA is preceded by constraints based LDA which reduces
#' the dimension to \code{useD}. \code{B} in this case is of rank \code{useD}.
#' 
#' @return list of the RCA results:
#' \item{B}{The RCA suggested Mahalanobis matrix.
#'          Distances between data points x1, x2 should be
#'          computed by (x2 - x1)' * B * (x2 - x1)}
#' \item{A}{The RCA suggested transformation of the data.
#'          The data should be transformed by A * data}
#' \item{newX}{The data after the RCA transformation (A).
#'             newData = A * data}
#'
#' The three returned argument are just different forms of the same output.
#' If one is interested in a Mahalanobis metric over the original data space,
#' the first argument is all she/he needs. If a transformation into another
#' space (where one can use the Euclidean metric) is preferred, the second
#' returned argument is sufficient. Using A and B is equivalent in the
#' following sense:
#'
#' if y1 = A * x1, y2 = A * y2  then
#' (x2 - x1)' * B * (x2 - x1) = (y2 - y1)' * (y2 - y1)
#'
#' @keywords rca transformation mahalanobis metric
#'
#' @aliases rca
#'
#' @note Note that any different sets of instances (chunklets),
#'       e.g. {1, 3, 7} and {4, 6}, might belong to the
#'       same class and might belong to different classes.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @seealso See \code{\link{dca}} for exploiting negative constrains.
#'
#' @export rca
#' @importFrom lfda %^%
#' @importFrom stats cov
#'
#' @references
#' Aharon Bar-Hillel, Tomer Hertz, Noam Shental, and Daphna Weinshall (2003).
#' Learning Distance Functions using Equivalence Relations.
#' \emph{Proceedings of 20th International Conference on
#' Machine Learning (ICML2003)}.
#'
#' @examples
#' \dontrun{
#' library("MASS") # generate synthetic multivariate normal data
#' set.seed(42)
#' k <- 100L # sample size of each class
#' n <- 3L # specify how many classes
#' N <- k * n # total sample size
#' x1 <- mvrnorm(k, mu = c(-16, 8), matrix(c(15, 1, 2, 10), ncol = 2))
#' x2 <- mvrnorm(k, mu = c(0, 0), matrix(c(15, 1, 2, 10), ncol = 2))
#' x3 <- mvrnorm(k, mu = c(16, -8), matrix(c(15, 1, 2, 10), ncol = 2))
#' x <- as.data.frame(rbind(x1, x2, x3)) # predictors
#' y <- gl(n, k) # response
#'
#' # fully labeled data set with 3 classes
#' # need to use a line in 2D to classify
#' plot(x[, 1L], x[, 2L],
#'   bg = c("#E41A1C", "#377EB8", "#4DAF4A")[y],
#'   pch = rep(c(22, 21, 25), each = k)
#' )
#' abline(a = -10, b = 1, lty = 2)
#' abline(a = 12, b = 1, lty = 2)
#'
#' # generate synthetic chunklets
#' chunks <- vector("list", 300)
#' for (i in 1:100) chunks[[i]] <- sample(1L:100L, 10L)
#' for (i in 101:200) chunks[[i]] <- sample(101L:200L, 10L)
#' for (i in 201:300) chunks[[i]] <- sample(201L:300L, 10L)
#'
#' chks <- x[unlist(chunks), ]
#'
#' # make "chunklet" vector to feed the chunks argument
#' chunksvec <- rep(-1L, nrow(x))
#' for (i in 1L:length(chunks)) {
#'   for (j in 1L:length(chunks[[i]])) {
#'     chunksvec[chunks[[i]][j]] <- i
#'   }
#' }
#'
#' # relevant component analysis
#' rcs <- rca(x, chunksvec)
#'
#' # learned transformation of the data
#' rcs$A
#'
#' # learned Mahalanobis distance metric
#' rcs$B
#'
#' # whitening transformation applied to the chunklets
#' chkTransformed <- as.matrix(chks) %*% rcs$A
#'
#' # original data after applying RCA transformation
#' # easier to classify - using only horizontal lines
#' xnew <- rcs$newX
#' plot(xnew[, 1L], xnew[, 2L],
#'   bg = c("#E41A1C", "#377EB8", "#4DAF4A")[gl(n, k)],
#'   pch = c(rep(22, k), rep(21, k), rep(25, k))
#' )
#' abline(a = -15, b = 0, lty = 2)
#' abline(a = 16, b = 0, lty = 2)
#' }

rca <- function(x, chunks, useD = NULL) {
  n <- nrow(x)
  d <- ncol(x)
  if(is.null(useD)) useD = d
  
  # subtract the mean
  TM <- colMeans(x)
  x <- x - (matrix(1L, n, 1L) %*% TM)
  
  # compute chunklet means and center data
  S <- max(chunks)
  
  Cdata <- matrix()
  AllInds <- vector()
  M <- matrix(NA, nrow = S, ncol = d)
  tmp <- vector('list', S)
  
  for (i in 1L:S) {
    inds <- which(chunks == i)
    M[i, ] <- colMeans(x[inds, ])
    tmp[[i]] <- x[inds, ] - matrix(1L, length(inds), 1L) %*% M[i, ]
    Cdata <- do.call(rbind, tmp)
    AllInds <- c(AllInds, inds)
  }
  
  # Compute inner covariance matrix
  InnerCov <- cov(Cdata) * ((nrow(Cdata) - 1L) / nrow(Cdata))
  
  # Optional cFLD: find optimal projection: min | A S_w A^t | / | A S_t A^t |
  
  if (useD < d) {
    # Compute total covariance using only chunkleted points
    TotalCov <- cov(x[AllInds, ])
    # TotalCov = cov(x)  # Compute total covariance using all the data.
    # More accurate but may lead to PSD problems
    tmp <- eigen(solve(TotalCov) %*% InnerCov)
    D <- tmp$values
    V <- tmp$vectors
    # reorder the vectors in descending order
    V <- V[ , ncol(V):1L]
    # A is the cFLD transformation
    # Acts on row data by multiplication from the right
    A <- V[ , 1L:useD]
    InnerCov <- t(A) %*% InnerCov %*% A
  } else {
    A <- diag(d)
  }
  
  # RCA: whiten the data w.r.t the inner covariance matrix
  tmp <- svd(InnerCov)
  U1 <- tmp$u
  S1 <- diag(tmp$d)
  V1 <- tmp$v
  
  # Operate from the right on row vectors
  A <- as.matrix(A %*% (U1 %*% S1^(0.5)))
  
  # The total mean subtracted is re-added before transformation.
  # This operation is not required for distance computations,
  # as it only adds a constant to all points
  newX <- as.matrix(x + matrix(1L, n, 1L) %*% TM) %*% A
  
  B <- A %*% t(A)

  out <- list("B" = B, "A" = A, "newX" = newX)

  class(out) <- "rca"
  return(out)
}
#' Print an rca object
#'
#' Print an rca object
#' @param x The result from rca function, which contains mahalanobis metric,
#' whitening transformation matrix, and transformed data
#' @param ... ignored
#' @export
#' @importFrom utils head
#' @method print rca
print.rca <- function(x, ...) {
  cat("Results for Relevant Component Analysis \n\n")
  cat("The Mahalanobis metric is: \n")
  print(head(x$B))

  cat("\n\n The whitening transformation matrix is:  \n")
  print(head(x$A))

  cat("\n\n The original data transformed is:  \n")
  print(head(x$newX))

  cat("\n")
  cat("Only partial output is shown above. Please see the model output for more details. \n")
  invisible(x)
}
