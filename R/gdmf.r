#' Global Distance Metric Learning
#'
#' Performs Global Distance Metric Learning (GDM) on the given data, learning a full matrix.
#'
#' Put GdmFull function details here.
#'
#' @param data \code{n * d} data matrix. \code{n} is the number of data points,
#'             \code{d} is the dimension of the data.
#'             Each data point is a row in the matrix.
#' @param simi \code{n * 2} matrix describing the similar constrains.
#'              Each row of matrix is serial number of a similar pair in the original data.
#'				For example, pair(1, 3) represents the first observation is similar the 3th observation in the original data.
#' @param dism \code{n * 2} matrix describing the dissimilar constrains as \code{simi}.
#'				Each row of matrix is serial number of a dissimilar pair in the original data.
#' @param maxiter numeric, the number of iteration.
#'
#' @return list of the GdmDiag results:
#' \item{newData}{GdmDiag transformed data}
#' \item{fullA}{suggested Mahalanobis matrix}
#' \item{dmlA}{matrix to transform data, square root of diagonalA }
#' \item{converged}{whether the iteration-projection optimization is converged or not}
#'
#' For every two original data points (x1, x2) in newData (y1, y2):
#'
#' \eqn{(x2 - x1)' * A * (x2 - x1) = || (x2 - x1) * B ||^2 = || y2 - y1 ||^2}
#'
#' @keywords GDM global distance metirc learning transformation mahalanobis metric
#'
#' @note Be sure to check whether the dimension of original data and constrains' format are valid for the function.
#'
#' @author Gao Tao <\url{http://www.gaotao.name}>
#'
#' @references
#' Steven C.H. Hoi, W. Liu, M.R. Lyu and W.Y. Ma (2003).
#' Distance metric learning, with application to clustering with side-information.
#  in \emph{Proc. NIPS}.
#'
#' @examples
#' \donotrun{
#' set.seed(123)
#' library(MASS)
#' library(scatterplot3d)
#'
#' # generate simulated Gaussian data
#' k = 100
#' m <- matrix(c(1, 0.5, 1, 0.5, 2, -1, 1, -1, 3), nrow =3, byrow = T)
#' x1 <- mvrnorm(k, mu = c(1, 1, 1), Sigma = m)
#' x2 <- mvrnorm(k, mu = c(-1, 0, 0), Sigma = m)
#' data <- rbind(x1, x2)
#'
#' # define similar constrains
#' simi <- rbind(t(combn(1:k, 2)), t(combn((k+1):(2*k), 2)))
#'
#' temp <-  as.data.frame(t(simi))
#' tol <- as.data.frame(combn(1:(2*k), 2))
#'
#' # define disimilar constrains
#' dism <- t(as.matrix(tol[!tol %in% simi]))
#'
#' # transform data using GdmFull
#' result <- GdmFull(data, simi, dism)
#' newData <- result$newData
#' # plot original data
#' color <- gl(2, k, labels = c("red", "blue"))
#' par(mfrow = c(2, 1), mar = rep(0, 4) + 0.1)
#' scatterplot3d(data, color = color, cex.symbols = 0.6,
#'			  xlim = range(data[, 1], newData[, 1]),
#'			  ylim = range(data[, 2], newData[, 2]),
#'			  zlim = range(data[, 3], newData[, 3]),
#'			  main = "Original Data")
#' # plot GdmFull transformed data
#' scatterplot3d(newData, color = color, cex.symbols = 0.6,
#'			  xlim = range(data[, 1], newData[, 1]),
#'			  ylim = range(data[, 2], newData[, 2]),
#'			  zlim = range(data[, 3], newData[, 3]),
#'			  main = "Transformed Data")
#' }
#'
GdmFull <- function(data, simi, dism, maxiter = 100) {
		data <- as.matrix(data)
		N <- dim(data)[1]
		d <- dim(data)[2]
		new.simi <- unique(t(apply(simi, 1, sort)))
		new.dism <- unique(t(apply(dism, 1, sort)))

		A <- diag(1, d) * 0.1
		W <- mat.or.vec(d, d)
		dij <- mat.or.vec(1, d)

		# sphereMult = cov(data)^(-0.5);
		# spheredata = data %*% sphereMult

		dist1.simi <- data[new.simi[, 1], ] - data[new.simi[, 2], ]
		dist2.ij <- t(apply(dist1.simi, 1, function(x) outer(x, x)))
		W <- matrix(apply(dist2.ij, 2, sum), ncol = d, byrow = TRUE)


		w <- matrix(W, ncol = 1)
		t0 <- as.numeric(crossprod(w, matrix(A, ncol = 1))/100)

		IterProjection <- function(data, simi, dism, A, w, t0 , maxiter = 100) {
							data <- as.matrix(data)
							N = dim(data)[1]     # number of examples
							d = dim(data)[2]     # dimensionality of examples
							# S1 <- mat.or.vec(N, N)
							# D1 <- mat.or.vec(N, N)
							# simi <- rbind(simi, simi[, c(2, 1)])
							# dism <- rbind(dism, dism[, c(2, 1)])
							# S1[simi] <- 1
							# D1[dism] <- 1
							new.simi <- unique(t(apply(simi, 1, sort)))
							new.dism <- unique(t(apply(dism, 1, sort)))

							# error1=1e5
							threshold2 <- 0.01  # error-bound of main A-update iteration
							epsilon <- 0.01   # error-bound of iterative projection on C1 and C2
							maxcount <- 200

							w1 <- w/norm(w, "F")    # make 'w' a unit vector
							t1 <- t0/norm(w, "F")

							count <- 1
							alpha <- 0.1    # initial step size along gradient

							GradProjection <- function(grad1, grad2, d) {
												g1 <- matrix(grad1, ncol = 1)
												g2 <- matrix(grad2, ncol = 1)

												g2 <- g2/norm(g2, "F")
												gtemp <- g1 - as.numeric(crossprod(g2, g1)) * g2
												gtemp <- gtemp/norm(gtemp, "F")
												grad.proj <- matrix(gtemp, d, d)
												return(grad.proj)
							}

							fS1 <- function(data, new.simi, A, N, d, fudge = 0.000001) {

									dist1.simi <- data[new.simi[, 1], ] - data[new.simi[, 2], ]
									dist2.ij <- t(apply(dist1.simi, 1, function(x) outer(x, x)))
									fs.1d <- matrix(apply(dist2.ij, 2, sum), ncol = d, byrow = TRUE)
									return(fs.1d)
							}

							fD1 <- function(data, new.simi, A, N, d, fudge = 0.000001) {
									dist1.dism <- data[new.dism[, 1], ] - data[new.dism[, 2], ]
									dist.ij <- numeric(dim(dist1.dism)[1])
									for (i in 1:dim(dist1.dism)[1]) {
										dist.ij[i] <- sqrt(t(dist1.dism[i, ]) %*% A %*% t(t(dist1.dism[i, ])))
									}
									sum.dist <- sum(dist.ij) + 0.000001
									Mij <- t(apply(dist1.dism, 1, function(x) outer(x, x)))
									temp <- cbind(Mij, t(t(dist.ij)))

									deri.ij <- 0.5 * temp[, 1:(d^2)]/(temp[, d^2 + 1] + (temp[, d^2 + 1] == 0) * fudge)
									sum.deri <- matrix(apply(deri.ij, 2, sum), ncol = d, byrow = TRUE)
									fd.1d <- sum.deri/sum.dist
									return(fd.1d)
							}

							fD <- function(data, new.dism, A, N, d) {
									dist1.dism <- data[new.dism[, 1], ] - data[new.dism[, 2], ]
									dist.ij <- numeric(dim(dist1.dism)[1])
									for (i in 1:dim(dist1.dism)[1]) {
										dist.ij[i] <- sqrt(t(dist1.dism[i, ]) %*% A %*% t(t(dist1.dism[i, ])))
									}
									fd <- sum(dist.ij) + 0.000001
									fd <- log(fd)
									return(fd)
							}

							grad1 <- fS1(data, new.simi, A, N, d);   # gradient of similarity constraint function
							grad2 <- fD1(data, new.dism, A, N, d);   # gradient of dissimilarity constraint func.
							M <- GradProjection(grad1, grad2, d); # gradient of fD1 orthognal to fS1


							A.last <- A        # initial A
							done <- 0
							delta <- 0
							converged <- 0
							while (done == 0) {
									projection.iters <- 0
									satisfy <- 0

									while (projection.iters < maxiter & satisfy == 0) {
										A0 <- A
										x0 <- matrix(A0, ncol = 1)
										if(crossprod(w, x0) <= t0)
											A <- A0
										else {
											x <- x0 + as.numeric(t1 - crossprod(w1, x0)) * w1
											A <- matrix(x, 3, 3)
										}

										A <- (A + t(A))/2
										vl <- eigen(A)
										vl[[1]][vl[[1]] < 0] = 0
										A <- vl[[2]] %*% diag(vl[[1]], d) %*% t(vl[[2]])

										fDC2 <- crossprod(w, matrix(A, ncol = 1))
										error1 <- as.numeric((fDC2 - t0)/t0)
										projection.iters <- projection.iters + 1
										satisfy <- as.numeric(ifelse(error1 > epsilon, 0, 1))
									}

									obj.previous <- fD(data, new.dism, A.last, N, d)
									obj <- fD(data, new.dism, A, N, d)

									if (obj > obj.previous & satisfy == 1) {
										alpha <-  alpha * 1.05
										A.last <- A
										grad2 <- fS1(data, new.simi, A, N, d)
										grad1 <- fD1(data, new.dism, A, N, d)
										M <- GradProjection(grad1, grad2, d)
										A <- A + alpha * M
									}
									else{
										alpha <- alpha/2
										A <- A.last + alpha * M
									}
									delta <- norm(alpha * M, "F")/norm(A.last, "F")
									count <- count + 1
									done <- ifelse(delta < threshold2 | count == maxcount, 1, 0)
							}
							converged <- ifelse(delta > threshold2, 0, 1)
							return(list("converged" = ifelse(converged == 1, "Yes", "No"), "fullA" = A))
		}

		iterproj <- IterProjection(data, simi, dism, A, w, t0)
		eigenvalue <- eigen(iterproj$fullA)
		dml <- eigenvalue[[2]] %*% sqrt(diag(eigenvalue[[1]], d))
		newData <- data %*% dml
		return(list("newData" = newData, "fullA" = iterproj[[2]], "dmlA" = dml, "converged" = iterproj[[1]]))
}
