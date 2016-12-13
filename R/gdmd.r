#' Global Distance Metric Learning
#'
#' Performs Global Distance Metric Learning (GDM) on the given data, learning a diagonal matrix.
#'
#' Put GdmDiag function details here.
#'
#' @param data \code{n * d} data matrix. \code{n} is the number of data points,
#'             \code{d} is the dimension of the data.
#'             Each data point is a row in the matrix.
#' @param simi \code{n * 2} matrix describing the similar constrains.
#'              Each row of matrix is serial number of a similar pair in the original data.
#'				For example, pair(1, 3) represents the first observation is similar the 3th observation in the original data.
#' @param dism \code{n * 2} matrix describing the dissimilar constrains as \code{simi}.
#'				Each row of matrix is serial number of a dissimilar pair in the original data.
#' @param C0 numeric, the bound of similar constrains.
#' @param threshold numeric, the threshold of stoping the learning iteration.
#'
#' @return list of the GdmDiag results:
#' \item{newData}{GdmDiag transformed data}
#' \item{diagonalA}{suggested Mahalanobis matrix}
#' \item{dmlA}{matrix to transform data, square root of diagonalA }
#' \item{error}{the precision of obtained distance metric by Newton-Raphson optimization }
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
#' @export GdmDiag
#' @import MASS
#'
#' @references
#' Steven C.H. Hoi, W. Liu, M.R. Lyu and W.Y. Ma (2003).
#' Distance metric learning, with application to clustering with side-information.
#  in \emph{Proc. NIPS}.
#'

#' @examples
#' \dontrun{
#' set.seed(602)
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
#' dism <- t(as.matrix(tol[!tol %in% temp]))
#'
#' # transform data using GdmDiag
#' result <- GdmDiag(data, simi, dism)
#' newData <- result$newData
#' # plot original data
#' color <- gl(2, k, labels = c("red", "blue"))
#' par(mfrow = c(2, 1), mar = rep(0, 4) + 0.1)
#' scatterplot3d(data, color = color, cex.symbols = 0.6,
#'			  xlim = range(data[, 1], newData[, 1]),
#'			  ylim = range(data[, 2], newData[, 2]),
#'			  zlim = range(data[, 3], newData[, 3]),
#'			  main = "Original Data")
#' # plot GdmDiag transformed data
#' scatterplot3d(newData, color = color, cex.symbols = 0.6,
#'			  xlim = range(data[, 1], newData[, 1]),
#'			  ylim = range(data[, 2], newData[, 2]),
#'			  zlim = range(data[, 3], newData[, 3]),
#'			  main = "Transformed Data")
#' }
#'
GdmDiag <- function(data, simi, dism, C0 = 1, threshold = 0.001) {
		fudge = 0.000001
		reduction = 2
		data <- as.matrix(data)
		
		# Check that simi and dism are k*2 matrices
		if (dim(simi)[2] != 2) {
		  stop(paste('simi needs to be of dimensions k*2 but has dimensions:', paste(dim(simi), collapse = ", ")))
		}
		if (dim(dism)[2] != 2) {
		  stop(paste('dism needs to be of dimensions k*2 but has dimensions:', paste(dim(dism), collapse = ", ")))
		}
		
		simi <- as.matrix(simi)
		dism <- as.matrix(dism)
		N <- dim(data)[1]
		d <- dim(data)[2]
		a <- matrix(rep(1, d), nrow = d)# initial diagonal A in the form of column vector
		# dij <- mat.or.vec(1, d)

		new.simi <- unique(t(apply(simi, 1, sort)))
		new.dism <- unique(t(apply(dism, 1, sort)))
    
		# Check that simi and dism do not overlap
		dup.pairs <- duplicated(rbind(new.simi, new.dism), MARGIN = 1)
		if (any(dup.pairs)) stop(paste('There are',sum(dup.pairs),'overlapping pairs in simi and dism.'))
		
		# Check that all indices in simi and dism are in the range 1:N
		if (any(new.simi < 1) || any(new.simi > N)) {
		  stop(paste('Some indices in simi are out of range of data.'))
		}
		if (any(new.dism < 1) || any(new.dism > N)) {
		  stop(paste('Some indices in dism are out of range of data.'))
		}
		
		######### contraints
		dist1.dism <- data[new.dism[, 1], ] - data[new.dism[, 2], ]
		dist.ij <- sqrt((dist1.dism^2) %*% a)
		sum.dist <- sum(dist.ij)
		temp <- cbind(dist1.dism^2, dist.ij)
		deri1.ij <-0.5 * temp[, 1:d]/(temp[, d + 1] + (temp[, d + 1] == 0) * fudge)
		sum.deri1 <- t(apply(deri1.ij, 2, sum))
		deri2.ij <- t(apply(dist1.dism, 1, function(x) outer(x, x)))
		temp1 <- cbind(deri2.ij, dist.ij^3)
		deri2.ij <- -0.25 * temp1[, 1:(d^2)]/(temp1[, d^2 + 1] + (temp1[, d^2 + 1] == 0) * fudge)
		sum.deri2 <- matrix(apply(deri2.ij, 2, sum), ncol = d, byrow = TRUE)

		fD <- log(sum.dist)
		fD.1d <- sum.deri1/sum.dist
		fD.2d <- sum.deri2/sum.dist - crossprod(sum.deri1, sum.deri1)/(sum.dist^2)

		####### objection is part of contraints
		# fD <- log(sum.dist)
		######################################
		dist1.dism <- data[new.dism[, 1], ] - data[new.dism[, 2], ]
		d.sum <- t(apply(dist1.dism^2, 2, sum))
		dist1.simi <- data[new.simi[, 1], ] - data[new.simi[, 2], ]
		s.sum <- t(apply(dist1.simi^2, 2, sum))

		# S1 <- mat.or.vec(N, N)
		# D1 <- mat.or.vec(N, N)

		# dism <- rbind(dism, dism[, c(2, 1)]) #
		# Dism <- rbind(Dism, Dism[, c(2, 1)])
		# S1[dism] <- 1
		# D1[Dism] <- 1

		error <- 1
		while (error > threshold) {
			obj.initial <- as.numeric(s.sum %*% a) + C0 * fD
			fS.1d <- s.sum

			gradient <- fS.1d - C0 * fD.1d
			hessian <- -C0 * fD.2d + fudge * diag(1, d)
			invhessian <- solve(hessian)
			cstep <- invhessian %*% t(gradient)

			lambda <- 1
			atemp <- a - lambda * cstep
			atemp[atemp < 0] <- 0

			fDo <- log(sum(sqrt((dist1.dism^2) %*% atemp)))
			obj <- as.numeric(s.sum %*% atemp) + C0 * fDo
			obj.previous <- obj * 1.1

			while (obj < obj.previous) {
				lambda.previous <- lambda
				obj.previous <- obj
				a.previous = atemp
				lambda <- lambda/reduction
				atemp <- a - lambda * cstep
				atemp[atemp < 0] <- 0
				fDo1 <- log(sum(sqrt((dist1.dism^2) %*% atemp)))
				obj <- as.numeric(s.sum %*% atemp) + C0 * fDo1
			}
		a <- a.previous
		error <- abs((obj.previous - obj.initial)/obj.previous)
		}
		diagonalA <- diag(as.numeric(a))
		dmlA <- sqrt(diagonalA)
		newData <- data %*% dmlA

		return(list("newData" = newData, "diagonalA" = diagonalA, "dmlA" = dmlA, "error" = error))
}
