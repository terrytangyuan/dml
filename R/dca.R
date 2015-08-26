#' Discriminative Component Analysis
#'
#' Performs discriminative component analysis on the given data.
#'
#' Put DCA function details here.
#'
#' @param data \code{n * d} data matrix. \code{n} is the number of data points,
#'             \code{d} is the dimension of the data.
#'             Each data point is a row in the matrix.
#' @param chunks length \code{n} vector describing the chunklets:
#'               \code{-1} in the \code{i} th place means point \code{i}
#'               doesn't belong to any chunklet;
#'               integer \code{j} in place \code{i} means point \code{i}
#'               belongs to chunklet j.
#'               The chunklets indexes should be 1:(number of chunklets).
#' @param neglinks \code{s * s} symmetric matrix describing the negative relationship
#'                 between all the \code{s} chunklets.
#'                 For the element \eqn{neglinks_{ij}}:
#'                 \eqn{neglinks_{ij} = 1} means chunklet \code{i} and chunklet {j}
#'                 have negative constraint(s);
#'                 \eqn{neglinks_{ij} = 0} means chunklet \code{i} and chunklet {j}
#'                 don't have negative constraints
#'                 or we don't have information about that.
#' @param useD Integer. Optional. When not given, DCA is done in the
#'             original dimension and B is full rank. When useD is given,
#'             DCA is preceded by constraints based LDA which reduces the
#'             dimension to useD. B in this case is of rank useD.
#'
#' @return list of the DCA results:
#' \item{B}{DCA suggested Mahalanobis matrix}
#' \item{DCA}{DCA suggested transformation of the data.
#'            The dimension is (original data dimension) * (useD)}
#' \item{newData}{DCA transformed data}
#'
#' For every two original data points (x1, x2) in newData (y1, y2):
#'
#' \eqn{(x2 - x1)' * B * (x2 - x1) = || (x2 - x1) * A ||^2 = || y2 - y1 ||^2}
#'
#' @keywords dca discriminant component transformation mahalanobis metric
#'
#' @aliases dca
#'
#' @note Put some note here.
#'
#' @author Xiao Nan <\url{http://www.road2stat.com}>
#'
#' @export dca
#'
#' @references
#' Steven C.H. Hoi, W. Liu, M.R. Lyu and W.Y. Ma (2006).
#' Learning Distance Metrics with Contextual Constraints for Image Retrieval.
#' \emph{Proceedings IEEE Conference on Computer Vision and Pattern Recognition
#' (CVPR2006)}.
#'
#' @examples
#' \donotrun{
#' set.seed(123)
#' require(MASS)  # generate synthetic Gaussian data
#' k = 100        # sample size of each class
#' n = 3          # specify how many class
#' N = k * n      # total sample number
#' x1 = mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
#' x2 = mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
#' x3 = mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
#' data = as.data.frame(rbind(x1, x2, x3))

#' # The fully labeled data set with 3 classes
#' plot(data$V1, data$V2, bg = c("#E41A1C", "#377EB8", "#4DAF4A")[gl(n, k)],
#'      pch = c(rep(22, k), rep(21, k), rep(25, k)))
#' Sys.sleep(3)

#' # Same data unlabeled; clearly the classes' structure is less evident
#' plot(x$V1, x$V2)
#' Sys.sleep(3)
#'
#' chunk1 = sample(1:100, 5)
#' chunk2 = sample(setdiff(1:100, chunk1), 5)
#' chunk3 = sample(101:200, 5)
#' chunk4 = sample(setdiff(101:200, chunk3), 5)
#' chunk5 = sample(201:300, 5)
#' chks = list(chunk1, chunk2, chunk3, chunk4, chunk5)

#' chunks = rep(-1, 300)

#' # positive samples in the chunks
#' for (i in 1:5) {
#'   for (j in chks[[i]]) {
#'     chunks[j] = i
#'   }
#' }
#'
#' # define the negative constrains between chunks
#' neglinks = matrix(c(
#' 		0, 0, 1, 1, 1,
#' 		0, 0, 1, 1, 1,
#' 		1, 1, 0, 0, 0,
#' 		1, 1, 0, 0, 1,
#' 		1, 1, 1, 1, 0),
#' 		ncol = 5, byrow = TRUE)
#'
#' dcaData = dca(data = data, chunks = chunks, neglinks = neglinks)$newData
#' # plot DCA transformed data
#' plot(dcaData[, 1], dcaData[, 2], bg = c("#E41A1C", "#377EB8", "#4DAF4A")[gl(n, k)],
#'      pch = c(rep(22, k), rep(21, k), rep(25, k)),
#'      xlim = c(-15, 15), ylim = c(-15, 15))
#' }
#'
dca <- function(data, chunks, neglinks, useD = NULL) {

	data     = t(as.matrix(data))
	chunks   = as.matrix(chunks)
	neglinks = as.matrix(neglinks)

	d = nrow(data)
	n = ncol(data)

	if(is.null(useD)) useD = d

	# 1. Compute means of chunks
	s = max(chunks)
	M = matrix(NA, ncol = s, nrow = d)

	for (i in 1:s) {
		inds = which(chunks == i)
		M[ , i] = as.matrix(rowMeans(data[ , inds]))
	}

	# 2. Compute Cb
	Cb = mat.or.vec(d, d)
	N_d = 0
	for (j in 1:s) {
		inds = which(neglinks[j, ] == 1)
		for (i in 1:length(inds)) {
			Cb = Cb + ((M[ , j] - M[ , inds[i]]) %*% t(M[ , j] - M[ , inds[i]]))
		}
		N_d = N_d + length(inds)
	}

	if (N_d == 0) {
		Cb = diag(d)
	} else {
		Cb = Cb/N_d
	}

	# 3. Compute Cw

	Cw = mat.or.vec(d, d)
	N_w = 0

	for (j in 1:s) {
		inds = which(chunks == j)
		for (i in 1:length(inds)) {
			Cw = Cw + ((data[ , inds[i]] - M[ , j]) %*% t(data[ , inds[i]] - M[ , j]))
		}
		N_w = N_w + length(inds)
	}

	Cw = Cw/N_w

	# 3. Diagonalize Cb

	eigTmp = eigen(Cb)
	eigVec = eigTmp$vectors
	eigVal = as.matrix(eigTmp$values)
	index = which(abs(eigVal) > 1e-9)  # find Non-Zero EigVals

	if (useD < d) {
	R = eigVec[ , index[1:useD]]  # R have already sorted eigenvalues
	} else {
	R = eigVec[ , index]  # Each col of D is an eigenvector
	}

	Db = t(R) %*% Cb %*% as.matrix(R)
	Z = as.matrix(R) %*% ((Db) %^% (-0.5))

	# Diagonalize t(Z) %*% Cw %*% Z
	Cz = t(Z) %*% Cw %*% as.matrix(Z)
	eigVal = eigen(Cz)$values
	if (length(eigVal) == 1) {
		Dw = as.matrix(eigVal)
	} else {
		Dw = diag(eigen(Cz)$values)
	}
	eigVec = eigen(Cz)$vectors

	DCA = (Dw %^% (-0.5)) %*% t(eigVec) %*% t(Z)
	B = t(DCA) %*% as.matrix(DCA)
	newData = t(DCA %*% data)

	return(list("B" = B, "DCA" = DCA, "newData" = newData))
}
