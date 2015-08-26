# repeat a matrix like the MATLAB grammar
repmat <- function(A, N, M) {
	kronecker(matrix(1, N, M), A)
}

# negative one half matrix power operator
"%^%" <- function(x, n) {
		with(eigen(as.matrix(x)), vectors %*% (values^n * t(vectors)))
	}

