context("test algorithms")

library("lfda")
library("MASS")

test_that("dca works", {
  k <- 100 # sample size of each class
  n <- 3 # specify how many class
  N <- k * n # total sample number
  x1 <- mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 <- mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 <- mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  data <- as.data.frame(rbind(x1, x2, x3))

  chunk1 <- sample(1:100, 5)
  chunk2 <- sample(setdiff(1:100, chunk1), 5)
  chunk3 <- sample(101:200, 5)
  chunk4 <- sample(setdiff(101:200, chunk3), 5)
  chunk5 <- sample(201:300, 5)
  chks <- list(chunk1, chunk2, chunk3, chunk4, chunk5)

  chunks <- rep(-1, 300)

  # positive samples in the chunks
  for (i in 1:5) {
    for (j in chks[[i]]) {
      chunks[j] <- i
    }
  }

  # define the negative constrains between chunks
  neglinks <- matrix(c(
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 1,
    1, 1, 1, 1, 0
  ),
  ncol = 5, byrow = TRUE
  )

  dca(data = data, chunks = chunks, neglinks = neglinks)
})

test_that("dca works with useD", {
  k <- 100 # sample size of each class
  n <- 3 # specify how many class
  N <- k * n # total sample number
  x1 <- mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 <- mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 <- mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  data <- as.data.frame(rbind(x1, x2, x3))

  chunk1 <- sample(1:100, 5)
  chunk2 <- sample(setdiff(1:100, chunk1), 5)
  chunk3 <- sample(101:200, 5)
  chunk4 <- sample(setdiff(101:200, chunk3), 5)
  chunk5 <- sample(201:300, 5)
  chks <- list(chunk1, chunk2, chunk3, chunk4, chunk5)

  chunks <- rep(-1, 300)

  # positive samples in the chunks
  for (i in 1:5) {
    for (j in chks[[i]]) {
      chunks[j] <- i
    }
  }

  # define the negative constrains between chunks
  neglinks <- matrix(c(
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 1,
    1, 1, 1, 1, 0
  ),
  ncol = 5, byrow = TRUE
  )

  useD <- 1

  expect_that(dca(data = data, chunks = chunks, neglinks = neglinks, useD = useD), not(throws_error()))
})

# generate necessary data set for gdmd and gdmf
k <- 100
m <- matrix(c(1, 0.5, 1, 0.5, 2, -1, 1, -1, 3), nrow = 3, byrow = T)
x1 <- mvrnorm(k, mu = c(1, 1, 1), Sigma = m)
x2 <- mvrnorm(k, mu = c(-1, 0, 0), Sigma = m)
data <- rbind(x1, x2)
# define similar constrains
simi <- rbind(t(combn(1:k, 2)), t(combn((k + 1):(2 * k), 2)))
temp <- as.data.frame(t(simi))
tol <- as.data.frame(combn(1:(2 * k), 2))
# define disimilar constrains
dism <- t(as.matrix(tol[!tol %in% simi]))
breakIT <- matrix(data = seq(100), ncol = 5, nrow = 20)

test_that("gdmd works", {
  GdmDiag(data, simi, dism)
  expect_error(GdmDiag(data = function(x) {
    4
  }, simi = simi, dism = dism),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmDiag(data = data, simi = function(x) {
    4
  }, dism = dism),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = function(x) {
    4
  }),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmDiag(data = data, simi = breakIT, dism = dism),
    regexp = "The object passed to simi should be an n x 2 matrix"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = breakIT),
    regexp = "The object passed to dism should be an n x 2 matrix"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = dism, C0 = seq(5)),
    regexp = "C0 of GdmDiag expects a number. An object of class"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = dism, C0 = "F"),
    regexp = "C0 of GdmDiag expects a number. An object of class"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = dism, threshold = seq(5)),
    regexp = "threshold of GdmDiag expects a number. An object of class"
  )
  expect_error(GdmDiag(data = data, simi = simi, dism = dism, threshold = "F"),
    regexp = "threshold of GdmDiag expects a number. An object of class"
  )
})

test_that("gdmf works", {
  GdmFull(data, simi, dism)
  expect_error(GdmFull(data = function(x) {
    4
  }, simi = simi, dism = dism),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmFull(data = data, simi = function(x) {
    4
  }, dism = dism),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmFull(data = data, simi = simi, dism = function(x) {
    4
  }),
  regexp = "expects an object coercible to a numeric matrix"
  )
  expect_error(GdmFull(data = data, simi = breakIT, dism = dism),
    regexp = "The object passed to simi should be an n x 2 matrix"
  )
  expect_error(GdmFull(data = data, simi = simi, dism = breakIT),
    regexp = "The object passed to dism should be an n x 2 matrix"
  )
  expect_error(GdmFull(data = data, simi = simi, dism = dism, maxiter = seq(5)),
    regexp = "maxiter of GdmFull expects a number. An object of class"
  )
  expect_error(GdmFull(data = data, simi = simi, dism = dism, maxiter = "F"),
    regexp = "maxiter of GdmFull expects a number. An object of class"
  )
})

test_that("rca works", {
  k <- 100 # sample size of each class
  n <- 3 # specify how many class
  N <- k * n # total sample number
  x1 <- mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 <- mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 <- mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  x <- as.data.frame(rbind(x1, x2, x3))
  x$V3 <- gl(n, k)

  chunk1 <- sample(1:100, 5)
  chunk2 <- sample(setdiff(1:100, chunk1), 5)
  chunk3 <- sample(101:200, 5)
  chunk4 <- sample(setdiff(101:200, chunk3), 5)
  chunk5 <- sample(201:300, 5)
  chks <- x[c(chunk1, chunk2, chunk3, chunk4, chunk5), ]
  chunks <- list(chunk1, chunk2, chunk3, chunk4, chunk5)

  rca(x[, 1:2], chunks)
})
