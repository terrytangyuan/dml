context('test algorithms')
library(lfda)
library(MASS)

test_that('dca works', {
  k = 100        # sample size of each class
  n = 3          # specify how many class
  N = k * n      # total sample number
  x1 = mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 = mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 = mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  data = as.data.frame(rbind(x1, x2, x3))

  chunk1 = sample(1:100, 5)
  chunk2 = sample(setdiff(1:100, chunk1), 5)
  chunk3 = sample(101:200, 5)
  chunk4 = sample(setdiff(101:200, chunk3), 5)
  chunk5 = sample(201:300, 5)
  chks = list(chunk1, chunk2, chunk3, chunk4, chunk5)

  chunks = rep(-1, 300)

  # positive samples in the chunks
  for (i in 1:5) {
    for (j in chks[[i]]) {
      chunks[j] = i
    }
  }

  # define the negative constrains between chunks
  neglinks = matrix(c(
  		0, 0, 1, 1, 1,
  		0, 0, 1, 1, 1,
  		1, 1, 0, 0, 0,
  		1, 1, 0, 0, 1,
  		1, 1, 1, 1, 0),
  		ncol = 5, byrow = TRUE)

  expect_that(dca(data = data, chunks = chunks, neglinks = neglinks), not(throws_error()))
})

test_that('dca works with useD', {
  k = 100        # sample size of each class
  n = 3          # specify how many class
  N = k * n      # total sample number
  x1 = mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 = mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 = mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  data = as.data.frame(rbind(x1, x2, x3))
  
  chunk1 = sample(1:100, 5)
  chunk2 = sample(setdiff(1:100, chunk1), 5)
  chunk3 = sample(101:200, 5)
  chunk4 = sample(setdiff(101:200, chunk3), 5)
  chunk5 = sample(201:300, 5)
  chks = list(chunk1, chunk2, chunk3, chunk4, chunk5)
  
  chunks = rep(-1, 300)
  
  # positive samples in the chunks
  for (i in 1:5) {
    for (j in chks[[i]]) {
      chunks[j] = i
    }
  }
  
  # define the negative constrains between chunks
  neglinks = matrix(c(
    0, 0, 1, 1, 1,
    0, 0, 1, 1, 1,
    1, 1, 0, 0, 0,
    1, 1, 0, 0, 1,
    1, 1, 1, 1, 0),
    ncol = 5, byrow = TRUE)
  
  useD = 1
  
  expect_that(dca(data = data, chunks = chunks, neglinks = neglinks, useD = useD), not(throws_error()))
})

# generate necessary data set for gdmd and gdmf
k = 100
m <- matrix(c(1, 0.5, 1, 0.5, 2, -1, 1, -1, 3), nrow =3, byrow = T)
x1 <- mvrnorm(k, mu = c(1, 1, 1), Sigma = m)
x2 <- mvrnorm(k, mu = c(-1, 0, 0), Sigma = m)
data <- rbind(x1, x2)
# define similar constrains
simi <- rbind(t(combn(1:k, 2)), t(combn((k+1):(2*k), 2)))
temp <-  as.data.frame(t(simi))
tol <- as.data.frame(combn(1:(2*k), 2))
# define disimilar constrains
dism <- t(as.matrix(tol[!tol %in% simi]))

test_that('gdmd works', {
  expect_that(GdmDiag(data, simi, dism), not(throws_error()))
})

test_that('gdmf works', {
  expect_that(GdmFull(data, simi, dism), not(throws_error()))
})

test_that('rca works', {
  k = 100        # sample size of each class
  n = 3          # specify how many class
  N = k * n      # total sample number
  x1 = mvrnorm(k, mu = c(-10, 6), matrix(c(10, 4, 4, 10), ncol = 2))
  x2 = mvrnorm(k, mu = c(0, 0), matrix(c(10, 4, 4, 10), ncol = 2))
  x3 = mvrnorm(k, mu = c(10, -6), matrix(c(10, 4, 4, 10), ncol = 2))
  x = as.data.frame(rbind(x1, x2, x3))
  x$V3 = gl(n, k)

  chunk1 = sample(1:100, 5)
  chunk2 = sample(setdiff(1:100, chunk1), 5)
  chunk3 = sample(101:200, 5)
  chunk4 = sample(setdiff(101:200, chunk3), 5)
  chunk5 = sample(201:300, 5)
  chks = x[c(chunk1, chunk2, chunk3, chunk4, chunk5), ]
  chunks = list(chunk1, chunk2, chunk3, chunk4, chunk5)

  expect_that(rca(x[ , 1:2], chunks), not(throws_error()))
})

