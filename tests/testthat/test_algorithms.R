context('test algorithms')
library(lfda)

test_that('lfda works', {
  k <- iris[,-5]
  y <- iris[,5]
  r <- 3
  model <- lfda(k, y, r, metric = "plain")
})
