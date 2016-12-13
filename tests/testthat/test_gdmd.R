set.seed(602)
library(MASS)
library(scatterplot3d)

# generate simulated Gaussian data
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

# transform data using GdmDiag
result <- GdmDiag(data, simi, dism)
newData <- result$newData
# plot original data
color <- gl(2, k, labels = c("red", "blue"))
par(mfrow = c(2, 1), mar = rep(0, 4) + 0.1)
scatterplot3d(data, color = color, cex.symbols = 0.6,
		  xlim = range(data[, 1], newData[, 1]),
		  ylim = range(data[, 2], newData[, 2]),
		  zlim = range(data[, 3], newData[, 3]),
		  main = "Original Data")
# plot GdmDiag transformed data
scatterplot3d(newData, color = color, cex.symbols = 0.6,
		  xlim = range(data[, 1], newData[, 1]),
		  ylim = range(data[, 2], newData[, 2]),
		  zlim = range(data[, 3], newData[, 3]),
		  main = "Transformed Data")
