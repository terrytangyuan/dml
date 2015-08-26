lfda_wp <- function(x, y, r = 3, metric = c("orthonormalized","plain","weighted"),knn = 5,
                    method = c('lfda', 'klfda', 'self'), ...){
  if(!requireNamespace('lfda')){stop("Please install lfda package to perform Local Fisher Discriminant Analysis. ")}
  method <- match.arg(method)
  metric <- match.arg(metric)
  modelArgs <- c(list(x, y, r, metric, knn, ...))
  out <- do.call(method, modelArgs)
  out$call <- NULL
  out
}
