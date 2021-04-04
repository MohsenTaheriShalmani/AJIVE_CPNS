library(minpack.lm)
mysphereFit<-function (x, initialCenter = NULL, geodesic = "great")
{
  if (is.null(initialCenter)) {
    initialCenter = apply(x, 1, mean)
  }
  op = nls.lm(par = initialCenter, fn = sphere.res, jac = sphere.jac,
              x = x, is.greatCircle = ifelse(geodesic == "great", TRUE, FALSE),
              control = nls.lm.control(maxiter = 10))

  center = coef(op)
  di = sqrt(apply((x - repmat(matrix(center, ncol = 1), 1,
                              ncol(x)))^2, 2, sum))
  if (geodesic == "great") {
    r = pi/2
  }
  else {
    r = mean(di)
  }
  list(center = center, r = r)
}
