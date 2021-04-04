mygetSubSphere<-function (x, geodesic = "great") 
{
  svd.x = svd(x)
  initialCenter = svd.x$u[, ncol(svd.x$u)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10
  while (err > TOL) {
    c0 = c0/norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = mysphereFit(x = TpX, initialCenter = rep(0, d - 1), geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi) {
      r = pi/2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi/2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30) {
      break
    }
  }
  i1save = list()
  i1save$Gnow = Gnow
  i1save$center = center
  i1save$r = r
  U = princomp(t(x))$loadings[, ]
  initialCenter = U[, ncol(U)]
  c0 = initialCenter
  TOL = 1e-10
  cnt = 0
  err = 1
  n = ncol(x)
  d = nrow(x)
  Gnow = 1e+10
  while (err > TOL) {
    c0 = c0/norm(c0, type = "2")
    rot = rotMat(c0)
    TpX = LogNPd(rot %*% x)
    fit = mysphereFit(x = TpX, initialCenter = rep(0, d - 1), geodesic = geodesic)
    newCenterTp = fit$center
    r = fit$r
    if (r > pi) {
      r = pi/2
      svd.TpX = svd(TpX)
      newCenterTp = svd.TpX$u[, ncol(svd.TpX$u)] * pi/2
    }
    newCenter = ExpNPd(newCenterTp)
    center = solve(rot, newCenter)
    Gnext = objfn(center, r, x)
    err = abs(Gnow - Gnext)
    Gnow = Gnext
    c0 = center
    cnt = cnt + 1
    if (cnt > 30) {
      break
    }
  }
  if (i1save$Gnow == min(Gnow, i1save$Gnow)) {
    center = i1save$center
    r = i1save$r
  }
  if (r > pi/2) {
    center = -center
    r = pi - r
  }
  return(list(center = c(center), r = r))
}
