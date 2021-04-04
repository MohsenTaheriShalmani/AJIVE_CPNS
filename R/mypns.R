mypns<-function (x, sphere.type = "seq.test", alpha = 0.1, R = 100,
          nlast.small.sphere = 0)
{
  n = ncol(x)
  k = nrow(x)
  if (abs(sum(apply(x^2, 2, sum)) - n) > 1e-08) {
    stop("Error: Each column of x should be a unit vector, ||x[ , i]|| = 1.")
  }
  svd.x = svd(x, nu = nrow(x))
  uu = svd.x$u
  maxd = which(svd.x$d < 1e-15)[1]
  if (is.na(maxd) | k > n) {
    maxd = min(k, n) + 1
  }
  nullspdim = k - maxd + 1
  d = k - 1
  cat("Message from pns() : dataset is on ", d, "-sphere. \n",
      sep = "")
  if (nullspdim > 0) {
    cat(" .. found null space of dimension ", nullspdim,
        ", to be trivially reduced. \n", sep = "")
  }
  resmat = matrix(NA, d, n)
  orthaxis = list()
  orthaxis[[d - 1]] = NA
  dist = rep(NA, d - 1)
  pvalues = matrix(NA, d - 1, 2)
  ratio = rep(NA, d - 1)
  currentSphere = x
  if (nullspdim > 0) {
    for (i in 1:nullspdim) {
      oaxis = uu[, ncol(uu) - i + 1]
      r = pi/2
      pvalues[i, ] = c(NaN, NaN)
      res = acos(t(oaxis) %*% currentSphere) - r
      orthaxis[[i]] = oaxis
      dist[i] = r
      resmat[i, ] = res
      NestedSphere = rotMat(oaxis) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                       k - i, 1)
      uu = rotMat(oaxis) %*% uu
      uu = uu[1:(k - i), ]/repmat(matrix(sqrt(1 - uu[nrow(uu),
                                                     ]^2), nrow = 1), k - i, 1)
      cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ",
          "NULL space \n", sep = "")
    }
  }
  if (sphere.type == "seq.test") {
    cat(" .. sequential tests with significance level ",
        alpha, "\n", sep = "")
    isIsotropic = FALSE
    for (i in (nullspdim + 1):(d - 1)) {
      if (!isIsotropic) {
        sp = mygetSubSphere(x = currentSphere, geodesic = "small")
        center.s = sp$center
        r.s = sp$r
        resSMALL = acos(t(center.s) %*% currentSphere) -
          r.s
        sp = mygetSubSphere(x = currentSphere, geodesic = "great")
        center.g = sp$center
        r.g = sp$r
        resGREAT = acos(t(center.g) %*% currentSphere) -
          r.g
        pval1 = LRTpval(resGREAT, resSMALL, n)
        pvalues[i, 1] = pval1
        if (pval1 > alpha) {
          center = center.g
          r = r.g
          pvalues[i, 2] = NA
          cat(d - i + 1, "-sphere to ", d - i,
              "-sphere, by GREAT sphere, p(LRT) = ",
              pval1, "\n", sep = "")
        }
        else {
          pval2 = vMFtest(currentSphere, R)
          pvalues[i, 2] = pval2
          if (pval2 > alpha) {
            center = center.g
            r = r.g
            cat(d - i + 1, "-sphere to ", d - i,
                "-sphere, by GREAT sphere, p(LRT) = ",
                pval1, ", p(vMF) = ", pval2, "\n",
                sep = "")
            isIsotropic = TRUE
          }
          else {
            center = center.s
            r = r.s
            cat(d - i + 1, "-sphere to ", d - i,
                "-sphere, by SMALL sphere, p(LRT) = ",
                pval1, ", p(vMF) = ", pval2, "\n",
                sep = "")
          }
        }
      }
      else if (isIsotropic) {
        sp = mygetSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        cat(d - i + 1, "-sphere to ", d - i, "-sphere, by GREAT sphere, restricted by testing vMF distn",
            "\n", sep = "")
        pvalues[i, 1] = NA
        pvalues[i, 2] = NA
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      cur.proj = project.subsphere(x = currentSphere, center = center,
                                   r = r)
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                       k - i, 1)
      if (nrow(currentSphere) == 3) {
        plot3d(x = t(currentSphere), xlab = "",
               ylab = "", zlab = "", xlim = c(-1,
                                              1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE,
               axes = FALSE, aspect = "iso")
        axis3d(edge = "x", labels = TRUE, tick = TRUE,
               at = c(-1, 1), pos = c(NA, 0, 0))
        axis3d(edge = "y", labels = TRUE, tick = TRUE,
               at = c(-1, 1), pos = c(0, NA, 0))
        axis3d(edge = "z", labels = TRUE, tick = TRUE,
               at = c(-1, 1), pos = c(0, 0, NA))
      }
      if (nrow(currentSphere) == 2) {
        points3d(t(cur.proj), col = "red", size = 2)
      }
    }
  }
  else if (sphere.type == "BIC") {
    cat(" .. with BIC \n")
    for (i in (nullspdim + 1):(d - 1)) {
      sp = mygetSubSphere(x = currentSphere, geodesic = "small")
      center.s = sp$center
      r.s = sp$r
      resSMALL = acos(t(center.s) %*% currentSphere) -
        r.s
      sp = mygetSubSphere(x = currentSphere, geodesic = "great")
      center.g = sp$center
      r.g = sp$r
      resGREAT = acos(t(center.g) %*% currentSphere) -
        r.g
      BICsmall = n * log(mean(resSMALL^2)) + (d - i + 1 +
                                                1) * log(n)
      BICgreat = n * log(mean(resGREAT^2)) + (d - i + 1) *
        log(n)
      cat("BICsm: ", BICsmall, ", BICgr: ",
          BICgreat, "\n", sep = "")
      if (BICsmall > BICgreat) {
        center = center.g
        r = r.g
        cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ",
            "GREAT sphere, BIC \n", sep = "")
      }
      else {
        center = center.s
        r = r.s
        cat(d - i + 1, "-sphere to ", d - i, "-sphere, by ",
            "SMALL sphere, BIC \n", sep = "")
      }
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                       k - i, 1)
    }
  }
  else if (sphere.type == "small" | sphere.type == "great") {
    pvalues = NaN
    for (i in (nullspdim + 1):(d - 1)) {
      sp = mygetSubSphere(x = currentSphere, geodesic = sphere.type)
      center = sp$center
      r = sp$r
      res = acos(t(center) %*% currentSphere) - r
      orthaxis[[i]] = center
      dist[i] = r
      resmat[i, ] = res
      NestedSphere = rotMat(center) %*% currentSphere
      currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                     NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                       k - i, 1)
    }
  }
  else if (sphere.type == "bi.sphere") {
    if (nlast.small.sphere < 0) {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be >= 0. \n")
      return(NULL)
    }
    mx = (d - 1) - nullspdim
    if (nlast.small.sphere > mx) {
      cat("!!! Error from pns(): \n")
      cat("!!! nlast.small.sphere should be <= ",
          mx, " for this data. \n", sep = "")
      return(NULL)
    }
    pvalues = NaN
    if (nlast.small.sphere != mx) {
      for (i in (nullspdim + 1):(d - 1 - nlast.small.sphere)) {
        sp = mygetSubSphere(x = currentSphere, geodesic = "great")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res
        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                       NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                         k - i, 1)
      }
    }
    if (nlast.small.sphere != 0) {
      for (i in (d - nlast.small.sphere):(d - 1)) {
        sp = mygetSubSphere(x = currentSphere, geodesic = "small")
        center = sp$center
        r = sp$r
        res = acos(t(center) %*% currentSphere) - r
        orthaxis[[i]] = center
        dist[i] = r
        resmat[i, ] = res
        NestedSphere = rotMat(center) %*% currentSphere
        currentSphere = NestedSphere[1:(k - i), ]/repmat(matrix(sqrt(1 -
                                                                       NestedSphere[nrow(NestedSphere), ]^2), nrow = 1),
                                                         k - i, 1)
      }
    }
  }
  else {
    print("!!! Error from pns():")
    print("!!! sphere.type must be 'seq.test', 'small', 'great', 'BIC', or 'bi.sphere'")
    print("!!!   Terminating execution     ")
    return(NULL)
  }
  S1toRadian = atan2(currentSphere[2, ], currentSphere[1, ])
  meantheta = geodmeanS1(S1toRadian)$geodmean
  orthaxis[[d]] = meantheta
  resmat[d, ] = mod(S1toRadian - meantheta + pi, 2 * pi) -
    pi
  par(mfrow = c(1, 1), mar = c(4, 4, 1, 1), mgp = c(2.5, 1,
                                                    0), cex = 0.8)
  plot(currentSphere[1, ], currentSphere[2, ], xlab = "",
       ylab = "", xlim = c(-1, 1), ylim = c(-1, 1), asp = 1)
  abline(h = 0, v = 0)
  points(cos(meantheta), sin(meantheta), pch = 1, cex = 3,
         col = "black", lwd = 5)
  abline(a = 0, b = sin(meantheta)/cos(meantheta), lty = 3)
  l = mod(S1toRadian - meantheta + pi, 2 * pi) - pi
  points(cos(S1toRadian[which.max(l)]), sin(S1toRadian[which.max(l)]),
         pch = 4, cex = 3, col = "blue")
  points(cos(S1toRadian[which.min(l)]), sin(S1toRadian[which.min(l)]),
         pch = 4, cex = 3, col = "red")
  legend("topright", legend = c("Geodesic mean",
                                "Max (+)ve from mean", "Min (-)ve from mean"),
         col = c("black", "blue", "red"), pch = c(1,
                                                  4, 4))
  {
    cat("\n")
    cat("length of BLUE from geodesic mean : ", max(l),
        " (", round(max(l) * 180/pi), " degree)",
        "\n", sep = "")
    cat("length of RED from geodesic mean : ", min(l),
        " (", round(min(l) * 180/pi), " degree)",
        "\n", sep = "")
    cat("\n")
  }
  radii = 1
  for (i in 1:(d - 1)) {
    radii = c(radii, prod(sin(dist[1:i])))
  }
  resmat = flipud0(repmat(matrix(radii, ncol = 1), 1, n) *
                     resmat)
  PNS = list()
  PNS$radii = radii
  PNS$orthaxis = orthaxis
  PNS$dist = dist
  PNS$pvalues = pvalues
  PNS$ratio = ratio
  PNS$basisu = NULL
  PNS$mean = c(PNSe2s(matrix(0, d, 1), PNS))
  if (sphere.type == "seq.test") {
    PNS$sphere.type = "seq.test"
  }
  else if (sphere.type == "small") {
    PNS$sphere.type = "small"
  }
  else if (sphere.type == "great") {
    PNS$sphere.type = "great"
  }
  else if (sphere.type == "BIC") {
    PNS$sphere.type = "BIC"
  }
  varPNS = apply(abs(resmat)^2, 1, sum)/n
  total = sum(varPNS)
  propPNS = varPNS/total * 100
  return(list(resmat = resmat, PNS = PNS, percent = propPNS))
}
