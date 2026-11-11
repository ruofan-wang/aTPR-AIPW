# The sign convention is J = -d E(score)/d parameter. The influence function
# equals the direct moment contribution + score %*% solve(t(J), gradient).
.atpr_influence <- function(nuisance, theta, alpha, G, groups, Q, V, y, H,
                            degree, prob_floor) {
  n <- length(groups)
  candidates <- Q == 1
  validated <- candidates & V == 1
  reference <- as.numeric(groups == 1L)
  ig <- ie <- it <- vector("list", 2L)
  cursor <- 0L
  for (s in seq_len(2L)) {
    ig[[s]] <- cursor + seq_along(nuisance[[s]]$gamma)
    cursor <- max(ig[[s]])
  }
  for (s in seq_len(2L)) {
    ie[[s]] <- cursor + seq_along(nuisance[[s]]$beta)
    cursor <- max(ie[[s]])
  }
  for (s in seq_len(2L)) {
    it[[s]] <- cursor + seq_len(2L)
    cursor <- max(it[[s]])
  }
  ia <- cursor + seq_len(degree + 1L)
  pdim <- max(ia)
  z <- numeric(pdim)
  for (s in seq_len(2L)) {
    z[ig[[s]]] <- nuisance[[s]]$gamma
    z[ie[[s]]] <- nuisance[[s]]$beta
    z[it[[s]]] <- theta[[s]]
  }
  z[ia] <- alpha

  components <- function(zz, with_scores = TRUE) {
    pp <- rep(1, n)
    mm <- eta <- numeric(n)
    for (s in seq_len(2L)) {
      i <- which(groups == s)
      c <- which(groups == s & candidates)
      pp[c] <- stats::plogis(drop(nuisance[[s]]$A[c, , drop = FALSE] %*% zz[ig[[s]]]))
      mm[c] <- stats::plogis(drop(nuisance[[s]]$W[c, , drop = FALSE] %*% zz[ie[[s]]]))
      eta[i] <- drop(G[i, , drop = FALSE] %*% zz[it[[s]]])
    }
    rr <- stats::plogis(eta)
    pb <- pmax(pp, prob_floor)
    pseudo <- Q * mm
    pseudo[validated] <- mm[validated] + (y[validated] - mm[validated]) / pb[validated]
    D <- .atpr_density_basis(eta, degree)
    lp <- drop(D %*% zz[ia])
    qq <- stats::plogis(lp)
    ww <- rep(1, n)
    ww[groups == 2L] <- exp(lp[groups == 2L])
    # The subgroup-size multiplier cancels from aTPR_2, including its derivative.
    B <- cbind(reference * rr * H, reference * rr,
               (1 - reference) * rr * ww * H, (1 - reference) * rr * ww)
    mu <- colMeans(B)
    at <- c(mu[1L] / mu[2L], mu[3L] / mu[4L])
    L <- NULL
    if (with_scores) {
      L <- matrix(0, n, pdim)
      for (s in seq_len(2L)) {
        i <- which(groups == s)
        c <- which(groups == s & candidates)
        v <- which(groups == s & validated)
        L[c, ig[[s]]] <- nuisance[[s]]$A[c, , drop = FALSE] * (V[c] - pp[c])
        L[v, ie[[s]]] <- nuisance[[s]]$W[v, , drop = FALSE] * (y[v] - mm[v])
        L[i, it[[s]]] <- G[i, , drop = FALSE] * (pseudo[i] - rr[i])
      }
      L[, ia] <- D * (reference - qq)
    }
    list(p = pp, pb = pb, m = mm, r = rr, eta = eta, D = D, q = qq,
         w = ww, B = B, mu = mu, at = at, delta = unname(at[2L] - at[1L]), L = L)
  }

  x <- components(z)
  J <- matrix(0, pdim, pdim)
  gradient <- numeric(pdim)
  observed_weight <- numeric(n)
  observed_weight[validated] <- 1 / x$pb[validated]
  slope <- if (degree == 2L) alpha[2L] + 2 * alpha[3L] * x$eta else rep(alpha[2L], n)
  for (s in seq_len(2L)) {
    i <- which(groups == s)
    c <- which(groups == s & candidates)
    v <- which(groups == s & validated)
    A <- nuisance[[s]]$A
    W <- nuisance[[s]]$W
    ac <- A[c, , drop = FALSE]
    wv <- W[v, , drop = FALSE]
    gi <- G[i, , drop = FALSE]
    J[ig[[s]], ig[[s]]] <- crossprod(ac, ac * (x$p[c] * (1 - x$p[c]))) / n
    J[ie[[s]], ie[[s]]] <- crossprod(wv, wv * (x$m[v] * (1 - x$m[v]))) / n
    J[it[[s]], it[[s]]] <- crossprod(gi, gi * (x$r[i] * (1 - x$r[i]))) / n
    pderiv <- (y[v] - x$m[v]) * x$p[v] * (1 - x$p[v]) /
      x$pb[v]^2 * (x$p[v] > prob_floor)
    J[it[[s]], ig[[s]]] <- crossprod(G[v, , drop = FALSE], A[v, , drop = FALSE] * pderiv) / n
    mderiv <- Q[i] * (1 - observed_weight[i]) * x$m[i] * (1 - x$m[i])
    J[it[[s]], ie[[s]]] <- -crossprod(gi, W[i, , drop = FALSE] * mderiv) / n
    derivative_basis <- if (degree == 2L) {
      cbind(0, 1, 2 * x$eta[i])
    } else cbind(rep(0, length(i)), rep(1, length(i)))
    derivative_score <- derivative_basis * (reference[i] - x$q[i]) -
      x$D[i, , drop = FALSE] * (x$q[i] * (1 - x$q[i]) * slope[i])
    J[ia, it[[s]]] <- -crossprod(derivative_score, gi) / n
    dmass <- if (s == 1L) {
      x$r[i] * (1 - x$r[i])
    } else x$w[i] * (x$r[i] * (1 - x$r[i]) + x$r[i] * slope[i])
    sign <- if (s == 1L) -1 else 1
    gradient[it[[s]]] <- sign * colSums(gi * ((H[i] - x$at[s]) * dmass)) /
      (n * x$mu[2L * s])
  }
  J[ia, ia] <- crossprod(x$D, x$D * (x$q * (1 - x$q))) / n
  i <- which(groups == 2L)
  gradient[ia] <- colSums(x$D[i, , drop = FALSE] *
                           ((H[i] - x$at[2L]) * x$r[i] * x$w[i])) / (n * x$mu[4L])
  direct <- (x$B[, 3L] - x$at[2L] * x$B[, 4L]) / x$mu[4L] -
    (x$B[, 1L] - x$at[1L] * x$B[, 2L]) / x$mu[2L]
  residual <- max(abs(colMeans(x$L)))
  scaled <- max(abs(colMeans(x$L)) / pmax(sqrt(colMeans(x$L^2)), 1e-12))
  h <- drop(solve(t(J), gradient))
  phi <- direct + drop(x$L %*% h)
  centered <- phi - stats::ave(phi, groups, FUN = mean)
  se_stratified <- sqrt(sum(centered^2)) / n
  se_iid <- sqrt(sum(phi^2)) / n
  list(z = z, J = J, gradient = gradient, components = components,
       phi = phi, se_stratified = se_stratified, se_iid = se_iid,
       score_residual = residual, scaled_score_residual = scaled,
       rcond = rcond(J), indices = list(validation = ig, outcome = ie,
                                       calibration = it, density = ia))
}
