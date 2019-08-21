KfilterYKP<-function (num, y, A, mu0, Sigma0, Phi, cQ, cR) 
{
  Q = t(cQ) %*% cQ
  R = t(cR) %*% cR
  Phi = as.matrix(Phi)
  pdim = nrow(Phi)
  y = as.matrix(y)
  qdim = ncol(y)
  xp = array(NA, dim = c(pdim, 1, num))
  Pp = array(NA, dim = c(pdim, pdim, num))
  xf = array(NA, dim = c(pdim, 1, num))
  Pf = array(NA, dim = c(pdim, pdim, num))
  innov = array(NA, dim = c(qdim, 1, num))
  sig = array(NA, dim = c(qdim, qdim, num))
  x00 = as.matrix(mu0, nrow = pdim, ncol = 1)
  P00 = as.matrix(Sigma0, nrow = pdim, ncol = pdim)
  xp[, , 1] = Phi %*% x00
  Pp[, , 1] = Phi %*% P00 %*% t(Phi) + Q
  sigtemp = A %*% Pp[, , 1] %*% t(A) + R
  sig[, , 1] = (t(sigtemp) + sigtemp)/2
  siginv = solve(sig[, , 1])
  K = Pp[, , 1] %*% t(A) %*% siginv
  innov[, , 1] = y[1, ] - A %*% xp[, , 1]
  xf[, , 1] = xp[, , 1] + K %*% innov[, , 1]
  Pf[, , 1] = Pp[, , 1] - K %*% A %*% Pp[, , 1]
  sigmat = as.matrix(sig[, , 1], nrow = qdim, ncol = qdim)
  like = log(det(sigmat)) + t(innov[, , 1]) %*% siginv %*% 
    innov[, , 1]
  for (i in 2:num) {
    if (num < 2) 
      break
    xp[, , i] = Phi %*% xf[, , i - 1]
    Pp[, , i] = Phi %*% Pf[, , i - 1] %*% t(Phi) + Q
    sigtemp = A %*% Pp[, , i] %*% t(A) + R
    sig[, , i] = (t(sigtemp) + sigtemp)/2
    siginv = solve(sig[, , i])
    K = Pp[, , i] %*% t(A) %*% siginv
    innov[, , i] = y[i, ] - A %*% xp[, , i]
    xf[, , i] = xp[, , i] + K %*% innov[, , i]
    Pf[, , i] = Pp[, , i] - K %*% A %*% Pp[, , i]
    sigmat = as.matrix(sig[, , i], nrow = qdim, ncol = qdim)
    like = like + log(det(sigmat)) + t(innov[, , i]) %*% 
      siginv %*% innov[, , i]
  }
  like = 0.5 * like
  list(xp = xp, Pp = Pp, xf = xf, Pf = Pf, like = like, innov = innov, 
       sig = sig, Kn = K)
}

KfilterYKP(num, y, A, mu0, Sigma0, Phi, cQ, cR)
