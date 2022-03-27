
FSAR <- function(X, p = 10, my_lags, n_basis = 10, percent_CPV = 0.95){
  Nt <- ncol(X)
  N <- nrow(X)
  basis <- create.bspline.basis(c(0, 1), nbasis = n_basis, norder = 4)
  # basis <- create.fourier.basis(c(0, 1), nbasis = basisfd)
  fdX <- Data2fd(argvals = s, t(X), basis)
  fdXpca <- pca.fd(fdX, nharm = p)
  p <- min(which(cumsum(fdXpca$values) >= percent_CPV * sum(fdXpca$values)))
  
  fdXpca <- pca.fd(fdX, nharm = p)
  eigenvalues <- fdXpca$values; scoresX <- fdXpca$scores
  # jth column of scoresX contains scores of the jth EFPC
  harmonicsX <- fdXpca$harmonics # extract the EFPC's
  varnceprop <- fdXpca$varprop # proportion of variance explained by the EFP's
  
  # Now using the least squares ---------------------------------------------
  k <- length(my_lags)
  X_LP <- matrix(0, N - max(my_lags), k*p)
  for(i in (max(my_lags)+1):N){
    X_LP_i <- rep(0, k*p)
    for(j in 1:k){
      X_LP_i[1:p + (j-1)*p] <- scoresX[i - my_lags[j], ]
    }
    X_LP[i - max(my_lags), ] <- X_LP_i
  }
  
  X_RP <- matrix(0, nrow = N - max(my_lags), ncol = p)
  for(i in (max(my_lags)+1):N){
    X_RP[i-max(my_lags), ] <- scoresX[i,]
  }
  
  # Least square estimation
  Phi_hat <- ginv(t(X_LP) %*% X_LP) %*% t(X_LP) %*% X_RP
  
  # Use it to estimate the kernel, there are k kernels
  Phi_hat <- lapply(1:length(my_lags), function(index){
    Phi_hat[1:p + (index-1)*p, ]
  })
  
  V <- sapply(1:p,  function(j){
    eval.fd(evalarg = s, harmonicsX[j])
  })
  
  kernel_est <- lapply(1:length(my_lags), function(index){
    kernel <- matrix(0, Nt, Nt)
    for(tk in 1:Nt){
      for(sj in 1:Nt){
        kernel[tk, sj] <- t(V[tk, ]) %*% Phi_hat[[index]] %*% V[sj, ]
      }
    }
    kernel
  })
  
  list(Phi_hat = Phi_hat, kernel_est = kernel_est)
}