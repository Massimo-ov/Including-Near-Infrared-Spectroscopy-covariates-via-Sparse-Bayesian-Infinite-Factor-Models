# pacchetti
  library(infinitefactor)   # CRAN/GitHub pack with jointRot, msf, lmean, etc
  
  # --- Input che ci aspettiamo ---
  # Lambda_samples : lista length T, ogni elemento p x k (matrici)
  # sigma_samples  : matrix T x p (o psi_samples), sigma_samples[t, ] = psi_t
  # (opzionale) eta_samples : list length T di factor scores (n x k)
  # labels (opzionale) : nomi delle variabili
  
  T <- length(Lambda_samples)
  p <- nrow(Lambda_samples[[1]])
  k <- ncol(Lambda_samples[[1]])
  
  # --------------------------
  # 1) Build pivot (reference) e allineamento
  # --------------------------
  # pivot = varimax(mean_loadings)
  pivot_mat <- as.matrix(varimax(lmean(Lambda_samples))[[1]])  # varimax returns loadings object
  
  # allineo campioni: varimax + match to pivot (msfOUT + aplr)
  aligned_Lambda <- vector("list", T)
  for(t in 1:T){
    va_t <- varimax(Lambda_samples[[t]])
    Lrot <- as.matrix(va_t[[1]])            # p x k
    # compute signed permutation to align Lrot -> pivot
    perm <- msfOUT(Lrot, pivot_mat)        # returns permutation+sign info
    Laligned <- aplr(Lrot, perm)           # apply permutation/sign
    aligned_Lambda[[t]] <- Laligned
  }
  
  # se hai eta_samples usa invece:
  # aligned <- jointRot(Lambda_samples, eta_samples)
  # aligned_Lambda <- aligned$lambda
  # aligned_eta <- aligned$eta
  
  # --------------------------
  # 2) Calcoli: per-sample explained variances, per-variable fractions, inclusion probs
  # --------------------------
  # prepara array p x k x T
  Larr <- array(NA, dim = c(p, k, T))
  for(t in 1:T) Larr[,,t] <- aligned_Lambda[[t]]
  
  # psi: prendi sigma_samples (T x p)
  Psi <- sigma_samples   # assumi dim T x p
  
  # Per-sample total variance and factor contributions
  factor_var <- matrix(NA, nrow = T, ncol = k)
  total_var  <- numeric(T)
  for(t in 1:T){
    Lt <- Larr[,,t]
    psi_t <- as.numeric(Psi[t, ])
    factor_var[t, ] <- colSums(Lt^2)                 # ||lambda_{.h}||^2
    total_var[t]    <- sum(rowSums(Lt^2) + psi_t)    # trace(Lambda Lambda' + Psi)
  }
  factor_prop <- sweep(factor_var, 1, total_var, "/")  # T x k
  
  # Per-variable fraction explained by factor h: p x k x T
  pervar_prop <- array(NA, dim = c(p, k, T))
  for(t in 1:T){
    Lt <- Larr[,,t]
    psi_t <- as.numeric(Psi[t, ])
    denom_j <- rowSums(Lt^2) + psi_t        # length p
    pervar_prop[,,t] <- sweep(Lt^2, 1, denom_j, "/")
  }
  
  # Inclusion probability P(|lambda| > eps)
  eps <- 1e-3
  incl_arr <- abs(Larr) > eps  # p x k x T logical
  incl_prob <- apply(incl_arr, c(1,2), mean)  # p x k
  
  # Riassunti posteriori:
  factor_prop_mean <- colMeans(factor_prop)     # mean % variance explained per factor
  factor_prop_ci   <- apply(factor_prop, 2, quantile, c(0.025, 0.975))
  
  # Per-variable summaries (mean + CI)
  pervar_mean <- apply(pervar_prop, c(1,2), mean)        # p x k
  pervar_ci_low <- apply(pervar_prop, c(1,2), quantile, 0.025)
  pervar_ci_high<- apply(pervar_prop, c(1,2), quantile, 0.975)

#PLOT FOR IDENTIFY THE INDEX OF H* ALIAS THE MOST PREDICTIVE FACTOR
 # posizione finta per Y
  x_y <- min(labels) - diff(range(labels))/20
  
  # limiti asse X includendo Y
  x_lim <- c(x_y, max(labels))
  
  # --- PLOT SOLO SPETTRO ---
  matplot(
    labels,
    pervar_mean[-1, ],
    type="l", lwd = 2,
    lty=1,
    col=1:k,
    xlab="Wavelength",
    ylab="Variance fraction",
    xlim = x_lim     # <-- QUESTO È IL FIX
  )
  
  # --- AGGIUNGO PUNTI PER Y ---
  for(h in 1:k){
    points(
      x_y,
      pervar_mean[1, h],
      col=h,
      pch=16,
      cex=2
    )
  }
  
  # etichetta Y leggermente sotto
  text(
    x_y,
    par("usr")[3],   # fondo asse Y
    "Y",
    pos=1,
    xpd=NA           # permette di uscire leggermente dal box
  )
  
  legend("topright", legend=paste0("F",1:k), col=1:k, lty=1, bty = "n",      # no box
         cex = 0.9,
         ncol = ifelse(k > 6, 2, 1))
