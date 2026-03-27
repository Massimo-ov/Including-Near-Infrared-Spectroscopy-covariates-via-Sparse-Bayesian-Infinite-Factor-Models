#Sparse Bayesian infinite factor models Bhattacharya,Dunson 2011, prior DELTA DURANTE VERSION 2017,  BY EWAN Poworoznek
  Y <- as.matrix(df_train)
  p <- ncol(Y)
  n <- nrow(Y)
  
  as = 1                          # gamma hyperparameters for residual precision
  bs = 0.3                        
  df = 3                          # gamma hyperparameters for t_{ij}
  ad1 = 2.1
  bd1 = 1                         # gamma hyperparameters for delta_1
  ad2 = 3.1
  bd2 = 1                         # gamma hyperparameters delta_h, h >= 2
  adf = 1
  bdf = 1 
  b0 = 1
  b1 = 0.0005
  prop = 1
  epsilon = 1e-3
  
  kinit = floor(log(p)*3)
  
  nrun <-30000
  thin <- 1
  burn <- 20000
  sp = floor((nrun - burn)/thin)        # number of posterior samples
  log_likelihood <- function(Y, Lambda, eta, ps){
    mu <- eta %*% t(Lambda)
    ll <- sum(dnorm(Y, mean = mu, sd = 1/sqrt(ps), log = TRUE))
    return(ll)
  }
 
  loglike <- numeric(nrun)
  k_star_history <- numeric(nrun)
  num = 0
  k=kinit                               # no. of factors to start with
  
  # --- Initial values --- #
  ps = rgamma(p, as, bs)
  Sigma = diag(1/ps)                             # Sigma = diagonal residual covariance
  Lambda = matrix(1, nrow = p, ncol = k)
  ta = matrix(rnorm(n*k), nrow = n, ncol = k)    # factor loadings & latent factors
  meta = matrix(0,nrow = n, ncol = k)
  veta = diag(k)                                 # latent factor distribution = standard normal
  
  psijh = matrix(rgamma(p*k, df/2, df/2), nrow = p, ncol = k)     # local shrinkage coefficients
  theta = c(1 / rgamma(1,ad1,bd1), 1 / rgamma(k-1,ad2,bd2))       # gobal shrinkage coefficients multilpliers
  tauh = 1 / cumprod(theta)                                       # global shrinkage coefficients
  Plam = t(t(psijh) * (tauh))                                     # precision of loadings rows
  
  # --- Allocate output object memory --- #
  LAMBDA = list()
  SIGMA =  matrix(NA, nrow = sp, ncol = p)
  K = rep(NA, sp)
  ind = 1
  
  #------start gibbs sampling-----
  
  
  for(i in 1:nrun) {
    # -- Update eta -- #
    Lmsg = Lambda * ps
    Veta1 = diag(k) + crossprod(Lmsg, Lambda)
    Tmat = chol(Veta1)
    S = backsolve(Tmat, diag(k))
    Veta = tcrossprod(S)                                              # Veta = inv(Veta1)
    Meta = Y %*% Lmsg %*% Veta                                      # n x k 
    eta = Meta + tcrossprod(matrix(rnorm(n*k), 
                                   nrow = n, 
                                   ncol = k), S)                    # update eta in a block
    
    # -- update Lambda (rue & held) -- #
    eta2 = crossprod(eta)    # prepare eta crossproduct before the loop
    zlams = rnorm(k*p)       # generate normal draws all at once 
    
    for(j in 1:p) {
      Llamt = chol(diag(Plam[j,]) + ps[j]*eta2)
      Lambda[j,] = t(backsolve(Llamt, zlams[1:k + (j-1)*k]) + 
                       backsolve(Llamt,
                                 forwardsolve(t(Llamt),
                                              ps[j] * crossprod(eta, Y[,j]))))
    }  
    
    
    
    
    #------Update psi_{jh}'s------#
    psijh = matrix(rgamma(p*k,
                          df/2 + 0.5,
                          df/2 + t(t(Lambda)^2 * (tauh))/2),
                   nrow = p, ncol = k)
    
    #------Update theta & tauh------#
    mat = psijh * Lambda^2
    ad = ad1 + 0.5*p*k
    bd = bd1 + 0.5 * theta[1] * sum(tauh*colSums(mat))
    theta[1] = 1 / rgamma(1,ad,bd)           
    tauh = 1 / cumprod(theta)
    
    
    for(h in 2:k) {
      ad = ad2 + 0.5*p*(k-h+1)
      bd = bd2 + 0.5 * theta[h] * sum(tauh[h:k]*colSums(mat[,h:k, drop = F]))
      theta[h] = 1 / rgamma(1,ad,bd)
      tauh = 1 / cumprod(theta)
    }
    
    # -- Update Sigma -- #
    Ytil = Y - tcrossprod(eta, Lambda)
    ps= rgamma(p, as + 0.5*n, bs+0.5*colSums(Ytil^2))
    Sigma = diag(1/ps)
    
    #---update precision parameters----#
    Plam = t(t(psijh) * tauh)
    
    #logpost[i] <- log_posterior(Y, Lambda, eta, ps,psijh, theta, tauh,as, bs, df, ad1, bd1, ad2, bd2)
    
    # Calcola e salva la log likelihood
    loglike[i] <- log_likelihood(
      Y, Lambda,eta,ps)
    
    #Kstar
     k_star_history[i] <- k
     
    # ----- make adaptations ----#
    prob = 1/exp(b0 + b1*i)                    # probability of adapting
    uu = runif(1)
    lind = colSums(abs(Lambda) < epsilon)/p    # proportion of elements in each column less than eps in magnitude
    vec = lind >= prop
    num = sum(vec)                             # number of redundant columns
    
    if((uu < prob)) {
      if((i > 20) & (num == 0) & all(lind < 0.995)) {
        k = k + 1
        Lambda = cbind(Lambda, rep(0,p))
        eta = cbind(eta,rnorm(n))
        psijh = cbind(psijh, rgamma(p,df/2,df/2))
        theta[k] = 1 / rgamma(1, ad2,bd2)
        tauh = 1 / cumprod(theta)
        Plam = t(t(psijh) * tauh)
      } else {
        if (num > 0) {
          k = max(k - num,1)
          Lambda = Lambda[,!vec, drop = F]
          psijh = psijh[,!vec, drop = F]
          eta = eta[,!vec, drop = F]
          theta = theta[!vec]
          tauh = 1 / cumprod(theta)
          Plam = t(t(psijh) * tauh)
        }
      }
    }
    # -- save sampled values (after thinning) -- #
    if((i %% thin == 0) & (i > burn)) {
      LAMBDA[[ind]] = Lambda
      SIGMA[ind,] = diag(Sigma)
      K[ind] = k
      ind = ind + 1
    }
    
    if(((i %% 1000) == 0)) {
      cat(i,"\n")
    }
  }
  
  Lambda_samples <- LAMBDA
  sigma_samples <- SIGMA


#CALCOLO MSE E BETA (added by the author Massimo Armano)
  y <- Y
  y_train <- y[,1]
  x_train <- y[,-1]
  y_test <- df_test[,1]
  x_test <- df_test[,-1]
  #ii. calcolo dei beta
  T     <- length(Lambda_samples)
  step  <- 25
  sel   <- seq(1, T, by = step) 
  Beta_samples <- matrix(NA, nrow = T/step, ncol = p-1)
  
  pb <- txtProgressBar(min = 0, max = T/step, style = 3)
  for (i in seq_along(sel)){
    Lambda_t   <- Lambda_samples[[i]]    # p × k_t
    psi_t <- sigma_samples[i, ]      # length p
    Omega_t  <- Lambda_t %*% t(Lambda_t) + diag(psi_t) 
    omega_xz_t <- Omega_t[2:p, 1]
    omega_xx_t <- Omega_t[-1,-1]
    Beta_samples[i, ] <- as.vector(solve(omega_xx_t) %*% omega_xz_t)
    setTxtProgressBar(pb, i) }
  
  close(pb)
  beta_map <- colMeans(Beta_samples) 
  
  #iii.calcolo predizioni e salvataggio mspe aape mape
  pred_test  <- Beta_samples %*% t(x_test)
  mean_pred_test  <- colMeans(pred_test)
  
  pred_train  <- Beta_samples %*% t(x_train)
  mean_pred_train  <- colMeans(pred_train)
  
  errs <- y_test-mean_pred_test
  mse_sc <- mean(errs^2)
  
  sd_y <- sd(labp[1,])
  mse_or <- mse_sc * sd_y^2
  
  cat("===== MODEL PERFORMANCE =====\n")
  cat(sprintf("MSE (scaled):    %.6f\n", mse_sc))
  cat(sprintf("MSE (original):  %.6f\n", mse_or))
  cat("=============================\n")

#PLOT LOGLIKELIHOOD
library(ggplot2)
library(gridExtra) 
  #loglike
  df_ll <- data.frame(
    Iterazione = 1:length(loglike),
    LogLik = loglike )
  p_ll <- ggplot(subset(df_ll,Iterazione >5000), aes(x = Iterazione, y = LogLik)) +
    geom_line(color = "steelblue", linewidth = 0.3) +
    theme_minimal(base_size = 13) +
    labs(x = "Iterazioni", y = "Log-verosimiglianza") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # andamento di k 
  df_k <- data.frame(
    Iterazione = 1:length(k_star_history),
    K_star = k_star_history)
  p_k <- ggplot(subset(df_k,Iterazione >5000), aes(x = Iterazione, y = K_star)) +
    geom_line(color = "darkgreen", linewidth = 0.3) +
    theme_minimal(base_size = 13) +
    labs(x = "Iterazioni", y = expression(k^"*")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  grid.arrange(p_ll, p_k, nrow=2)
