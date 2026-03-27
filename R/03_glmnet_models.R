# RIDGE/LASSO REGRESSION MULTI-OUTPUT, NB LASSO alpha = 1 , RIDGE alpha = 0
  ########################
  
  # install.packages("glmnet")
  library(glmnet)
  library(fds)
  library(MASS)
  data(labp)
  data(labc)
  data(nirp)
  data(nirc)
  
  # DATA
  y_train_df <- as.data.frame(t(labc))
  y_test_df  <- as.data.frame(t(labp))
  X_train_df <- as.data.frame(t(nirc$y))
  X_test_df  <- as.data.frame(t(nirp$y))
  
  # REMOVE obs 23 (outlier)
  y_train_df <- y_train_df[-23, ]
  X_train_df <- X_train_df[-23, ]
  
  #nirs p reduction 256 BROWN 1380-2400 , NB REMOVE THIS REDUCTION BLOCK TO OBTAIN THE MSPE ON P 700
  start_wl <- 1380
  end_wl   <- 2400
  i0 <- (start_wl - 1100) / 2 + 1
  i1 <- (end_wl   - 1100) / 2 + 1
  step_cols <- 4 / 2
  keep_idx <- seq(from = i0, to = i1, by = step_cols)
  X_train_df <-  X_train_df[, keep_idx]
  X_test_df <- X_test_df[, keep_idx]
  labels <- nirc$x
  labels <- labels[keep_idx]
  
  # standardizzazione
  X_train <- scale(X_train_df)
  X_test  <- scale(X_test_df)
  y_train <- scale(y_train_df)
  y_test  <- scale(y_test_df)
  
  # n component to predict
  n_resp <- ncol(y_train)
  mse_vec <- numeric(n_resp)
  
  ########################
  # 2. LOOP ON EVERY COMPONENT
  ########################
  for (j in 1:n_resp) {
    
    y_tr <- y_train[, j]
    y_te <- y_test[, j]
    
    # GRID LAMBDA
    lambda_grid <- 10^seq(5, -2, length = 100)
    
    # CROSS-VALIDATION 
    cv_ridge <- cv.glmnet(
      x = X_train,
      y = y_tr,
      alpha = 1,   # alpha = 1 LASSO / alpha = 0 RIDGE       
      lambda = lambda_grid,
      nfolds = 5
    )
    
    best_lambda <- cv_ridge$lambda.min
    
    # FINAL MODEL
    ridge_final <- glmnet(
      x = X_train,
      y = y_tr,
      alpha = 1,  # alpha = 1 LASSO / alpha = 0 RIDGE  
      lambda = best_lambda
    )
    
    # PREDICTION
    pred <- predict(ridge_final, newx = X_test, s = best_lambda)
    
    # MSE
    mse_vec[j] <- mean((y_te - pred)^2)
    
    cat("Costituente", j, "- lambda:", round(best_lambda, 5),
        "- MSE:", round(mse_vec[j], 6), "\n")
  }
  
  ########################
  # 3. STAMPA MSE COMPLESSIVI
  ########################
  cat("\nMSE per tutte le componenti:\n")
  print(mse_vec)
}
