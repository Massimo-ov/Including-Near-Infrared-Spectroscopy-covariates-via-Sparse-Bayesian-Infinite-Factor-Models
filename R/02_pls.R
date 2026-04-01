# PLS MODEL, RUN AFTER 01_preprocess_data.R
library(pls)
pls_model <- plsr(y_train ~ X_train,
                   ncomp = 30,
                   validation = "CV")
  
  # automatic choice of the component
  ncomp_pls <- which.min(RMSEP(pls_model)$val[1,1,])
  
  # Prediction
  pls_pred <- predict(pls_model,
                      newdata = X_test,
                      ncomp = ncomp_pls)
  loadings <- pls_model$loading.weights
  
  
  # MSE test
  mse_pls <- mean((y_test - pls_pred)^2)
