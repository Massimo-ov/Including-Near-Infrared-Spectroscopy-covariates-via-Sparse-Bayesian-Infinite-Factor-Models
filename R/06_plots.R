#DOTPLOT FOR MSPE
#MSPE VALOUR
vars <- c("Fat", "Sugar", "Water", "Flour")
  
  PLSp700   <- c(0.094, 0.056, 0.015, 0.071)
  RIDGEp700 <- c(0.046, 0.050, 0.019, 0.074)
  LASSOp700 <- c(0.0132, 0.078, 0.098, 0.068)
  SBIFMp700 <- c(0.0127, 0.056, 0.042, 0.076)
  
  PLSp256   <- c(0.023, 0.060, 0.050, 0.063)
  LASSOp256 <- c(0.0105, 0.077, 0.055, 0.065)
  RIDGEp256 <- c(0.020, 0.055, 0.066, 0.054)
  SBIFMp256 <- c(0.0133, 0.055, 0.042, 0.074)
  
  mse_mat <- rbind(
    PLSp700, LASSOp700, RIDGEp700, SBIFMp700,
    PLSp256, RIDGEp256, LASSOp256, SBIFMp256
  )

library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  df <- data.frame(
    Variable = vars,
    PLSp700, RIDGEp700, LASSOp700, SBIFMp700,
    PLSp256, RIDGEp256, LASSOp256, SBIFMp256
  )
  
  df_long <- df %>%
    pivot_longer(-Variable,
                 names_to = "Model",
                 values_to = "MSE") %>%
    mutate(
      Variable  = factor(Variable, levels = vars),
      Method    = gsub("p[0-9]+", "", Model),
      Dimension = ifelse(grepl("700", Model), "p = 700", "p = 256")
    )
  
  # Coordinate
  x_right <- length(vars)
  y_min   <- min(df_long$MSE)
  y_max   <- max(df_long$MSE)
  y_range <- y_max - y_min
  
  # Spostiamo più a destra
  x_leg <- x_right +0.25
  y_leg <- y_min + 0.18 * y_range
  
  p <- ggplot(df_long,
              aes(x = Variable,
                  y = MSE,
                  color = Method,
                  shape = Dimension,
                  group = interaction(Method, Dimension))) +
    
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8) +
    
    scale_shape_manual(values = c(
      "p = 700" = 16,
      "p = 256" = 17
    )) +
    
    scale_color_manual(values = c(
      PLS   = "orange",
      RIDGE = "red",
      LASSO = "green",
      SBIFM = "black"
    )) +
    
    labs(
      x = "Response Variable",
      y = "MSE"
    ) +
    
    # ---- LINEETTE METODI ----
  annotate("segment", x = x_leg - 0.55, xend = x_leg - 0.35,
           y = y_leg + 0.15*y_range,
           yend = y_leg + 0.15*y_range,
           color = "orange", linewidth = 0.9) +
    
    annotate("segment", x = x_leg - 0.55, xend = x_leg - 0.35,
             y = y_leg + 0.10*y_range,
             yend = y_leg + 0.10*y_range,
             color = "red", linewidth = 0.9) +
    
    annotate("segment", x = x_leg - 0.55, xend = x_leg - 0.35,
             y = y_leg + 0.05*y_range,
             yend = y_leg + 0.05*y_range,
             color = "green", linewidth = 0.9) +
    
    annotate("segment", x = x_leg - 0.55, xend = x_leg - 0.35,
             y = y_leg,
             yend = y_leg,
             color = "black", linewidth = 0.9) +
    
    # ---- TESTI METODI ----
  annotate("text", x = x_leg - 0.30, y = y_leg + 0.15*y_range,
           label = "PLS", hjust = 0, size = 3.2) +
    annotate("text", x = x_leg - 0.30, y = y_leg + 0.10*y_range,
             label = "Ridge", hjust = 0, size = 3.2) +
    annotate("text", x = x_leg - 0.30, y = y_leg + 0.05*y_range,
             label = "Lasso", hjust = 0, size = 3.2) +
    annotate("text", x = x_leg - 0.30, y = y_leg,
             label = "SBIFM", hjust = 0, size = 3.2) +
    
    # ---- SIMBOLI DIMENSIONE ----
  annotate("point", x = x_leg - 0.48,
           y = y_leg - 0.05*y_range,
           shape = 16, size = 2.5) +
    
    annotate("point", x = x_leg - 0.48,
             y = y_leg - 0.10*y_range,
             shape = 17, size = 2.5) +
    
    annotate("text", x = x_leg - 0.30,
             y = y_leg - 0.05*y_range,
             label = "p = 700", hjust = 0, size = 3.2) +
    
    annotate("text", x = x_leg - 0.30,
             y = y_leg - 0.10*y_range,
             label = "p = 256", hjust = 0, size = 3.2) +
    
    theme_classic() +
    
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text  = element_text(size = 10),
      #axis.text.y = element_text(angle = 90, hjust = 0.5),
      axis.line = element_line(linewidth = 0.6)
    )
print(p)


#PLOT BETA SBIFM, BETA RIDGE, MOST PREDICTIVE FACTOR 
library(ggplot2)
  
  ############################
  # 1. FATTORE SELEZIONATO
  ############################
  beta_mean  <- colMeans(Beta_samples)
  beta_low   <- apply(Beta_samples, 2, quantile, 0.025)
  beta_high  <- apply(Beta_samples, 2, quantile, 0.975)
  
  h_star <- 3  # <-- metti quello giusto
  
  factor_mean <- pervar_mean[-1, h_star]
  factor_low  <- pervar_ci_low[-1, h_star]
  factor_high <- pervar_ci_high[-1, h_star]
  
  ############################
  # 2. DATAFRAME (SENZA SCALING)
  ############################
  
  df_plot <- data.frame(
    wavelength  = labels,
    beta_mean   = beta_mean,
    beta_low    = beta_low,
    beta_high   = beta_high,
    beta_ridge  = beta_ridge,
    factor_mean = factor_mean,
    factor_low  = factor_low,
    factor_high = factor_high
  )
  
  ############################
  # 3. PLOT
  ############################
p <- ggplot(df_plot, aes(x = wavelength)) +
    
    geom_ribbon(aes(ymin = beta_low, ymax = beta_high),
                fill = "grey70", alpha = 0.4) +
    
    geom_ribbon(aes(ymin = factor_low, ymax = factor_high),
                fill = "purple", alpha = 0.15) +
    
    geom_line(aes(y = beta_mean, color = "Beta Sbifm"),
              linewidth = 0.9) +
    
    geom_line(aes(y = beta_ridge, color = "Beta Ridge"),
              linewidth = 0.9) +
    
    geom_line(aes(y = factor_mean, color = "MIF"),
              linewidth = 0.9) +
    
    scale_color_manual(
      values = c("Beta Sbifm" = "black",
                 "Beta Ridge" = "red",
                 "MIF" = "purple"),
      name = NULL
    ) +
    
    labs(
      x = "Wavelength",
      y = NULL
    ) +
    
    theme_classic() +
    
    theme(
      #legend.position = c(0.82, 0.88),
      legend.position = "none",
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 12),
      legend.key = element_blank()
    )
  print(p)
