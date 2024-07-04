 

opt_irt_PSE_see <- function(data, lrt_pvalue_threshold = 0.05) {
  models <- list()
  
  # Fit all possible models: 1-PL, 2-PL, 3-PL
  models[["1PL"]] <- ltm::ltm(data ~ z1) # 1-PL (Rasch model)
  models[["2PL"]] <- ltm::ltm(data ~ z1, IRT.param = TRUE) # 2-PL
  models[["3PL"]] <- ltm::tpm(data, type = c("latent.trait"), IRT.param = TRUE, start.val = "random") # 3-PL
  
  # Extract AIC, BIC, and log-likelihood
  model_info <- sapply(models, function(x) {
    c(AIC = AIC(x),
      BIC = BIC(x),
      LogLik = logLik(x))
  })
  
  # Perform LRT for nested models using log-likelihood values
  lrt_2PL_vs_1PL <- if ("2PL" %in% names(models) && "1PL" %in% names(models)) {
    chisq <- -2 * (model_info["LogLik", "1PL"] - model_info["LogLik", "2PL"])
    pchisq(chisq, df = 1, lower.tail = FALSE)
  } else {
    NA
  }
  
  lrt_3PL_vs_2PL <- if ("3PL" %in% names(models) && "2PL" %in% names(models)) {
    chisq <- -2 * (model_info["LogLik", "2PL"] - model_info["LogLik", "3PL"])
    pchisq(chisq, df = 1, lower.tail = FALSE)
  } else {
    NA
  }
  
  # Find best models according to AIC, BIC, and LRT
  best_aic_model_name <- names(which.min(model_info["AIC", ]))
  best_bic_model_name <- names(which.min(model_info["BIC", ]))
  best_lrt_model_name <- if (!is.na(lrt_2PL_vs_1PL) && lrt_2PL_vs_1PL < lrt_pvalue_threshold) {
    if (!is.na(lrt_3PL_vs_2PL) && lrt_3PL_vs_2PL < lrt_pvalue_threshold) "3PL" else "2PL"
  } else {
    "1PL"
  }
  
  return(list(models = models, model_info = model_info, 
              aic_mod = best_aic_model_name, bic_mod = best_bic_model_name, 
              lrt_mod = best_lrt_model_name, pvalue_lrt = c(lrt_2PL_vs_1PL, lrt_3PL_vs_2PL),
              best_aic_mod = models[[best_aic_model_name]],
              best_bic_mod = models[[best_bic_model_name]],
              best_lrt_mod = models[[best_lrt_model_name]]))
}

