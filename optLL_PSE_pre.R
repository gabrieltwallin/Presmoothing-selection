optLL_PSE_pre <- function(XAfreq, YAfreq, X_range, Y_range, A_range, max_power = 5, lrt_pvalue_threshold = 0.05) {
  # Fit Poisson regression models for XAfreq and YAfreq
  XAresult <- presmooth_select(data = XAfreq, 
                               max_power = 5,
                           #    max_interaction_power = 3, 
                               lrt_pvalue_threshold = 0.05, 
                               x_var = "X", 
                               a_var = "A")
  
  YAresult <- presmooth_select(data = YAfreq, 
                               max_power = 5,
                           #    max_interaction_power = 3, 
                               lrt_pvalue_threshold = 0.05, 
                               x_var = "Y", 
                               a_var = "A")
  
  # Get the best XA and YA models based on AIC, BIC, and LRT
  best_XAmodels_names <- c(XAresult$aic_mod_interaction, 
                           XAresult$bic_mod_interaction, 
                           XAresult$lrt_mod_interaction)
  
  best_YAmodels_names <- c(YAresult$aic_mod_interaction, 
                           YAresult$bic_mod_interaction, 
                           YAresult$lrt_mod_interaction)
  
  best_XAmodels <- XAresult$interaction_models[best_XAmodels_names]
  best_YAmodels <- YAresult$interaction_models[best_YAmodels_names]
  
  
  
  # Initialize variables to store best models and SEE values
  best_see <- Inf
  best_XAmodel <- NULL
  best_YAmodel <- NULL
  best_XAmodel_name <- NULL
  best_YAmodel_name <- NULL
  
  # Loop through the best XA and YA models
  for (XAmodel_name in names(best_XAmodels)) {
    XAmodel <- best_XAmodels[[XAmodel_name]]
    for (YAmodel_name in names(best_YAmodels)) {
      YAmodel <- best_YAmodels[[YAmodel_name]]
      
      # Run kequate with the current XA and YA models
      keq_result <- kequate("NEAT_PSE", X_range, Y_range, XAmodel, YAmodel)
      
      # Calculate average SEE
      avg_see <- mean(keq_result@PRE$PREYx)  # PREAx
      
      # Update best models and SEE if current models have a lower average SEE
      if (avg_see < best_see) {
        best_see <- avg_see
        best_XAmodel <- XAmodel
        best_YAmodel <- YAmodel
        best_XAmodel_name <- XAmodel_name
        best_YAmodel_name <- YAmodel_name
        best_XAmodel_criteria <- ifelse(XAmodel_name == XAresult$aic_mod_interaction, "AIC", ifelse(XAmodel_name == XAresult$bic_mod_interaction, "BIC", "LRT"))
        best_YAmodel_criteria <- ifelse(YAmodel_name == YAresult$aic_mod_interaction, "AIC", ifelse(YAmodel_name == YAresult$bic_mod_interaction, "BIC", "LRT"))
      }
    }
  }
  
  # Return the best models, their names, and the lowest average SEE value
  return(list(best_XAmodel = best_XAmodel, 
              best_YAmodel = best_YAmodel, 
              best_XAmodel_name = best_XAmodel_name, 
              best_YAmodel_name = best_YAmodel_name,
              best_XAmodel_criteria = best_XAmodel_criteria,
              best_YAmodel_criteria = best_YAmodel_criteria,
              best_see = best_see))
}
