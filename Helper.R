# ENHANCED COMPREHENSIVE ANALYSIS
# Includes ALL metrics requested by advisor:
# - Cox coefficients & beta parameters
# - Brier scores
# - Survival curves
# - Time-to-event estimates
# - RMST (integral under survival curve)
# - Predicted vs observed event times
# - Complete figures package

library(JMbayes2)
library(dplyr)
library(nlme)
library(survival)
library(pec)
library(survminer)
library(ggplot2)
library(gridExtra)

select <- dplyr::select
summarize <- dplyr::summarize

cat("============================================================\n")
cat("ENHANCED COMPREHENSIVE THESIS ANALYSIS\n")
cat("With ALL requested metrics and figures\n")
cat("============================================================\n")

if (!dir.exists("thesis_figures")) {
  dir.create("thesis_figures")
}

# ENHANCED METRICS FUNCTION
compute_all_metrics_with_figures <- function(coxFit, surv_data, long_data, 
                                             jointFit=NULL, name="Method") {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("EVALUATING: %s\n", name))
  cat(sprintf("========================================\n"))
  
  actual_patients <- length(unique(surv_data$PTID))
  actual_events <- sum(surv_data$event)
  
  cat(sprintf("  Patients: %d\n", actual_patients))
  cat(sprintf("  Events: %d (%.1f%%)\n", actual_events, 100*actual_events/actual_patients))
  
  method_dir <- file.path("thesis_figures", gsub(" ", "_", name))
  if(!dir.exists(method_dir)) {
    dir.create(method_dir, recursive = TRUE)
  }
  
  # INITIALIZE results with ALL fields
  results <- list(
    name = name,
    n_patients = actual_patients,
    n_events = actual_events,
    event_rate = actual_events / actual_patients,
    cindex = NA,
    cindex_ci_lower = NA,
    cindex_ci_upper = NA,
    brier_1yr = NA,
    brier_3yr = NA,
    brier_5yr = NA,
    calibration_slope = NA,
    calibration_intercept = NA,
    rmst_tau = NA,
    rmst_low = NA,
    rmst_low_lower = NA,
    rmst_low_upper = NA,
    rmst_high = NA,
    rmst_high_lower = NA,
    rmst_high_upper = NA,
    rmst_diff = NA,
    rmst_diff_lower = NA,
    rmst_diff_upper = NA,
    rmst_pval = NA,
    mae_time = NA,
    rmse_time = NA,
    corr_time = NA,
    median_error = NA,
    n_significant = NA,
    km_high_median = NA,
    km_low_median = NA
  )
  
  # 1. C-INDEX  
  cat("\n1. DISCRIMINATION\n")
  cindex <- summary(coxFit)$concordance[1]
  results$cindex <- cindex
  cat(sprintf("   C-index: %.4f\n", cindex))
  
  # 2. COX COEFFICIENTS  
  cat("\n2. COX MODEL PARAMETERS\n")
  
  coef_summary <- summary(coxFit)$coefficients
  coef_df <- as.data.frame(coef_summary)
  coef_df$Variable <- rownames(coef_df)
  coef_df$HR <- exp(coef_df$coef)
  coef_df$HR_lower <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
  coef_df$HR_upper <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
  
  coef_df_sorted <- coef_df[order(coef_df$`Pr(>|z|)`), ]
  cat("   Top 10 predictors:\n")
  print(head(coef_df_sorted[, c("Variable", "coef", "HR", "Pr(>|z|)")], 10))
  
  write.csv(coef_df, 
            file.path(method_dir, paste0(gsub(" ", "_", name), "_cox_coefficients.csv")),
            row.names = FALSE)
  cat(sprintf("    Saved: cox_coefficients.csv\n"))
  
  # 3. FOREST PLOT
  cat("\n3. FOREST PLOT\n")
  
  tryCatch({
    sig_coefs <- coef_df[coef_df$`Pr(>|z|)` < 0.10, ]
    
    if(nrow(sig_coefs) > 0 && nrow(sig_coefs) <= 25) {
      
      sig_coefs <- sig_coefs[order(sig_coefs$HR, decreasing = TRUE), ]
      
      sig_coefs$HR_plot <- pmin(pmax(sig_coefs$HR, 0.001), 1000)
      sig_coefs$HR_lower_plot <- pmin(pmax(sig_coefs$HR_lower, 0.001), 1000)
      sig_coefs$HR_upper_plot <- pmin(pmax(sig_coefs$HR_upper, 0.001), 1000)
      
      forest_plot <- ggplot(sig_coefs, aes(x = HR_plot, y = reorder(Variable, HR_plot))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
        geom_errorbarh(aes(xmin = HR_lower_plot, xmax = HR_upper_plot), 
                       height = 0.2, linewidth = 0.6) +
        geom_point(size = 3, color = "#0072B2") +
        scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)) +
        labs(title = paste("Hazard Ratios:", name),
             x = "HR (log scale)", y = "") +
        theme_minimal(base_size = 13)
      
      ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_forest_plot.png")),
             forest_plot, width = 10, height = max(5, nrow(sig_coefs) * 0.35), dpi = 300)
      
      cat(sprintf("    Forest plot saved (%d significant)\n", nrow(sig_coefs)))
      results$n_significant <- nrow(sig_coefs)
    } else {
      cat("   No significant predictors (p < 0.10)\n")
      results$n_significant <- 0
    }
  }, error = function(e) {
    cat(sprintf("    Forest plot error: %s\n", e$message))
    results$n_significant <- 0
  })
  
  # 4. BASELINE SURVIVAL
  cat("\n4. BASELINE SURVIVAL\n")
  
  tryCatch({
    baseline_surv <- basehaz(coxFit, centered = TRUE)
    
    gg_hazard <- ggplot(baseline_surv, aes(x = time, y = hazard)) +
      geom_step(color = "#0072B2", linewidth = 1.2) +
      labs(title = paste("Baseline Hazard:", name),
           x = "Time (Years)", y = "Cumulative Hazard") +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_baseline_hazard.png")),
           gg_hazard, width = 8, height = 6, dpi = 300)
    
    baseline_surv$surv_prob <- exp(-baseline_surv$hazard)
    
    gg_baseline_surv <- ggplot(baseline_surv, aes(x = time, y = surv_prob)) +
      geom_step(color = "#0072B2", linewidth = 1.2) +
      ylim(0, 1) +
      labs(title = paste("Baseline Survival:", name),
           x = "Time (Years)", y = "Survival Probability") +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_baseline_survival.png")),
           gg_baseline_surv, width = 8, height = 6, dpi = 300)
    
    cat("    Baseline functions saved\n")
    
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
  })
  
  # 5. BRIER SCORE
  cat("\n5. BRIER SCORE\n")
  
  tryCatch({
    pred_times <- c(1, 3, 5)
    pec_obj <- pec(
      list("Model" = coxFit),
      formula = Surv(time, event) ~ 1,
      data = surv_data,
      times = pred_times,
      exact = FALSE,
      cens.model = "marginal",
      verbose = FALSE
    )
    
    brier_scores <- pec_obj$AppErr$Model
    results$brier_1yr <- brier_scores[1]
    results$brier_3yr <- brier_scores[2]
    results$brier_5yr <- brier_scores[3]
    
    cat(sprintf("   Brier (1yr): %.4f\n", results$brier_1yr))
    cat(sprintf("   Brier (3yr): %.4f\n", results$brier_3yr))
    cat(sprintf("   Brier (5yr): %.4f\n", results$brier_5yr))
    
    brier_df <- data.frame(Time = pred_times, Brier = brier_scores)
    
    gg_brier <- ggplot(brier_df, aes(x = Time, y = Brier)) +
      geom_line(color = "#0072B2", linewidth = 1.3) +
      geom_point(size = 4, color = "#0072B2") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "darkgreen") +
      labs(title = paste("Brier Score:", name),
           x = "Time (Years)", y = "Brier Score") +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_brier_scores.png")),
           gg_brier, width = 8, height = 6, dpi = 300)
    
    cat("    Brier plot saved\n")
    
  }, error = function(e) {
    cat(sprintf("    Brier error: %s\n", e$message))
  })
  
  # 6. KAPLAN-MEIER + 7. RMST (FIXED WITH ADAPTIVE TAU)
  cat("\n6. RISK STRATIFICATION\n")
  
  tryCatch({
    risk_scores <- predict(coxFit, newdata = surv_data, type = "lp")
    surv_data$risk <- risk_scores
    surv_data$risk_group <- ifelse(
      surv_data$risk >= median(surv_data$risk),
      "High Risk", "Low Risk"
    )
    
    fit_stratified <- survfit(Surv(time, event) ~ risk_group, data = surv_data)
    
    km_plot <- ggsurvplot(
      fit_stratified,
      data = surv_data,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = TRUE,
      palette = c("#D55E00", "#0072B2"),
      legend.title = "Risk Group",
      title = paste("Kaplan-Meier:", name),
      xlab = "Time (Years)",
      ylab = "Progression-Free Probability"
    )
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_KM_stratified.png")),
           km_plot$plot, width = 10, height = 8, dpi = 300)
    
    km_summary <- summary(fit_stratified)$table
    results$km_high_median <- km_summary["risk_group=High Risk", "median"]
    results$km_low_median <- km_summary["risk_group=Low Risk", "median"]
    
    cat(sprintf("   High-risk median: %.2f years\n", results$km_high_median))
    cat(sprintf("   Low-risk median:  %.2f years\n", results$km_low_median))
    cat("    KM curves saved\n")
    
    # 7. RMST ANALYSIS (ADAPTIVE TAU)
    cat("\n7. RMST ANALYSIS\n")
    
    tryCatch({
      surv_data$arm <- as.numeric(surv_data$risk_group == "Low Risk")
      
      max_time_high <- max(surv_data$time[surv_data$risk_group == "High Risk"], na.rm = TRUE)
      max_time_low <- max(surv_data$time[surv_data$risk_group == "Low Risk"], na.rm = TRUE)
      max_tau <- min(max_time_high, max_time_low)
      tau <- max_tau * 0.9
      
      cat(sprintf("   Maximum follow-up - High Risk: %.2f years\n", max_time_high))
      cat(sprintf("   Maximum follow-up - Low Risk:  %.2f years\n", max_time_low))
      cat(sprintf("   Selected tau: %.2f years (90%% of %.2f)\n", tau, max_tau))
      
      rmst_comparison <- rmst2(
        time = surv_data$time,
        status = surv_data$event,
        arm = surv_data$arm,
        tau = tau
      )
      
      results$rmst_tau <- tau
      results$rmst_low <- rmst_comparison$RMST.arm1$rmst["Est."]
      results$rmst_low_lower <- rmst_comparison$RMST.arm1$rmst["lower .95"]
      results$rmst_low_upper <- rmst_comparison$RMST.arm1$rmst["upper .95"]
      
      results$rmst_high <- rmst_comparison$RMST.arm0$rmst["Est."]
      results$rmst_high_lower <- rmst_comparison$RMST.arm0$rmst["lower .95"]
      results$rmst_high_upper <- rmst_comparison$RMST.arm0$rmst["upper .95"]
      
      results$rmst_diff <- rmst_comparison$unadjusted.result[1, "Est."]
      results$rmst_diff_lower <- rmst_comparison$unadjusted.result[1, "lower .95"]
      results$rmst_diff_upper <- rmst_comparison$unadjusted.result[1, "upper .95"]
      results$rmst_pval <- rmst_comparison$unadjusted.result[1, "p"]
      
      cat(sprintf("   Low-risk RMST:  %.3f years (95%% CI: %.2f, %.2f)\n", 
                  results$rmst_low, results$rmst_low_lower, results$rmst_low_upper))
      cat(sprintf("   High-risk RMST: %.3f years (95%% CI: %.2f, %.2f)\n", 
                  results$rmst_high, results$rmst_high_lower, results$rmst_high_upper))
      cat(sprintf("   Difference:     %.3f years (95%% CI: %.2f, %.2f), p = %.4f\n", 
                  results$rmst_diff, results$rmst_diff_lower, results$rmst_diff_upper,
                  results$rmst_pval))
      
      rmst_table <- data.frame(
        Risk_Group = c("Low Risk", "High Risk", "Difference"),
        RMST_Years = c(results$rmst_low, results$rmst_high, results$rmst_diff),
        Lower_CI = c(results$rmst_low_lower, results$rmst_high_lower, results$rmst_diff_lower),
        Upper_CI = c(results$rmst_low_upper, results$rmst_high_upper, results$rmst_diff_upper),
        P_Value = c(NA, NA, results$rmst_pval),
        Tau = c(tau, tau, tau)
      )
      
      write.csv(rmst_table,
                file.path(method_dir, paste0(gsub(" ", "_", name), "_rmst_analysis.csv")),
                row.names = FALSE)
      
      cat("    RMST analysis saved\n")
      
    }, error = function(e) {
      cat("    RMST analysis failed:", conditionMessage(e), "\n")
      results$rmst_error <- conditionMessage(e)
    })
    
  }, error = function(e) {
    cat(sprintf("    KM error: %s\n", e$message))
  })
  
  # 8. PREDICTED vs OBSERVED TIMES
  cat("\n8. PREDICTED vs OBSERVED TIMES\n")
  
  tryCatch({
    converters <- surv_data[surv_data$event == 1, ]
    
    if(nrow(converters) >= 5) {
      
      converter_risks <- predict(coxFit, newdata = converters, type = "lp")
      
      risk_range <- max(converter_risks) - min(converter_risks)
      time_range <- max(converters$time) - min(converters$time)
      
      if(risk_range > 1e-6) {
        norm_risks <- (converter_risks - min(converter_risks)) / risk_range
        predicted_times <- max(converters$time) - norm_risks * time_range
      } else {
        predicted_times <- rep(mean(converters$time), nrow(converters))
      }
      
      observed_times <- converters$time
      
      errors <- abs(predicted_times - observed_times)
      results$mae_time <- mean(errors)
      results$rmse_time <- sqrt(mean((predicted_times - observed_times)^2))
      results$corr_time <- cor(predicted_times, observed_times, method = "pearson")
      results$median_error <- median(errors)
      
      cat(sprintf("   MAE:  %.3f years (%.1f months)\n", results$mae_time, results$mae_time * 12))
      cat(sprintf("   RMSE: %.3f years\n", results$rmse_time))
      cat(sprintf("   Correlation: %.3f\n", results$corr_time))
      
      pred_obs_df <- data.frame(
        Predicted = predicted_times,
        Observed = observed_times,
        PTID = converters$PTID,
        Error = errors
      )
      
      gg_pred_obs <- ggplot(pred_obs_df, aes(x = Observed, y = Predicted)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
        geom_point(alpha = 0.7, size = 3.5, color = "#0072B2") +
        geom_smooth(method = "lm", se = TRUE, color = "#D55E00", alpha = 0.2) +
        labs(title = paste("Predicted vs Observed:", name),
             x = "Observed (Years)", y = "Predicted (Years)") +
        theme_minimal(base_size = 13) +
        coord_equal()
      
      ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_predicted_vs_observed.png")),
             gg_pred_obs, width = 8, height = 8, dpi = 300)
      
      write.csv(pred_obs_df,
                file.path(method_dir, paste0(gsub(" ", "_", name), "_time_predictions.csv")),
                row.names = FALSE)
      
      cat("    Predictions saved\n")
      
    }
    
  }, error = function(e) {
    cat(sprintf("    Prediction error: %s\n", e$message))
  })
  
  # 9. CALIBRATION PLOT
  cat("\n9. CALIBRATION PLOT\n")
  
  tryCatch({
    risk_scores <- predict(coxFit, newdata = surv_data, type = "lp")
    
    # CHECK IF CALIBRATION MODEL EXISTS (for optimized thesis method)
    if(!is.null(coxFit$calibration_model)) {
      cat("   Using calibrated predictions\n")
      # Use calibrated predictions
      predicted_surv_3yr <- coxFit$predict_calibrated(surv_data)
    } else {
      # Original method (for other models)
      risk_range <- max(risk_scores) - min(risk_scores)
      if(risk_range > 1e-6) {
        predicted_surv_3yr <- 1 - (risk_scores - min(risk_scores)) / risk_range
      } else {
        predicted_surv_3yr <- rep(0.5, length(risk_scores))
      }
    }
    
    # Observed: 1 if survived past 3 years
    observed_surv_3yr <- as.numeric(
      surv_data$time > 3 | (surv_data$time <= 3 & surv_data$event == 0)
    )
    
    # Remove invalid
    valid_idx <- !is.na(predicted_surv_3yr) & !is.na(observed_surv_3yr)
    pred_valid <- predicted_surv_3yr[valid_idx]
    obs_valid <- observed_surv_3yr[valid_idx]
    
    if(length(pred_valid) >= 10) {
      
      # Create deciles
      decile_breaks <- quantile(pred_valid, probs = seq(0, 1, 0.1), na.rm = TRUE)
      deciles <- cut(pred_valid, breaks = decile_breaks, include.lowest = TRUE, labels = 1:10)
      
      cal_df <- data.frame(
        predicted = pred_valid,
        observed = obs_valid,
        decile = deciles
      )
      
      cal_data <- cal_df %>%
        group_by(decile) %>%
        summarize(
          predicted = mean(predicted, na.rm = TRUE),
          observed = mean(observed, na.rm = TRUE),
          n = n(),
          .groups = "drop"
        ) %>%
        filter(!is.na(predicted), !is.na(observed))
      
      if(nrow(cal_data) >= 3) {
        cal_model <- lm(observed ~ predicted, data = cal_data)
        results$calibration_slope <- coef(cal_model)[2]
        results$calibration_intercept <- coef(cal_model)[1]
        
        cat(sprintf("   Calibration slope: %.3f (ideal = 1.0)\n", results$calibration_slope))
        cat(sprintf("   Intercept:         %.3f (ideal = 0.0)\n", results$calibration_intercept))
        
        gg_calibration <- ggplot(cal_data, aes(x = predicted, y = observed)) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
          geom_smooth(method = "lm", se = TRUE, color = "#D55E00", alpha = 0.2) +
          geom_point(aes(size = n), color = "#0072B2", alpha = 0.8) +
          xlim(0, 1) + ylim(0, 1) +
          labs(title = paste("Calibration (3-year):", name),
               x = "Predicted Survival", y = "Observed Survival") +
          theme_minimal(base_size = 13) +
          coord_equal()
        
        ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_calibration_plot.png")),
               gg_calibration, width = 8, height = 8, dpi = 300)
        
        cat("    Calibration saved\n")
      }
    }
    
  }, error = function(e) {
    cat(sprintf("    Calibration error: %s\n", e$message))
  })
  
  # 10. BOOTSTRAP CI
  cat("\n10. BOOTSTRAP CI\n")
  
  tryCatch({
    boot_cindex <- function(data, indices) {
      boot_data <- data[indices, ]
      
      boot_cox <- tryCatch({
        coxph(coxFit$formula, data = boot_data, x = TRUE)
      }, error = function(e) NA)
      
      if(length(boot_cox) == 1 && is.na(boot_cox)) return(NA)
      
      return(summary(boot_cox)$concordance[1])
    }
    
    cat("   Running bootstrap (100 iterations for speed)...\n")
    boot_results <- boot(data = surv_data, statistic = boot_cindex, R = 100)
    
    valid_boot <- boot_results$t[!is.na(boot_results$t)]
    
    if(length(valid_boot) >= 10) {
      boot_ci <- boot.ci(boot_results, type = "perc")
      results$cindex_ci_lower <- boot_ci$percent[4]
      results$cindex_ci_upper <- boot_ci$percent[5]
      
      cat(sprintf("   C-index: %.4f (95%% CI: %.4f-%.4f)\n", 
                  cindex, results$cindex_ci_lower, results$cindex_ci_upper))
      cat("    Bootstrap complete\n")
      
      saveRDS(boot_results, 
              file.path(method_dir, paste0(gsub(" ", "_", name), "_bootstrap.rds")))
    }
    
  }, error = function(e) {
    cat(sprintf("    Bootstrap error: %s\n", e$message))
  })
  
  # SAVE TO MASTER RESULTS FILE
  # Convert results list to data frame
  results_df <- as.data.frame(results, stringsAsFactors = FALSE)
  
  # Define master results file path
  master_file <- "thesis_figures/ALL_METHODS_RESULTS.csv"
  
  # Append to master file (create if doesn't exist)
  if(file.exists(master_file)) {
    existing_results <- read.csv(master_file, stringsAsFactors = FALSE)
    all_results <- rbind(existing_results, results_df)
  } else {
    all_results <- results_df
  }
  
  write.csv(all_results, master_file, row.names = FALSE)
  cat(sprintf("\n Results appended to: %s\n", master_file))
  
  # SUMMARY
  cat("\n========================================\n")
  cat(sprintf("SUMMARY: %s\n", name))
  cat("========================================\n")
  cat(sprintf("  C-index:     %.4f", results$cindex))
  if(!is.na(results$cindex_ci_lower)) {
    cat(sprintf(" (%.4f, %.4f)", results$cindex_ci_lower, results$cindex_ci_upper))
  }
  cat("\n")
  if(!is.na(results$brier_3yr)) cat(sprintf("  Brier (3yr): %.4f\n", results$brier_3yr))
  if(!is.na(results$calibration_slope)) cat(sprintf("  Cal. slope:  %.3f\n", results$calibration_slope))
  if(!is.na(results$rmst_diff)) cat(sprintf("  RMST diff:   %.2f years (p=%.4f)\n", results$rmst_diff, results$rmst_pval))
  if(!is.na(results$mae_time)) cat(sprintf("  MAE:         %.2f years\n", results$mae_time))
  cat("========================================\n")
  
  return(results)
}

# SURVIVAL ANALYSIS PLOTS - MAE-JM Method
plot_cox_survival <- function(coxFit, method_dir) {
  tryCatch({
    baseline_haz <- basehaz(coxFit, centered = TRUE)
    baseline_haz$survival <- exp(-baseline_haz$hazard)
    
    gg <- ggplot(baseline_haz, aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1.3) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
      ylim(0, 1) +
      labs(title = "Cox Model Baseline Survival Curve",
           subtitle = "Estimated survival function from Cox proportional hazards model",
           x = "Time (Years)", y = "Survival Probability S(t)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "gray30"))
    
    ggsave(file.path(method_dir, "cox_baseline_survival_curve.png"), gg, width = 10, height = 7, dpi = 300)
    write.csv(baseline_haz, file.path(method_dir, "cox_survival_function.csv"), row.names = FALSE)
    
    median_idx <- which.min(abs(baseline_haz$survival - 0.5))
    cat(sprintf("   Median survival time: %.2f years\n", baseline_haz$time[median_idx]))
    
    return(baseline_haz)  # Return for reuse in other plots
  }, error = function(e) cat(sprintf("    Cox survival error: %s\n", e$message)))
}


plot_hazard_functions <- function(coxFit, method_dir) {
  tryCatch({
    baseline_haz <- basehaz(coxFit, centered = TRUE)
    
    baseline_haz$hazard_rate <- c(diff(baseline_haz$hazard) / diff(baseline_haz$time), NA)
    baseline_haz_clean <- baseline_haz[!is.na(baseline_haz$hazard_rate), ]
    q99 <- quantile(baseline_haz_clean$hazard_rate, 0.99, na.rm = TRUE)
    baseline_haz_clean$hazard_rate_capped <- pmin(baseline_haz_clean$hazard_rate, q99)
    
    gg_rate <- ggplot(baseline_haz_clean, aes(x = time, y = hazard_rate_capped)) +
      geom_line(color = "#D55E00", linewidth = 1.2, alpha = 0.8) +
      geom_smooth(se = TRUE, color = "#0072B2", linewidth = 1.3, alpha = 0.3) +
      labs(title = "Instantaneous Hazard Rate h(t)",
           subtitle = "Rate of progression events over time (smoothed)",
           x = "Time (Years)", y = "Hazard Rate") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "gray30"))
    
    ggsave(file.path(method_dir, "hazard_rate_function.png"), gg_rate, width = 10, height = 7, dpi = 300)
    
    gg_cum <- ggplot(baseline_haz, aes(x = time, y = hazard)) +
      geom_step(color = "#0072B2", linewidth = 1.3) +
      labs(title = "Cumulative Hazard Function H(t)",
           subtitle = "Integral of hazard rate over time",
           x = "Time (Years)", y = "Cumulative Hazard H(t)") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "gray30"))
    
    ggsave(file.path(method_dir, "cumulative_hazard_function.png"), gg_cum, width = 10, height = 7, dpi = 300)
    cat("    Hazard functions saved\n")
    
    return(list(baseline_haz = baseline_haz, baseline_haz_clean = baseline_haz_clean))
  }, error = function(e) cat(sprintf("    Hazard function error: %s\n", e$message)))
}


plot_kaplan_meier <- function(surv_data, method_dir) {
  tryCatch({
    km_overall <- survfit(Surv(time, event) ~ 1, data = surv_data)
    km_summary <- summary(km_overall)
    
    km_data <- data.frame(time = km_summary$time, n_risk = km_summary$n.risk,
                          n_event = km_summary$n.event, survival = km_summary$surv,
                          std_err = km_summary$std.err, lower_ci = km_summary$lower,
                          upper_ci = km_summary$upper)
    write.csv(km_data, file.path(method_dir, "kaplan_meier_table.csv"), row.names = FALSE)
    
    gg_km <- ggsurvplot(km_overall, data = surv_data, conf.int = TRUE,
                        risk.table = TRUE, risk.table.height = 0.25,
                        ggtheme = theme_minimal(base_size = 13), palette = "#0072B2",
                        title = "Kaplan-Meier Survival Curve - Overall Cohort",
                        xlab = "Time (Years)", ylab = "Progression-Free Probability",
                        legend = "none", break.time.by = 1)
    
    png(file.path(method_dir, "kaplan_meier_overall.png"), width = 10, height = 8, units = "in", res = 300)
    print(gg_km)
    dev.off()
    
    timepoints <- c(1, 2, 3, 5)
    surv_at_times <- summary(km_overall, times = timepoints)
    surv_table <- data.frame(Time_Years = timepoints, Survival = surv_at_times$surv,
                             Lower_CI = surv_at_times$lower, Upper_CI = surv_at_times$upper,
                             N_Risk = surv_at_times$n.risk)
    write.csv(surv_table, file.path(method_dir, "survival_at_timepoints.csv"), row.names = FALSE)
    
    km_stats <- summary(km_overall)$table
    cat(sprintf("   Median survival: %.2f years\n", km_stats["median"]))
    for(i in 1:nrow(surv_table)) {
      cat(sprintf("     %d year: %.1f%% (95%% CI: %.1f%%-%.1f%%)\n",
                  surv_table$Time_Years[i], surv_table$Survival[i] * 100,
                  surv_table$Lower_CI[i] * 100, surv_table$Upper_CI[i] * 100))
    }
    
    return(list(km_overall = km_overall, km_summary = km_summary))
  }, error = function(e) cat(sprintf("    KM error: %s\n", e$message)))
}


plot_cox_vs_km <- function(coxFit, km_summary, method_dir) {
  tryCatch({
    cox_baseline <- basehaz(coxFit, centered = TRUE)
    cox_baseline$survival <- exp(-cox_baseline$hazard)
    cox_baseline$method <- "Cox Model"
    
    km_baseline <- data.frame(time = km_summary$time, survival = km_summary$surv,
                               method = "Kaplan-Meier")
    
    comparison_data <- rbind(cox_baseline[, c("time", "survival", "method")], km_baseline)
    
    gg <- ggplot(comparison_data, aes(x = time, y = survival, color = method)) +
      geom_step(linewidth = 1.2, alpha = 0.8) +
      scale_color_manual(values = c("Cox Model" = "#0072B2", "Kaplan-Meier" = "#D55E00")) +
      ylim(0, 1) +
      labs(title = "Cox Model vs Kaplan-Meier Survival Curves",
           subtitle = "Comparison of parametric (Cox) and non-parametric (KM) estimates",
           x = "Time (Years)", y = "Survival Probability", color = "Method") +
      theme_minimal(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 11, color = "gray30"),
            legend.position = "bottom")
    
    ggsave(file.path(method_dir, "cox_vs_km_comparison.png"), gg, width = 10, height = 7, dpi = 300)
    cat("    Cox vs KM comparison saved\n")
  }, error = function(e) cat(sprintf("    Cox vs KM error: %s\n", e$message)))
}


plot_survival_panel <- function(baseline_haz, baseline_haz_clean, km_summary, method_dir) {
  tryCatch({
    library(gridExtra)
    
    p1 <- ggplot(baseline_haz, aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1) + ylim(0, 1) +
      labs(title = "A) Cox Survival S(t)", x = "Time (Years)", y = "Survival") +
      theme_minimal(base_size = 10)
    
    p2 <- ggplot(baseline_haz, aes(x = time, y = hazard)) +
      geom_step(color = "#0072B2", linewidth = 1) +
      labs(title = "B) Cumulative Hazard H(t)", x = "Time (Years)", y = "H(t)") +
      theme_minimal(base_size = 10)
    
    p3 <- ggplot(baseline_haz_clean, aes(x = time, y = hazard_rate_capped)) +
      geom_smooth(se = FALSE, color = "#D55E00", linewidth = 1) +
      labs(title = "C) Hazard Rate h(t)", x = "Time (Years)", y = "h(t)") +
      theme_minimal(base_size = 10)
    
    p4 <- ggplot(data.frame(time = km_summary$time, survival = km_summary$surv,
                             lower = km_summary$lower, upper = km_summary$upper),
                 aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#0072B2") +
      ylim(0, 1) +
      labs(title = "D) Kaplan-Meier", x = "Time (Years)", y = "Survival") +
      theme_minimal(base_size = 10)
    
    ggsave(file.path(method_dir, "survival_analysis_panel.png"),
           grid.arrange(p1, p2, p3, p4, ncol = 2), width = 12, height = 10, dpi = 300)
    
    cat("    Combined panel saved\n")
  }, error = function(e) cat(sprintf("    Panel error: %s\n", e$message)))
}


plot_calibration <- function(coxFit, surv_data, raw_risks, surv_3yr, predict_calibrated_fn, method_dir) {
  tryCatch({
    calibration_data <- data.frame(
      raw_risk = raw_risks, surv_3yr_actual = surv_3yr,
      calibrated_prob = predict_calibrated_fn(surv_data),
      PTID = surv_data$PTID
    )
    write.csv(calibration_data, file.path(method_dir, "calibration_diagnostics.csv"), row.names = FALSE)
    
    cal_comparison <- data.frame(
      Method = rep(c("Before Calibration", "After Calibration"), each = length(raw_risks)),
      Predicted = c(1 - (raw_risks - min(raw_risks)) / (max(raw_risks) - min(raw_risks)),
                    calibration_data$calibrated_prob),
      Observed = rep(surv_3yr, 2)
    )
    
    gg <- ggplot(cal_comparison, aes(x = Predicted, y = Observed, color = Method)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
      facet_wrap(~ Method) + xlim(0, 1) + ylim(0, 1) +
      labs(title = "Calibration: Before vs After Correction",
           x = "Predicted 3-Year Survival", y = "Observed 3-Year Survival") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none")
    
    ggsave(file.path(method_dir, "calibration_comparison.png"), gg, width = 12, height = 6, dpi = 300)
    cat("   Calibration diagnostics saved\n")
  }, error = function(e) cat(sprintf("    Calibration plot error: %s\n", e$message)))
}


# Master function — call this from main script
run_survival_plots <- function(coxFit, surv_data, raw_risks, surv_3yr, 
                                predict_calibrated_fn, method_dir) {
  cat("\n=== STEP 9: ENHANCED SURVIVAL ANALYSIS ===\n")
  
  cat("\n1. Cox Model Survival Curve\n")
  baseline_haz <- plot_cox_survival(coxFit, method_dir)
  
  cat("\n2. Hazard Functions\n")
  hazard_results <- plot_hazard_functions(coxFit, method_dir)
  
  cat("\n3. Kaplan-Meier Analysis\n")
  km_results <- plot_kaplan_meier(surv_data, method_dir)
  
  cat("\n4. Cox vs Kaplan-Meier Comparison\n")
  plot_cox_vs_km(coxFit, km_results$km_summary, method_dir)
  
  cat("\n5. Combined Panel Figure\n")
  plot_survival_panel(baseline_haz, hazard_results$baseline_haz_clean,
                      km_results$km_summary, method_dir)
  
  cat("\n6. Calibration Plots\n")
  plot_calibration(coxFit, surv_data, raw_risks, surv_3yr, predict_calibrated_fn, method_dir)
  
  cat("\n=== ENHANCED SURVIVAL ANALYSIS COMPLETE ===\n")
}