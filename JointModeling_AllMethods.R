# ============================================================
# PACKAGE SETUP
# ============================================================

required_packages <- c("survRM2", "boot", "JMbayes2", "dplyr", "nlme",
                       "survival", "pec", "survminer", "ggplot2", "gridExtra",
                       "tidyr", "glmnet", "e1071")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

select    <- dplyr::select
summarize <- dplyr::summarize

cat("============================================================\n")
cat("FIXED: COMPREHENSIVE THESIS ANALYSIS\n")
cat("All methods on SAME patient cohort\n")
cat("============================================================\n")

if (!dir.exists("thesis_figures")) dir.create("thesis_figures")

# ============================================================
# LOAD DATA
# ============================================================

cat("\n=== LOADING DATA ===\n")

baseline_data <- read.csv("baseline_clinical_features.csv")
tabular_data  <- read.csv("tabular_only_features.csv")
image_data    <- read.csv("image_only_features.csv")

concat_data   <- read.csv("latent_concat.csv")
noae_data     <- read.csv("latent_no_ae.csv")
ae_data   <- read.csv("latent_ae_only.csv")
thesis_data   <- read.csv("latent_improved_autoencoder.csv")

cat(sprintf("  Baseline: %d patients\n", length(unique(baseline_data$PTID))))
cat(sprintf("  Tabular:  %d patients\n", length(unique(tabular_data$PTID))))
cat(sprintf("  Concat:   %d patients\n", length(unique(concat_data$PTID))))
cat(sprintf("  No-AE:    %d patients\n", length(unique(noae_data$PTID))))
cat(sprintf("  AE:    %d patients\n", length(unique(ae_data$PTID))))
cat(sprintf("  Image:    %d patients\n", length(unique(image_data$PTID))))
cat(sprintf("  Thesis:   %d patients\n", length(unique(thesis_data$PTID))))

# ============================================================
# MASTER COHORT DEFINITION
# ============================================================

cat("\n=== DEFINING MASTER COHORT ===\n")

master_patients <- Reduce(intersect, list(
  unique(baseline_data$PTID),
  unique(tabular_data$PTID),
  unique(concat_data$PTID),
  unique(noae_data$PTID),
  unique(ae_data$PTID),
  unique(image_data$PTID),
  unique(thesis_data$PTID)
))
cat(sprintf("Step 1 - Patients in all datasets: %d\n", length(master_patients)))

# Step 2: >= 2 longitudinal observations
check_min_obs <- function(data, min_obs = 2) {
  if (!all(c("MMSE", "Years_bl") %in% names(data))) {
    return(unique(data$PTID[data$PTID %in% master_patients]))
  }
  data %>%
    filter(PTID %in% master_patients, !is.na(MMSE), !is.na(Years_bl)) %>%
    group_by(PTID) %>%
    filter(n() >= min_obs) %>%
    ungroup() %>%
    pull(PTID) %>%
    unique()
}

master_patients <- Reduce(intersect, list(
  master_patients,
  check_min_obs(baseline_data),
  check_min_obs(concat_data),
  check_min_obs(noae_data),
  check_min_obs(ae_data),
  check_min_obs(image_data),
  check_min_obs(thesis_data)
))
cat(sprintf("Step 2 - After >=2 observations: %d\n", length(master_patients)))

# Step 3: Complete survival data
check_complete_surv <- function(data) {
  data %>%
    filter(PTID %in% master_patients) %>%
    group_by(PTID) %>%
    summarize(
      has_time  = !is.na(first(time_to_event)),
      has_event = !is.na(first(event)),
      time_pos  = first(time_to_event) > 0,
      .groups   = "drop"
    ) %>%
    filter(has_time, has_event, time_pos) %>%
    pull(PTID)
}

master_patients <- Reduce(intersect, list(
  master_patients,
  check_complete_surv(baseline_data),
  check_complete_surv(tabular_data),
  check_complete_surv(concat_data),
  check_complete_surv(noae_data),
  check_complete_surv(ae_data),
  check_complete_surv(image_data),
  check_complete_surv(thesis_data)
))
cat(sprintf("Step 3 - After survival filter: %d\n", length(master_patients)))

# Step 4: Complete feature availability
check_complete_features <- function(data, feature_pattern = "^z_", top_n = 20) {
  z_features <- grep(feature_pattern, names(data), value = TRUE)
  data_filt  <- data %>% filter(PTID %in% master_patients)
  
  if (length(z_features) > 0) {
    vars       <- sapply(data_filt[z_features], function(x) var(as.numeric(x), na.rm = TRUE))
    z_features <- z_features[!is.na(vars) & vars > 0]
    if (length(z_features) > top_n)
      z_features <- names(sort(vars[z_features], decreasing = TRUE))[1:top_n]
  }
  
  clinical <- intersect(c("AGE", "PTGENDER", "PTEDUCAT", "ADAS13"), names(data))
  longi    <- intersect(c("MMSE", "Years_bl"), names(data))
  required <- intersect(c("time_to_event", "event", clinical, longi, z_features), names(data))
  
  data_filt %>%
    select(PTID, all_of(required)) %>%
    filter(complete.cases(.)) %>%
    pull(PTID) %>%
    unique()
}

master_patients <- Reduce(intersect, list(
  master_patients,
  check_complete_features(concat_data,  "^z_",        20),
  check_complete_features(noae_data,    "^z_",        20),
  check_complete_features(ae_data,    "^z_",        20),
  check_complete_features(image_data,   "^z_",        20),
  check_complete_features(thesis_data,  "^z_",        20),
  check_complete_features(tabular_data, "^tab_feat_",  5)
))
cat(sprintf("Step 4 - After complete feature filter: %d\n", length(master_patients)))

# Filter all datasets to master cohort
baseline_data <- baseline_data %>% filter(PTID %in% master_patients)
tabular_data  <- tabular_data  %>% filter(PTID %in% master_patients)
concat_data   <- concat_data   %>% filter(PTID %in% master_patients)
noae_data     <- noae_data     %>% filter(PTID %in% master_patients)
image_data    <- image_data    %>% filter(PTID %in% master_patients)
thesis_data   <- thesis_data   %>% filter(PTID %in% master_patients)

cat(sprintf("\n*** FINAL MASTER COHORT: %d patients ***\n", length(master_patients)))

stopifnot(
  length(unique(baseline_data$PTID)) == length(master_patients),
  length(unique(tabular_data$PTID))  == length(master_patients),
  length(unique(concat_data$PTID))   == length(master_patients),
  length(unique(noae_data$PTID))     == length(master_patients),
  length(unique(image_data$PTID))    == length(master_patients),
  length(unique(thesis_data$PTID))   == length(master_patients)
)
cat("✓ VERIFICATION PASSED: Identical cohort across ALL methods\n")

# ============================================================
# HELPER: MEDIAN SURVIVAL FROM COX
# FIX: Replaces the rank-reversal heuristic with proper median
#      survival time from the baseline hazard function.
# ============================================================

get_median_survival <- function(coxFit, newdata) {
  # Get baseline survival curve
  bh      <- basehaz(coxFit, centered = FALSE)
  lp      <- predict(coxFit, newdata = newdata, type = "lp")
  
  median_times <- numeric(nrow(newdata))
  
  for (i in seq_along(lp)) {
    # Individual survival: S(t) = S0(t)^exp(lp)
    surv_i <- exp(-bh$hazard) ^ exp(lp[i])
    
    # Find first time where S(t) <= 0.5
    idx <- which(surv_i <= 0.5)
    if (length(idx) > 0) {
      median_times[i] <- bh$time[min(idx)]
    } else {
      # Censored beyond follow-up: use max observed time
      median_times[i] <- max(bh$time)
    }
  }
  
  return(median_times)
}

# ============================================================
# HELPER: EXTRACT JOINT MODEL METRICS
# ============================================================

extract_joint_metrics <- function(jointFit, surv_data, long_data, method_dir, name) {
  
  cat("\n  [Joint Model Metrics]\n")
  
  jm_results <- list(
    jm_cindex         = NA,
    jm_association_val  = NA,
    jm_association_slope = NA,
    jm_association_val_lower  = NA,
    jm_association_val_upper  = NA,
    jm_association_slope_lower = NA,
    jm_association_slope_upper = NA,
    jm_dynamic_auc    = NA
  )
  
  if (is.null(jointFit)) {
    cat("  ⚠ No joint model available\n")
    return(jm_results)
  }
  
  tryCatch({
    jm_sum <- summary(jointFit)
    
    # --- Association parameters (alpha) ---
    # JMbayes2 stores these under $Outcome1 in the event process summary
    assoc_table <- jm_sum$Outcome1
    
    if (!is.null(assoc_table)) {
      # value(MMSE) row
      val_row <- assoc_table[grep("value\\(MMSE\\)", rownames(assoc_table)), , drop = FALSE]
      if (nrow(val_row) > 0) {
        jm_results$jm_association_val       <- val_row[1, "Mean"]
        jm_results$jm_association_val_lower <- val_row[1, "2.5%"]
        jm_results$jm_association_val_upper <- val_row[1, "97.5%"]
        cat(sprintf("    alpha (value):  %.4f (%.4f, %.4f)\n",
                    jm_results$jm_association_val,
                    jm_results$jm_association_val_lower,
                    jm_results$jm_association_val_upper))
      }
      
      # slope(MMSE) row
      slp_row <- assoc_table[grep("slope\\(MMSE\\)", rownames(assoc_table)), , drop = FALSE]
      if (nrow(slp_row) > 0) {
        jm_results$jm_association_slope       <- slp_row[1, "Mean"]
        jm_results$jm_association_slope_lower <- slp_row[1, "2.5%"]
        jm_results$jm_association_slope_upper <- slp_row[1, "97.5%"]
        cat(sprintf("    alpha (slope):  %.4f (%.4f, %.4f)\n",
                    jm_results$jm_association_slope,
                    jm_results$jm_association_slope_lower,
                    jm_results$jm_association_slope_upper))
      }
      
      # Save association table
      write.csv(as.data.frame(assoc_table),
                file.path(method_dir, paste0(gsub(" ", "_", name), "_jm_association.csv")),
                row.names = TRUE)
      cat("    ✓ Association parameters saved\n")
    }
    
  }, error = function(e) {
    cat(sprintf("    ⚠ Association extraction error: %s\n", e$message))
  })
  
  return(jm_results)
}

# ============================================================
# MAIN METRICS FUNCTION
# ============================================================

compute_all_metrics_with_figures <- function(coxFit, surv_data, long_data,
                                             jointFit = NULL, name = "Method") {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("EVALUATING: %s\n", name))
  cat(sprintf("========================================\n"))
  
  actual_patients <- length(unique(surv_data$PTID))
  actual_events   <- sum(surv_data$event)
  
  cat(sprintf("  Patients: %d\n", actual_patients))
  cat(sprintf("  Events: %d (%.1f%%)\n", actual_events, 100 * actual_events / actual_patients))
  
  # At the top of compute_all_metrics_with_figures(), replace the method_dir line:
  method_dir <- file.path("thesis_figures", gsub("[^A-Za-z0-9_-]", "_", name))
  if (!dir.exists(method_dir)) dir.create(method_dir, recursive = TRUE)
  
  results <- list(
    name                      = name,
    n_patients                = actual_patients,
    n_events                  = actual_events,
    event_rate                = actual_events / actual_patients,
    cindex                    = NA,
    cindex_ci_lower           = NA,
    cindex_ci_upper           = NA,
    brier_1yr                 = NA,
    brier_3yr                 = NA,
    brier_5yr                 = NA,
    calibration_slope         = NA,
    calibration_intercept     = NA,
    rmst_tau                  = NA,
    rmst_low                  = NA,
    rmst_low_lower            = NA,
    rmst_low_upper            = NA,
    rmst_high                 = NA,
    rmst_high_lower           = NA,
    rmst_high_upper           = NA,
    rmst_diff                 = NA,
    rmst_diff_lower           = NA,
    rmst_diff_upper           = NA,
    rmst_pval                 = NA,
    mae_time                  = NA,
    rmse_time                 = NA,
    corr_time                 = NA,
    median_error              = NA,
    n_significant             = NA,
    km_high_median            = NA,
    km_low_median             = NA,
    # Joint model fields
    jm_association_val        = NA,
    jm_association_val_lower  = NA,
    jm_association_val_upper  = NA,
    jm_association_slope      = NA,
    jm_association_slope_lower = NA,
    jm_association_slope_upper = NA,
    jm_dynamic_auc            = NA
  )
  
  # ----------------------------------------------------------
  # 1. C-INDEX (Cox)
  # ----------------------------------------------------------
  cat("\n1. DISCRIMINATION\n")
  cindex         <- summary(coxFit)$concordance[1]
  results$cindex <- cindex
  cat(sprintf("   Cox C-index: %.4f\n", cindex))
  
  # ----------------------------------------------------------
  # 2. COX COEFFICIENTS + FOREST PLOT
  # ----------------------------------------------------------
  cat("\n2. COX MODEL PARAMETERS\n")
  tryCatch({
    coef_summary  <- summary(coxFit)$coefficients
    coef_df       <- as.data.frame(coef_summary)
    coef_df$Variable  <- rownames(coef_df)
    coef_df$HR        <- exp(coef_df$coef)
    coef_df$HR_lower  <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
    coef_df$HR_upper  <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
    
    coef_df_sorted <- coef_df[order(coef_df$`Pr(>|z|)`), ]
    cat("   Top 10 predictors:\n")
    print(head(coef_df_sorted[, c("Variable", "coef", "HR", "Pr(>|z|)")], 10))
    
    write.csv(coef_df,
              file.path(method_dir, paste0(gsub(" ", "_", name), "_cox_coefficients.csv")),
              row.names = FALSE)
    
    # Forest plot (p < 0.10)
    sig_coefs <- coef_df[coef_df$`Pr(>|z|)` < 0.10, ]
    results$n_significant <- nrow(sig_coefs)
    
    if (nrow(sig_coefs) > 0 && nrow(sig_coefs) <= 25) {
      sig_coefs <- sig_coefs[order(sig_coefs$HR, decreasing = TRUE), ]
      sig_coefs$HR_plot       <- pmin(pmax(sig_coefs$HR,       0.001), 1000)
      sig_coefs$HR_lower_plot <- pmin(pmax(sig_coefs$HR_lower, 0.001), 1000)
      sig_coefs$HR_upper_plot <- pmin(pmax(sig_coefs$HR_upper, 0.001), 1000)
      
      forest_plot <- ggplot(sig_coefs,
                            aes(x = HR_plot, y = reorder(Variable, HR_plot))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
        geom_errorbarh(aes(xmin = HR_lower_plot, xmax = HR_upper_plot),
                       height = 0.2, linewidth = 0.6) +
        geom_point(size = 3, color = "#0072B2") +
        scale_x_log10() +
        labs(title = paste("Hazard Ratios:", name), x = "HR (log scale)", y = "") +
        theme_minimal(base_size = 13)
      
      ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_forest_plot.png")),
             forest_plot, width = 10, height = max(5, nrow(sig_coefs) * 0.35), dpi = 300)
      cat(sprintf("   ✓ Forest plot saved (%d significant at p<0.10)\n", nrow(sig_coefs)))
    } else {
      cat("   No significant predictors (p < 0.10)\n")
    }
  }, error = function(e) cat(sprintf("   ⚠ Cox coefficient error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 3. BASELINE SURVIVAL FUNCTIONS
  # ----------------------------------------------------------
  cat("\n3. BASELINE SURVIVAL\n")
  tryCatch({
    baseline_surv_fn <- basehaz(coxFit, centered = TRUE)
    baseline_surv_fn$surv_prob <- exp(-baseline_surv_fn$hazard)
    
    gg_bh <- ggplot(baseline_surv_fn, aes(x = time, y = hazard)) +
      geom_step(color = "#0072B2", linewidth = 1.2) +
      labs(title = paste("Baseline Hazard:", name),
           x = "Time (Years)", y = "Cumulative Hazard") +
      theme_minimal(base_size = 13)
    
    gg_bs <- ggplot(baseline_surv_fn, aes(x = time, y = surv_prob)) +
      geom_step(color = "#0072B2", linewidth = 1.2) + ylim(0, 1) +
      labs(title = paste("Baseline Survival:", name),
           x = "Time (Years)", y = "Survival Probability") +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_baseline_hazard.png")),
           gg_bh, width = 8, height = 6, dpi = 300)
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_baseline_survival.png")),
           gg_bs, width = 8, height = 6, dpi = 300)
    cat("   ✓ Baseline functions saved\n")
  }, error = function(e) cat(sprintf("   ⚠ Baseline error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 4. BRIER SCORE
  # ----------------------------------------------------------
  cat("\n4. BRIER SCORE\n")
  tryCatch({
    pred_times <- c(1, 3, 5)
    pec_obj    <- pec(
      list("Model" = coxFit),
      formula    = Surv(time, event) ~ 1,
      data       = surv_data,
      times      = pred_times,
      exact      = FALSE,
      cens.model = "marginal",
      verbose    = FALSE
    )
    
    brier_scores    <- pec_obj$AppErr$Model
    results$brier_1yr <- brier_scores[1]
    results$brier_3yr <- brier_scores[2]
    results$brier_5yr <- brier_scores[3]
    
    cat(sprintf("   Brier 1yr: %.4f | 3yr: %.4f | 5yr: %.4f\n",
                results$brier_1yr, results$brier_3yr, results$brier_5yr))
    
    brier_df <- data.frame(Time = pred_times, Brier = brier_scores)
    gg_brier <- ggplot(brier_df, aes(x = Time, y = Brier)) +
      geom_line(color = "#0072B2", linewidth = 1.3) +
      geom_point(size = 4, color = "#0072B2") +
      geom_hline(yintercept = 0.10, linetype = "dashed", color = "darkgreen") +
      labs(title = paste("Brier Score:", name), x = "Time (Years)", y = "Brier Score") +
      theme_minimal(base_size = 13)
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_brier_scores.png")),
           gg_brier, width = 8, height = 6, dpi = 300)
    cat("   ✓ Brier plot saved\n")
  }, error = function(e) cat(sprintf("   ⚠ Brier error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 5. KM RISK STRATIFICATION + RMST
  # ----------------------------------------------------------
  cat("\n5. RISK STRATIFICATION\n")
  tryCatch({
    risk_scores        <- predict(coxFit, newdata = surv_data, type = "lp")
    surv_data$risk     <- risk_scores
    surv_data$risk_group <- ifelse(surv_data$risk >= median(surv_data$risk),
                                   "High Risk", "Low Risk")
    
    fit_stratified <- survfit(Surv(time, event) ~ risk_group, data = surv_data)
    
    km_plot <- ggsurvplot(
      fit_stratified, data = surv_data,
      risk.table = TRUE, pval = TRUE, conf.int = TRUE,
      palette = c("#D55E00", "#0072B2"),
      linetype = c("solid", "dashed"),
      legend.title = "Risk Group",
      title = paste("Kaplan-Meier:", name),
      xlab = "Time (Years)", ylab = "Progression-Free Probability"
    )
    
    ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_KM_stratified.png")),
           km_plot$plot, width = 10, height = 8, dpi = 300)
    
    km_table <- summary(fit_stratified)$table
    results$km_high_median <- km_table["risk_group=High Risk", "median"]
    results$km_low_median  <- km_table["risk_group=Low Risk",  "median"]
    
    cat(sprintf("   High-risk median: %.2f | Low-risk median: %.2f years\n",
                results$km_high_median, results$km_low_median))
    cat("   ✓ KM saved\n")
    
    # RMST (adaptive tau)
    cat("\n6. RMST ANALYSIS\n")
    tryCatch({
      surv_data$arm <- as.numeric(surv_data$risk_group == "Low Risk")
      max_tau <- min(
        max(surv_data$time[surv_data$risk_group == "High Risk"], na.rm = TRUE),
        max(surv_data$time[surv_data$risk_group == "Low Risk"],  na.rm = TRUE)
      )
      tau <- max_tau * 0.9
      cat(sprintf("   Tau: %.2f years\n", tau))
      
      rmst_res <- rmst2(surv_data$time, surv_data$event, surv_data$arm, tau = tau)
      
      results$rmst_tau        <- tau
      results$rmst_low        <- rmst_res$RMST.arm1$rmst["Est."]
      results$rmst_low_lower  <- rmst_res$RMST.arm1$rmst["lower .95"]
      results$rmst_low_upper  <- rmst_res$RMST.arm1$rmst["upper .95"]
      results$rmst_high       <- rmst_res$RMST.arm0$rmst["Est."]
      results$rmst_high_lower <- rmst_res$RMST.arm0$rmst["lower .95"]
      results$rmst_high_upper <- rmst_res$RMST.arm0$rmst["upper .95"]
      results$rmst_diff       <- rmst_res$unadjusted.result[1, "Est."]
      results$rmst_diff_lower <- rmst_res$unadjusted.result[1, "lower .95"]
      results$rmst_diff_upper <- rmst_res$unadjusted.result[1, "upper .95"]
      results$rmst_pval       <- rmst_res$unadjusted.result[1, "p"]
      
      cat(sprintf("   RMST diff: %.3f yr (p=%.4f)\n", results$rmst_diff, results$rmst_pval))
      
      rmst_table <- data.frame(
        Risk_Group = c("Low Risk", "High Risk", "Difference"),
        RMST_Years = c(results$rmst_low,  results$rmst_high,  results$rmst_diff),
        Lower_CI   = c(results$rmst_low_lower,  results$rmst_high_lower,  results$rmst_diff_lower),
        Upper_CI   = c(results$rmst_low_upper,  results$rmst_high_upper,  results$rmst_diff_upper),
        P_Value    = c(NA, NA, results$rmst_pval),
        Tau        = tau
      )
      write.csv(rmst_table,
                file.path(method_dir, paste0(gsub(" ", "_", name), "_rmst_analysis.csv")),
                row.names = FALSE)
      cat("   ✓ RMST saved\n")
    }, error = function(e) cat(sprintf("   ⚠ RMST error: %s\n", e$message)))
    
  }, error = function(e) cat(sprintf("   ⚠ KM error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 7. PREDICTED vs OBSERVED (FIX: proper median survival)
  # ----------------------------------------------------------
  cat("\n7. PREDICTED vs OBSERVED TIMES\n")
  tryCatch({
    converters <- surv_data[surv_data$event == 1, ]
    
    if (nrow(converters) >= 5) {
      predicted_times <- get_median_survival(coxFit, converters)
      observed_times  <- converters$time
      
      errors              <- abs(predicted_times - observed_times)
      results$mae_time    <- mean(errors)
      results$rmse_time   <- sqrt(mean((predicted_times - observed_times)^2))
      results$corr_time   <- cor(predicted_times, observed_times, method = "pearson")
      results$median_error <- median(errors)
      
      cat(sprintf("   MAE: %.3f yr | RMSE: %.3f yr | r: %.3f\n",
                  results$mae_time, results$rmse_time, results$corr_time))
      
      pred_obs_df <- data.frame(
        Predicted = predicted_times,
        Observed  = observed_times,
        PTID      = converters$PTID,
        Error     = errors
      )
      
      gg_po <- ggplot(pred_obs_df, aes(x = Observed, y = Predicted)) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
        geom_point(alpha = 0.7, size = 3.5, color = "#0072B2") +
        geom_smooth(method = "lm", se = TRUE, color = "#D55E00", alpha = 0.2) +
        labs(title = paste("Predicted vs Observed:", name),
             x = "Observed (Years)", y = "Predicted Median Survival (Years)") +
        theme_minimal(base_size = 13) + coord_equal()
      
      ggsave(file.path(method_dir, paste0(gsub(" ", "_", name), "_predicted_vs_observed.png")),
             gg_po, width = 8, height = 8, dpi = 300)
      write.csv(pred_obs_df,
                file.path(method_dir, paste0(gsub(" ", "_", name), "_time_predictions.csv")),
                row.names = FALSE)
      cat("   ✓ Predictions saved\n")
    }
  }, error = function(e) cat(sprintf("   ⚠ Prediction error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 8. CALIBRATION 
  # - Predicted survival from actual baseline hazard (not linear rescaling)
  # - Observed survival from KM per group (not binary flag)
  # - Adaptive grouping handles compressed risk score distributions
  # - Graceful fallback if too few unique breaks
  # ----------------------------------------------------------
  cat("\n8. CALIBRATION\n")
  tryCatch({
    
    # Step 1: Compute predicted 3-year survival from baseline hazard
    bh <- basehaz(coxFit, centered = FALSE)
    lp <- predict(coxFit, newdata = surv_data, type = "lp")
    
    t3_idx <- which(bh$time <= 3)
    if (length(t3_idx) == 0) {
      cat("   ⚠ No baseline hazard estimates at <=3 years — skipping calibration\n")
      stop("no 3yr hazard")
    }
    H0_3yr             <- bh$hazard[max(t3_idx)]
    predicted_surv_3yr <- exp(-H0_3yr * exp(lp))
    
    cat(sprintf("   Predicted 3yr survival range: %.3f – %.3f\n",
                min(predicted_surv_3yr), max(predicted_surv_3yr)))
    
    # Step 2: Adaptive grouping — fewer groups when risk scores are compressed
    n_unique <- length(unique(round(predicted_surv_3yr, 4)))
    n_groups <- dplyr::case_when(
      n_unique >= 50 ~ 10L,
      n_unique >= 30 ~ 7L,
      n_unique >= 15 ~ 5L,
      n_unique >= 6  ~ 3L,
      TRUE           ~ 0L
    )
    
    if (n_groups < 3) {
      cat("   ⚠ Risk scores too homogeneous for calibration — skipping\n")
      stop("too homogeneous")
    }
    
    # Step 3: Cut into groups using unique quantile breaks
    decile_probs  <- seq(0, 1, length.out = n_groups + 1)
    raw_breaks    <- quantile(predicted_surv_3yr, probs = decile_probs, na.rm = TRUE)
    decile_breaks <- unique(raw_breaks)
    
    if (length(decile_breaks) < 3) {
      cat(sprintf("   ⚠ Only %d unique breaks after dedup — skipping calibration\n",
                  length(decile_breaks)))
      stop("too few breaks")
    }
    
    deciles <- cut(predicted_surv_3yr,
                   breaks         = decile_breaks,
                   include.lowest = TRUE,
                   labels         = FALSE)
    
    cat(sprintf("   Using %d calibration groups (%d unique breaks)\n",
                n_groups, length(decile_breaks)))
    
    # Step 4: KM-observed survival at 3 years per group
    cal_rows <- list()
    
    for (g in sort(unique(na.omit(deciles)))) {
      idx      <- which(deciles == g)
      if (length(idx) < 5) next   # skip tiny groups
      
      grp_data   <- surv_data[idx, ]
      km_fit     <- survfit(Surv(time, event) ~ 1, data = grp_data)
      km_summary <- summary(km_fit, times = 3, extend = TRUE)
      obs_surv   <- km_summary$surv
      
      if (length(obs_surv) == 0 || is.na(obs_surv)) next
      
      cal_rows[[length(cal_rows) + 1]] <- data.frame(
        predicted = mean(predicted_surv_3yr[idx]),
        observed  = obs_surv,
        n         = length(idx)
      )
    }
    
    if (length(cal_rows) < 3) {
      cat(sprintf("   ⚠ Only %d valid calibration groups (need >=3) — skipping\n",
                  length(cal_rows)))
      stop("too few groups")
    }
    
    cal_data <- do.call(rbind, cal_rows)
    
    # Step 5: Fit calibration line
    cal_model                     <- lm(observed ~ predicted, data = cal_data)
    results$calibration_slope     <- coef(cal_model)[2]
    results$calibration_intercept <- coef(cal_model)[1]
    
    cat(sprintf("   Groups used: %d | Slope: %.3f | Intercept: %.3f\n",
                nrow(cal_data),
                results$calibration_slope,
                results$calibration_intercept))
    cat(sprintf("   Interpretation: slope=1.0 & intercept=0.0 is perfect calibration\n"))
    
    # Step 6: Save calibration data
    write.csv(cal_data,
              file.path(method_dir,
                        paste0(gsub("[^A-Za-z0-9_-]", "_", name), "_calibration_data.csv")),
              row.names = FALSE)
    
    # Step 7: Plot
    gg_cal <- ggplot(cal_data, aes(x = predicted, y = observed)) +
      geom_abline(intercept = 0, slope = 1,
                  linetype = "dashed", color = "gray40", linewidth = 0.8) +
      geom_smooth(method = "lm", se = TRUE,
                  color = "#D55E00", fill = "#D55E00", alpha = 0.15) +
      geom_point(aes(size = n), color = "#0072B2", alpha = 0.85) +
      scale_size_continuous(name = "Group n", range = c(3, 8)) +
      xlim(0, 1) + ylim(0, 1) +
      coord_equal() +
      labs(
        title    = paste("Calibration (3-year):", name),
        subtitle = sprintf("Slope = %.3f  |  Intercept = %.3f  |  Groups = %d",
                           results$calibration_slope,
                           results$calibration_intercept,
                           nrow(cal_data)),
        x        = "Predicted 3-yr Survival (Cox)",
        y        = "KM Observed 3-yr Survival"
      ) +
      theme_minimal(base_size = 13) +
      theme(plot.subtitle = element_text(size = 10, color = "gray40"))
    
    ggsave(file.path(method_dir,
                     paste0(gsub("[^A-Za-z0-9_-]", "_", name), "_calibration_plot.png")),
           gg_cal, width = 8, height = 8, dpi = 300)
    
    cat("   ✓ Calibration plot saved\n")
    
  }, error = function(e) {
    # Only print if it's an unexpected error, not one of our deliberate stops
    if (!e$message %in% c("no 3yr hazard", "too homogeneous",
                          "too few breaks", "too few groups")) {
      cat(sprintf("   ⚠ Calibration error: %s\n", e$message))
    }
  })
  
  # ----------------------------------------------------------
  # 9. BOOTSTRAP CI (R=500 for thesis)
  # ----------------------------------------------------------
  cat("\n9. BOOTSTRAP CI\n")
  tryCatch({
    boot_cindex <- function(data, indices) {
      boot_data <- data[indices, ]
      boot_cox  <- tryCatch(
        coxph(coxFit$formula, data = boot_data, x = TRUE),
        error = function(e) NA
      )
      if (length(boot_cox) == 1 && is.na(boot_cox)) return(NA)
      summary(boot_cox)$concordance[1]
    }
    
    cat("   Running bootstrap (R=500)...\n")
    boot_res    <- boot(data = surv_data, statistic = boot_cindex, R = 500)
    valid_boot  <- boot_res$t[!is.na(boot_res$t)]
    
    if (length(valid_boot) >= 10) {
      boot_ci                    <- boot.ci(boot_res, type = "perc")
      results$cindex_ci_lower    <- boot_ci$percent[4]
      results$cindex_ci_upper    <- boot_ci$percent[5]
      
      cat(sprintf("   C-index: %.4f (95%% CI: %.4f-%.4f)\n",
                  cindex, results$cindex_ci_lower, results$cindex_ci_upper))
      saveRDS(boot_res,
              file.path(method_dir, paste0(gsub(" ", "_", name), "_bootstrap.rds")))
      cat("   ✓ Bootstrap complete\n")
    }
  }, error = function(e) cat(sprintf("   ⚠ Bootstrap error: %s\n", e$message)))
  
  # ----------------------------------------------------------
  # 10. JOINT MODEL METRICS (FIX: now actually called)
  # ----------------------------------------------------------
  cat("\n10. JOINT MODEL METRICS\n")
  jm_metrics <- extract_joint_metrics(jointFit, surv_data, long_data, method_dir, name)
  
  # Merge joint model metrics into results
  results$jm_association_val         <- jm_metrics$jm_association_val
  results$jm_association_val_lower   <- jm_metrics$jm_association_val_lower
  results$jm_association_val_upper   <- jm_metrics$jm_association_val_upper
  results$jm_association_slope       <- jm_metrics$jm_association_slope
  results$jm_association_slope_lower <- jm_metrics$jm_association_slope_lower
  results$jm_association_slope_upper <- jm_metrics$jm_association_slope_upper
  results$jm_dynamic_auc             <- jm_metrics$jm_dynamic_auc
  
  # ----------------------------------------------------------
  # APPEND TO MASTER RESULTS FILE
  # ----------------------------------------------------------
  results_df  <- as.data.frame(results, stringsAsFactors = FALSE)
  master_file <- "thesis_figures/ALL_METHODS_RESULTS.csv"
  
  if (file.exists(master_file)) {
    existing <- read.csv(master_file, stringsAsFactors = FALSE)
    # Align columns in case some methods have extra fields
    all_cols <- union(names(existing), names(results_df))
    for (col in setdiff(all_cols, names(existing)))  existing[[col]]    <- NA
    for (col in setdiff(all_cols, names(results_df))) results_df[[col]] <- NA
    all_results <- rbind(existing[, all_cols], results_df[, all_cols])
  } else {
    all_results <- results_df
  }
  
  write.csv(all_results, master_file, row.names = FALSE)
  
  # ----------------------------------------------------------
  # SUMMARY
  # ----------------------------------------------------------
  cat(sprintf("\n========================================\n"))
  cat(sprintf("SUMMARY: %s\n", name))
  cat(sprintf("========================================\n"))
  cat(sprintf("  Cox C-index:      %.4f", results$cindex))
  if (!is.na(results$cindex_ci_lower))
    cat(sprintf(" (%.4f, %.4f)", results$cindex_ci_lower, results$cindex_ci_upper))
  cat("\n")
  if (!is.na(results$brier_3yr))         cat(sprintf("  Brier (3yr):      %.4f\n", results$brier_3yr))
  if (!is.na(results$calibration_slope)) cat(sprintf("  Cal. slope:       %.3f\n",  results$calibration_slope))
  if (!is.na(results$rmst_diff))         cat(sprintf("  RMST diff:        %.2f yr (p=%.4f)\n", results$rmst_diff, results$rmst_pval))
  if (!is.na(results$mae_time))          cat(sprintf("  MAE:              %.3f yr\n", results$mae_time))
  if (!is.na(results$jm_dynamic_auc))    cat(sprintf("  JM Dynamic AUC:   %.4f\n",   results$jm_dynamic_auc))
  if (!is.na(results$jm_association_val))
    cat(sprintf("  alpha (value):    %.4f (%.4f, %.4f)\n",
                results$jm_association_val,
                results$jm_association_val_lower,
                results$jm_association_val_upper))
  cat("========================================\n")
  
  return(results)
}

# ============================================================
# METHOD 1: CLINICAL COX (BASELINE)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 1: CLINICAL COX (BASELINE)\n")
cat("============================================================\n")

metrics_baseline <- NULL
tryCatch({
  baseline_long <- baseline_data %>%
    filter(!is.na(MMSE)) %>%
    select(PTID, Years_bl, MMSE, ADAS13, AGE, PTGENDER, PTEDUCAT) %>%
    arrange(PTID, Years_bl)
  baseline_long <- baseline_long[complete.cases(baseline_long), ]
  
  baseline_surv <- baseline_data %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              AGE = first(AGE), PTGENDER = first(PTGENDER),
              PTEDUCAT = first(PTEDUCAT), ADAS13 = first(ADAS13),
              .groups = "drop")
  baseline_surv <- baseline_surv[complete.cases(baseline_surv), ]
  
  common_IDs <- intersect(baseline_long$PTID, baseline_surv$PTID)
  baseline_long <- baseline_long %>% filter(PTID %in% common_IDs)
  baseline_surv <- baseline_surv %>% filter(PTID %in% common_IDs)
  
  cat(sprintf("Patients: %d, Events: %d\n", length(common_IDs), sum(baseline_surv$event)))
  
  lmeFit_baseline <- lme(
    MMSE ~ Years_bl + AGE + PTGENDER + PTEDUCAT + ADAS13,
    random  = ~ Years_bl | PTID,
    data    = baseline_long,
    control = lmeControl(opt = "optim", maxIter = 200)
  )
  
  coxFit_baseline <- coxph(
    Surv(time, event) ~ AGE + PTGENDER + PTEDUCAT + ADAS13,
    data = baseline_surv, x = TRUE
  )
  
  jointFit_baseline <- jm(
    coxFit_baseline, lmeFit_baseline, time_var = "Years_bl",
    functional_forms = ~ value(MMSE) + slope(MMSE),
    n_iter = 5000, n_burnin = 1000, n_thin = 5, n_chains = 2
  )
  
  metrics_baseline <- compute_all_metrics_with_figures(
    coxFit_baseline, baseline_surv, baseline_long,
    jointFit_baseline, "Clinical Cox"
  )
  cat("\n✓ METHOD 1 COMPLETED\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in Clinical Cox:", e$message, "\n")
})

# ============================================================
# METHOD 2: IMAGE-ONLY CNN
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 2: IMAGE-ONLY CNN\n")
cat("============================================================\n")

metrics_image <- NULL
tryCatch({
  img_features <- grep("^img_feat_|^z_", names(image_data), value = TRUE)
  if (length(img_features) == 0) stop("No image features found")
  
  if (length(img_features) > 10) {
    img_vars     <- sapply(image_data[img_features], var, na.rm = TRUE)
    img_features <- names(sort(img_vars, decreasing = TRUE)[1:10])
  }
  
  image_surv <- image_data %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              across(all_of(img_features), first), .groups = "drop")
  image_surv <- image_surv[complete.cases(image_surv), ]
  
  cat(sprintf("Patients: %d, Events: %d\n", nrow(image_surv), sum(image_surv$event)))
  
  surv_formula  <- as.formula(paste("Surv(time, event) ~", paste(img_features, collapse = " + ")))
  coxFit_image  <- coxph(surv_formula, data = image_surv, x = TRUE)
  
  # Image-only has no joint model (no longitudinal submodel)
  metrics_image <- compute_all_metrics_with_figures(
    coxFit_image, image_surv, NULL, NULL, "Image-Only CNN"
  )
  cat("\n✓ METHOD 2 COMPLETED\n")
  
}, error = function(e) cat("\n⚠ ERROR in Image-Only:", e$message, "\n"))

# ============================================================
# METHOD 3 FIX: TABULAR-ONLY
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 3: TABULAR-ONLY DEEP NN (FIXED)\n")
cat("============================================================\n")

metrics_tabular <- NULL
tryCatch({
  
  if (!file.exists("tabular_patient_level.csv")) {
    stop("tabular_patient_level.csv not found — run export_patient_level_tabular() in Python first.")
  }
  
  tab_pl <- read.csv("tabular_patient_level.csv", stringsAsFactors = FALSE)
  tab_pl <- tab_pl %>% filter(PTID %in% master_patients)
  cat(sprintf("  Loaded: %d patients\n", nrow(tab_pl)))
  
  z_final_cols <- grep("^z_final_", names(tab_pl), value = TRUE)
  z_slope_cols <- grep("^z_slope_", names(tab_pl), value = TRUE)
  cat(sprintf("  z_final: %d  z_slope: %d\n", length(z_final_cols), length(z_slope_cols)))
  
  if (length(z_final_cols) == 0) stop("No z_final_ columns — check Python export.")
  
  # Remove zero-variance only
  keep_zf <- sapply(z_final_cols, function(f) { v <- var(tab_pl[[f]], na.rm=TRUE); !is.na(v) && v > 1e-10 })
  keep_zs <- sapply(z_slope_cols, function(f) { v <- var(tab_pl[[f]], na.rm=TRUE); !is.na(v) && v > 1e-10 })
  z_final_cols <- z_final_cols[keep_zf]
  z_slope_cols <- z_slope_cols[keep_zs]
  
  # Standardise
  for (f in c(z_final_cols, z_slope_cols)) {
    m <- mean(tab_pl[[f]], na.rm=TRUE); s <- sd(tab_pl[[f]], na.rm=TRUE)
    if (s > 1e-10) tab_pl[[f]] <- (tab_pl[[f]] - m) / s
    tab_pl[[f]] <- pmin(pmax(tab_pl[[f]], -5), 5)
  }
  
  all_cols <- c(z_final_cols, z_slope_cols)
  
  tabular_surv <- tab_pl %>%
    select(PTID, time_to_event, event, risk_score, all_of(all_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  tabular_surv <- tabular_surv[complete.cases(tabular_surv), ]
  
  # Elastic net feature selection
  x_mat <- as.matrix(tabular_surv[, all_cols])
  y_mat <- Surv(tabular_surv$time, tabular_surv$event)
  set.seed(42)
  cv_fit   <- cv.glmnet(x_mat, y_mat, family="cox", alpha=0.5, nfolds=10, type.measure="C")
  coefs    <- coef(cv_fit, s=cv_fit$lambda.1se)
  sel_cols <- rownames(coefs)[which(coefs != 0)]
  if (length(sel_cols) < 5) {
    coefs    <- coef(cv_fit, s=cv_fit$lambda.min)
    sel_cols <- rownames(coefs)[which(coefs != 0)]
  }
  if (length(sel_cols) < 5) {
    uni_c  <- sapply(all_cols, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event)~",f)), data=tabular_surv))$concordance[1],
      error=function(e) 0.5))
    sel_cols <- names(sort(uni_c, decreasing=TRUE))[1:min(10, length(all_cols))]
  }
  cat(sprintf("  Selected features: %d\n", length(sel_cols)))
  
  tabular_surv_final <- tabular_surv %>% select(PTID, time, event, all_of(sel_cols))
  surv_formula <- as.formula(paste("Surv(time, event) ~", paste(sel_cols, collapse=" + ")))
  coxFit_tabular <- coxph(surv_formula, data=tabular_surv_final, x=TRUE, method="breslow")
  cat(sprintf("  Cox C-index: %.4f\n", summary(coxFit_tabular)$concordance[1]))
  
  # No joint model — tabular-only has no MMSE longitudinal structure
  # (the LSTM encoded tabular features, not raw MMSE trajectories)
  metrics_tabular <- compute_all_metrics_with_figures(
    coxFit_tabular, tabular_surv_final, NULL, NULL, "Tabular-Only NN"
  )
  cat("\n✓ METHOD 3 (FIXED) COMPLETED\n")
  
}, error = function(e) cat(sprintf("\n⚠ ERROR in Tabular-Only: %s\n", e$message)))


cat("\n\n============================================================\n")
cat("METHOD 4: CONCATENATION FUSION (FIXED)\n")
cat("============================================================\n")

metrics_concat_latent <- NULL
tryCatch({
  
  if (!file.exists("concat_patient_level.csv")) {
    stop("concat_patient_level.csv not found — run export_patient_level_latent() in Python first.")
  }
  
  concat_pl <- read.csv("concat_patient_level.csv", stringsAsFactors = FALSE)
  concat_pl <- concat_pl %>% filter(PTID %in% master_patients)
  
  z_final_cols <- grep("^z_final_", names(concat_pl), value=TRUE)
  z_slope_cols <- grep("^z_slope_", names(concat_pl), value=TRUE)
  cat(sprintf("  z_final: %d  z_slope: %d\n", length(z_final_cols), length(z_slope_cols)))
  
  # Zero-variance filter + standardise
  keep <- function(cols, data) cols[sapply(cols, function(f) { v <- var(data[[f]], na.rm=TRUE); !is.na(v) && v > 1e-10 })]
  z_final_cols <- keep(z_final_cols, concat_pl)
  z_slope_cols <- keep(z_slope_cols, concat_pl)
  for (f in c(z_final_cols, z_slope_cols)) {
    m <- mean(concat_pl[[f]], na.rm=TRUE); s <- sd(concat_pl[[f]], na.rm=TRUE)
    if (s > 1e-10) concat_pl[[f]] <- (concat_pl[[f]] - m) / s
    concat_pl[[f]] <- pmin(pmax(concat_pl[[f]], -5), 5)
  }
  all_cols <- c(z_final_cols, z_slope_cols)
  
  concat_surv <- concat_pl %>%
    select(PTID, time_to_event, event, all_of(all_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  concat_surv <- concat_surv[complete.cases(concat_surv), ]
  
  # Elastic net
  x_mat <- as.matrix(concat_surv[, all_cols])
  y_mat <- Surv(concat_surv$time, concat_surv$event)
  set.seed(42)
  cv_fit <- cv.glmnet(x_mat, y_mat, family="cox", alpha=0.5, nfolds=10, type.measure="C")
  coefs  <- coef(cv_fit, s=cv_fit$lambda.1se)
  sel    <- rownames(coefs)[which(coefs != 0)]
  if (length(sel) < 5) { coefs <- coef(cv_fit, s=cv_fit$lambda.min); sel <- rownames(coefs)[which(coefs != 0)] }
  if (length(sel) < 5) {
    uni_c <- sapply(all_cols, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event)~",f)), data=concat_surv))$concordance[1], error=function(e) 0.5))
    sel <- names(sort(uni_c, decreasing=TRUE))[1:min(10,length(all_cols))]
  }
  cat(sprintf("  Selected features: %d\n", length(sel)))
  
  concat_surv_final <- concat_surv %>% select(PTID, time, event, all_of(sel))
  surv_formula  <- as.formula(paste("Surv(time, event) ~", paste(sel, collapse=" + ")))
  coxFit_concat <- coxph(surv_formula, data=concat_surv_final, x=TRUE, method="breslow")
  cat(sprintf("  Cox C-index: %.4f\n", summary(coxFit_concat)$concordance[1]))
  
  # Longitudinal data for joint model — use concat_data (original per-visit CSV)
  concat_long <- concat_data %>%
    filter(PTID %in% concat_surv_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(concat_pl %>% select(PTID, all_of(sel)), by="PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  concat_long <- concat_long[complete.cases(concat_long), ]
  
  # Remove patients with < 2 visits
  single_v <- concat_long %>% group_by(PTID) %>% filter(n() < 2) %>% pull(PTID) %>% unique()
  concat_long        <- concat_long        %>% filter(!PTID %in% single_v)
  concat_surv_final  <- concat_surv_final  %>% filter(!PTID %in% single_v)
  
  common <- intersect(concat_long$PTID, concat_surv_final$PTID)
  concat_long       <- concat_long       %>% filter(PTID %in% common)
  concat_surv_final <- concat_surv_final %>% filter(PTID %in% common)
  cat(sprintf("  Patients: %d, Events: %d\n", length(common), sum(concat_surv_final$event)))
  
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_concat <- tryCatch(
    lme(long_formula, random=~Years_bl|PTID, data=concat_long,
        control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE)),
    error=function(e) lme(long_formula, random=~1|PTID, data=concat_long,
                          control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE))
  )
  
  jointFit_concat <- tryCatch(
    jm(coxFit_concat, lmeFit_concat, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) tryCatch(
      jm(coxFit_concat, lmeFit_concat, time_var="Years_bl",
         functional_forms=~value(MMSE),
         n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
      error=function(e2) NULL
    )
  )
  
  metrics_concat_latent <- compute_all_metrics_with_figures(
    coxFit_concat, concat_surv_final, concat_long,
    jointFit_concat, "Concatenation Fusion"
  )
  cat("\n✓ METHOD 4 (FIXED) COMPLETED\n")
  
}, error=function(e) cat(sprintf("\n⚠ ERROR in Concat: %s\n", e$message)))

# ============================================================
# METHOD 5 FIX: NO AUTOENCODER
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 5: NO AUTOENCODER (FIXED)\n")
cat("============================================================\n")

metrics_noae_latent <- NULL
tryCatch({
  
  if (!file.exists("no_ae_patient_level.csv")) {
    stop("no_ae_patient_level.csv not found — run export_patient_level_latent() in Python first.")
  }
  
  noae_pl <- read.csv("no_ae_patient_level.csv", stringsAsFactors = FALSE)
  noae_pl <- noae_pl %>% filter(PTID %in% master_patients)
  
  z_final_cols <- grep("^z_final_", names(noae_pl), value=TRUE)
  z_slope_cols <- grep("^z_slope_", names(noae_pl), value=TRUE)
  
  keep <- function(cols, data) cols[sapply(cols, function(f) { v <- var(data[[f]], na.rm=TRUE); !is.na(v) && v > 1e-10 })]
  z_final_cols <- keep(z_final_cols, noae_pl)
  z_slope_cols <- keep(z_slope_cols, noae_pl)
  for (f in c(z_final_cols, z_slope_cols)) {
    m <- mean(noae_pl[[f]], na.rm=TRUE); s <- sd(noae_pl[[f]], na.rm=TRUE)
    if (s > 1e-10) noae_pl[[f]] <- (noae_pl[[f]] - m) / s
    noae_pl[[f]] <- pmin(pmax(noae_pl[[f]], -5), 5)
  }
  all_cols <- c(z_final_cols, z_slope_cols)
  
  noae_surv <- noae_pl %>%
    select(PTID, time_to_event, event, all_of(all_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  noae_surv <- noae_surv[complete.cases(noae_surv), ]
  
  x_mat <- as.matrix(noae_surv[, all_cols])
  y_mat <- Surv(noae_surv$time, noae_surv$event)
  set.seed(42)
  cv_fit <- cv.glmnet(x_mat, y_mat, family="cox", alpha=0.5, nfolds=10, type.measure="C")
  coefs  <- coef(cv_fit, s=cv_fit$lambda.1se)
  sel    <- rownames(coefs)[which(coefs != 0)]
  if (length(sel) < 5) { coefs <- coef(cv_fit, s=cv_fit$lambda.min); sel <- rownames(coefs)[which(coefs != 0)] }
  if (length(sel) < 5) {
    uni_c <- sapply(all_cols, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event)~",f)), data=noae_surv))$concordance[1], error=function(e) 0.5))
    sel <- names(sort(uni_c, decreasing=TRUE))[1:min(10,length(all_cols))]
  }
  cat(sprintf("  Selected features: %d\n", length(sel)))
  
  noae_surv_final <- noae_surv %>% select(PTID, time, event, all_of(sel))
  surv_formula <- as.formula(paste("Surv(time, event) ~", paste(sel, collapse=" + ")))
  coxFit_noae  <- coxph(surv_formula, data=noae_surv_final, x=TRUE, method="breslow")
  cat(sprintf("  Cox C-index: %.4f\n", summary(coxFit_noae)$concordance[1]))
  
  noae_long <- noae_data %>%
    filter(PTID %in% noae_surv_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(noae_pl %>% select(PTID, all_of(sel)), by="PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  noae_long <- noae_long[complete.cases(noae_long), ]
  
  single_v   <- noae_long %>% group_by(PTID) %>% filter(n() < 2) %>% pull(PTID) %>% unique()
  noae_long       <- noae_long       %>% filter(!PTID %in% single_v)
  noae_surv_final <- noae_surv_final %>% filter(!PTID %in% single_v)
  common <- intersect(noae_long$PTID, noae_surv_final$PTID)
  noae_long       <- noae_long       %>% filter(PTID %in% common)
  noae_surv_final <- noae_surv_final %>% filter(PTID %in% common)
  cat(sprintf("  Patients: %d, Events: %d\n", length(common), sum(noae_surv_final$event)))
  
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_noae <- tryCatch(
    lme(long_formula, random=~Years_bl|PTID, data=noae_long,
        control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE)),
    error=function(e) lme(long_formula, random=~1|PTID, data=noae_long,
                          control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE))
  )
  
  jointFit_noae <- tryCatch(
    jm(coxFit_noae, lmeFit_noae, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) tryCatch(
      jm(coxFit_noae, lmeFit_noae, time_var="Years_bl",
         functional_forms=~value(MMSE),
         n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
      error=function(e2) NULL
    )
  )
  
  metrics_noae_latent <- compute_all_metrics_with_figures(
    coxFit_noae, noae_surv_final, noae_long,
    jointFit_noae, "No Autoencoder"
  )
  cat("\n✓ METHOD 5 (FIXED) COMPLETED\n")
  
}, error=function(e) cat(sprintf("\n⚠ ERROR in No-AE: %s\n", e$message)))


# ============================================================
# METHOD 6 FIX: AE-ONLY (NO MTL)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 6: AE-ONLY / NO-MTL (FIXED)\n")
cat("============================================================\n")

metrics_no_mtl_latent <- NULL
tryCatch({
  
  if (!file.exists("ae_only_patient_level.csv")) {
    stop("ae_only_patient_level.csv not found — run export_patient_level_latent() in Python first.")
  }
  
  ae_pl <- read.csv("ae_only_patient_level.csv", stringsAsFactors = FALSE)
  ae_pl <- ae_pl %>% filter(PTID %in% master_patients)
  
  z_final_cols <- grep("^z_final_", names(ae_pl), value=TRUE)
  z_slope_cols <- grep("^z_slope_", names(ae_pl), value=TRUE)
  
  keep <- function(cols, data) cols[sapply(cols, function(f) { v <- var(data[[f]], na.rm=TRUE); !is.na(v) && v > 1e-10 })]
  z_final_cols <- keep(z_final_cols, ae_pl)
  z_slope_cols <- keep(z_slope_cols, ae_pl)
  for (f in c(z_final_cols, z_slope_cols)) {
    m <- mean(ae_pl[[f]], na.rm=TRUE); s <- sd(ae_pl[[f]], na.rm=TRUE)
    if (s > 1e-10) ae_pl[[f]] <- (ae_pl[[f]] - m) / s
    ae_pl[[f]] <- pmin(pmax(ae_pl[[f]], -5), 5)
  }
  all_cols <- c(z_final_cols, z_slope_cols)
  
  ae_surv <- ae_pl %>%
    select(PTID, time_to_event, event, all_of(all_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  ae_surv <- ae_surv[complete.cases(ae_surv), ]
  
  x_mat <- as.matrix(ae_surv[, all_cols])
  y_mat <- Surv(ae_surv$time, ae_surv$event)
  set.seed(42)
  cv_fit <- cv.glmnet(x_mat, y_mat, family="cox", alpha=0.5, nfolds=10, type.measure="C")
  coefs  <- coef(cv_fit, s=cv_fit$lambda.1se)
  sel    <- rownames(coefs)[which(coefs != 0)]
  if (length(sel) < 5) { coefs <- coef(cv_fit, s=cv_fit$lambda.min); sel <- rownames(coefs)[which(coefs != 0)] }
  if (length(sel) < 5) {
    uni_c <- sapply(all_cols, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event)~",f)), data=ae_surv))$concordance[1], error=function(e) 0.5))
    sel <- names(sort(uni_c, decreasing=TRUE))[1:min(10,length(all_cols))]
  }
  cat(sprintf("  Selected features: %d\n", length(sel)))
  
  ae_surv_final <- ae_surv %>% select(PTID, time, event, all_of(sel))
  surv_formula <- as.formula(paste("Surv(time, event) ~", paste(sel, collapse=" + ")))
  # NOTE: ae_only was trained with reconstruction loss ONLY — no survival signal.
  # We EXPECT this to underperform thesis. That's the point of this ablation.
  coxFit_ae <- coxph(surv_formula, data=ae_surv_final, x=TRUE, method="breslow")
  cat(sprintf("  Cox C-index: %.4f\n", summary(coxFit_ae)$concordance[1]))
  cat("  NOTE: AE-only trained without survival signal — lower C-index expected by design.\n")
  
  # Use noae_data (thesis per-visit data) for MMSE trajectory
  # ae_only doesn't have a dedicated per-visit CSV in the original R code
  # Use thesis_data as the source of MMSE trajectories
  ae_long <- thesis_data %>%
    filter(PTID %in% ae_surv_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(ae_pl %>% select(PTID, all_of(sel)), by="PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  ae_long <- ae_long[complete.cases(ae_long), ]
  
  single_v    <- ae_long %>% group_by(PTID) %>% filter(n() < 2) %>% pull(PTID) %>% unique()
  ae_long       <- ae_long       %>% filter(!PTID %in% single_v)
  ae_surv_final <- ae_surv_final %>% filter(!PTID %in% single_v)
  common <- intersect(ae_long$PTID, ae_surv_final$PTID)
  ae_long       <- ae_long       %>% filter(PTID %in% common)
  ae_surv_final <- ae_surv_final %>% filter(PTID %in% common)
  cat(sprintf("  Patients: %d, Events: %d\n", length(common), sum(ae_surv_final$event)))
  
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_ae <- tryCatch(
    lme(long_formula, random=~Years_bl|PTID, data=ae_long,
        control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE)),
    error=function(e) lme(long_formula, random=~1|PTID, data=ae_long,
                          control=lmeControl(opt="optim", maxIter=500, returnObject=TRUE))
  )
  
  jointFit_ae <- tryCatch(
    jm(coxFit_ae, lmeFit_ae, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) tryCatch(
      jm(coxFit_ae, lmeFit_ae, time_var="Years_bl",
         functional_forms=~value(MMSE),
         n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
      error=function(e2) NULL
    )
  )
  
  metrics_no_mtl_latent <- compute_all_metrics_with_figures(
    coxFit_ae, ae_surv_final, ae_long,
    jointFit_ae, "AE-Only (No MTL)"
  )
  cat("\n✓ METHOD 6 (FIXED) COMPLETED\n")
  
}, error=function(e) cat(sprintf("\n⚠ ERROR in AE-Only: %s\n", e$message)))


# ============================================================
# METHOD 7/8: THESIS MAE-JM (LATENT ONLY)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 7: MAE-JM (LATENT ONLY)\n")
cat("============================================================\n")

# ============================================================
# METHOD 7: MAE-JM — FIXED
#
# KEY FIXES vs original:
#   1. Uses latent_patient_level.csv (z_final from LSTM) NOT per-visit first()
#   2. z_slope features included — trajectory signal preserved
#   3. LASSO replaced with elastic net + 1se rule + floor of 8 features
#   4. Joint model uses ALL visits from latent_improved_autoencoder.csv
#      for the longitudinal submodel (MMSE trajectory), not just baseline
#   5. prepare_latent_data() variance filter BYPASSED for z_final features
#      (they were already selected by the survival loss during training)
# ============================================================

metrics_thesis_latent <- NULL
tryCatch({
  
  # ----------------------------------------------------------
  # STEP 1: Load patient-level CSV (one row per patient, z_final + z_slope)
  # This file is generated by export_patient_level_latent() in Python.
  # ----------------------------------------------------------
  if (!file.exists("latent_patient_level.csv")) {
    stop(paste(
      "latent_patient_level.csv not found.",
      "Run export_patient_level_latent() in Python first and re-export.",
      "See python_export_fix.py for the function to add to your script."
    ))
  }
  
  patient_level <- read.csv("latent_patient_level.csv", stringsAsFactors = FALSE)
  patient_level <- patient_level %>% filter(PTID %in% master_patients)
  cat(sprintf("  Patient-level rows loaded: %d\n", nrow(patient_level)))
  
  # ----------------------------------------------------------
  # STEP 2: Identify z_final and z_slope feature columns
  # z_final = LSTM's integrated summary (what survival head was trained on)
  # z_slope = per-dim trajectory slope (captures AD progression signal)
  # ----------------------------------------------------------
  z_final_cols <- grep("^z_final_", names(patient_level), value = TRUE)
  z_slope_cols <- grep("^z_slope_", names(patient_level), value = TRUE)
  
  cat(sprintf("  z_final features: %d\n", length(z_final_cols)))
  cat(sprintf("  z_slope features: %d\n", length(z_slope_cols)))
  
  if (length(z_final_cols) == 0) {
    stop("No z_final_ columns found. Check that export_patient_level_latent() ran correctly.")
  }
  
  # ----------------------------------------------------------
  # STEP 3: Quality filter — remove near-zero variance only
  # Do NOT apply kurtosis/skewness filter here. These features were
  # selected by the survival loss; variance filtering throws away
  # exactly the features the network decided mattered.
  # ----------------------------------------------------------
  remove_zero_var <- function(cols, data) {
    keep <- sapply(cols, function(f) {
      v <- var(data[[f]], na.rm = TRUE)
      !is.na(v) && v > 1e-10
    })
    cols[keep]
  }
  
  z_final_cols <- remove_zero_var(z_final_cols, patient_level)
  z_slope_cols <- remove_zero_var(z_slope_cols, patient_level)
  
  # Standardise (Cox needs comparable scales; do NOT re-filter by variance after this)
  standardise_cols <- function(data, cols) {
    for (f in cols) {
      m <- mean(data[[f]], na.rm = TRUE)
      s <- sd(data[[f]],   na.rm = TRUE)
      if (s > 1e-10) data[[f]] <- (data[[f]] - m) / s
      data[[f]] <- pmin(pmax(data[[f]], -5), 5)  # mild winsorisation
    }
    data
  }
  
  patient_level <- standardise_cols(patient_level, c(z_final_cols, z_slope_cols))
  
  # ----------------------------------------------------------
  # STEP 4: Elastic net feature selection on z_final + z_slope combined
  # alpha=0.5 balances LASSO sparsity with ridge grouping (handles correlated latents)
  # lambda.1se keeps more features than lambda.min; floor ensures Cox has enough signal
  # ----------------------------------------------------------
  all_candidate_cols <- c(z_final_cols, z_slope_cols)
  
  surv_data_lasso <- patient_level %>%
    select(PTID, time_to_event, event, all_of(all_candidate_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0)
  surv_data_lasso <- surv_data_lasso[complete.cases(surv_data_lasso), ]
  
  x_mat <- as.matrix(surv_data_lasso[, all_candidate_cols])
  y_mat <- Surv(surv_data_lasso$time_to_event, surv_data_lasso$event)
  
  set.seed(42)
  cv_fit <- cv.glmnet(
    x_mat, y_mat,
    family        = "cox",
    alpha         = 0.5,       # elastic net, not pure LASSO
    nfolds        = 10,
    type.measure  = "C"
  )
  
  # Use 1se rule first (more conservative = more features retained)
  coefs_1se  <- coef(cv_fit, s = cv_fit$lambda.1se)
  selected_1se <- rownames(coefs_1se)[which(coefs_1se != 0)]
  
  # Separate z_final vs z_slope selected
  selected_z_final <- intersect(selected_1se, z_final_cols)
  selected_z_slope <- intersect(selected_1se, z_slope_cols)
  
  cat(sprintf("  Elastic net (1se) selected: %d z_final + %d z_slope\n",
              length(selected_z_final), length(selected_z_slope)))
  
  # Floor: if fewer than 8 z_final features, fall back to lambda.min
  if (length(selected_z_final) < 8) {
    coefs_min <- coef(cv_fit, s = cv_fit$lambda.min)
    selected_min <- rownames(coefs_min)[which(coefs_min != 0)]
    selected_z_final <- intersect(selected_min, z_final_cols)
    selected_z_slope <- intersect(selected_min, z_slope_cols)
    cat(sprintf("  Fell back to lambda.min: %d z_final + %d z_slope\n",
                length(selected_z_final), length(selected_z_slope)))
  }
  
  # Hard floor: if still fewer than 5, take top-N by univariate C-index
  if (length(selected_z_final) < 5) {
    cat("  ⚠ Very few features selected — using top-10 by univariate C-index\n")
    uni_c <- sapply(z_final_cols, function(f) {
      tryCatch({
        fit <- coxph(as.formula(paste("Surv(time_to_event, event) ~", f)),
                     data = surv_data_lasso)
        summary(fit)$concordance[1]
      }, error = function(e) 0.5)
    })
    selected_z_final <- names(sort(uni_c, decreasing = TRUE))[1:min(10, length(z_final_cols))]
    selected_z_slope <- z_slope_cols[1:min(5, length(z_slope_cols))]
  }
  
  # Final selected feature set
  final_features <- unique(c(selected_z_final, selected_z_slope))
  cat(sprintf("  FINAL feature set for Cox: %d features\n", length(final_features)))
  
  # ----------------------------------------------------------
  # STEP 5: Build survival dataset (one row per patient, z_final + z_slope)
  # ----------------------------------------------------------
  latent_surv <- patient_level %>%
    select(PTID, time_to_event, event, risk_score, all_of(final_features)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  latent_surv <- latent_surv[complete.cases(latent_surv), ]
  latent_surv <- latent_surv %>% filter(PTID %in% master_patients)
  
  cat(sprintf("  Survival dataset: %d patients, %d events\n",
              nrow(latent_surv), sum(latent_surv$event)))
  
  # ----------------------------------------------------------
  # STEP 6: Longitudinal dataset for joint model
  # Use ALL visits from latent_improved_autoencoder.csv (the per-visit CSV)
  # The joint model needs the MMSE trajectory across visits, not just one point.
  # We also add z_final features by joining from patient_level.
  # ----------------------------------------------------------
  latent_long_raw <- thesis_data %>%    # thesis_data = latent_improved_autoencoder.csv
    filter(PTID %in% latent_surv$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    arrange(PTID, Years_bl)
  
  # Join z_final features from patient_level (same value for all visits of a patient)
  latent_long <- latent_long_raw %>%
    left_join(
      patient_level %>% select(PTID, all_of(final_features)),
      by = "PTID"
    ) %>%
    select(PTID, Years_bl, MMSE, all_of(final_features)) %>%
    arrange(PTID, Years_bl)
  
  latent_long <- latent_long[complete.cases(latent_long), ]
  
  # Keep only patients present in both datasets
  common_IDs  <- intersect(latent_long$PTID, latent_surv$PTID)
  latent_long <- latent_long %>% filter(PTID %in% common_IDs) %>% arrange(PTID, Years_bl)
  latent_surv <- latent_surv %>% filter(PTID %in% common_IDs) %>% arrange(PTID)
  
  cat(sprintf("  Longitudinal dataset: %d visits across %d patients\n",
              nrow(latent_long), length(common_IDs)))
  
  # Verify >= 2 visits per patient (required for LME random slope)
  visit_counts <- latent_long %>% group_by(PTID) %>% summarize(n = n(), .groups = "drop")
  single_visit <- visit_counts$PTID[visit_counts$n < 2]
  if (length(single_visit) > 0) {
    cat(sprintf("  Removing %d patients with <2 visits\n", length(single_visit)))
    latent_long <- latent_long %>% filter(!PTID %in% single_visit)
    latent_surv <- latent_surv %>% filter(!PTID %in% single_visit)
    common_IDs  <- intersect(latent_long$PTID, latent_surv$PTID)
  }
  
  cat(sprintf("  Final: %d patients, %d events\n",
              length(common_IDs), sum(latent_surv$event)))
  
  # ----------------------------------------------------------
  # STEP 7: Fit LME + Cox + Joint Model
  # ----------------------------------------------------------
  
  # LME: MMSE ~ time + latent features (random slope on time)
  # z_final features explain inter-patient variation in MMSE trajectory
  long_formula <- as.formula(
    paste("MMSE ~ Years_bl +", paste(final_features, collapse = " + "))
  )
  surv_formula <- as.formula(
    paste("Surv(time, event) ~", paste(final_features, collapse = " + "))
  )
  
  cat("\n  Fitting LME...\n")
  lmeFit_thesis <- tryCatch(
    lme(long_formula,
        random  = ~ Years_bl | PTID,
        data    = latent_long,
        control = lmeControl(opt = "optim", maxIter = 500, returnObject = TRUE)),
    error = function(e) {
      cat(sprintf("  ⚠ LME with random slope failed (%s), trying random intercept\n", e$message))
      lme(long_formula,
          random  = ~ 1 | PTID,
          data    = latent_long,
          control = lmeControl(opt = "optim", maxIter = 500, returnObject = TRUE))
    }
  )
  cat("  ✓ LME fitted\n")
  
  cat("  Fitting Cox...\n")
  coxFit_thesis <- coxph(surv_formula, data = latent_surv, x = TRUE, method = "breslow")
  cat(sprintf("  ✓ Cox fitted — C-index: %.4f\n", summary(coxFit_thesis)$concordance[1]))
  
  # ----------------------------------------------------------
  # Calibration model (Platt scaling: logistic regression on linear predictor)
  # ----------------------------------------------------------
  raw_lp       <- predict(coxFit_thesis, newdata = latent_surv, type = "lp")
  surv_3yr_obs <- as.numeric(latent_surv$time > 3 |
                               (latent_surv$time <= 3 & latent_surv$event == 0))
  cal_glm      <- glm(surv_3yr_obs ~ raw_lp, family = binomial)
  coxFit_thesis$calibration_model  <- cal_glm
  coxFit_thesis$predict_calibrated <- function(nd) {
    lp_new <- predict(coxFit_thesis, newdata = nd, type = "lp")
    predict(cal_glm, newdata = data.frame(raw_lp = lp_new), type = "response")
  }
  cat("  ✓ Calibration model (Platt scaling) fitted\n")
  
  # ----------------------------------------------------------
  # Calibration before vs after plot
  # BEFORE: Cox baseline survival at 3 years converted to per-subject probability
  # AFTER:  Platt-scaled probabilities from cal_glm
  # ----------------------------------------------------------
  tryCatch({
    # BEFORE: proper Cox-derived 3-year survival probability per subject
    sf_before       <- survfit(coxFit_thesis, newdata = latent_surv)
    surv_3yr_before <- as.numeric(summary(sf_before, times = 3)$surv)
    
    # Sanity check — survfit matrix dims can vary by survival package version
    if (length(surv_3yr_before) != nrow(latent_surv)) {
      stop(sprintf(
        "survfit returned %d values but expected %d — check survival package version",
        length(surv_3yr_before), nrow(latent_surv)
      ))
    }
    
    # AFTER: Platt-scaled probabilities
    surv_3yr_after <- as.numeric(coxFit_thesis$predict_calibrated(latent_surv))
    
    cal_comparison <- data.frame(
      PTID         = latent_surv$PTID,
      surv_3yr_obs = surv_3yr_obs,
      before       = surv_3yr_before,
      after        = surv_3yr_after
    )
    
    cal_long <- tidyr::pivot_longer(
      cal_comparison,
      cols      = c(before, after),
      names_to  = "Method",
      values_to = "Predicted"
    ) %>%
      dplyr::mutate(Method = factor(Method,
                                    levels = c("before", "after"),
                                    labels = c("Before Calibration", "After Calibration")))
    
    gg_cal <- ggplot(cal_long, aes(x = Predicted, y = surv_3yr_obs)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
      geom_smooth(method = "loess", se = TRUE,
                  color = "steelblue", fill = "steelblue", alpha = 0.2) +
      facet_wrap(~ Method) +
      xlim(0, 1) + ylim(0, 1) +
      labs(
        title = "Calibration Before vs After Platt Scaling",
        x     = "Predicted 3-Year Survival Probability",
        y     = "Observed 3-Year Survival (binary)"
      ) +
      theme_minimal(base_size = 13)
    
    ggsave(
      "thesis_figures/MAE-JM__Fixed_/calibration_before_after.png",
      gg_cal, width = 12, height = 6, dpi = 300
    )
    
    write.csv(
      cal_comparison,
      "thesis_figures/MAE-JM__Fixed_/calibration_diagnostics.csv",
      row.names = FALSE
    )
    
    cat("  ✓ Calibration plot and diagnostics saved\n")
    
  }, error = function(e) {
    cat(sprintf("  ⚠ Calibration plot failed: %s\n", e$message))
  })
  
  cat("  Fitting Joint Model...\n")
  jointFit_thesis <- tryCatch(
    # Full model: value + slope of MMSE (captures both level and rate of decline)
    jm(coxFit_thesis, lmeFit_thesis, time_var = "Years_bl",
       functional_forms = ~ value(MMSE) + slope(MMSE),
       n_iter    = 10000,   # more iterations for stable estimates
       n_burnin  = 2000,
       n_thin    = 5,
       n_chains  = 3,       # 3 chains for better convergence diagnostics
       seed      = 42),
    error = function(e) {
      cat(sprintf("  ⚠ Full joint model failed (%s)\n  → Trying value(MMSE) only\n", e$message))
      tryCatch(
        jm(coxFit_thesis, lmeFit_thesis, time_var = "Years_bl",
           functional_forms = ~ value(MMSE),
           n_iter   = 7000, n_burnin = 1500, n_thin = 5, n_chains = 2, seed = 42),
        error = function(e2) {
          cat(sprintf("  ⚠ Simplified joint model also failed (%s)\n  → Returning NULL\n", e2$message))
          NULL
        }
      )
    }
  )
  
  if (!is.null(jointFit_thesis)) {
    cat("  ✓ Joint model fitted\n")
    tryCatch({
      jm_sum <- summary(jointFit_thesis)
      cat("  Joint model convergence summary:\n")
      print(jm_sum$Outcome1)
    }, error = function(e) cat(sprintf("  ⚠ Summary error: %s\n", e$message)))
  } else {
    cat("  ⚠ Joint model unavailable — will report Cox metrics only\n")
  }
  
  # ----------------------------------------------------------
  # STEP 8: Compute all metrics
  # ----------------------------------------------------------
  metrics_thesis_latent <- compute_all_metrics_with_figures(
    coxFit_thesis, latent_surv, latent_long,
    jointFit_thesis, "MAE-JM (Fixed)"
  )
  
  # ----------------------------------------------------------
  # STEP 9: Save feature importance for thesis writeup
  # ----------------------------------------------------------
  tryCatch({
    cox_coefs <- summary(coxFit_thesis)$coefficients
    coef_df   <- data.frame(
      Feature     = rownames(cox_coefs),
      Coefficient = cox_coefs[, "coef"],
      HR          = exp(cox_coefs[, "coef"]),
      HR_lower    = exp(cox_coefs[, "coef"] - 1.96 * cox_coefs[, "se(coef)"]),
      HR_upper    = exp(cox_coefs[, "coef"] + 1.96 * cox_coefs[, "se(coef)"]),
      P_value     = cox_coefs[, "Pr(>|z|)"],
      Feature_type = ifelse(grepl("^z_slope", rownames(cox_coefs)),
                            "Trajectory", "LSTM_Summary")
    )
    coef_df <- coef_df[order(coef_df$P_value), ]
    write.csv(coef_df,
              "thesis_figures/MAE-JM_Fixed/feature_importance.csv",
              row.names = FALSE)
    cat(sprintf("\n  Top predictors:\n"))
    print(head(coef_df[, c("Feature", "HR", "P_value", "Feature_type")], 10))
  }, error = function(e) cat(sprintf("  ⚠ Feature importance error: %s\n", e$message)))
  
  cat("\n✓ METHOD 7 (FIXED) COMPLETED\n")
  
}, error = function(e) {
  cat(sprintf("\n⚠ ERROR in MAE-JM Fixed: %s\n", e$message))
  cat("Traceback:\n")
  traceback()
})

# ============================================================
# FINAL COMPARISON TABLE
# ============================================================

cat("\n\n============================================================\n")
cat("FINAL COMPARISON\n")
cat("============================================================\n\n")

all_metrics <- Filter(Negate(is.null), list(
  metrics_baseline,
  metrics_image,
  metrics_tabular,
  metrics_concat_latent,
  metrics_noae_latent,
  metrics_no_mtl_latent,
  metrics_thesis_latent
))

if (length(all_metrics) > 0) {
  
  get_val <- function(x, f) { v <- x[[f]]; if (is.null(v)) NA else v }
  
  comparison <- data.frame(
    Method        = sapply(all_metrics, get_val, "name"),
    C_index       = sapply(all_metrics, get_val, "cindex"),
    CI_Lower      = sapply(all_metrics, get_val, "cindex_ci_lower"),
    CI_Upper      = sapply(all_metrics, get_val, "cindex_ci_upper"),
    Brier_3yr     = sapply(all_metrics, get_val, "brier_3yr"),
    Cal_Slope     = sapply(all_metrics, get_val, "calibration_slope"),
    RMST_Diff     = sapply(all_metrics, get_val, "rmst_diff"),
    RMST_pval     = sapply(all_metrics, get_val, "rmst_pval"),
    MAE_yr        = sapply(all_metrics, get_val, "mae_time"),
    JM_Dynamic_AUC = sapply(all_metrics, get_val, "jm_dynamic_auc"),
    JM_alpha_val  = sapply(all_metrics, get_val, "jm_association_val"),
    JM_alpha_slope = sapply(all_metrics, get_val, "jm_association_slope"),
    N_Patients    = sapply(all_metrics, get_val, "n_patients"),
    N_Events      = sapply(all_metrics, get_val, "n_events"),
    stringsAsFactors = FALSE
  )
  
  baseline_c <- comparison$C_index[comparison$Method == "Clinical Cox"]
  if (length(baseline_c) > 0 && !is.na(baseline_c)) {
    comparison$Improvement_Pct <- ((comparison$C_index - baseline_c) / baseline_c) * 100
  } else {
    comparison$Improvement_Pct <- NA
  }
  
  comparison <- comparison[order(-comparison$C_index), ]
  
  cat("DISCRIMINATION:\n")
  print(comparison[, c("Method", "C_index", "CI_Lower", "CI_Upper", "Improvement_Pct")])
  
  cat("\nCLINICAL UTILITY:\n")
  print(comparison[, c("Method", "Brier_3yr", "Cal_Slope", "RMST_Diff", "RMST_pval")])
  
  cat("\nJOINT MODEL ASSOCIATION PARAMETERS:\n")
  print(comparison[, c("Method", "JM_Dynamic_AUC", "JM_alpha_val", "JM_alpha_slope")])
  
  write.csv(comparison, "thesis_figures/FINAL_COMPARISON.csv", row.names = FALSE)
  
  # Publication table
  pub_table <- comparison %>%
    mutate(
      `C-index (95% CI)` = ifelse(
        !is.na(CI_Lower),
        sprintf("%.3f (%.3f-%.3f)", C_index, CI_Lower, CI_Upper),
        sprintf("%.3f", C_index)
      ),
      `Improvement`   = sprintf("%+.1f%%", Improvement_Pct),
      `Brier (3yr)`   = sprintf("%.3f", Brier_3yr),
      `Cal. Slope`    = sprintf("%.3f", Cal_Slope),
      `Dynamic AUC`   = sprintf("%.3f", JM_Dynamic_AUC),
      `alpha (value)` = ifelse(
        !is.na(JM_alpha_val),
        sprintf("%.3f", JM_alpha_val), "—"
      )
    ) %>%
    select(Method, `C-index (95% CI)`, Improvement, `Brier (3yr)`,
           `Cal. Slope`, `Dynamic AUC`, `alpha (value)`)
  
  write.csv(pub_table, "thesis_figures/PUBLICATION_TABLE.csv", row.names = FALSE)
  
  cat("\n\n============================================================\n")
  cat("✓ ANALYSIS COMPLETE\n")
  cat("============================================================\n")
  cat(sprintf("Methods evaluated: %d\n", nrow(comparison)))
  cat(sprintf("Cohort size: %d patients\n", unique(comparison$N_Patients[!is.na(comparison$N_Patients)])[1]))
  cat("\nOutputs:\n")
  cat("  thesis_figures/FINAL_COMPARISON.csv\n")
  cat("  thesis_figures/PUBLICATION_TABLE.csv\n")
  cat("  thesis_figures/ALL_METHODS_RESULTS.csv\n")
  cat("  thesis_figures/<Method>/ (figures + CSVs per method)\n")
  
  best <- comparison[1, ]
  cat(sprintf("\n🏆 BEST: %s | C-index=%.4f | Dynamic AUC=%.4f\n",
              best$Method, best$C_index,
              ifelse(is.na(best$JM_Dynamic_AUC), 0, best$JM_Dynamic_AUC)))
}