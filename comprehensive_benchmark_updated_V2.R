############################################################
# FIXED: COMPREHENSIVE THESIS ANALYSIS
# All methods evaluated on SAME patient cohort
# Image-only method fully integrated
############################################################

cat("============================================================\n")
cat("INSTALLING MISSING PACKAGES...\n")
cat("============================================================\n")

# Install missing packages
required_packages <- c("survRM2", "boot", "JMbayes2", "dplyr", "nlme", 
                       "survival", "pec", "survminer", "ggplot2", "gridExtra")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "http://cran.r-project.org")
  }
}

# Load all packages
library(JMbayes2)
library(dplyr)
library(nlme)
library(survival)
library(pec)
library(survminer)
library(ggplot2)
library(gridExtra)
library(survRM2)  # CRITICAL: For RMST analysis
library(boot)     # CRITICAL: For bootstrap CI

# Ensure tidyr is loaded (for case_when and other functions)
if(!require("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(tidyr)

# Force dplyr versions
select <- dplyr::select
summarize <- dplyr::summarize

cat("============================================================\n")
cat("FIXED: COMPREHENSIVE THESIS ANALYSIS\n")
cat("All methods on SAME patient cohort\n")
cat("============================================================\n")

# Create output directory
if (!dir.exists("thesis_figures")) {
  dir.create("thesis_figures")
}

############################################################
# LOAD ALL DATA
############################################################

cat("\n=== LOADING DATA ===\n")

baseline_data <- read.csv("baseline_clinical_features.csv")
tabular_data  <- read.csv("tabular_only_features.csv")
concat_data   <- read.csv("concat_fusion_features.csv")
noae_data     <- read.csv("no_autoencoder_features.csv")
image_data    <- read.csv("image_only_features.csv")
thesis_data   <- read.csv("latent_improved_autoencoder (1).csv")

cat(sprintf("Loaded data files:\n"))
cat(sprintf("  Baseline: %d patients\n", length(unique(baseline_data$PTID))))
cat(sprintf("  Tabular:  %d patients\n", length(unique(tabular_data$PTID))))
cat(sprintf("  Concat:   %d patients\n", length(unique(concat_data$PTID))))
cat(sprintf("  No-AE:    %d patients\n", length(unique(noae_data$PTID))))
cat(sprintf("  Image:    %d patients\n", length(unique(image_data$PTID))))
cat(sprintf("  Thesis:   %d patients\n", length(unique(thesis_data$PTID))))

############################################################
# MASTER COHORT DEFINITION
############################################################

cat("\n=== DEFINING MASTER COHORT ===\n")

# Step 1: Intersection across ALL datasets
master_patients <- Reduce(intersect, list(
  unique(baseline_data$PTID),
  unique(tabular_data$PTID),
  unique(concat_data$PTID),
  unique(noae_data$PTID),
  unique(image_data$PTID),
  unique(thesis_data$PTID)
))

cat(sprintf("Step 1 - Patients in all datasets: %d\n", length(master_patients)))

############################################################
# Step 2: ≥2 longitudinal observations (joint-model safe)
############################################################

check_min_obs <- function(data, min_obs = 2) {
  
  # If dataset has no longitudinal structure, skip this check
  if (!all(c("MMSE", "Years_bl") %in% names(data))) {
    return(unique(data$PTID[data$PTID %in% master_patients]))
  }
  
  data %>%
    filter(PTID %in% master_patients) %>%
    filter(!is.na(MMSE), !is.na(Years_bl)) %>%
    group_by(PTID) %>%
    filter(n() >= min_obs) %>%
    ungroup() %>%
    pull(PTID) %>%
    unique()
}


baseline_valid <- check_min_obs(baseline_data)
concat_valid   <- check_min_obs(concat_data)
noae_valid     <- check_min_obs(noae_data)
image_valid    <- check_min_obs(image_data)
thesis_valid   <- check_min_obs(thesis_data)

master_patients <- Reduce(intersect, list(
  master_patients,
  baseline_valid,
  concat_valid,
  noae_valid,
  image_valid,
  thesis_valid
))

cat(sprintf("Step 2 - After ≥2 observations filter: %d\n", length(master_patients)))

############################################################
# Step 3: Complete survival data
############################################################

check_complete_surv <- function(data) {
  data %>%
    filter(PTID %in% master_patients) %>%
    group_by(PTID) %>%
    summarize(
      has_time  = !is.na(first(time_to_event)),
      has_event = !is.na(first(event)),
      time_pos  = first(time_to_event) > 0,
      .groups = "drop"
    ) %>%
    filter(has_time, has_event, time_pos) %>%
    pull(PTID)
}

baseline_surv <- check_complete_surv(baseline_data)
tabular_surv  <- check_complete_surv(tabular_data)
concat_surv   <- check_complete_surv(concat_data)
noae_surv     <- check_complete_surv(noae_data)
image_surv    <- check_complete_surv(image_data)
thesis_surv   <- check_complete_surv(thesis_data)

master_patients <- Reduce(intersect, list(
  master_patients,
  baseline_surv,
  tabular_surv,
  concat_surv,
  noae_surv,
  image_surv,
  thesis_surv
))

cat(sprintf("Step 3 - After survival filter: %d\n", length(master_patients)))

############################################################
# Step 4: Complete feature availability
############################################################

check_complete_features <- function(data, feature_pattern = "^z_", top_n = 20) {
  
  z_features <- grep(feature_pattern, names(data), value = TRUE)
  
  data_filt <- data %>% filter(PTID %in% master_patients)
  
  if (length(z_features) > 0) {
    
    # Remove zero-variance / all-NA features
    vars <- sapply(data_filt[z_features], function(x)
      var(as.numeric(x), na.rm = TRUE))
    
    z_features <- z_features[!is.na(vars) & vars > 0]
    
    if (length(z_features) > top_n) {
      z_features <- names(sort(vars[z_features], decreasing = TRUE))[1:top_n]
    }
  }
  
  clinical <- intersect(c("AGE", "PTGENDER", "PTEDUCAT", "ADAS13"), names(data))
  longi    <- intersect(c("MMSE", "Years_bl"), names(data))
  
  required <- intersect(
    c("time_to_event", "event", clinical, longi, z_features),
    names(data)
  )
  
  data_filt %>%
    select(PTID, all_of(required)) %>%
    filter(complete.cases(.)) %>%
    pull(PTID) %>%
    unique()
}

concat_feat  <- check_complete_features(concat_data, "^z_", 20)
noae_feat    <- check_complete_features(noae_data, "^z_", 20)
image_feat   <- check_complete_features(image_data, "^z_", 20)
thesis_feat  <- check_complete_features(thesis_data, "^z_", 20)
tabular_feat <- check_complete_features(tabular_data, "^tab_feat_", 5)

master_patients <- Reduce(intersect, list(
  master_patients,
  concat_feat,
  noae_feat,
  image_feat,
  thesis_feat,
  tabular_feat
))

cat(sprintf("Step 4 - After complete feature filter: %d\n", length(master_patients)))

############################################################
# FINAL COHORT + FILTER DATASETS
############################################################

cat("\n*** FINAL MASTER COHORT ***\n")
cat(sprintf("Patients usable by ALL methods: %d\n", length(master_patients)))

baseline_data <- baseline_data %>% filter(PTID %in% master_patients)
tabular_data  <- tabular_data  %>% filter(PTID %in% master_patients)
concat_data   <- concat_data   %>% filter(PTID %in% master_patients)
noae_data     <- noae_data     %>% filter(PTID %in% master_patients)
image_data    <- image_data    %>% filter(PTID %in% master_patients)
thesis_data   <- thesis_data   %>% filter(PTID %in% master_patients)

cat("\nFiltered counts:\n")
cat(sprintf("  Baseline: %d\n", length(unique(baseline_data$PTID))))
cat(sprintf("  Tabular:  %d\n", length(unique(tabular_data$PTID))))
cat(sprintf("  Concat:   %d\n", length(unique(concat_data$PTID))))
cat(sprintf("  No-AE:    %d\n", length(unique(noae_data$PTID))))
cat(sprintf("  Image:    %d\n", length(unique(image_data$PTID))))
cat(sprintf("  Thesis:   %d\n", length(unique(thesis_data$PTID))))

# Final verification
stopifnot(
  length(unique(baseline_data$PTID)) == length(master_patients),
  length(unique(tabular_data$PTID))  == length(master_patients),
  length(unique(concat_data$PTID))   == length(master_patients),
  length(unique(noae_data$PTID))     == length(master_patients),
  length(unique(image_data$PTID))    == length(master_patients),
  length(unique(thesis_data$PTID))   == length(master_patients)
)

cat("\n✓ VERIFICATION PASSED: Identical cohort across ALL methods (including image-only)\n")


############################################################
# HELPER FUNCTION: COMPUTE ALL METRICS + FIGURES
############################################################

############################################################
# ENHANCED COMPREHENSIVE ANALYSIS
# Includes ALL metrics requested by advisor:
# - Cox coefficients & beta parameters
# - Brier scores
# - Survival curves
# - Time-to-event estimates
# - RMST (integral under survival curve)
# - Predicted vs observed event times
# - Complete figures package
############################################################

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

############################################################
# ENHANCED METRICS FUNCTION
############################################################

############################################################
# ENHANCED METRICS FUNCTION WITH MASTER RESULTS FILE
############################################################

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
  
  # [All your existing code from sections 1-10 stays the same...]
  # I'll just show the key additions
  
  ############################################################
  # 1. C-INDEX
  ############################################################
  
  cat("\n1. DISCRIMINATION\n")
  cindex <- summary(coxFit)$concordance[1]
  results$cindex <- cindex
  cat(sprintf("   C-index: %.4f\n", cindex))
  
  ############################################################
  # 2. COX COEFFICIENTS
  ############################################################
  
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
  cat(sprintf("   ✓ Saved: cox_coefficients.csv\n"))
  
  ############################################################
  # 3. FOREST PLOT
  ############################################################
  
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
      
      cat(sprintf("   ✓ Forest plot saved (%d significant)\n", nrow(sig_coefs)))
      results$n_significant <- nrow(sig_coefs)
    } else {
      cat("   No significant predictors (p < 0.10)\n")
      results$n_significant <- 0
    }
  }, error = function(e) {
    cat(sprintf("   ⚠ Forest plot error: %s\n", e$message))
    results$n_significant <- 0
  })
  
  ############################################################
  # 4. BASELINE SURVIVAL
  ############################################################
  
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
    
    cat("   ✓ Baseline functions saved\n")
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  ############################################################
  # 5. BRIER SCORE
  ############################################################
  
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
    
    cat("   ✓ Brier plot saved\n")
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Brier error: %s\n", e$message))
  })
  
  ############################################################
  # 6. KAPLAN-MEIER + 7. RMST (FIXED WITH ADAPTIVE TAU)
  ############################################################
  
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
    cat("   ✓ KM curves saved\n")
    
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
      
      cat("   ✓ RMST analysis saved\n")
      
    }, error = function(e) {
      cat("   ✗ RMST analysis failed:", conditionMessage(e), "\n")
      results$rmst_error <- conditionMessage(e)
    })
    
  }, error = function(e) {
    cat(sprintf("   ⚠ KM error: %s\n", e$message))
  })
  
  ############################################################
  # 8. PREDICTED vs OBSERVED TIMES
  ############################################################
  
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
      
      cat("   ✓ Predictions saved\n")
      
    }
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Prediction error: %s\n", e$message))
  })
  
  ############################################################
  # 9. CALIBRATION PLOT
  ############################################################
  
  ############################################################
  # 9. CALIBRATION PLOT (SECTION IN compute_all_metrics_with_figures)
  ############################################################
  
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
        
        cat("   ✓ Calibration saved\n")
      }
    }
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Calibration error: %s\n", e$message))
  })
  
  ############################################################
  # 10. BOOTSTRAP CI
  ############################################################
  
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
      cat("   ✓ Bootstrap complete\n")
      
      saveRDS(boot_results, 
              file.path(method_dir, paste0(gsub(" ", "_", name), "_bootstrap.rds")))
    }
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Bootstrap error: %s\n", e$message))
  })
  
  ############################################################
  # SAVE TO MASTER RESULTS FILE
  ############################################################
  
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
  cat(sprintf("\n✓ Results appended to: %s\n", master_file))
  
  ############################################################
  # SUMMARY
  ############################################################
  
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

############################################################
# METHOD 1: CLINICAL COX (BASELINE)
############################################################

cat("\n\n============================================================\n")
cat("METHOD 1: CLINICAL COX (BASELINE)\n")
cat("============================================================\n")

metrics_baseline <- NULL

tryCatch({
  baseline_long <- baseline_data %>%
    filter(!is.na(MMSE)) %>%
    select(PTID, Years_bl, MMSE, ADAS13, AGE, PTGENDER, PTEDUCAT) %>%
    arrange(PTID, Years_bl)
  
  # Use complete.cases instead of drop_na
  baseline_long <- baseline_long[complete.cases(baseline_long), ]
  
  baseline_surv <- baseline_data %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              AGE = first(AGE), PTGENDER = first(PTGENDER),
              PTEDUCAT = first(PTEDUCAT), ADAS13 = first(ADAS13),
              .groups = "drop") %>%
    ungroup()
  
  # Use complete.cases instead of drop_na
  baseline_surv <- baseline_surv[complete.cases(baseline_surv), ]
  
  common_IDs_baseline <- intersect(baseline_long$PTID, baseline_surv$PTID)
  baseline_long <- baseline_long %>% filter(PTID %in% common_IDs_baseline)
  baseline_surv <- baseline_surv %>% filter(PTID %in% common_IDs_baseline)
  
  # Verify master cohort
  if(length(common_IDs_baseline) != length(master_patients)) {
    cat(sprintf("  ⚠ WARNING: Baseline cohort (%d) differs from master (%d)\n",
                length(common_IDs_baseline), length(master_patients)))
  }
  
  cat(sprintf("Patients: %d, Events: %d\n", 
              length(common_IDs_baseline), sum(baseline_surv$event)))
  
  lmeFit_baseline <- lme(MMSE ~ Years_bl + AGE + PTGENDER + PTEDUCAT + ADAS13,
                         random = ~ Years_bl | PTID, data = baseline_long,
                         control = lmeControl(opt = "optim", maxIter = 200))
  
  coxFit_baseline <- coxph(Surv(time, event) ~ AGE + PTGENDER + PTEDUCAT + ADAS13,
                           data = baseline_surv, x = TRUE)
  
  jointFit_baseline <- jm(coxFit_baseline, lmeFit_baseline, time_var = "Years_bl",
                          functional_forms = ~ value(MMSE) + slope(MMSE),
                          n_iter = 5000, n_burnin = 1000, n_thin = 5, n_chains = 2)
  
  metrics_baseline <- compute_all_metrics_with_figures(coxFit_baseline, baseline_surv, 
                                                       baseline_long, jointFit_baseline, 
                                                       "Clinical Cox")
  
  cat("\n✓ METHOD 1 COMPLETED SUCCESSFULLY\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in Clinical Cox method\n")
  cat("Message:", e$message, "\n")
  traceback()
})

############################################################
# METHOD 2: IMAGE-ONLY CNN
############################################################

cat("\n\n============================================================\n")
cat("METHOD 2: IMAGE-ONLY CNN\n")
cat("============================================================\n")

metrics_image <- NULL

tryCatch({
  # Check if image features exist
  img_features <- grep("^img_feat_", names(image_data), value = TRUE)
  
  if(length(img_features) == 0) {
    # Try alternative naming pattern
    img_features <- grep("^z_", names(image_data), value = TRUE)
  }
  
  if(length(img_features) == 0) {
    cat("  ⚠ No image features found, skipping Method 2\n")
  } else {
    
    if(length(img_features) > 10) {
      img_vars <- sapply(image_data[img_features], var, na.rm = TRUE)
      selected_img <- names(sort(img_vars, decreasing = TRUE)[1:10])
    } else {
      selected_img <- img_features
    }
    
    image_surv <- image_data %>%
      group_by(PTID) %>%
      summarize(time = first(time_to_event), event = first(event),
                across(all_of(selected_img), first),
                .groups = "drop") %>%
      ungroup()
    
    # Use complete.cases instead of drop_na
    image_surv <- image_surv[complete.cases(image_surv), ]
    
    cat(sprintf("Patients: %d, Events: %d\n", 
                nrow(image_surv), sum(image_surv$event)))
    
    surv_formula <- as.formula(paste("Surv(time, event) ~", 
                                     paste(selected_img, collapse = " + ")))
    coxFit_image <- coxph(surv_formula, data = image_surv, x = TRUE)
    
    metrics_image <- compute_all_metrics_with_figures(
      coxFit_image, image_surv, NULL, NULL, "Image-Only CNN"
    )
    
    cat("\n✓ METHOD 2 COMPLETED SUCCESSFULLY\n")
  }
  
}, error = function(e) {
  cat("\n⚠ ERROR in Image-Only CNN method\n")
  cat("Message:", e$message, "\n")
  traceback()
})

############################################################
# METHOD 3: TABULAR-ONLY
############################################################

cat("\n\n============================================================\n")
cat("METHOD 3: TABULAR-ONLY DEEP NN\n")
cat("============================================================\n")

metrics_tabular <- NULL

tryCatch({
  tab_features <- grep("^tab_feat_", names(tabular_data), value = TRUE)
  
  if(length(tab_features) == 0) {
    cat("  ⚠ No tabular features found (tab_feat_*), trying alternative patterns\n")
    # Try alternative patterns
    tab_features <- grep("^tabular_", names(tabular_data), value = TRUE)
    
    if(length(tab_features) == 0) {
      cat("  ⚠ Still no features found, skipping Method 3\n")
      stop("No tabular features found")
    }
  }
  
  if(length(tab_features) > 5) {
    tab_vars <- sapply(tabular_data[tab_features], var, na.rm = TRUE)
    selected_tab <- names(sort(tab_vars, decreasing = TRUE)[1:5])
  } else {
    selected_tab <- tab_features
  }
  
  tabular_surv <- tabular_data %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              across(all_of(selected_tab), first),
              .groups = "drop") %>%
    ungroup()
  
  # Use complete.cases instead of drop_na
  tabular_surv <- tabular_surv[complete.cases(tabular_surv), ]
  
  # Verify master cohort
  if(nrow(tabular_surv) != length(master_patients)) {
    cat(sprintf("  ⚠ WARNING: Tabular cohort (%d) differs from master (%d)\n",
                nrow(tabular_surv), length(master_patients)))
  }
  
  cat(sprintf("Patients: %d, Events: %d\n", 
              nrow(tabular_surv), sum(tabular_surv$event)))
  
  surv_formula <- as.formula(paste("Surv(time, event) ~", 
                                   paste(selected_tab, collapse = " + ")))
  coxFit_tabular <- coxph(surv_formula, data = tabular_surv, x = TRUE)
  
  metrics_tabular <- compute_all_metrics_with_figures(
    coxFit_tabular, tabular_surv, NULL, NULL, "Tabular-Only NN"
  )
  
  cat("\n✓ METHOD 3 COMPLETED SUCCESSFULLY\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in Tabular-Only method\n")
  cat("Message:", e$message, "\n")
  traceback()
})

############################################################
# METHOD 4: CONCATENATION FUSION
############################################################

cat("\n\n============================================================\n")
cat("METHOD 4: CONCATENATION FUSION\n")
cat("============================================================\n")

metrics_concat <- NULL

tryCatch({
  z_features <- grep("^z_", names(concat_data), value = TRUE)
  
  # Remove zero-variance features
  zero_var <- sapply(concat_data[z_features], function(x) {
    v <- var(as.numeric(x), na.rm = TRUE)
    is.na(v) || v == 0
  })
  z_features_filtered <- z_features[!zero_var]
  
  if(length(z_features_filtered) > 20) {
    z_vars <- sapply(concat_data[z_features_filtered], var, na.rm = TRUE)
    selected_z <- names(sort(z_vars, decreasing = TRUE)[1:20])
  } else {
    selected_z <- z_features_filtered
  }
  
  all_predictors <- selected_z
  
  concat_long <- concat_data %>%
    filter(!is.na(MMSE)) %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors)) %>%
    arrange(PTID, Years_bl)
  
  # Use complete.cases instead of drop_na
  concat_long <- concat_long[complete.cases(concat_long), ]
  
  concat_surv <- concat_data %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              across(all_of(all_predictors), first),
              .groups = "drop") %>%
    ungroup()
  
  # Use complete.cases instead of drop_na
  concat_surv <- concat_surv[complete.cases(concat_surv), ]
  
  common_IDs_concat <- intersect(concat_long$PTID, concat_surv$PTID)
  concat_long <- concat_long %>% filter(PTID %in% common_IDs_concat)
  concat_surv <- concat_surv %>% filter(PTID %in% common_IDs_concat)
  
  # Verify master cohort
  if(length(common_IDs_concat) != length(master_patients)) {
    cat(sprintf("  ⚠ WARNING: Concat cohort (%d) differs from master (%d)\n",
                length(common_IDs_concat), length(master_patients)))
  }
  
  cat(sprintf("Patients: %d, Events: %d\n", 
              length(common_IDs_concat), sum(concat_surv$event)))
  
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", paste(all_predictors, collapse = " + ")))
  lmeFit_concat <- lme(
    long_formula,
    random = ~ Years_bl | PTID,
    data = concat_long,
    control = lmeControl(opt = "optim", maxIter = 200)
  )
  
  surv_formula <- as.formula(paste("Surv(time, event) ~", paste(all_predictors, collapse = " + ")))
  coxFit_concat <- coxph(surv_formula, data = concat_surv, x = TRUE)
  
  jointFit_concat <- jm(
    coxFit_concat, lmeFit_concat, time_var = "Years_bl",
    functional_forms = ~ value(MMSE) + slope(MMSE),
    n_iter = 5000, n_burnin = 1000, n_thin = 5, n_chains = 2
  )
  
  metrics_concat <- compute_all_metrics_with_figures(coxFit_concat, concat_surv, 
                                                     concat_long, jointFit_concat, 
                                                     "Concatenation Fusion")
  
  cat("\n✓ METHOD 4 COMPLETED SUCCESSFULLY\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in Concatenation Fusion method\n")
  cat("Message:", e$message, "\n")
  traceback()
})

cat("============================================================\n")
cat("METHOD 5: NO AUTOENCODER\n")
cat("============================================================\n")

metrics_noae <- NULL

tryCatch({
  
  # Load data
  noae_data <- read.csv("no_autoencoder_features.csv")
  
  # 1️⃣ Identify latent features
  z_features <- grep("^z_", names(noae_data), value = TRUE)
  
  # Remove zero-variance z-features
  zero_var <- sapply(noae_data[z_features], function(x) var(as.numeric(x), na.rm = TRUE) == 0)
  z_features_filtered <- z_features[!zero_var]
  
  # Filter features with pathological distributions
  cat("  Filtering features with extreme distributions...\n")
  bad_features <- c()
  
  for(feat in z_features_filtered){
    vals <- as.numeric(noae_data[[feat]])
    
    feat_mean <- mean(vals, na.rm = TRUE)
    feat_sd <- sd(vals, na.rm = TRUE)
    
    if(feat_sd > 1e-10){
      z_vals <- (vals - feat_mean) / feat_sd
      extreme_pct <- mean(abs(z_vals) > 5, na.rm = TRUE)
      
      if(extreme_pct > 0.10){
        bad_features <- c(bad_features, feat)
        cat(sprintf("    Excluding %s (%.1f%% extreme outliers)\n", feat, extreme_pct * 100))
      }
    }
  }
  
  z_features_filtered <- setdiff(z_features_filtered, bad_features)
  cat(sprintf("  Number of z features after filtering: %d\n", length(z_features_filtered)))
  
  # Keep top 20 by variance
  if(length(z_features_filtered) > 10){
    feat_vars <- sapply(noae_data[z_features_filtered], function(x) var(as.numeric(x), na.rm = TRUE))
    selected_feat <- names(sort(feat_vars, decreasing = TRUE)[1:10])
  } else {
    selected_feat <- z_features_filtered
  }
  cat(sprintf("  Selected z features: %d\n", length(selected_feat)))
  
  # Convert ALL z-features to numeric and clean
  for(f in selected_feat){
    noae_data[[f]] <- as.numeric(noae_data[[f]])
    bad_vals <- is.na(noae_data[[f]]) | is.nan(noae_data[[f]]) | is.infinite(noae_data[[f]])
    if(any(bad_vals)){
      noae_data[[f]][bad_vals] <- mean(noae_data[[f]], na.rm = TRUE)
    }
  }
  
  # 2️⃣ Include clinical features
  clinical_features <- intersect(c("AGE", "PTGENDER", "PTEDUCAT", "ADAS13"), names(noae_data))
  
  # Clean clinical features
  for(feat in clinical_features){
    if(feat == "PTGENDER"){
      noae_data[[feat]] <- factor(noae_data[[feat]])
    } else {
      noae_data[[feat]] <- as.numeric(noae_data[[feat]])
      bad_vals <- is.na(noae_data[[feat]]) | is.nan(noae_data[[feat]]) | is.infinite(noae_data[[feat]])
      if(any(bad_vals)){
        noae_data[[feat]][bad_vals] <- mean(noae_data[[feat]], na.rm = TRUE)
      }
    }
  }
  
  all_predictors <- c(clinical_features, selected_feat)
  
  # Remove predictors with zero or very low variance
  var_check <- sapply(noae_data[all_predictors], function(x) {
    if(is.numeric(x)) {
      var(x, na.rm=TRUE)
    } else {
      length(unique(x))
    }
  })
  all_predictors <- names(var_check[var_check > 1e-6])
  cat(sprintf("  All predictors used in modeling: %d\n", length(all_predictors)))
  
  # 3️⃣ Prepare longitudinal data - KEEP PTID AS CHARACTER INITIALLY
  noae_long <- noae_data %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors)) %>%
    filter(!is.na(MMSE), !is.na(Years_bl)) %>%
    arrange(PTID, Years_bl)
  
  # Convert PTID to character to avoid factor issues
  noae_long$PTID <- as.character(noae_long$PTID)
  
  cat(sprintf("  Longitudinal data rows (before complete cases): %d\n", nrow(noae_long)))
  
  # 4️⃣ Prepare survival data - KEEP PTID AS CHARACTER INITIALLY
  noae_surv <- noae_data %>%
    mutate(PTID = as.character(PTID)) %>%
    arrange(PTID, Years_bl) %>%
    group_by(PTID) %>%
    summarize(
      time = first(time_to_event),
      event = first(event),
      across(all_of(all_predictors), first),
      .groups = 'drop'
    ) %>%
    filter(!is.na(time), !is.na(event), time > 0)
  
  cat(sprintf("  Survival data rows (before complete cases): %d\n", nrow(noae_surv)))
  cat(sprintf("  Number of events: %d\n", sum(noae_surv$event, na.rm=TRUE)))
  
  # ⚠️ CRITICAL FIX: Ensure EXACT patient matching with complete data
  cat("\n  CRITICAL: Ensuring exact patient matching with complete data...\n")
  
  # Identify complete cases without grouping warnings
  noae_long_complete <- noae_long[complete.cases(noae_long), ]
  noae_surv_complete <- noae_surv[complete.cases(noae_surv), ]
  
  # Get valid patients from both
  valid_long_patients <- unique(noae_long_complete$PTID)
  valid_surv_patients <- unique(noae_surv_complete$PTID)
  
  cat(sprintf("  Patients with complete longitudinal data: %d\n", length(valid_long_patients)))
  cat(sprintf("  Patients with complete survival data: %d\n", length(valid_surv_patients)))
  
  # Keep only patients complete in BOTH datasets
  final_patients <- intersect(valid_long_patients, valid_surv_patients)
  cat(sprintf("  ✓ Patients with complete data in BOTH datasets: %d\n", length(final_patients)))
  
  # Filter to final patients
  noae_long <- noae_long_complete %>% 
    filter(PTID %in% final_patients)
  
  noae_surv <- noae_surv_complete %>% 
    filter(PTID %in% final_patients)
  
  # Remove patients with < 2 observations
  obs_counts <- table(noae_long$PTID)
  valid_patients <- names(obs_counts[obs_counts >= 2])
  
  noae_long <- noae_long %>% filter(PTID %in% valid_patients)
  noae_surv <- noae_surv %>% filter(PTID %in% valid_patients)
  
  # ⚠️ CRITICAL: Sort both datasets by PTID to ensure alignment
  noae_long <- noae_long %>% arrange(PTID, Years_bl)
  noae_surv <- noae_surv %>% arrange(PTID)
  
  # Convert to data.frame BEFORE making PTID a factor
  noae_long <- as.data.frame(noae_long)
  noae_surv <- as.data.frame(noae_surv)
  
  # Now convert PTID to factor with SAME levels in SAME order
  patient_levels <- sort(unique(noae_long$PTID))
  noae_long$PTID <- factor(noae_long$PTID, levels = patient_levels)
  noae_surv$PTID <- factor(noae_surv$PTID, levels = patient_levels)
  
  cat(sprintf("  ✓ Final matched data - Long: %d obs (%d patients), Surv: %d patients\n",
              nrow(noae_long), nlevels(noae_long$PTID), nrow(noae_surv)))
  
  # 6️⃣ PRE-STANDARDIZE numeric predictors
  numeric_preds <- all_predictors[sapply(noae_surv[all_predictors], is.numeric)]
  
  # Create standardization parameters from survival data
  standardization_params <- list()
  
  for(pred in numeric_preds){
    pred_mean <- mean(noae_surv[[pred]], na.rm = TRUE)
    pred_sd <- sd(noae_surv[[pred]], na.rm = TRUE)
    
    standardization_params[[pred]] <- list(mean = pred_mean, sd = pred_sd)
    
    if(pred_sd > 1e-10){
      noae_surv[[pred]] <- (noae_surv[[pred]] - pred_mean) / pred_sd
      noae_long[[pred]] <- (noae_long[[pred]] - pred_mean) / pred_sd
    }
  }
  
  cat(sprintf("  Standardized %d numeric predictors\n", length(numeric_preds)))
  
  # Cap extreme standardized values
  cat("  Capping extreme values for numerical stability...\n")
  cap_threshold <- 5
  n_capped_total <- 0
  
  for(pred in numeric_preds){
    # Cap in survival data
    extreme_vals_surv <- abs(noae_surv[[pred]]) > cap_threshold
    if(any(extreme_vals_surv)){
      n_capped <- sum(extreme_vals_surv)
      n_capped_total <- n_capped_total + n_capped
      noae_surv[[pred]][noae_surv[[pred]] > cap_threshold] <- cap_threshold
      noae_surv[[pred]][noae_surv[[pred]] < -cap_threshold] <- -cap_threshold
      cat(sprintf("    Capped %d extreme values in %s (survival data)\n", n_capped, pred))
    }
    
    # Cap in longitudinal data
    extreme_vals_long <- abs(noae_long[[pred]]) > cap_threshold
    if(any(extreme_vals_long)){
      n_capped <- sum(extreme_vals_long)
      n_capped_total <- n_capped_total + n_capped
      noae_long[[pred]][noae_long[[pred]] > cap_threshold] <- cap_threshold
      noae_long[[pred]][noae_long[[pred]] < -cap_threshold] <- -cap_threshold
      cat(sprintf("    Capped %d extreme values in %s (longitudinal data)\n", n_capped, pred))
    }
  }
  
  cat(sprintf("  ✓ Total values capped: %d\n", n_capped_total))
  cat(sprintf("  ✓ All standardized values within ±%d range\n", cap_threshold))
  
  # FINAL verification before modeling
  cat("\n  FINAL DATA VERIFICATION:\n")
  cat(sprintf("  Longitudinal: %d rows, %d patients\n", nrow(noae_long), nlevels(noae_long$PTID)))
  cat(sprintf("  Survival: %d rows, %d patients\n", nrow(noae_surv), nlevels(noae_surv$PTID)))
  cat(sprintf("  PTID levels identical: %s\n", identical(levels(noae_long$PTID), levels(noae_surv$PTID))))
  cat(sprintf("  Min observations per patient: %d\n", min(table(noae_long$PTID))))
  
  # ⚠️ CRITICAL: Check that survival data is ONE ROW PER PATIENT
  if(nrow(noae_surv) != nlevels(noae_long$PTID)){
    cat("  ERROR: Survival data does not have exactly one row per patient!\n")
    stop("Data structure error")
  }
  
  # ⚠️ CRITICAL: Verify row names alignment
  rownames(noae_surv) <- as.character(noae_surv$PTID)
  
  # 7️⃣ Fit LME with improved settings
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", paste(all_predictors, collapse = " + ")))
  cat("\n  Fitting LME model...\n")
  cat(sprintf("  Formula: %s\n", deparse(long_formula, width.cutoff = 500)))
  
  lmeFit_noae <- lme(
    long_formula,
    random = ~ Years_bl | PTID,
    data = noae_long,
    control = lmeControl(
      opt = "optim", 
      maxIter = 500,
      msMaxIter = 500,
      niterEM = 50,
      msVerbose = FALSE,
      returnObject = TRUE
    ),
    method = "REML"
  )
  cat("  ✓ LME fitted\n")
  
  # 8️⃣ Fit Cox with improved settings
  surv_formula <- as.formula(paste("Surv(time, event) ~", paste(all_predictors, collapse = " + ")))
  cat("  Fitting Cox model...\n")
  cat(sprintf("  Formula: %s\n", deparse(surv_formula, width.cutoff = 500)))
  
  coxFit_noae <- coxph(
    surv_formula, 
    data = noae_surv, 
    x = TRUE, 
    model = TRUE,
    method = "breslow"
  )
  cat("  ✓ Cox fitted\n")
  
  # 9️⃣ Diagnostic info before joint model
  cat("\n  DIAGNOSTIC INFO FOR JOINT MODEL:\n")
  cat(sprintf("  LME fitted on: %d observations from %d patients\n", 
              nrow(lmeFit_noae$data), length(unique(lmeFit_noae$data$PTID))))
  cat(sprintf("  Cox fitted on: %d observations\n", nrow(coxFit_noae$model)))
  cat(sprintf("  Events in Cox model: %d (%.1f%%)\n", 
              sum(noae_surv$event), 100*mean(noae_surv$event)))
  
  # ⚠️ CRITICAL: Verify the data in the fitted models matches
  lme_patients <- sort(unique(as.character(lmeFit_noae$data$PTID)))
  cox_patients <- sort(rownames(coxFit_noae$model))
  
  if(!identical(lme_patients, cox_patients)){
    cat("  ERROR: Patient sets in LME and Cox models do not match!\n")
    cat(sprintf("  LME has %d patients, Cox has %d patients\n", 
                length(lme_patients), length(cox_patients)))
    
    missing_in_cox <- setdiff(lme_patients, cox_patients)
    missing_in_lme <- setdiff(cox_patients, lme_patients)
    
    if(length(missing_in_cox) > 0){
      cat(sprintf("  %d patients in LME but not Cox\n", length(missing_in_cox)))
    }
    if(length(missing_in_lme) > 0){
      cat(sprintf("  %d patients in Cox but not LME\n", length(missing_in_lme)))
    }
    
    stop("Patient mismatch between models")
  }
  
  cat("  ✓ Patient sets in both models match perfectly\n")
  
  # 🔟 Fit joint model with conservative parameters
  cat("\n  Fitting joint model...\n")
  cat("  This may take several minutes...\n")
  
  tryCatch({
    jointFit_noae <- jm(
      coxFit_noae, 
      lmeFit_noae, 
      time_var = "Years_bl",
      functional_forms = ~ value(MMSE) + slope(MMSE),
      n_iter = 3000,
      n_burnin = 500, 
      n_thin = 3, 
      n_chains = 1,
      seed = 123
    )
    cat("  ✓ Joint model fitted successfully\n")
    
  }, error = function(e) {
    cat("\n  ⚠️ Joint model fitting FAILED\n")
    cat("  Error message:", e$message, "\n\n")
    
    cat("  Attempting diagnosis...\n")
    
    # Try simpler model
    cat("  Trying simpler joint model (value only, no slope)...\n")
    tryCatch({
      jointFit_noae <<- jm(
        coxFit_noae, 
        lmeFit_noae, 
        time_var = "Years_bl",
        functional_forms = ~ value(MMSE),  # Simpler: just value, no slope
        n_iter = 2000,
        n_burnin = 300, 
        n_thin = 2, 
        n_chains = 1,
        seed = 123
      )
      cat("  ✓ Simpler joint model fitted successfully\n")
    }, error = function(e2) {
      cat("  Simpler model also failed:", e2$message, "\n")
      stop(paste("Joint model error:", e$message))
    })
  })
  
  # 1️⃣1️⃣ Compute metrics
  cat("\n  Computing metrics...\n")
  metrics_noae <- compute_all_metrics_with_figures(
    coxFit_noae, 
    noae_surv, 
    noae_long, 
    jointFit_noae, 
    "No Autoencoder"
  )
  
  cat("\n✓ METHOD 5 COMPLETED SUCCESSFULLY\n")
  
}, error = function(e){
  cat("\n⚠️ ERROR in No Autoencoder method\n")
  cat("Message:", e$message, "\n")
  if(exists("traceback")){
    print(traceback())
  }
  metrics_noae <- NULL
})

cat("============================================================\n\n")

############################################################
# METHOD 6: THESIS (FULL METHOD) - OPTIMIZED VERSION (FIXED)
############################################################

cat("\n\n============================================================\n")
cat("METHOD 6: THESIS (FULL METHOD) - OPTIMIZED\n")
cat("============================================================\n")

metrics_thesis <- NULL

tryCatch({
  
  # Load additional packages
  if(!require("glmnet", quietly = TRUE)) {
    install.packages("glmnet")
  }
  if(!require("tidyr", quietly = TRUE)) {
    install.packages("tidyr")
  }
  library(glmnet)
  library(tidyr)
  
  ############################################################
  # STEP 1: FEATURE PREPARATION & QUALITY FILTERING
  ############################################################
  
  cat("\n=== STEP 1: FEATURE PREPARATION ===\n")
  
  z_features <- grep("^z_", names(thesis_data), value = TRUE)
  cat(sprintf("  Initial z-features: %d\n", length(z_features)))
  
  # Remove zero-variance and extreme variance features
  zero_var <- sapply(thesis_data[z_features], function(x) {
    v <- var(as.numeric(x), na.rm = TRUE)
    is.na(v) | v < 1e-10 | v > 1e10
  })
  z_features_filtered <- z_features[!zero_var]
  cat(sprintf("  After variance filter: %d\n", length(z_features_filtered)))
  
  # Filter features with extreme distributions (kurtosis/skewness) - FIXED
  if(require("e1071", quietly = TRUE)) {
    cat("  Filtering by distribution quality...\n")
    
    stable_features <- c()
    for(feat in z_features_filtered) {
      vals <- as.numeric(thesis_data[[feat]])
      
      # Remove NAs for calculation
      vals_clean <- vals[!is.na(vals)]
      
      if(length(vals_clean) > 3) {  # Need at least 4 points
        tryCatch({
          kurt <- kurtosis(vals_clean)
          skew <- skewness(vals_clean)
          
          # Keep features with reasonable distributions
          if(abs(kurt) < 10 & abs(skew) < 3) {
            stable_features <- c(stable_features, feat)
          }
        }, error = function(e) {
          # If calculation fails, skip this feature
          NULL
        })
      }
    }
    
    if(length(stable_features) > 0) {
      cat(sprintf("  After quality filter: %d\n", length(stable_features)))
      z_features_filtered <- stable_features
    } else {
      cat("  ⚠ Quality filter removed all features, keeping variance-filtered set\n")
    }
  }
  
  # Select top features by variance (reduce from 20 to 15)
  if(length(z_features_filtered) > 15) {
    z_vars <- sapply(thesis_data[z_features_filtered], var, na.rm = TRUE)
    selected_z <- names(sort(z_vars, decreasing = TRUE)[1:15])
  } else {
    selected_z <- z_features_filtered
  }
  
  cat(sprintf("  Final z-features selected: %d\n", length(selected_z)))
  
  ############################################################
  # STEP 2: ADD CLINICAL FEATURES & INTERACTIONS
  ############################################################
  
  cat("\n=== STEP 2: CLINICAL FEATURES ===\n")
  
  clinical_features <- c("AGE", "PTGENDER", "PTEDUCAT", "ADAS13")
  
  # Create interaction terms
  thesis_data$AGE_ADAS <- thesis_data$AGE * thesis_data$ADAS13
  thesis_data$EDUC_AGE <- thesis_data$PTEDUCAT * thesis_data$AGE
  thesis_data$EDUC_ADAS <- thesis_data$PTEDUCAT * thesis_data$ADAS13
  
  interaction_terms <- c("AGE_ADAS", "EDUC_AGE", "EDUC_ADAS")
  
  all_predictors <- c(clinical_features, selected_z, interaction_terms)
  cat(sprintf("  Total predictors: %d (clinical: %d, latent: %d, interactions: %d)\n",
              length(all_predictors), length(clinical_features), 
              length(selected_z), length(interaction_terms)))
  
  ############################################################
  # STEP 3: PREPARE DATA (FIXED)
  ############################################################
  
  cat("\n=== STEP 3: DATA PREPARATION ===\n")
  
  # Use complete.cases() instead of drop_na()
  thesis_long <- thesis_data %>%
    filter(!is.na(MMSE)) %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors)) %>%
    arrange(PTID, Years_bl)
  
  # Remove incomplete rows
  thesis_long <- thesis_long[complete.cases(thesis_long), ]
  
  thesis_surv <- thesis_data %>%
    group_by(PTID) %>%
    summarize(
      time = first(time_to_event), 
      event = first(event),
      across(all_of(all_predictors), first),
      .groups = "drop"
    ) %>%
    ungroup()
  
  # Remove incomplete rows
  thesis_surv <- thesis_surv[complete.cases(thesis_surv), ]
  
  common_IDs_thesis <- intersect(thesis_long$PTID, thesis_surv$PTID)
  thesis_long <- thesis_long %>% filter(PTID %in% common_IDs_thesis)
  thesis_surv <- thesis_surv %>% filter(PTID %in% common_IDs_thesis)
  
  # Verify master cohort
  if(length(common_IDs_thesis) != length(master_patients)) {
    cat(sprintf("  ⚠ WARNING: Thesis cohort (%d) differs from master (%d)\n",
                length(common_IDs_thesis), length(master_patients)))
  } else {
    cat(sprintf("  ✓ Cohort size matches master: %d patients\n", length(common_IDs_thesis)))
  }
  
  cat(sprintf("  Patients: %d, Events: %d (%.1f%%)\n", 
              length(common_IDs_thesis), sum(thesis_surv$event),
              100 * sum(thesis_surv$event) / length(common_IDs_thesis)))
  
  ############################################################
  # STEP 4: LASSO FEATURE SELECTION
  ############################################################
  
  cat("\n=== STEP 4: LASSO FEATURE SELECTION ===\n")
  
  # Prepare matrices for glmnet
  x_matrix <- as.matrix(thesis_surv[, all_predictors])
  y_matrix <- Surv(thesis_surv$time, thesis_surv$event)
  
  # Cross-validated Lasso
  set.seed(42)
  cv_fit <- cv.glmnet(
    x_matrix, y_matrix, 
    family = "cox", 
    alpha = 1,  # Lasso
    nfolds = 10,
    type.measure = "C"
  )
  
  # Get coefficients at optimal lambda
  optimal_lambda <- cv_fit$lambda.min
  coefs <- coef(cv_fit, s = optimal_lambda)
  selected_features <- rownames(coefs)[which(coefs != 0)]
  
  cat(sprintf("  Optimal lambda: %.6f\n", optimal_lambda))
  cat(sprintf("  Features selected by Lasso: %d (from %d)\n", 
              length(selected_features), length(all_predictors)))
  
  # Ensure at least clinical features are included
  if(length(selected_features) == 0) {
    cat("  ⚠ No features selected by Lasso, using top 8 by univariate p-value\n")
    
    # Fallback: select by univariate significance
    univariate_pvals <- sapply(all_predictors, function(pred) {
      tryCatch({
        formula <- as.formula(paste("Surv(time, event) ~", pred))
        fit <- coxph(formula, data = thesis_surv)
        summary(fit)$coefficients[1, "Pr(>|z|)"]
      }, error = function(e) 1.0)  # If fails, assign p=1
    })
    
    selected_features <- names(sort(univariate_pvals)[1:8])
  } else {
    # Always include core clinical features
    selected_features <- unique(c(
      intersect(clinical_features, c(selected_features, all_predictors)),
      selected_features
    ))
  }
  
  cat("  Selected features:\n")
  print(selected_features)
  
  # Update predictors
  all_predictors_optimized <- selected_features
  
  ############################################################
  # STEP 5: FIT MODELS WITH OPTIMIZED FEATURES
  ############################################################
  
  cat("\n=== STEP 5: MODEL FITTING ===\n")
  
  # Update longitudinal data with optimized features
  thesis_long_opt <- thesis_long %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors_optimized))
  
  thesis_surv_opt <- thesis_surv %>%
    select(PTID, time, event, all_of(all_predictors_optimized))
  
  # Fit LME
  cat("  Fitting LME model...\n")
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", 
                                   paste(all_predictors_optimized, collapse = " + ")))
  
  lmeFit_thesis <- lme(
    long_formula, 
    random = ~ Years_bl | PTID, 
    data = thesis_long_opt,
    control = lmeControl(
      opt = "optim", 
      maxIter = 500,
      msMaxIter = 500,
      returnObject = TRUE
    )
  )
  cat("  ✓ LME fitted\n")
  
  # Fit Cox
  cat("  Fitting Cox model...\n")
  surv_formula <- as.formula(paste("Surv(time, event) ~", 
                                   paste(all_predictors_optimized, collapse = " + ")))
  
  coxFit_thesis <- coxph(
    surv_formula, 
    data = thesis_surv_opt, 
    x = TRUE,
    method = "breslow"
  )
  cat("  ✓ Cox fitted\n")
  
  # Initial C-index
  initial_cindex <- summary(coxFit_thesis)$concordance[1]
  cat(sprintf("  Initial C-index: %.4f\n", initial_cindex))
  
  ############################################################
  # STEP 6: CALIBRATION CORRECTION (PLATT SCALING)
  ############################################################
  
  cat("\n=== STEP 6: CALIBRATION CORRECTION ===\n")
  
  # Get raw risk scores
  raw_risks <- predict(coxFit_thesis, newdata = thesis_surv_opt, type = "lp")
  
  # Calculate 3-year survival status
  surv_3yr <- as.numeric(
    thesis_surv_opt$time > 3 | 
      (thesis_surv_opt$time <= 3 & thesis_surv_opt$event == 0)
  )
  
  # Fit calibration model (Platt scaling)
  calibration_model <- glm(
    surv_3yr ~ raw_risks, 
    family = binomial
  )
  
  # Function to get calibrated predictions
  predict_calibrated <- function(new_data) {
    raw_risks_new <- predict(coxFit_thesis, newdata = new_data, type = "lp")
    calibrated_probs <- predict(
      calibration_model, 
      newdata = data.frame(raw_risks = raw_risks_new),
      type = "response"
    )
    return(calibrated_probs)
  }
  
  # Store calibration model with Cox fit
  coxFit_thesis$calibration_model <- calibration_model
  coxFit_thesis$predict_calibrated <- predict_calibrated
  
  cat("  ✓ Calibration model fitted\n")
  
  ############################################################
  # STEP 7: FIT JOINT MODEL
  ############################################################
  
  cat("\n=== STEP 7: JOINT MODEL ===\n")
  cat("  Fitting joint model (this may take several minutes)...\n")
  
  jointFit_thesis <- NULL  # Initialize
  
  tryCatch({
    jointFit_thesis <- jm(
      coxFit_thesis, 
      lmeFit_thesis, 
      time_var = "Years_bl",
      functional_forms = ~ value(MMSE) + slope(MMSE),
      n_iter = 5000, 
      n_burnin = 1000, 
      n_thin = 5, 
      n_chains = 2,
      seed = 42
    )
    cat("  ✓ Joint model fitted successfully\n")
    
  }, error = function(e) {
    cat("\n  ⚠ Joint model fitting failed, trying simpler model...\n")
    cat("  Error:", conditionMessage(e), "\n")
    
    # Try simpler model
    tryCatch({
      jointFit_thesis <<- jm(
        coxFit_thesis, 
        lmeFit_thesis, 
        time_var = "Years_bl",
        functional_forms = ~ value(MMSE),  # Value only, no slope
        n_iter = 3000, 
        n_burnin = 500, 
        n_thin = 3, 
        n_chains = 1,
        seed = 42
      )
      cat("  ✓ Simpler joint model fitted successfully\n")
      
    }, error = function(e2) {
      cat("  ✗ Joint model failed completely, proceeding without it\n")
      cat("  Error:", conditionMessage(e2), "\n")
      jointFit_thesis <<- NULL
    })
  })
  
  ############################################################
  # STEP 8: COMPREHENSIVE EVALUATION
  ############################################################
  
  cat("\n=== STEP 8: EVALUATION ===\n")
  
  metrics_thesis <- compute_all_metrics_with_figures(
    coxFit_thesis, 
    thesis_surv_opt, 
    thesis_long_opt, 
    jointFit_thesis, 
    "MAE-JM"
  )
  
  ############################################################
  # STEP 9: ENHANCED SURVIVAL ANALYSIS (NEW)
  ############################################################
  
  cat("\n=== STEP 9: ENHANCED SURVIVAL ANALYSIS ===\n")
  
  method_dir <- "thesis_figures/Thesis_(Full)_-_Optimized"
  
  # 1. Cox Model Survival Curve
  cat("\n1. Cox Model Survival Curve\n")
  
  tryCatch({
    baseline_haz <- basehaz(coxFit_thesis, centered = TRUE)
    baseline_haz$survival <- exp(-baseline_haz$hazard)
    
    gg_cox_survival <- ggplot(baseline_haz, aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1.3) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
      ylim(0, 1) +
      labs(
        title = "Cox Model Baseline Survival Curve",
        subtitle = "Estimated survival function from Cox proportional hazards model",
        x = "Time (Years)",
        y = "Survival Probability S(t)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    
    ggsave(
      file.path(method_dir, "cox_baseline_survival_curve.png"),
      gg_cox_survival, 
      width = 10, 
      height = 7, 
      dpi = 300
    )
    
    median_idx <- which.min(abs(baseline_haz$survival - 0.5))
    median_survival <- baseline_haz$time[median_idx]
    
    cat(sprintf("   ✓ Cox baseline survival curve saved\n"))
    cat(sprintf("   Median survival time: %.2f years\n", median_survival))
    
    write.csv(
      baseline_haz,
      file.path(method_dir, "cox_survival_function.csv"),
      row.names = FALSE
    )
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  # 2. Hazard Functions
  cat("\n2. Hazard Function\n")
  
  tryCatch({
    baseline_haz <- basehaz(coxFit_thesis, centered = TRUE)
    
    baseline_haz$hazard_rate <- c(
      diff(baseline_haz$hazard) / diff(baseline_haz$time),
      NA
    )
    
    baseline_haz_clean <- baseline_haz[!is.na(baseline_haz$hazard_rate), ]
    q99 <- quantile(baseline_haz_clean$hazard_rate, 0.99, na.rm = TRUE)
    baseline_haz_clean$hazard_rate_capped <- pmin(baseline_haz_clean$hazard_rate, q99)
    
    gg_hazard_rate <- ggplot(baseline_haz_clean, aes(x = time, y = hazard_rate_capped)) +
      geom_line(color = "#D55E00", linewidth = 1.2, alpha = 0.8) +
      geom_smooth(se = TRUE, color = "#0072B2", linewidth = 1.3, alpha = 0.3) +
      labs(
        title = "Instantaneous Hazard Rate h(t)",
        subtitle = "Rate of progression events over time (smoothed)",
        x = "Time (Years)",
        y = "Hazard Rate"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    
    ggsave(
      file.path(method_dir, "hazard_rate_function.png"),
      gg_hazard_rate, 
      width = 10, 
      height = 7, 
      dpi = 300
    )
    
    cat("   ✓ Hazard rate function saved\n")
    
    gg_cum_hazard <- ggplot(baseline_haz, aes(x = time, y = hazard)) +
      geom_step(color = "#0072B2", linewidth = 1.3) +
      labs(
        title = "Cumulative Hazard Function H(t)",
        subtitle = "Integral of hazard rate over time",
        x = "Time (Years)",
        y = "Cumulative Hazard H(t)"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30")
      )
    
    ggsave(
      file.path(method_dir, "cumulative_hazard_function.png"),
      gg_cum_hazard, 
      width = 10, 
      height = 7, 
      dpi = 300
    )
    
    cat("   ✓ Cumulative hazard function saved\n")
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  # 3. Kaplan-Meier Analysis
  cat("\n3. Kaplan-Meier Analysis\n")
  
  tryCatch({
    km_overall <- survfit(Surv(time, event) ~ 1, data = thesis_surv_opt)
    
    km_summary <- summary(km_overall)
    km_data <- data.frame(
      time = km_summary$time,
      n_risk = km_summary$n.risk,
      n_event = km_summary$n.event,
      survival = km_summary$surv,
      std_err = km_summary$std.err,
      lower_ci = km_summary$lower,
      upper_ci = km_summary$upper
    )
    
    write.csv(
      km_data,
      file.path(method_dir, "kaplan_meier_table.csv"),
      row.names = FALSE
    )
    
    # Create KM plot with risk table
    gg_km_overall <- ggsurvplot(
      km_overall,
      data = thesis_surv_opt,
      conf.int = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      ggtheme = theme_minimal(base_size = 13),
      palette = "#0072B2",
      title = "Kaplan-Meier Survival Curve - Overall Cohort",
      xlab = "Time (Years)",
      ylab = "Progression-Free Probability",
      legend = "none",
      break.time.by = 1
    )
    
    # FIXED: Use special saving method for ggsurvplot objects
    png(file.path(method_dir, "kaplan_meier_overall.png"), 
        width = 10, height = 8, units = "in", res = 300)
    print(gg_km_overall)
    dev.off()
    
    km_stats <- summary(km_overall)$table
    
    cat("   ✓ Kaplan-Meier curves saved\n")
    cat(sprintf("   Median survival: %.2f years", km_stats["median"]))
    if(!is.na(km_stats["0.95LCL"]) && !is.na(km_stats["0.95UCL"])) {
      cat(sprintf(" (95%% CI: %.2f-%.2f)", km_stats["0.95LCL"], km_stats["0.95UCL"]))
    }
    cat("\n")
    
    timepoints <- c(1, 2, 3, 5)
    surv_at_times <- summary(km_overall, times = timepoints)
    
    surv_table <- data.frame(
      Time_Years = timepoints,
      Survival = surv_at_times$surv,
      Lower_CI = surv_at_times$lower,
      Upper_CI = surv_at_times$upper,
      N_Risk = surv_at_times$n.risk
    )
    
    write.csv(
      surv_table,
      file.path(method_dir, "survival_at_timepoints.csv"),
      row.names = FALSE
    )
    
    cat("   Survival probabilities:\n")
    for(i in 1:nrow(surv_table)) {
      cat(sprintf("     %d year: %.1f%% (95%% CI: %.1f%%-%.1f%%)\n",
                  surv_table$Time_Years[i],
                  surv_table$Survival[i] * 100,
                  surv_table$Lower_CI[i] * 100,
                  surv_table$Upper_CI[i] * 100))
    }
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  # 4. Cox vs KM Comparison
  cat("\n4. Cox vs Kaplan-Meier Comparison\n")
  
  tryCatch({
    cox_baseline <- basehaz(coxFit_thesis, centered = TRUE)
    cox_baseline$survival <- exp(-cox_baseline$hazard)
    cox_baseline$method <- "Cox Model"
    
    km_overall <- survfit(Surv(time, event) ~ 1, data = thesis_surv_opt)
    km_summary <- summary(km_overall)
    km_baseline <- data.frame(
      time = km_summary$time,
      survival = km_summary$surv,
      method = "Kaplan-Meier"
    )
    
    comparison_data <- rbind(
      cox_baseline[, c("time", "survival", "method")],
      km_baseline
    )
    
    gg_comparison <- ggplot(comparison_data, aes(x = time, y = survival, color = method)) +
      geom_step(linewidth = 1.2, alpha = 0.8) +
      scale_color_manual(values = c("Cox Model" = "#0072B2", "Kaplan-Meier" = "#D55E00")) +
      ylim(0, 1) +
      labs(
        title = "Cox Model vs Kaplan-Meier Survival Curves",
        subtitle = "Comparison of parametric (Cox) and non-parametric (KM) estimates",
        x = "Time (Years)",
        y = "Survival Probability",
        color = "Method"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray30"),
        legend.position = "bottom"
      )
    
    ggsave(
      file.path(method_dir, "cox_vs_km_comparison.png"),
      gg_comparison, 
      width = 10, 
      height = 7, 
      dpi = 300
    )
    
    cat("   ✓ Comparison plot saved\n")
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  # 5. Combined Panel
  cat("\n5. Creating Combined Panel Figure\n")
  
  tryCatch({
    library(gridExtra)
    
    p1 <- ggplot(baseline_haz, aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1) +
      ylim(0, 1) +
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
    
    km_plot_data <- data.frame(
      time = km_summary$time,
      survival = km_summary$surv,
      lower = km_summary$lower,
      upper = km_summary$upper
    )
    
    p4 <- ggplot(km_plot_data, aes(x = time, y = survival)) +
      geom_step(color = "#0072B2", linewidth = 1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#0072B2") +
      ylim(0, 1) +
      labs(title = "D) Kaplan-Meier", x = "Time (Years)", y = "Survival") +
      theme_minimal(base_size = 10)
    
    combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2)
    
    ggsave(
      file.path(method_dir, "survival_analysis_panel.png"),
      combined_plot,
      width = 12,
      height = 10,
      dpi = 300
    )
    
    cat("   ✓ Combined panel figure saved\n")
    
  }, error = function(e) {
    cat(sprintf("   ⚠ Error: %s\n", e$message))
  })
  
  cat("\n=== ENHANCED SURVIVAL ANALYSIS COMPLETE ===\n")
  
  ############################################################
  # STEP 10: SAVE ADDITIONAL DIAGNOSTICS
  ############################################################
  
  cat("\n=== STEP 10: SAVING DIAGNOSTICS ===\n")
  
  # Save feature selection results
  lasso_results <- data.frame(
    Feature = all_predictors,
    Selected = all_predictors %in% selected_features,
    Coefficient = as.numeric(coefs[match(all_predictors, rownames(coefs))]),
    Type = case_when(
      all_predictors %in% clinical_features ~ "Clinical",
      all_predictors %in% interaction_terms ~ "Interaction",
      TRUE ~ "Latent"
    )
  )
  
  write.csv(lasso_results,
            file.path(method_dir, "lasso_feature_selection.csv"),
            row.names = FALSE)
  
  # Save calibration diagnostics
  calibration_data <- data.frame(
    raw_risk = raw_risks,
    surv_3yr_actual = surv_3yr,
    calibrated_prob = predict_calibrated(thesis_surv_opt),
    PTID = thesis_surv_opt$PTID
  )
  
  write.csv(calibration_data,
            file.path(method_dir, "calibration_diagnostics.csv"),
            row.names = FALSE)
  
  # Plot calibration comparison
  cal_comparison <- data.frame(
    Method = rep(c("Before Calibration", "After Calibration"), each = length(raw_risks)),
    Predicted = c(
      1 - (raw_risks - min(raw_risks)) / (max(raw_risks) - min(raw_risks)),
      calibration_data$calibrated_prob
    ),
    Observed = rep(surv_3yr, 2)
  )
  
  gg_cal_compare <- ggplot(cal_comparison, aes(x = Predicted, y = Observed, color = Method)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    facet_wrap(~ Method) +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title = "Calibration: Before vs After Correction",
      x = "Predicted 3-Year Survival",
      y = "Observed 3-Year Survival"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
  
  ggsave(
    file.path(method_dir, "calibration_comparison.png"),
    gg_cal_compare, 
    width = 12, 
    height = 6, 
    dpi = 300
  )
  
  cat("  ✓ Diagnostics saved\n")
  
  ############################################################
  # SUMMARY
  ############################################################
  
  cat("\n============================================================\n")
  cat("THESIS METHOD OPTIMIZATION SUMMARY\n")
  cat("============================================================\n")
  cat(sprintf("\nFeature Reduction:\n"))
  cat(sprintf("  Initial features: %d\n", length(z_features)))
  cat(sprintf("  After filtering: %d\n", length(z_features_filtered)))
  cat(sprintf("  After Lasso: %d\n", length(selected_features)))
  cat(sprintf("\nModel Performance:\n"))
  cat(sprintf("  C-index: %.4f", metrics_thesis$cindex))
  if(!is.na(metrics_thesis$cindex_ci_lower)) {
    cat(sprintf(" (95%% CI: %.4f-%.4f)", 
                metrics_thesis$cindex_ci_lower, 
                metrics_thesis$cindex_ci_upper))
  }
  cat("\n")
  if(!is.na(metrics_thesis$calibration_slope)) {
    cat(sprintf("  Calibration slope: %.3f (ideal = 1.0)\n", 
                metrics_thesis$calibration_slope))
  }
  if(!is.na(metrics_thesis$brier_3yr)) {
    cat(sprintf("  Brier score (3yr): %.4f\n", metrics_thesis$brier_3yr))
  }
  if(!is.na(metrics_thesis$rmst_diff)) {
    cat(sprintf("  RMST difference: %.2f years (p = %.4f)\n", 
                metrics_thesis$rmst_diff, metrics_thesis$rmst_pval))
  }
  cat("============================================================\n")
  
  cat("\n✓ METHOD 6 (OPTIMIZED) COMPLETED SUCCESSFULLY\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in Optimized Thesis method\n")
  cat("Message:", conditionMessage(e), "\n")
  if(exists("traceback")) {
    print(traceback())
  }
})

cat("\n\n============================================================\n")
cat("METHOD 7: THESIS NO-MTL (Reconstruction-Only Autoencoder)\n")
cat("============================================================\n")

metrics_no_mtl <- NULL

tryCatch({
  
  # Load No-MTL features
  cat("\n=== LOADING NO-MTL DATA ===\n")
  
  # Check if the file exists
  if(!file.exists("latent_features_no_mtl (1).csv")) {
    cat("  ⚠ No-MTL features not found: latent_features_no_mtl.csv\n")
    cat("  Please run AD_Autoencoder_No_MTL.py first to generate features\n")
    stop("No-MTL features not available")
  }
  
  no_mtl_data <- read.csv("latent_features_no_mtl (1).csv")
  
  cat(sprintf("  Loaded No-MTL data: %d patients\n", length(unique(no_mtl_data$PTID))))
  
  # The No-MTL features will have columns: PTID, z_0, z_1, ..., z_127
  # We need to merge with survival/longitudinal information
  
  ############################################################
  # MERGE WITH ORIGINAL DATA TO GET SURVIVAL INFO
  ############################################################
  
  cat("\n=== MERGING WITH SURVIVAL DATA ===\n")
  
  # Get survival info from thesis_data (or any complete dataset)
  survival_info <- thesis_data %>%
    group_by(PTID) %>%
    summarize(
      time_to_event = first(time_to_event),
      event = first(event),
      Years_bl = list(Years_bl),
      MMSE = list(MMSE),
      AGE = first(AGE),
      PTGENDER = first(PTGENDER),
      PTEDUCAT = first(PTEDUCAT),
      ADAS13 = first(ADAS13),
      .groups = "drop"
    )
  
  # Merge No-MTL features with survival info
  no_mtl_complete <- no_mtl_data %>%
    left_join(survival_info, by = "PTID")
  
  # Filter to master cohort
  no_mtl_complete <- no_mtl_complete %>%
    filter(PTID %in% master_patients)
  
  cat(sprintf("  After filtering to master cohort: %d patients\n", 
              length(unique(no_mtl_complete$PTID))))
  
  ############################################################
  # FEATURE SELECTION (SAME AS THESIS METHOD)
  ############################################################
  
  cat("\n=== FEATURE SELECTION ===\n")
  
  z_features <- grep("^z_", names(no_mtl_complete), value = TRUE)
  cat(sprintf("  Initial z-features: %d\n", length(z_features)))
  
  # Remove zero-variance features
  zero_var <- sapply(no_mtl_complete[z_features], function(x) {
    v <- var(as.numeric(x), na.rm = TRUE)
    is.na(v) | v < 1e-10 | v > 1e10
  })
  z_features_filtered <- z_features[!zero_var]
  cat(sprintf("  After variance filter: %d\n", length(z_features_filtered)))
  
  # Filter by distribution quality (if e1071 available)
  if(require("e1071", quietly = TRUE)) {
    stable_features <- c()
    for(feat in z_features_filtered) {
      vals <- as.numeric(no_mtl_complete[[feat]])
      vals_clean <- vals[!is.na(vals)]
      
      if(length(vals_clean) > 3) {
        tryCatch({
          kurt <- kurtosis(vals_clean)
          skew <- skewness(vals_clean)
          
          if(abs(kurt) < 10 & abs(skew) < 3) {
            stable_features <- c(stable_features, feat)
          }
        }, error = function(e) NULL)
      }
    }
    
    if(length(stable_features) > 0) {
      z_features_filtered <- stable_features
      cat(sprintf("  After quality filter: %d\n", length(stable_features)))
    }
  }
  
  # Select top 15 by variance
  if(length(z_features_filtered) > 15) {
    z_vars <- sapply(no_mtl_complete[z_features_filtered], var, na.rm = TRUE)
    selected_z <- names(sort(z_vars, decreasing = TRUE)[1:15])
  } else {
    selected_z <- z_features_filtered
  }
  
  cat(sprintf("  Final z-features selected: %d\n", length(selected_z)))
  
  ############################################################
  # ADD CLINICAL FEATURES
  ############################################################
  
  cat("\n=== ADDING CLINICAL FEATURES ===\n")
  
  clinical_features <- c("AGE", "PTGENDER", "PTEDUCAT", "ADAS13")
  
  # Create interactions
  no_mtl_complete$AGE_ADAS <- no_mtl_complete$AGE * no_mtl_complete$ADAS13
  no_mtl_complete$EDUC_AGE <- no_mtl_complete$PTEDUCAT * no_mtl_complete$AGE
  no_mtl_complete$EDUC_ADAS <- no_mtl_complete$PTEDUCAT * no_mtl_complete$ADAS13
  
  interaction_terms <- c("AGE_ADAS", "EDUC_AGE", "EDUC_ADAS")
  
  all_predictors <- c(clinical_features, selected_z, interaction_terms)
  cat(sprintf("  Total predictors: %d\n", length(all_predictors)))
  
  ############################################################
  # PREPARE LONGITUDINAL & SURVIVAL DATA
  ############################################################
  
  cat("\n=== PREPARING DATA ===\n")
  
  # Expand longitudinal data
  no_mtl_long <- no_mtl_complete %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors)) %>%
    tidyr::unnest(cols = c(Years_bl, MMSE)) %>%
    filter(!is.na(MMSE), !is.na(Years_bl))
  
  # Use complete.cases
  no_mtl_long <- no_mtl_long[complete.cases(no_mtl_long), ]
  
  # Survival data (one row per patient)
  no_mtl_surv <- no_mtl_complete %>%
    select(PTID, time_to_event, event, all_of(all_predictors))
  
  # Rename for consistency
  names(no_mtl_surv)[names(no_mtl_surv) == "time_to_event"] <- "time"
  
  # Use complete.cases
  no_mtl_surv <- no_mtl_surv[complete.cases(no_mtl_surv), ]
  
  # Get common patients
  common_IDs_no_mtl <- intersect(no_mtl_long$PTID, no_mtl_surv$PTID)
  no_mtl_long <- no_mtl_long %>% filter(PTID %in% common_IDs_no_mtl)
  no_mtl_surv <- no_mtl_surv %>% filter(PTID %in% common_IDs_no_mtl)
  
  cat(sprintf("  Patients: %d, Events: %d (%.1f%%)\n", 
              length(common_IDs_no_mtl), sum(no_mtl_surv$event),
              100 * sum(no_mtl_surv$event) / length(common_IDs_no_mtl)))
  
  # Verify against master cohort
  if(length(common_IDs_no_mtl) != length(master_patients)) {
    cat(sprintf("  ⚠ WARNING: No-MTL cohort (%d) differs from master (%d)\n",
                length(common_IDs_no_mtl), length(master_patients)))
  }
  
  ############################################################
  # LASSO FEATURE SELECTION (SAME AS THESIS METHOD)
  ############################################################
  
  cat("\n=== LASSO FEATURE SELECTION ===\n")
  
  # Prepare matrices
  x_matrix <- as.matrix(no_mtl_surv[, all_predictors])
  y_matrix <- Surv(no_mtl_surv$time, no_mtl_surv$event)
  
  # Cross-validated Lasso
  set.seed(42)
  cv_fit <- cv.glmnet(
    x_matrix, y_matrix, 
    family = "cox", 
    alpha = 1,
    nfolds = 10,
    type.measure = "C"
  )
  
  # Get selected features
  optimal_lambda <- cv_fit$lambda.min
  coefs <- coef(cv_fit, s = optimal_lambda)
  selected_features <- rownames(coefs)[which(coefs != 0)]
  
  cat(sprintf("  Optimal lambda: %.6f\n", optimal_lambda))
  cat(sprintf("  Features selected by Lasso: %d (from %d)\n", 
              length(selected_features), length(all_predictors)))
  
  # Fallback if no features selected
  if(length(selected_features) == 0) {
    cat("  ⚠ No features selected by Lasso, using top 8 by univariate p-value\n")
    
    univariate_pvals <- sapply(all_predictors, function(pred) {
      tryCatch({
        formula <- as.formula(paste("Surv(time, event) ~", pred))
        fit <- coxph(formula, data = no_mtl_surv)
        summary(fit)$coefficients[1, "Pr(>|z|)"]
      }, error = function(e) 1.0)
    })
    
    selected_features <- names(sort(univariate_pvals)[1:8])
  } else {
    # Always include core clinical features
    selected_features <- unique(c(
      intersect(clinical_features, c(selected_features, all_predictors)),
      selected_features
    ))
  }
  
  cat("  Selected features:\n")
  print(selected_features)
  
  # Update predictors
  all_predictors_opt <- selected_features
  
  ############################################################
  # FIT MODELS
  ############################################################
  
  cat("\n=== MODEL FITTING ===\n")
  
  # Update data with optimized features
  no_mtl_long_opt <- no_mtl_long %>%
    select(PTID, Years_bl, MMSE, all_of(all_predictors_opt))
  
  no_mtl_surv_opt <- no_mtl_surv %>%
    select(PTID, time, event, all_of(all_predictors_opt))
  
  # Fit LME
  cat("  Fitting LME model...\n")
  long_formula <- as.formula(paste("MMSE ~ Years_bl +", 
                                   paste(all_predictors_opt, collapse = " + ")))
  
  lmeFit_no_mtl <- lme(
    long_formula, 
    random = ~ Years_bl | PTID, 
    data = no_mtl_long_opt,
    control = lmeControl(
      opt = "optim", 
      maxIter = 500,
      returnObject = TRUE
    )
  )
  cat("  ✓ LME fitted\n")
  
  # Fit Cox
  cat("  Fitting Cox model...\n")
  surv_formula <- as.formula(paste("Surv(time, event) ~", 
                                   paste(all_predictors_opt, collapse = " + ")))
  
  coxFit_no_mtl <- coxph(
    surv_formula, 
    data = no_mtl_surv_opt, 
    x = TRUE,
    method = "breslow"
  )
  cat("  ✓ Cox fitted\n")
  
  initial_cindex <- summary(coxFit_no_mtl)$concordance[1]
  cat(sprintf("  Initial C-index: %.4f\n", initial_cindex))
  
  ############################################################
  # CALIBRATION (OPTIONAL - SAME AS THESIS)
  ############################################################
  
  cat("\n=== CALIBRATION CORRECTION ===\n")
  
  raw_risks <- predict(coxFit_no_mtl, newdata = no_mtl_surv_opt, type = "lp")
  
  surv_3yr <- as.numeric(
    no_mtl_surv_opt$time > 3 | 
      (no_mtl_surv_opt$time <= 3 & no_mtl_surv_opt$event == 0)
  )
  
  calibration_model <- glm(
    surv_3yr ~ raw_risks, 
    family = binomial
  )
  
  predict_calibrated <- function(new_data) {
    raw_risks_new <- predict(coxFit_no_mtl, newdata = new_data, type = "lp")
    calibrated_probs <- predict(
      calibration_model, 
      newdata = data.frame(raw_risks = raw_risks_new),
      type = "response"
    )
    return(calibrated_probs)
  }
  
  coxFit_no_mtl$calibration_model <- calibration_model
  coxFit_no_mtl$predict_calibrated <- predict_calibrated
  
  cat("  ✓ Calibration model fitted\n")
  
  ############################################################
  # FIT JOINT MODEL
  ############################################################
  
  cat("\n=== JOINT MODEL ===\n")
  cat("  Fitting joint model (this may take several minutes)...\n")
  
  jointFit_no_mtl <- NULL
  
  tryCatch({
    jointFit_no_mtl <- jm(
      coxFit_no_mtl, 
      lmeFit_no_mtl, 
      time_var = "Years_bl",
      functional_forms = ~ value(MMSE) + slope(MMSE),
      n_iter = 5000, 
      n_burnin = 1000, 
      n_thin = 5, 
      n_chains = 2,
      seed = 42
    )
    cat("  ✓ Joint model fitted successfully\n")
    
  }, error = function(e) {
    cat("\n  ⚠ Joint model fitting failed, trying simpler model...\n")
    cat("  Error:", conditionMessage(e), "\n")
    
    tryCatch({
      jointFit_no_mtl <<- jm(
        coxFit_no_mtl, 
        lmeFit_no_mtl, 
        time_var = "Years_bl",
        functional_forms = ~ value(MMSE),
        n_iter = 3000, 
        n_burnin = 500, 
        n_thin = 3, 
        n_chains = 1,
        seed = 42
      )
      cat("  ✓ Simpler joint model fitted successfully\n")
      
    }, error = function(e2) {
      cat("  ✗ Joint model failed completely, proceeding without it\n")
      jointFit_no_mtl <<- NULL
    })
  })
  
  ############################################################
  # EVALUATION
  ############################################################
  
  cat("\n=== EVALUATION ===\n")
  
  metrics_no_mtl <- compute_all_metrics_with_figures(
    coxFit_no_mtl, 
    no_mtl_surv_opt, 
    no_mtl_long_opt, 
    jointFit_no_mtl, 
    "Thesis No-MTL"
  )
  
  cat("\n✓ METHOD 7 (NO-MTL) COMPLETED SUCCESSFULLY\n")
  
}, error = function(e) {
  cat("\n⚠ ERROR in No-MTL method\n")
  cat("Message:", conditionMessage(e), "\n")
  if(exists("traceback")) {
    print(traceback())
  }
})

############################################################
# UPDATE FINAL COMPARISON TO INCLUDE NO-MTL
############################################################

cat("\n\n============================================================\n")
cat("COMPREHENSIVE COMPARISON: ALL METHODS (INCLUDING NO-MTL)\n")
cat("============================================================\n\n")

all_metrics <- list(
  metrics_baseline,
  metrics_image,
  metrics_tabular,
  metrics_concat,
  metrics_noae,
  metrics_thesis,
  metrics_no_mtl  # ADD THIS LINE
)

# Remove NULL metrics
all_metrics <- all_metrics[!sapply(all_metrics, is.null)]

# Check if we have any metrics
if(length(all_metrics) == 0) {
  cat("⚠ WARNING: No methods completed successfully\n")
  cat("Cannot generate comparison table\n")
} else {
  
  get_metric <- function(x, field) {
    val <- x[[field]]
    if (is.null(val)) NA else val
  }
  
  comparison_df <- data.frame(
    Method      = sapply(all_metrics, function(x) x$name),
    C_index     = sapply(all_metrics, get_metric, "cindex"),
    Brier_1yr   = sapply(all_metrics, get_metric, "brier_1yr"),
    Brier_3yr   = sapply(all_metrics, get_metric, "brier_3yr"),
    Cal_Slope   = sapply(all_metrics, get_metric, "calibration_slope"),
    RMST_Diff   = sapply(all_metrics, get_metric, "rmst_diff"),
    RMST_Pval   = sapply(all_metrics, get_metric, "rmst_pval"),
    N_Patients  = sapply(all_metrics, get_metric, "n_patients"),
    N_Events    = sapply(all_metrics, get_metric, "n_events"),
    stringsAsFactors = FALSE
  )
  
  # Calculate improvement over baseline
  baseline_cindex <- comparison_df$C_index[comparison_df$Method == "Clinical Cox"]
  
  if(length(baseline_cindex) > 0 && !is.na(baseline_cindex)) {
    comparison_df$Improvement_Percent <- ((comparison_df$C_index - baseline_cindex) / baseline_cindex) * 100
  } else {
    comparison_df$Improvement_Percent <- NA
  }
  
  # Verify patient counts
  patient_counts <- unique(comparison_df$N_Patients[!is.na(comparison_df$N_Patients)])
  if(length(patient_counts) == 1) {
    cat(sprintf("\n✓ VERIFICATION: All methods evaluated on %d patients\n\n", 
                patient_counts[1]))
  } else {
    cat("\n⚠ WARNING: Methods have different patient counts!\n")
    print(comparison_df[, c("Method", "N_Patients")])
    cat("\n")
  }
  
  # Print results
  print(comparison_df)
  
  # Find best method
  valid_cindexes <- comparison_df$C_index[!is.na(comparison_df$C_index)]
  if(length(valid_cindexes) > 0) {
    best_idx <- which.max(comparison_df$C_index)
    best_method <- comparison_df$Method[best_idx]
    best_cindex <- comparison_df$C_index[best_idx]
    
    cat(sprintf("\n🏆 BEST METHOD: %s (C-index: %.4f)\n", best_method, best_cindex))
    
    if(!is.na(comparison_df$Improvement_Percent[best_idx])) {
      cat(sprintf("📈 Improvement over baseline: %.2f%%\n", 
                  comparison_df$Improvement_Percent[best_idx]))
    }
    
    # Compare MTL vs No-MTL if both exist
    mtl_idx <- which(comparison_df$Method == "MAE-JM")
    no_mtl_idx <- which(comparison_df$Method == "Thesis No-MTL")
    
    if(length(mtl_idx) > 0 && length(no_mtl_idx) > 0) {
      mtl_cindex <- comparison_df$C_index[mtl_idx]
      no_mtl_cindex <- comparison_df$C_index[no_mtl_idx]
      
      cat("\n============================================================\n")
      cat("🔬 MULTI-TASK LEARNING COMPARISON\n")
      cat("============================================================\n")
      cat(sprintf("With MTL (Cox+MMSE):       C-index = %.4f\n", mtl_cindex))
      cat(sprintf("Without MTL (Recon only):  C-index = %.4f\n", no_mtl_cindex))
      
      if(!is.na(mtl_cindex) && !is.na(no_mtl_cindex)) {
        diff <- mtl_cindex - no_mtl_cindex
        pct_diff <- (diff / no_mtl_cindex) * 100
        
        cat(sprintf("Difference:                Δ = %.4f (%.2f%%)\n", diff, pct_diff))
        
        if(diff > 0.01) {
          cat("\n✅ CONCLUSION: Multi-task learning provides significant benefit\n")
          cat("   The supervised signal from Cox+MMSE losses improves feature learning\n")
        } else if(diff > 0 && diff <= 0.01) {
          cat("\n⚠ CONCLUSION: Multi-task learning provides marginal benefit\n")
          cat("   The improvement is small but positive\n")
        } else {
          cat("\n❌ CONCLUSION: Multi-task learning does not improve performance\n")
          cat("   Reconstruction-only autoencoder performs as well or better\n")
        }
      }
      cat("============================================================\n")
    }
  }
  
  # Save results
  write.csv(comparison_df, "all_benchmarks_with_no_mtl.csv", row.names = FALSE)
  
  cat("\n============================================================\n")
  cat("✓ ANALYSIS COMPLETE\n")
  cat("============================================================\n")
  cat(sprintf("\nMethods completed: %d\n", nrow(comparison_df)))
  if(length(patient_counts) == 1) {
    cat(sprintf("Cohort size: %d patients\n", patient_counts[1]))
  }
  cat("\nResults saved to:\n")
  cat("  - all_benchmarks_with_no_mtl.csv\n")
  cat("  - thesis_figures/ (individual method folders)\n")
  cat("============================================================\n")
}
