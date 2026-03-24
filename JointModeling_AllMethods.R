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
cat("THESIS ANALYSIS v5\n")
cat("Improvements vs v4:\n")
cat("  A. MAE-JM JM uses single best Cox feature (not all 3) -> better dynamic AUC\n")
cat("  B. MAE-JM MCMC: 20000 iter / 5000 burnin for tighter posteriors\n")
cat("  C. tvAUC window extended to Dt=1.5 for more events per landmark\n")
cat("============================================================\n")

if (!dir.exists("thesis_figures")) dir.create("thesis_figures")

# ============================================================
# HELPERS
# ============================================================

safe_cv_glmnet <- function(x_mat, y_mat, alpha = 0.5) {
  n <- nrow(x_mat)
  if (n < 5) {
    cat(sprintf("  ⚠ Only %d patients — skipping glmnet\n", n))
    return(NULL)
  }
  nfolds <- min(10L, n)
  if (nfolds < 10)
    cat(sprintf("  ⚠ Small N (%d): nfolds=%d\n", n, nfolds))
  tryCatch(
    cv.glmnet(x_mat, y_mat, family = "cox", alpha = alpha,
              nfolds = nfolds, type.measure = "C"),
    error = function(e) {
      cat(sprintf("  ⚠ cv.glmnet failed (%s)\n", e$message)); NULL
    }
  )
}

cap_features_by_epe <- function(sel, events, epe = 10) {
  sel    <- as.character(sel)
  events <- as.integer(events)
  max_f  <- max(1L, floor(events / epe))
  if (length(sel) > max_f) {
    cat(sprintf("  EPE cap: %d events → max %d features (was %d)\n",
                events, max_f, length(sel)))
    sel <- sel[seq_len(max_f)]
  }
  sel
}

standardise_features <- function(train_df, test_df, feat_cols) {
  feat_cols <- intersect(feat_cols, intersect(names(train_df), names(test_df)))
  for (f in feat_cols) {
    train_df[[f]] <- as.numeric(train_df[[f]])
    test_df[[f]]  <- as.numeric(test_df[[f]])
    m <- mean(train_df[[f]], na.rm = TRUE)
    s <- sd(train_df[[f]],   na.rm = TRUE)
    if (!is.na(s) && s > 1e-10) {
      train_df[[f]] <- pmin(pmax((train_df[[f]] - m) / s, -5), 5)
      test_df[[f]]  <- pmin(pmax((test_df[[f]]  - m) / s, -5), 5)
    } else {
      train_df[[f]] <- 0.0
      test_df[[f]]  <- 0.0
    }
  }
  list(train = train_df, test = test_df)
}

select_features <- function(train_df, all_cols, n_events, method_name = "") {
  all_cols <- intersect(all_cols, names(train_df))
  ok_var   <- sapply(all_cols, function(f) {
    v <- var(as.numeric(train_df[[f]]), na.rm = TRUE)
    !is.na(v) && v > 1e-10
  })
  all_cols <- all_cols[ok_var]
  
  if (length(all_cols) == 0) {
    cat(sprintf("  ⚠ %s: no valid columns\n", method_name))
    return(character(0))
  }
  
  x_mat  <- as.matrix(train_df[, all_cols, drop = FALSE])
  y_mat  <- Surv(train_df$time, train_df$event)
  cv_fit <- safe_cv_glmnet(x_mat, y_mat)
  
  sel <- character(0)
  if (!is.null(cv_fit)) {
    coefs <- coef(cv_fit, s = cv_fit$lambda.1se)
    sel   <- all_cols[which(coefs != 0)]
    if (length(sel) < 5) {
      coefs <- coef(cv_fit, s = cv_fit$lambda.min)
      sel   <- all_cols[which(coefs != 0)]
    }
  }
  
  if (length(sel) < 5) {
    cat(sprintf("  ⚠ %s: glmnet gave <5 features — top-10 univariate\n", method_name))
    uni_c <- sapply(all_cols, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                    data = train_df))$concordance[1],
      error = function(e) 0.5
    ))
    top_idx <- order(uni_c, decreasing = TRUE)[seq_len(min(10L, length(all_cols)))]
    sel     <- all_cols[top_idx]
  }
  
  sel <- cap_features_by_epe(sel, n_events)
  cat(sprintf("  Selected %d features\n", length(sel)))
  
  if (length(sel) > 1) {
    test_cox <- tryCatch(
      coxph(as.formula(paste("Surv(time,event) ~", paste(sel, collapse=" + "))),
            data = train_df, x = FALSE,
            control = coxph.control(iter.max = 30)),
      error = function(e) NULL
    )
    if (!is.null(test_cox)) {
      cv <- tryCatch(coef(test_cox), error = function(e) NULL)
      if (!is.null(cv) && any(!is.finite(cv))) {
        cat(sprintf("  ⚠ Infinite coefficients with %d features — falling back to top-1\n",
                    length(sel)))
        uni_c2  <- sapply(sel, function(f) tryCatch(
          summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                        data = train_df))$concordance[1],
          error = function(e) 0.5
        ))
        sel <- sel[which.max(uni_c2)]
        cat(sprintf("  Fallback feature: %s\n", sel))
      }
    }
  }
  
  sel
}

fit_ridge_cox <- function(train_df, sel, method_name = "") {
  x_mat  <- as.matrix(train_df[, sel, drop = FALSE])
  y_mat  <- Surv(train_df$time, train_df$event)
  nfolds <- min(10L, nrow(x_mat))
  lambda_grid <- exp(seq(log(0.001), log(10), length.out = 100))
  
  cv_ridge <- tryCatch(
    cv.glmnet(x_mat, y_mat, family = "cox", alpha = 0,
              lambda = lambda_grid,
              nfolds = nfolds, type.measure = "C"),
    error = function(e) {
      cat(sprintf("  ⚠ %s: ridge cv failed — unpenalised fallback\n", method_name))
      NULL
    }
  )
  
  if (is.null(cv_ridge)) {
    return(coxph(as.formula(paste("Surv(time,event) ~", paste(sel, collapse=" + "))),
                 data = train_df, x = TRUE, method = "breslow"))
  }
  
  ridge_coefs <- as.numeric(coef(cv_ridge, s = cv_ridge$lambda.min))
  names(ridge_coefs) <- sel
  
  if (max(abs(ridge_coefs)) < 0.001) {
    cat(sprintf("  ⚠ %s: ridge coefs near-zero — using lambda at max CV C\n", method_name))
    best_lambda <- cv_ridge$lambda[which.max(cv_ridge$cvm)]
    ridge_coefs <- as.numeric(coef(cv_ridge, s = best_lambda))
    names(ridge_coefs) <- sel
    cat(sprintf("  Fallback lambda=%.5f\n", best_lambda))
  }
  
  cv_c <- max(cv_ridge$cvm, na.rm = TRUE)
  cat(sprintf("  Ridge CV C-index: %.4f  lambda=%.5f\n", cv_c, cv_ridge$lambda.min))
  cat(sprintf("  Coef range: [%.4f, %.4f]\n", min(ridge_coefs), max(ridge_coefs)))
  
  tryCatch(
    coxph(as.formula(paste("Surv(time,event) ~", paste(sel, collapse=" + "))),
          data = train_df, x = TRUE, method = "breslow",
          init = ridge_coefs, control = coxph.control(iter.max = 0)),
    error = function(e) {
      cat(sprintf("  ⚠ Ridge refit failed (%s) — unpenalised\n", e$message))
      coxph(as.formula(paste("Surv(time,event) ~", paste(sel, collapse=" + "))),
            data = train_df, x = TRUE, method = "breslow")
    }
  )
}

get_median_survival <- function(coxFit, newdata) {
  bh <- basehaz(coxFit, centered = FALSE)
  lp <- predict(coxFit, newdata = newdata, type = "lp")
  sapply(lp, function(l) {
    surv_i <- exp(-bh$hazard) ^ exp(l)
    idx    <- which(surv_i <= 0.5)
    if (length(idx) > 0) bh$time[min(idx)] else max(bh$time)
  })
}

# ============================================================
# FIX 3: JM association — read from $Survival (not $Outcome1)
# $Outcome1 contains LME fixed effects; $Survival contains the
# alpha (association) parameters we actually want.
# ============================================================

extract_jm_association <- function(jointFit, method_dir, name) {
  res <- list(
    jm_association_val         = NA, jm_association_val_lower   = NA,
    jm_association_val_upper   = NA, jm_association_slope       = NA,
    jm_association_slope_lower = NA, jm_association_slope_upper = NA
  )
  if (is.null(jointFit)) return(res)
  
  sm <- tryCatch(summary(jointFit), error = function(e) NULL)
  if (is.null(sm)) { cat("    ⚠ summary(jointFit) failed\n"); return(res) }
  
  # $Survival holds the association (alpha) parameters in JMbayes2
  tbl <- tryCatch(sm$Survival, error = function(e) NULL)
  if (is.null(tbl)) {
    cat("    ⚠ sm$Survival is NULL — names:", paste(names(sm), collapse=", "), "\n")
    return(res)
  }
  tbl <- as.data.frame(tbl)
  cn  <- colnames(tbl)
  mean_col  <- cn[grep("^[Mm]ean$|^[Ee]stimate$", cn)][1]
  lower_col <- cn[grep("2\\.5|[Ll]ower", cn)][1]
  upper_col <- cn[grep("97\\.5|[Uu]pper", cn)][1]
  rn <- rownames(tbl)
  
  cat("    Survival table rownames:", paste(rn, collapse=", "), "\n")
  
  vi <- grep("value|Assoct$", rn, ignore.case = TRUE)
  if (length(vi) > 0) {
    res$jm_association_val       <- as.numeric(tbl[vi[1], mean_col])
    res$jm_association_val_lower <- if (!is.na(lower_col)) as.numeric(tbl[vi[1], lower_col]) else NA
    res$jm_association_val_upper <- if (!is.na(upper_col)) as.numeric(tbl[vi[1], upper_col]) else NA
    cat(sprintf("    alpha(value): %.4f (%.4f, %.4f)\n",
                res$jm_association_val,
                res$jm_association_val_lower,
                res$jm_association_val_upper))
  } else {
    cat("    ⚠ No value(MMSE) row. Rownames:", paste(rn, collapse=", "), "\n")
  }
  
  si <- grep("slope|Assoct_slope", rn, ignore.case = TRUE)
  if (length(si) > 0) {
    res$jm_association_slope       <- as.numeric(tbl[si[1], mean_col])
    res$jm_association_slope_lower <- if (!is.na(lower_col)) as.numeric(tbl[si[1], lower_col]) else NA
    res$jm_association_slope_upper <- if (!is.na(upper_col)) as.numeric(tbl[si[1], upper_col]) else NA
    cat(sprintf("    alpha(slope): %.4f (%.4f, %.4f)\n",
                res$jm_association_slope,
                res$jm_association_slope_lower,
                res$jm_association_slope_upper))
  } else {
    cat("    ⚠ No slope(MMSE) row. Rownames:", paste(rn, collapse=", "), "\n")
  }
  
  write.csv(tbl,
            file.path(method_dir, paste0(gsub(" ", "_", name), "_jm_association.csv")),
            row.names = TRUE)
  res
}

# ============================================================
# FIX 1: tvAUC — pass real survival times, not landmark-censored
# JMbayes2 tvAUC needs actual time/event in newdata to identify
# who converts after Tstart. The previous mutate(time=t0, event=0)
# caused the "no data on subjects" error because tvAUC could not
# find anyone with an event time beyond Tstart.
# ============================================================

compute_dynamic_auc <- function(jointFit, long_data, surv_data,
                                patient_level_data, cox_vars,
                                method_dir, name) {
  res_auc <- NA
  if (is.null(jointFit) || is.null(long_data) || nrow(long_data) == 0) {
    cat("    ⚠ Skipping dynamic AUC — no joint model or no longitudinal data\n")
    return(res_auc)
  }
  
  long_data$PTID     <- as.character(long_data$PTID)
  long_data$Years_bl <- as.numeric(long_data$Years_bl)
  long_data$MMSE     <- as.numeric(long_data$MMSE)
  surv_data$PTID     <- as.character(surv_data$PTID)
  
  cat(sprintf("    Years_bl range: [%.3f, %.3f]\n",
              min(long_data$Years_bl, na.rm = TRUE),
              max(long_data$Years_bl, na.rm = TRUE)))
  
  max_obs   <- max(surv_data$time, na.rm = TRUE)
  landmarks <- c(0.5, 1, 1.5, 2)
  landmarks <- landmarks[landmarks < max_obs * 0.80]
  
  # Auto-detect Cox covariates from joint model
  cox_vars_detected <- tryCatch({
    f <- jointFit$model_info$terms_Surv_noResp
    if (!is.null(f)) all.vars(f) else character(0)
  }, error = function(e) character(0))
  cox_vars_detected <- setdiff(cox_vars_detected,
                               c("time", "event", "PTID", "Years_bl", "MMSE"))
  if (length(cox_vars_detected) == 0 && length(cox_vars) > 0) {
    cox_vars_detected <- cox_vars
    cat(sprintf("    Cox covariate hint: %s\n", paste(cox_vars_detected, collapse = ", ")))
  } else {
    cat(sprintf("    Cox covariates detected: %s\n",
                if (length(cox_vars_detected) == 0) "(none)"
                else paste(cox_vars_detected, collapse = ", ")))
  }
  
  surv_outcome <- surv_data %>% select(PTID, time, event)
  aucs <- c()
  
  for (t0 in landmarks) {
    tryCatch({
      at_risk_ids <- surv_outcome$PTID[surv_outcome$time > t0]
      if (length(at_risk_ids) < 5) {
        cat(sprintf("    ⚠ t=%.1f: only %d at risk\n", t0, length(at_risk_ids))); next
      }
      
      eps     <- 1e-6
      long_t0 <- long_data %>%
        filter(PTID %in% at_risk_ids, Years_bl <= t0 + eps)
      has_visits <- unique(long_t0$PTID)
      
      cat(sprintf("    t=%.1f: %d at-risk, %d have visits ≤ t0\n",
                  t0, length(at_risk_ids), length(has_visits)))
      
      if (length(has_visits) < 5) {
        cat(sprintf("    ⚠ t=%.1f: too few patients with visits\n", t0)); next
      }
      
      # KEY FIX: join REAL survival times — tvAUC needs actual time/event
      # (not landmark-censored values) to identify converters after t0.
      long_t0 <- long_t0 %>%
        left_join(surv_outcome, by = "PTID")
      
      # Join Cox submodel covariates if needed
      if (length(cox_vars_detected) > 0 && !is.null(patient_level_data)) {
        missing_cox <- setdiff(cox_vars_detected, names(long_t0))
        if (length(missing_cox) > 0) {
          avail <- intersect(c("PTID", missing_cox), names(patient_level_data))
          if (length(avail) > 1) {
            long_t0 <- long_t0 %>%
              left_join(
                patient_level_data %>%
                  select(all_of(avail)) %>%
                  filter(PTID %in% has_visits) %>%
                  distinct(PTID, .keep_all = TRUE),
                by = "PTID"
              )
          }
        }
      }
      
      still_missing <- setdiff(cox_vars_detected, names(long_t0))
      if (length(still_missing) > 0) {
        cat(sprintf("    ⚠ t=%.1f: missing %s\n", t0,
                    paste(still_missing, collapse = ", "))); next
      }
      
      long_t0 <- long_t0[complete.cases(
        long_t0[, intersect(c("PTID", "Years_bl", "MMSE", "time", "event"),
                            names(long_t0)), drop = FALSE]), ]
      
      n_pts_t0 <- length(unique(long_t0$PTID))
      n_evt_t0 <- sum(surv_outcome$event == 1 &
                        surv_outcome$time > t0 &
                        surv_outcome$time <= (t0 + 1.5) &
                        surv_outcome$PTID %in% has_visits)
      cat(sprintf("    t=%.1f: %d patients, %d events in (t0, t0+1.5]\n",
                  t0, n_pts_t0, n_evt_t0))
      
      if (n_pts_t0 < 5 || n_evt_t0 < 2) {
        cat(sprintf("    ⚠ t=%.1f: insufficient data\n", t0)); next
      }
      
      auc_obj <- tvAUC(jointFit, newdata = long_t0, Tstart = t0, Dt = 1.5)
      aucs <- c(aucs, auc_obj$auc)
      cat(sprintf("    Dynamic AUC at t=%.1f: %.4f\n", t0, auc_obj$auc))
      
    }, error = function(e)
      cat(sprintf("    ⚠ tvAUC at t=%.1f failed: %s\n", t0,
                  gsub("\n", " ", trimws(e$message)))))
  }
  
  if (length(aucs) > 0) {
    res_auc <- mean(aucs, na.rm = TRUE)
    cat(sprintf("    Mean dynamic AUC: %.4f\n", res_auc))
    write.csv(data.frame(Landmark = landmarks[seq_along(aucs)], AUC = aucs),
              file.path(method_dir, paste0(gsub(" ", "_", name), "_dynamic_auc.csv")),
              row.names = FALSE)
  } else {
    cat("    ⚠ No valid landmark AUCs computed\n")
  }
  res_auc
}

# ============================================================
# MAIN EVALUATION FUNCTION
# ============================================================

compute_all_metrics_with_figures <- function(coxFit, surv_data, long_data,
                                             jointFit = NULL, name = "Method",
                                             patient_level_data = NULL,
                                             cox_vars_hint = character(0)) {
  cat(sprintf("\n========================================\nEVALUATING: %s\n", name))
  
  n_pts    <- length(unique(surv_data$PTID))
  n_events <- sum(surv_data$event)
  cat(sprintf("  Patients: %d | Events: %d (%.1f%%)\n",
              n_pts, n_events, 100 * n_events / n_pts))
  
  method_dir <- file.path("thesis_figures", gsub("[^A-Za-z0-9_-]", "_", name))
  if (!dir.exists(method_dir)) dir.create(method_dir, recursive = TRUE)
  
  res <- list(
    name = name, n_patients = n_pts, n_events = n_events,
    event_rate = n_events / n_pts,
    cindex = NA, cindex_ci_lower = NA, cindex_ci_upper = NA,
    brier_1yr = NA, brier_3yr = NA, brier_5yr = NA,
    calibration_slope = NA, calibration_intercept = NA,
    rmst_tau = NA, rmst_low = NA, rmst_low_lower = NA, rmst_low_upper = NA,
    rmst_high = NA, rmst_high_lower = NA, rmst_high_upper = NA,
    rmst_diff = NA, rmst_diff_lower = NA, rmst_diff_upper = NA, rmst_pval = NA,
    mae_time = NA, rmse_time = NA, corr_time = NA, median_error = NA,
    n_significant = NA, km_high_median = NA, km_low_median = NA,
    jm_association_val = NA, jm_association_val_lower = NA,
    jm_association_val_upper = NA, jm_association_slope = NA,
    jm_association_slope_lower = NA, jm_association_slope_upper = NA,
    jm_dynamic_auc = NA
  )
  
  # 1. C-index
  cat("\n1. DISCRIMINATION\n")
  res$cindex <- summary(coxFit)$concordance[1]
  cat(sprintf("   C-index: %.4f\n", res$cindex))
  
  coef_vals <- tryCatch(coef(coxFit), error = function(e) NULL)
  if (!is.null(coef_vals) && any(!is.finite(coef_vals))) {
    cat("   ⚠ WARNING: Infinite/NA coefficients — complete separation.\n")
    cat("     C-index is UNRELIABLE.\n")
  }
  
  # 2. Cox coefficients + forest plot
  cat("\n2. COX COEFFICIENTS\n")
  tryCatch({
    cf <- as.data.frame(summary(coxFit)$coefficients)
    cf$Variable  <- rownames(cf)
    cf$HR        <- exp(cf$coef)
    cf$HR_lower  <- exp(cf$coef - 1.96 * cf$`se(coef)`)
    cf$HR_upper  <- exp(cf$coef + 1.96 * cf$`se(coef)`)
    print(head(cf[order(cf$`Pr(>|z|)`), c("Variable","coef","HR","Pr(>|z|)")], 10))
    write.csv(cf, file.path(method_dir, paste0(gsub(" ","_",name),"_cox_coefficients.csv")),
              row.names = FALSE)
    res$n_significant <- sum(cf$`Pr(>|z|)` < 0.10)
    sig <- cf[cf$`Pr(>|z|)` < 0.10, ]
    if (nrow(sig) > 0 && nrow(sig) <= 25) {
      sig <- sig[order(sig$HR, decreasing = TRUE), ]
      sig$HR_plot       <- pmin(pmax(sig$HR,       0.001), 1000)
      sig$HR_lower_plot <- pmin(pmax(sig$HR_lower, 0.001), 1000)
      sig$HR_upper_plot <- pmin(pmax(sig$HR_upper, 0.001), 1000)
      fp <- ggplot(sig, aes(x = HR_plot, y = reorder(Variable, HR_plot))) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
        geom_errorbarh(aes(xmin = HR_lower_plot, xmax = HR_upper_plot),
                       height = 0.2, linewidth = 0.6) +
        geom_point(size = 3, color = "#0072B2") + scale_x_log10() +
        labs(title = paste("Hazard Ratios:", name), x = "HR (log scale)", y = "") +
        theme_minimal(base_size = 13)
      ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_forest_plot.png")),
             fp, width = 10, height = max(5, nrow(sig) * 0.35), dpi = 300)
      cat(sprintf("   ✓ Forest plot (%d features at p<0.10)\n", nrow(sig)))
    }
  }, error = function(e) cat(sprintf("   ⚠ Coef error: %s\n", e$message)))
  
  # 3. Baseline survival
  cat("\n3. BASELINE SURVIVAL\n")
  tryCatch({
    bsf <- basehaz(coxFit, centered = TRUE)
    bsf$surv_prob <- exp(-bsf$hazard)
    ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_baseline_hazard.png")),
           ggplot(bsf, aes(time, hazard)) + geom_step(color="#0072B2", linewidth=1.2) +
             labs(title=paste("Baseline Hazard:",name), x="Time (Years)", y="Cum. Hazard") +
             theme_minimal(base_size=13),
           width=8, height=6, dpi=300)
    ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_baseline_survival.png")),
           ggplot(bsf, aes(time, surv_prob)) + geom_step(color="#0072B2", linewidth=1.2) +
             ylim(0,1) + labs(title=paste("Baseline Survival:",name),
                              x="Time (Years)", y="Survival Prob.") +
             theme_minimal(base_size=13),
           width=8, height=6, dpi=300)
    cat("   ✓ Saved\n")
  }, error = function(e) cat(sprintf("   ⚠ %s\n", e$message)))
  
  # 4. Brier score
  cat("\n4. BRIER SCORE\n")
  tryCatch({
    if (n_pts < 10) stop("too few patients")
    spec <- surv_data
    if (!"time" %in% names(spec) && "time_to_event" %in% names(spec))
      spec$time <- spec$time_to_event
    pt  <- c(1, 3, 5)
    pec_obj <- pec(list("Model" = coxFit), formula = Surv(time, event) ~ 1,
                   data = spec, times = pt, exact = FALSE,
                   cens.model = "marginal", verbose = FALSE)
    bs <- pec_obj$AppErr$Model
    res$brier_1yr <- bs[1]; res$brier_3yr <- bs[2]; res$brier_5yr <- bs[3]
    cat(sprintf("   1yr: %.4f | 3yr: %.4f | 5yr: %.4f\n", bs[1], bs[2], bs[3]))
    ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_brier_scores.png")),
           ggplot(data.frame(Time=pt, Brier=bs), aes(Time, Brier)) +
             geom_line(color="#0072B2",linewidth=1.3) + geom_point(size=4,color="#0072B2") +
             geom_hline(yintercept=0.10, linetype="dashed", color="darkgreen") +
             labs(title=paste("Brier Score:",name), x="Time (Years)", y="Brier Score") +
             theme_minimal(base_size=13),
           width=8, height=6, dpi=300)
    cat("   ✓ Saved\n")
  }, error = function(e) {
    if (e$message != "too few patients") cat(sprintf("   ⚠ %s\n", e$message))
    else cat("   ⚠ Skipped (<10 patients)\n")
  })
  
  # 5. KM + RMST
  cat("\n5. KM STRATIFICATION\n")
  tryCatch({
    surv_data$risk       <- predict(coxFit, newdata = surv_data, type = "lp")
    surv_data$risk_group <- ifelse(surv_data$risk >= median(surv_data$risk),
                                   "High Risk", "Low Risk")
    fs  <- survfit(Surv(time, event) ~ risk_group, data = surv_data)
    kmp <- ggsurvplot(fs, data = surv_data, risk.table = TRUE, pval = TRUE,
                      conf.int = TRUE, palette = c("#D55E00","#0072B2"),
                      linetype = c("solid","dashed"), legend.title = "Risk Group",
                      title = paste("Kaplan-Meier:", name),
                      xlab = "Time (Years)", ylab = "Progression-Free Probability")
    ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_KM_stratified.png")),
           kmp$plot, width=10, height=8, dpi=300)
    kt <- summary(fs)$table
    res$km_high_median <- kt["risk_group=High Risk", "median"]
    res$km_low_median  <- kt["risk_group=Low Risk",  "median"]
    cat(sprintf("   High: %.2f yr | Low: %.2f yr\n",
                res$km_high_median, res$km_low_median))
    cat("   ✓ KM saved\n")
    
    cat("\n6. RMST\n")
    tryCatch({
      surv_data$arm <- as.numeric(surv_data$risk_group == "Low Risk")
      tau <- min(
        max(surv_data$time[surv_data$risk_group == "High Risk"], na.rm = TRUE),
        max(surv_data$time[surv_data$risk_group == "Low Risk"],  na.rm = TRUE)
      ) * 0.9
      cat(sprintf("   Tau: %.2f yr\n", tau))
      rr <- rmst2(surv_data$time, surv_data$event, surv_data$arm, tau = tau)
      res$rmst_tau        <- tau
      res$rmst_low        <- rr$RMST.arm1$rmst["Est."]
      res$rmst_low_lower  <- rr$RMST.arm1$rmst["lower .95"]
      res$rmst_low_upper  <- rr$RMST.arm1$rmst["upper .95"]
      res$rmst_high       <- rr$RMST.arm0$rmst["Est."]
      res$rmst_high_lower <- rr$RMST.arm0$rmst["lower .95"]
      res$rmst_high_upper <- rr$RMST.arm0$rmst["upper .95"]
      res$rmst_diff       <- rr$unadjusted.result[1, "Est."]
      res$rmst_diff_lower <- rr$unadjusted.result[1, "lower .95"]
      res$rmst_diff_upper <- rr$unadjusted.result[1, "upper .95"]
      res$rmst_pval       <- rr$unadjusted.result[1, "p"]
      cat(sprintf("   Diff: %.3f yr (p=%.4f)\n", res$rmst_diff, res$rmst_pval))
      write.csv(data.frame(
        Risk_Group = c("Low Risk","High Risk","Difference"),
        RMST_Years = c(res$rmst_low, res$rmst_high, res$rmst_diff),
        Lower_CI   = c(res$rmst_low_lower, res$rmst_high_lower, res$rmst_diff_lower),
        Upper_CI   = c(res$rmst_low_upper, res$rmst_high_upper, res$rmst_diff_upper),
        P_Value    = c(NA, NA, res$rmst_pval), Tau = tau),
        file.path(method_dir, paste0(gsub(" ","_",name),"_rmst_analysis.csv")),
        row.names = FALSE)
      cat("   ✓ RMST saved\n")
    }, error = function(e) cat(sprintf("   ⚠ RMST: %s\n", e$message)))
  }, error = function(e) cat(sprintf("   ⚠ KM: %s\n", e$message)))
  
  # 7. Predicted vs observed
  cat("\n7. PREDICTED vs OBSERVED\n")
  tryCatch({
    conv <- surv_data[surv_data$event == 1, ]
    if (nrow(conv) < 5) stop("too few converters")
    pred <- get_median_survival(coxFit, conv)
    obs  <- conv$time
    err  <- abs(pred - obs)
    res$mae_time     <- mean(err)
    res$rmse_time    <- sqrt(mean((pred - obs)^2))
    res$corr_time    <- cor(pred, obs, method = "pearson")
    res$median_error <- median(err)
    cat(sprintf("   MAE: %.3f yr | RMSE: %.3f yr | r: %.3f\n",
                res$mae_time, res$rmse_time, res$corr_time))
    pod <- data.frame(Predicted=pred, Observed=obs, PTID=conv$PTID, Error=err)
    ggsave(file.path(method_dir, paste0(gsub(" ","_",name),"_predicted_vs_observed.png")),
           ggplot(pod, aes(Observed, Predicted)) +
             geom_abline(intercept=0, slope=1, linetype="dashed", color="gray40") +
             geom_point(alpha=0.7, size=3.5, color="#0072B2") +
             geom_smooth(method="lm", se=TRUE, color="#D55E00", alpha=0.2) +
             labs(title=paste("Predicted vs Observed:",name),
                  x="Observed (Years)", y="Predicted (Years)") +
             theme_minimal(base_size=13) + coord_equal(),
           width=8, height=8, dpi=300)
    write.csv(pod, file.path(method_dir, paste0(gsub(" ","_",name),"_time_predictions.csv")),
              row.names=FALSE)
    cat("   ✓ Saved\n")
  }, error = function(e) cat(sprintf("   ⚠ %s\n", e$message)))
  
  # 8. Calibration
  cat("\n8. CALIBRATION\n")
  tryCatch({
    bh  <- basehaz(coxFit, centered = FALSE)
    lp  <- predict(coxFit, newdata = surv_data, type = "lp")
    t3  <- which(bh$time <= 3)
    if (length(t3) == 0) t3 <- 1L
    H0  <- bh$hazard[max(t3)]
    ps  <- exp(-H0 * exp(lp))
    if (diff(range(ps, na.rm = TRUE)) < 0.05) stop("range too narrow")
    nu  <- length(unique(round(ps, 4)))
    ng  <- if (nu >= 50) 10L else if (nu >= 30) 7L else if (nu >= 15) 5L else
      if (nu >= 6) 3L else 0L
    if (ng < 3) stop("too homogeneous")
    brk <- unique(quantile(ps, seq(0, 1, length.out = ng + 1), na.rm = TRUE))
    if (length(brk) < 3) stop("too few breaks")
    dec <- cut(ps, breaks = brk, include.lowest = TRUE, labels = FALSE)
    cal_rows <- lapply(sort(unique(na.omit(dec))), function(g) {
      idx <- which(dec == g)
      if (length(idx) < 5) return(NULL)
      km  <- summary(survfit(Surv(time,event)~1, data=surv_data[idx,]),
                     times=3, extend=TRUE)
      if (length(km$surv) == 0 || is.na(km$surv)) return(NULL)
      data.frame(predicted=mean(ps[idx]), observed=km$surv, n=length(idx))
    })
    cal_rows <- Filter(Negate(is.null), cal_rows)
    if (length(cal_rows) < 3) stop("too few groups")
    cd  <- do.call(rbind, cal_rows)
    cm  <- lm(observed ~ predicted, data = cd)
    res$calibration_slope     <- coef(cm)[2]
    res$calibration_intercept <- coef(cm)[1]
    cat(sprintf("   Groups: %d | Slope: %.3f | Intercept: %.3f\n",
                nrow(cd), res$calibration_slope, res$calibration_intercept))
    write.csv(cd, file.path(method_dir, paste0(gsub("[^A-Za-z0-9_-]","_",name),
                                               "_calibration_data.csv")), row.names=FALSE)
    ggsave(file.path(method_dir, paste0(gsub("[^A-Za-z0-9_-]","_",name),
                                        "_calibration_plot.png")),
           ggplot(cd, aes(predicted, observed)) +
             geom_abline(intercept=0,slope=1,linetype="dashed",color="gray40",linewidth=0.8) +
             geom_smooth(method="lm",se=TRUE,color="#D55E00",fill="#D55E00",alpha=0.15) +
             geom_point(aes(size=n),color="#0072B2",alpha=0.85) +
             scale_size_continuous(name="Group n",range=c(3,8)) +
             xlim(0,1)+ylim(0,1)+coord_equal() +
             labs(title=paste("Calibration (3-year):",name),
                  subtitle=sprintf("Slope=%.3f | Intercept=%.3f | Groups=%d",
                                   res$calibration_slope,res$calibration_intercept,nrow(cd)),
                  x="Predicted 3-yr Survival",y="Observed 3-yr Survival") +
             theme_minimal(base_size=13),
           width=8, height=8, dpi=300)
    cat("   ✓ Saved\n")
  }, error = function(e) {
    known <- c("range too narrow","too homogeneous","too few breaks","too few groups")
    if (!e$message %in% known) cat(sprintf("   ⚠ Calibration: %s\n", e$message))
    else cat(sprintf("   ⚠ Skipped (%s)\n", e$message))
  })
  
  # 9. Bootstrap CI
  # FIX 2: filter out C=1.0 replicates (complete separation in resample);
  # fall back from "perc" to "norm" if too few valid replicates remain.
  cat("\n9. BOOTSTRAP CI\n")
  tryCatch({
    if (n_pts < 20) stop("too few patients")
    # Safely recover predictors — ridge Cox with iter.max=0 may not store $formula
    boot_pred_cols <- names(coef(coxFit))
    boot_formula   <- as.formula(
      paste("Surv(time, event) ~", paste(boot_pred_cols, collapse = " + "))
    )
    boot_fn <- function(data, idx) {
      bd  <- data[idx, ]
      bfx <- tryCatch(coxph(boot_formula, data=bd, x=TRUE,
                            control=coxph.control(iter.max=100)),
                      error=function(e) NULL)
      if (is.null(bfx)) return(NA_real_)
      cv <- tryCatch(coef(bfx), error = function(e) NULL)
      # Return NA if any coefficient is non-finite (separation in resample)
      if (!is.null(cv) && any(!is.finite(cv))) return(NA_real_)
      summary(bfx)$concordance[1]
    }
    cat("   Running bootstrap R=500...\n")
    br   <- boot(data=surv_data, statistic=boot_fn, R=500)
    
    # Filter replicates: remove NA and values == 1.0 (degenerate)
    vld  <- br$t[!is.na(br$t) & br$t < 1.0]
    pct_valid <- 100 * length(vld) / 500
    cat(sprintf("   Valid replicates: %d/500 (%.0f%%)\n", length(vld), pct_valid))
    
    if (length(vld) >= 50) {
      # Percentile CI on filtered replicates
      ci_lo <- quantile(vld, 0.025)
      ci_hi <- quantile(vld, 0.975)
      res$cindex_ci_lower <- ci_lo
      res$cindex_ci_upper <- ci_hi
      cat(sprintf("   C-index: %.4f (%.4f, %.4f) [percentile, filtered]\n",
                  res$cindex, ci_lo, ci_hi))
    } else if (length(vld) >= 10) {
      # Too few for percentile — use normal approximation on filtered set
      se_boot <- sd(vld)
      res$cindex_ci_lower <- res$cindex - 1.96 * se_boot
      res$cindex_ci_upper <- res$cindex + 1.96 * se_boot
      cat(sprintf("   C-index: %.4f (%.4f, %.4f) [normal approx, %d valid replicates]\n",
                  res$cindex, res$cindex_ci_lower, res$cindex_ci_upper, length(vld)))
      cat("   ⚠ Few valid replicates — interpret CI cautiously\n")
    } else {
      cat(sprintf("   ⚠ Only %d valid replicates — CI not reported\n", length(vld)))
    }
    saveRDS(br, file.path(method_dir, paste0(gsub(" ","_",name),"_bootstrap.rds")))
    cat("   ✓ Done\n")
  }, error = function(e) {
    if (e$message != "too few patients") cat(sprintf("   ⚠ Bootstrap: %s\n", e$message))
    else cat("   ⚠ Skipped (<20 patients)\n")
  })
  
  # 10. Joint model: association parameters + dynamic AUC
  cat("\n10. JOINT MODEL METRICS\n")
  if (!is.null(jointFit)) {
    jm_assoc <- extract_jm_association(jointFit, method_dir, name)
    res$jm_association_val         <- jm_assoc$jm_association_val
    res$jm_association_val_lower   <- jm_assoc$jm_association_val_lower
    res$jm_association_val_upper   <- jm_assoc$jm_association_val_upper
    res$jm_association_slope       <- jm_assoc$jm_association_slope
    res$jm_association_slope_lower <- jm_assoc$jm_association_slope_lower
    res$jm_association_slope_upper <- jm_assoc$jm_association_slope_upper
    
    res$jm_dynamic_auc <- compute_dynamic_auc(
      jointFit, long_data, surv_data,
      patient_level_data, cox_vars_hint,
      method_dir, name
    )
  } else {
    cat("    ⚠ No joint model — skipping\n")
  }
  
  # Append to master results file
  rdf <- as.data.frame(res, stringsAsFactors = FALSE)
  mf  <- "thesis_figures/ALL_METHODS_RESULTS.csv"
  if (file.exists(mf)) {
    ex  <- read.csv(mf, stringsAsFactors=FALSE)
    ac  <- union(names(ex), names(rdf))
    for (col in setdiff(ac, names(ex)))  ex[[col]]  <- NA
    for (col in setdiff(ac, names(rdf))) rdf[[col]] <- NA
    rdf <- rbind(ex[,ac], rdf[,ac])
  }
  write.csv(rdf, mf, row.names=FALSE)
  
  cat(sprintf("\n========================================\nSUMMARY: %s\n", name))
  cat(sprintf("  C-index: %.4f", res$cindex))
  if (!is.na(res$cindex_ci_lower))
    cat(sprintf(" (%.4f, %.4f)", res$cindex_ci_lower, res$cindex_ci_upper))
  cat("\n")
  if (!is.na(res$brier_3yr))         cat(sprintf("  Brier 3yr:   %.4f\n", res$brier_3yr))
  if (!is.na(res$calibration_slope)) cat(sprintf("  Cal slope:   %.3f\n",  res$calibration_slope))
  if (!is.na(res$rmst_diff))         cat(sprintf("  RMST diff:   %.3f yr (p=%.4f)\n",
                                                 res$rmst_diff, res$rmst_pval))
  if (!is.na(res$mae_time))          cat(sprintf("  MAE:         %.3f yr\n", res$mae_time))
  cat("========================================\n")
  return(res)
}

# ============================================================
# HELPER: build patient-level survival frame from a
# longitudinal CSV that already has a split column.
# ============================================================

build_split_surv <- function(pl_csv, z_final_pat = "^z_final_",
                             z_slope_pat = "^z_slope_",
                             extra_surv_cols = character(0),
                             method_name = "") {
  pl <- read.csv(pl_csv, stringsAsFactors = FALSE)
  pl$PTID <- as.character(pl$PTID)
  if (!"split" %in% names(pl))
    stop(sprintf("No 'split' column in %s — re-export from Python.", pl_csv))
  
  cat(sprintf("  %s — loaded %d patients\n", method_name, nrow(pl)))
  cat("  Split:\n"); print(table(pl$split, useNA="always"))
  
  zf <- grep(z_final_pat, names(pl), value = TRUE)
  zs <- grep(z_slope_pat, names(pl), value = TRUE)
  all_cols <- c(zf, zs, extra_surv_cols)
  all_cols <- intersect(all_cols, names(pl))
  
  surv_full <- pl %>%
    select(PTID, split, time_to_event, event, all_of(all_cols)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  surv_full <- surv_full[complete.cases(surv_full), ]
  
  train_raw <- surv_full %>% filter(split == "train")
  test_raw  <- surv_full %>% filter(split == "val")
  
  std <- standardise_features(train_raw, test_raw, all_cols)
  list(train = std$train, test = std$test,
       all_cols = all_cols, zf = zf, zs = zs, pl = pl)
}

# ============================================================
# METHOD 1: CLINICAL COX (BASELINE)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 1: CLINICAL COX (BASELINE)\n")
cat("============================================================\n")

metrics_baseline <- NULL
tryCatch({
  dat <- read.csv("baseline_clinical_features.csv", stringsAsFactors = FALSE)
  dat$PTID <- as.character(dat$PTID)
  if (!"split" %in% names(dat)) stop("No 'split' column in baseline_clinical_features.csv")
  cat("  Split:\n"); print(table(dat$split, useNA="always"))
  
  long_full <- dat %>%
    filter(!is.na(MMSE), !is.na(Years_bl)) %>%
    select(PTID, split, Years_bl, MMSE, ADAS13, AGE, PTGENDER, PTEDUCAT) %>%
    arrange(PTID, Years_bl)
  long_full <- long_full[complete.cases(long_full), ]
  
  surv_full <- dat %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              AGE = first(AGE), PTGENDER = first(PTGENDER),
              PTEDUCAT = first(PTEDUCAT), ADAS13 = first(ADAS13),
              split = first(split), .groups = "drop")
  surv_full <- surv_full[complete.cases(surv_full), ]
  
  surv_train <- surv_full %>% filter(split == "train")
  surv_test  <- surv_full %>% filter(split == "val")
  long_train <- long_full %>% filter(PTID %in% surv_train$PTID)
  long_test  <- long_full %>% filter(PTID %in% surv_test$PTID)
  
  single_v <- long_train %>% group_by(PTID) %>% filter(n() < 2) %>% pull(PTID) %>% unique()
  long_train  <- long_train  %>% filter(!PTID %in% single_v)
  surv_train  <- surv_train  %>% filter(!PTID %in% single_v)
  common_tr   <- intersect(surv_train$PTID, long_train$PTID)
  surv_train  <- surv_train %>% filter(PTID %in% common_tr)
  long_train  <- long_train %>% filter(PTID %in% common_tr)
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(surv_train), sum(surv_train$event)))
  cat(sprintf("  Test:  %d patients, %d events\n", nrow(surv_test),  sum(surv_test$event)))
  
  lmeFit_baseline <- lme(
    MMSE ~ Years_bl + AGE + PTGENDER + PTEDUCAT + ADAS13,
    random  = ~ Years_bl | PTID,
    data    = long_train,
    control = lmeControl(opt = "optim", maxIter = 200)
  )
  coxFit_baseline <- coxph(
    Surv(time, event) ~ AGE + PTGENDER + PTEDUCAT + ADAS13,
    data = surv_train, x = TRUE, ties = "breslow"
  )
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_baseline)$concordance[1]))
  
  jointFit_baseline <- tryCatch(
    jm(coxFit_baseline, lmeFit_baseline, time_var = "Years_bl",
       functional_forms = ~value(MMSE) + slope(MMSE),
       n_iter = 5000, n_burnin = 1000, n_thin = 5, n_chains = 2),
    error = function(e) { cat(sprintf("  ⚠ JM failed: %s\n", e$message)); NULL }
  )
  
  metrics_baseline <- compute_all_metrics_with_figures(
    coxFit_baseline, surv_test, long_test, jointFit_baseline,
    "Clinical Cox", patient_level_data = surv_test
  )
  cat("\n✓ METHOD 1 COMPLETED\n")
}, error = function(e) cat(sprintf("\n⚠ ERROR in Clinical Cox: %s\n", e$message)))

# ============================================================
# METHOD 2: IMAGE-ONLY CNN
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 2: IMAGE-ONLY CNN\n")
cat("============================================================\n")

metrics_image <- NULL
tryCatch({
  dat <- read.csv("image_only_features.csv", stringsAsFactors = FALSE)
  dat$PTID <- as.character(dat$PTID)
  if (!"split" %in% names(dat)) stop("No 'split' column in image_only_features.csv")
  
  img_features <- grep("^img_feat_|^z_", names(dat), value = TRUE)
  if (length(img_features) == 0) stop("No image features found")
  
  surv_full <- dat %>%
    group_by(PTID) %>%
    summarize(time = first(time_to_event), event = first(event),
              split = first(split),
              across(all_of(img_features), first), .groups = "drop")
  surv_full <- surv_full[complete.cases(surv_full), ]
  
  train_raw <- surv_full %>% filter(split == "train")
  test_raw  <- surv_full %>% filter(split == "val")
  std       <- standardise_features(train_raw, test_raw, img_features)
  surv_train <- std$train
  surv_test  <- std$test
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(surv_train), sum(surv_train$event)))
  cat(sprintf("  Test:  %d patients, %d events\n", nrow(surv_test),  sum(surv_test$event)))
  
  if (length(img_features) > 10) {
    img_vars     <- sapply(surv_train[img_features], var, na.rm = TRUE)
    top_idx      <- order(img_vars, decreasing = TRUE)[1:10]
    img_features <- img_features[top_idx]
  }
  img_features <- cap_features_by_epe(img_features, sum(surv_train$event))
  
  coxFit_image <- fit_ridge_cox(surv_train, img_features, "Image-Only")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_image)$concordance[1]))
  
  surv_test_final <- surv_test %>% select(PTID, time, event, all_of(img_features))
  metrics_image <- compute_all_metrics_with_figures(
    coxFit_image, surv_test_final, NULL, NULL, "Image-Only CNN"
  )
  cat("\n✓ METHOD 2 COMPLETED\n")
}, error = function(e) cat(sprintf("\n⚠ ERROR in Image-Only: %s\n", e$message)))

# ============================================================
# METHOD 3: TABULAR-ONLY DEEP NN
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 3: TABULAR-ONLY DEEP NN\n")
cat("============================================================\n")

metrics_tabular <- NULL
tryCatch({
  if (!file.exists("tabular_patient_level.csv")) stop("tabular_patient_level.csv not found")
  
  out <- build_split_surv("tabular_patient_level.csv", method_name = "Tabular")
  tab_train <- out$train
  tab_test  <- out$test
  all_cols  <- out$all_cols
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(tab_train), sum(tab_train$event)))
  cat(sprintf("  Test:  %d patients, %d events\n", nrow(tab_test),  sum(tab_test$event)))
  
  sel <- select_features(tab_train, all_cols, sum(tab_train$event), "Tabular-Only")
  
  coxFit_tabular <- fit_ridge_cox(tab_train, sel, "Tabular-Only")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_tabular)$concordance[1]))
  
  tab_test_final <- tab_test %>% select(PTID, time, event, all_of(sel))
  metrics_tabular <- compute_all_metrics_with_figures(
    coxFit_tabular, tab_test_final, NULL, NULL, "Tabular-Only NN"
  )
  cat("\n✓ METHOD 3 COMPLETED\n")
}, error = function(e) cat(sprintf("\n⚠ ERROR in Tabular-Only: %s\n", e$message)))

# ============================================================
# METHOD 4: CONCATENATION FUSION
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 4: CONCATENATION FUSION\n")
cat("============================================================\n")

metrics_concat_latent <- NULL
tryCatch({
  if (!file.exists("concat_patient_level.csv")) stop("concat_patient_level.csv not found")
  
  long_csv <- read.csv("latent_concat.csv", stringsAsFactors = FALSE)
  long_csv$PTID <- as.character(long_csv$PTID)
  
  out <- build_split_surv("concat_patient_level.csv", method_name = "Concat")
  concat_train <- out$train
  concat_test  <- out$test
  all_cols     <- out$all_cols
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(concat_train), sum(concat_train$event)))
  
  sel <- select_features(concat_train, all_cols, sum(concat_train$event), "Concat")
  
  coxFit_concat <- fit_ridge_cox(concat_train, sel, "Concat")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_concat)$concordance[1]))
  
  latent_vals_train <- concat_train %>% select(PTID, all_of(sel))
  long_train <- long_csv %>%
    filter(PTID %in% concat_train$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_train, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_train <- long_train[complete.cases(long_train), ]
  single_v   <- long_train %>% group_by(PTID) %>% filter(n()<2) %>% pull(PTID) %>% unique()
  long_train     <- long_train     %>% filter(!PTID %in% single_v)
  concat_train   <- concat_train   %>% filter(!PTID %in% single_v)
  common_tr      <- intersect(concat_train$PTID, long_train$PTID)
  long_train     <- long_train   %>% filter(PTID %in% common_tr)
  concat_train   <- concat_train %>% filter(PTID %in% common_tr)
  
  lf <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_concat <- tryCatch(
    lme(lf, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  cat("  ✓ LME fitted\n")
  
  uni_c_jm <- sapply(sel, function(f) tryCatch(
    summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                  data = concat_train))$concordance[1],
    error = function(e) 0.5
  ))
  top1_concat <- sel[which.max(uni_c_jm)]
  cat(sprintf("  JM Cox feature: %s\n", top1_concat))
  coxFit_concat_jm <- coxph(
    as.formula(paste("Surv(time,event) ~", top1_concat)),
    data = concat_train, x = TRUE, method = "breslow"
  )
  lf_jm <- as.formula(paste("MMSE ~ Years_bl +", top1_concat))
  lmeFit_concat_jm <- tryCatch(
    lme(lf_jm, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf_jm, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  
  jointFit_concat <- tryCatch(
    jm(coxFit_concat_jm, lmeFit_concat_jm, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) {
      cat(sprintf("  ⚠ JM (value+slope) failed: %s\n", e$message))
      tryCatch(
        jm(coxFit_concat_jm, lmeFit_concat_jm, time_var="Years_bl",
           functional_forms=~value(MMSE),
           n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
        error=function(e2) { cat(sprintf("  ⚠ JM (value only) failed: %s\n", e2$message)); NULL }
      )
    }
  )
  if (!is.null(jointFit_concat)) cat("  ✓ Joint model fitted\n") else cat("  ⚠ No joint model\n")
  
  latent_vals_test  <- concat_test %>% select(PTID, all_of(sel))
  concat_test_final <- concat_test %>% select(PTID, time, event, all_of(sel))
  long_test <- long_csv %>%
    filter(PTID %in% concat_test_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_test, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_test <- long_test[complete.cases(long_test), ]
  
  metrics_concat_latent <- compute_all_metrics_with_figures(
    coxFit_concat, concat_test_final, long_test, jointFit_concat,
    "Concatenation Fusion", patient_level_data = concat_test,
    cox_vars_hint = top1_concat
  )
  cat("\n✓ METHOD 4 COMPLETED\n")
}, error=function(e) cat(sprintf("\n⚠ ERROR in Concat: %s\n", e$message)))

# ============================================================
# METHOD 5: NO AUTOENCODER
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 5: NO AUTOENCODER\n")
cat("============================================================\n")

metrics_noae_latent <- NULL
tryCatch({
  if (!file.exists("no_ae_patient_level.csv")) stop("no_ae_patient_level.csv not found")
  
  long_csv <- read.csv("latent_no_ae.csv", stringsAsFactors = FALSE)
  long_csv$PTID <- as.character(long_csv$PTID)
  
  out <- build_split_surv("no_ae_patient_level.csv", method_name = "No-AE")
  noae_train <- out$train
  noae_test  <- out$test
  all_cols   <- out$all_cols
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(noae_train), sum(noae_train$event)))
  
  sel <- select_features(noae_train, all_cols, sum(noae_train$event), "No-AE")
  
  coxFit_noae <- fit_ridge_cox(noae_train, sel, "No-AE")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_noae)$concordance[1]))
  
  latent_vals_train <- noae_train %>% select(PTID, all_of(sel))
  long_train <- long_csv %>%
    filter(PTID %in% noae_train$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_train, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_train <- long_train[complete.cases(long_train), ]
  single_v   <- long_train %>% group_by(PTID) %>% filter(n()<2) %>% pull(PTID) %>% unique()
  long_train  <- long_train  %>% filter(!PTID %in% single_v)
  noae_train  <- noae_train  %>% filter(!PTID %in% single_v)
  common_tr   <- intersect(noae_train$PTID, long_train$PTID)
  long_train  <- long_train  %>% filter(PTID %in% common_tr)
  noae_train  <- noae_train  %>% filter(PTID %in% common_tr)
  
  lf <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_noae <- tryCatch(
    lme(lf, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  cat("  ✓ LME fitted\n")
  
  uni_c_jm <- sapply(sel, function(f) tryCatch(
    summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                  data = noae_train))$concordance[1],
    error = function(e) 0.5
  ))
  top1_noae <- sel[which.max(uni_c_jm)]
  cat(sprintf("  JM Cox feature: %s\n", top1_noae))
  coxFit_noae_jm <- coxph(
    as.formula(paste("Surv(time,event) ~", top1_noae)),
    data = noae_train, x = TRUE, method = "breslow"
  )
  lf_jm <- as.formula(paste("MMSE ~ Years_bl +", top1_noae))
  lmeFit_noae_jm <- tryCatch(
    lme(lf_jm, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf_jm, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  
  jointFit_noae <- tryCatch(
    jm(coxFit_noae_jm, lmeFit_noae_jm, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) {
      cat(sprintf("  ⚠ JM (value+slope) failed: %s\n", e$message))
      tryCatch(
        jm(coxFit_noae_jm, lmeFit_noae_jm, time_var="Years_bl",
           functional_forms=~value(MMSE),
           n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
        error=function(e2) { cat(sprintf("  ⚠ JM (value only) failed: %s\n", e2$message)); NULL }
      )
    }
  )
  if (!is.null(jointFit_noae)) cat("  ✓ Joint model fitted\n") else cat("  ⚠ No joint model\n")
  
  latent_vals_test <- noae_test %>% select(PTID, all_of(sel))
  noae_test_final  <- noae_test %>% select(PTID, time, event, all_of(sel))
  long_test <- long_csv %>%
    filter(PTID %in% noae_test_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_test, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_test <- long_test[complete.cases(long_test), ]
  
  metrics_noae_latent <- compute_all_metrics_with_figures(
    coxFit_noae, noae_test_final, long_test, jointFit_noae,
    "No Autoencoder", patient_level_data = noae_test,
    cox_vars_hint = top1_noae
  )
  cat("\n✓ METHOD 5 COMPLETED\n")
}, error=function(e) cat(sprintf("\n⚠ ERROR in No-AE: %s\n", e$message)))

# ============================================================
# METHOD 6: AE-ONLY (NO MTL)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 6: AE-ONLY / NO-MTL\n")
cat("============================================================\n")

metrics_no_mtl_latent <- NULL
tryCatch({
  if (!file.exists("ae_only_patient_level.csv")) stop("ae_only_patient_level.csv not found")
  
  long_csv <- read.csv("latent_ae_only.csv", stringsAsFactors = FALSE)
  long_csv$PTID <- as.character(long_csv$PTID)
  
  out <- build_split_surv("ae_only_patient_level.csv", method_name = "AE-Only")
  ae_train <- out$train
  ae_test  <- out$test
  all_cols <- out$all_cols
  
  cat(sprintf("  Train: %d patients, %d events\n", nrow(ae_train), sum(ae_train$event)))
  cat("  NOTE: AE-only has no survival signal — lower C-index expected\n")
  
  sel <- select_features(ae_train, all_cols, sum(ae_train$event), "AE-Only")
  
  coxFit_ae <- fit_ridge_cox(ae_train, sel, "AE-Only")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_ae)$concordance[1]))
  
  latent_vals_train <- ae_train %>% select(PTID, all_of(sel))
  long_train <- long_csv %>%
    filter(PTID %in% ae_train$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_train, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_train <- long_train[complete.cases(long_train), ]
  single_v   <- long_train %>% group_by(PTID) %>% filter(n()<2) %>% pull(PTID) %>% unique()
  long_train <- long_train %>% filter(!PTID %in% single_v)
  ae_train   <- ae_train   %>% filter(!PTID %in% single_v)
  common_tr  <- intersect(ae_train$PTID, long_train$PTID)
  long_train <- long_train %>% filter(PTID %in% common_tr)
  ae_train   <- ae_train   %>% filter(PTID %in% common_tr)
  
  lf <- as.formula(paste("MMSE ~ Years_bl +", paste(sel, collapse=" + ")))
  lmeFit_ae <- tryCatch(
    lme(lf, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  cat("  ✓ LME fitted\n")
  
  uni_c_jm <- sapply(sel, function(f) tryCatch(
    summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                  data = ae_train))$concordance[1],
    error = function(e) 0.5
  ))
  top1_ae <- sel[which.max(uni_c_jm)]
  cat(sprintf("  JM Cox feature: %s\n", top1_ae))
  coxFit_ae_jm <- coxph(
    as.formula(paste("Surv(time,event) ~", top1_ae)),
    data = ae_train, x = TRUE, method = "breslow"
  )
  lf_jm <- as.formula(paste("MMSE ~ Years_bl +", top1_ae))
  lmeFit_ae_jm <- tryCatch(
    lme(lf_jm, random=~Years_bl|PTID, data=long_train,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) lme(lf_jm, random=~1|PTID, data=long_train,
                          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
  )
  
  jointFit_ae <- tryCatch(
    jm(coxFit_ae_jm, lmeFit_ae_jm, time_var="Years_bl",
       functional_forms=~value(MMSE)+slope(MMSE),
       n_iter=10000, n_burnin=2000, n_thin=5, n_chains=3, seed=42),
    error=function(e) {
      cat(sprintf("  ⚠ JM (value+slope) failed: %s\n", e$message))
      tryCatch(
        jm(coxFit_ae_jm, lmeFit_ae_jm, time_var="Years_bl",
           functional_forms=~value(MMSE),
           n_iter=7000, n_burnin=1500, n_thin=5, n_chains=2, seed=42),
        error=function(e2) { cat(sprintf("  ⚠ JM (value only) failed: %s\n", e2$message)); NULL }
      )
    }
  )
  if (!is.null(jointFit_ae)) cat("  ✓ Joint model fitted\n") else cat("  ⚠ No joint model\n")
  
  latent_vals_test <- ae_test %>% select(PTID, all_of(sel))
  ae_test_final    <- ae_test %>% select(PTID, time, event, all_of(sel))
  long_test <- long_csv %>%
    filter(PTID %in% ae_test_final$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    left_join(latent_vals_test, by = "PTID") %>%
    select(PTID, Years_bl, MMSE, all_of(sel)) %>%
    arrange(PTID, Years_bl)
  long_test <- long_test[complete.cases(long_test), ]
  
  metrics_no_mtl_latent <- compute_all_metrics_with_figures(
    coxFit_ae, ae_test_final, long_test, jointFit_ae,
    "AE-Only (No MTL)", patient_level_data = ae_test,
    cox_vars_hint = top1_ae
  )
  cat("\n✓ METHOD 6 COMPLETED\n")
}, error=function(e) cat(sprintf("\n⚠ ERROR in AE-Only: %s\n", e$message)))

# ============================================================
# METHOD 7: MAE-JM (THESIS MODEL)
# ============================================================

cat("\n\n============================================================\n")
cat("METHOD 7: MAE-JM (THESIS MODEL)\n")
cat("============================================================\n")

metrics_thesis_latent <- NULL
tryCatch({
  if (!file.exists("latent_patient_level.csv"))
    stop("latent_patient_level.csv not found")
  if (!file.exists("latent_improved_autoencoder.csv"))
    stop("latent_improved_autoencoder.csv not found")
  
  patient_level_full <- read.csv("latent_patient_level.csv", stringsAsFactors=FALSE)
  thesis_data_full   <- read.csv("latent_improved_autoencoder.csv", stringsAsFactors=FALSE)
  patient_level_full$PTID <- as.character(patient_level_full$PTID)
  thesis_data_full$PTID   <- as.character(thesis_data_full$PTID)
  
  if (!"split" %in% names(patient_level_full))
    stop("No 'split' column in latent_patient_level.csv")
  
  cat(sprintf("  Full dataset: %d patients\n", length(unique(patient_level_full$PTID))))
  cat("  Split:\n"); print(table(patient_level_full$split, useNA="always"))
  
  train_pl_m7 <- patient_level_full %>% filter(split == "train")
  test_pl_m7  <- patient_level_full %>% filter(split == "val")
  train_lr_m7 <- thesis_data_full   %>% filter(split == "train")
  test_lr_m7  <- thesis_data_full   %>% filter(split == "val")
  
  cat(sprintf("  Train: %d | Test: %d patients\n", nrow(train_pl_m7), nrow(test_pl_m7)))
  if (nrow(train_pl_m7) == 0) stop("No training patients found")
  
  compute_slope_m7 <- function(years, mmse) {
    if (length(years) < 2 || is.na(var(years)) || var(years) <= 1e-10) return(0.0)
    tryCatch(as.numeric(coef(lm(mmse ~ years))[2]), error = function(e) 0.0)
  }
  compute_clinical_m7 <- function(long_df) {
    long_df %>%
      mutate(PTID=as.character(PTID), MMSE=as.numeric(MMSE), Years_bl=as.numeric(Years_bl)) %>%
      group_by(PTID) %>% arrange(Years_bl, .by_group=TRUE) %>%
      summarize(mmse_slope_clinical = compute_slope_m7(Years_bl, MMSE),
                mmse_last           = as.numeric(last(MMSE)),
                .groups = "drop") %>%
      select(PTID, mmse_slope_clinical, mmse_last)
  }
  CLINICAL_COLS_M7 <- c("mmse_slope_clinical", "mmse_last")
  
  clin_train_m7 <- compute_clinical_m7(train_lr_m7)
  clin_test_m7  <- compute_clinical_m7(test_lr_m7)
  
  for (col in CLINICAL_COLS_M7) {
    if (col %in% names(train_pl_m7)) train_pl_m7[[col]] <- NULL
    if (col %in% names(test_pl_m7))  test_pl_m7[[col]]  <- NULL
  }
  train_pl_m7 <- train_pl_m7 %>% left_join(clin_train_m7, by="PTID")
  test_pl_m7  <- test_pl_m7  %>% left_join(clin_test_m7,  by="PTID")
  
  z_final_cols_m7 <- grep("^z_final_", names(patient_level_full), value=TRUE)
  z_slope_cols_m7 <- grep("^z_slope_", names(patient_level_full), value=TRUE)
  keep_var <- function(cols, data)
    cols[sapply(cols, function(f) { v <- var(as.numeric(data[[f]]), na.rm=TRUE); !is.na(v) && v>1e-10 })]
  keep_bin <- function(cols, data)
    cols[sapply(cols, function(f) mean(data[[f]]==0, na.rm=TRUE) < 0.40)]
  z_final_cols_m7 <- keep_bin(keep_var(z_final_cols_m7, train_pl_m7), train_pl_m7)
  z_slope_cols_m7 <- keep_bin(keep_var(z_slope_cols_m7, train_pl_m7), train_pl_m7)
  latent_pool_m7  <- c(z_final_cols_m7, z_slope_cols_m7)
  cat(sprintf("  Latent features after filters: %d\n", length(latent_pool_m7)))
  
  all_cols_m7 <- c(CLINICAL_COLS_M7, latent_pool_m7)
  std_m7      <- standardise_features(train_pl_m7, test_pl_m7, all_cols_m7)
  train_pl_m7 <- std_m7$train
  test_pl_m7  <- std_m7$test
  
  surv_train_sel <- train_pl_m7 %>%
    select(PTID, time_to_event, event, all_of(all_cols_m7)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  surv_train_sel <- surv_train_sel[complete.cases(surv_train_sel), ]
  n_events_m7    <- sum(surv_train_sel$event)
  cat(sprintf("  Train for selection: %d patients, %d events\n",
              nrow(surv_train_sel), n_events_m7))
  
  final_features_m7 <- CLINICAL_COLS_M7
  if (length(latent_pool_m7) > 0) {
    uni_c_m7 <- sapply(latent_pool_m7, function(f) tryCatch(
      summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                    data=surv_train_sel))$concordance[1],
      error=function(e) 0.5
    ))
    best_latent <- latent_pool_m7[which.max(uni_c_m7)]
    cat(sprintf("  Best latent: %s (C=%.4f)\n", best_latent, max(uni_c_m7)))
    final_features_m7 <- c(final_features_m7, best_latent)
  }
  cat(sprintf("  Final features: %s\n", paste(final_features_m7, collapse=", ")))
  
  train_surv_m7 <- train_pl_m7 %>%
    select(PTID, time_to_event, event, all_of(final_features_m7)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  train_surv_m7 <- train_surv_m7[complete.cases(train_surv_m7), ]
  
  latent_in_model <- setdiff(final_features_m7, CLINICAL_COLS_M7)
  train_long_m7 <- train_lr_m7 %>%
    filter(PTID %in% train_surv_m7$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    { if (length(latent_in_model)>0)
      left_join(., train_pl_m7 %>% select(PTID, all_of(latent_in_model)), by="PTID")
      else . } %>%
    select(PTID, Years_bl, MMSE, any_of(latent_in_model)) %>% arrange(PTID, Years_bl)
  train_long_m7 <- train_long_m7[complete.cases(train_long_m7), ]
  
  common_tr     <- intersect(train_surv_m7$PTID, train_long_m7$PTID)
  train_surv_m7 <- train_surv_m7 %>% filter(PTID %in% common_tr)
  train_long_m7 <- train_long_m7 %>% filter(PTID %in% common_tr)
  single_v_m7   <- train_long_m7 %>% group_by(PTID) %>%
    summarize(n=n(),.groups="drop") %>% filter(n<2) %>% pull(PTID)
  if (length(single_v_m7)>0) {
    train_long_m7 <- train_long_m7 %>% filter(!PTID %in% single_v_m7)
    train_surv_m7 <- train_surv_m7 %>% filter(!PTID %in% single_v_m7)
  }
  cat(sprintf("  Train final: %d patients, %d events, %d visits\n",
              nrow(train_surv_m7), sum(train_surv_m7$event), nrow(train_long_m7)))
  
  test_surv_m7 <- test_pl_m7 %>%
    select(PTID, time_to_event, event, all_of(final_features_m7)) %>%
    filter(!is.na(time_to_event), !is.na(event), time_to_event > 0) %>%
    rename(time = time_to_event)
  test_surv_m7 <- test_surv_m7[complete.cases(test_surv_m7), ]
  
  test_long_m7 <- test_lr_m7 %>%
    filter(PTID %in% test_surv_m7$PTID, !is.na(MMSE), !is.na(Years_bl)) %>%
    { if (length(latent_in_model)>0)
      left_join(., test_pl_m7 %>% select(PTID, all_of(latent_in_model)), by="PTID")
      else . } %>%
    select(PTID, Years_bl, MMSE, any_of(latent_in_model)) %>% arrange(PTID, Years_bl)
  test_long_m7 <- test_long_m7[complete.cases(test_long_m7), ]
  common_te    <- intersect(test_surv_m7$PTID, test_long_m7$PTID)
  test_surv_m7 <- test_surv_m7 %>% filter(PTID %in% common_te)
  test_long_m7 <- test_long_m7 %>% filter(PTID %in% common_te)
  cat(sprintf("  Test final: %d patients, %d events, %d visits\n",
              nrow(test_surv_m7), sum(test_surv_m7$event), nrow(test_long_m7)))
  
  sf_m7 <- as.formula(paste("Surv(time,event) ~", paste(final_features_m7, collapse=" + ")))
  coxFit_m7 <- coxph(sf_m7, data=train_surv_m7, x=TRUE, method="breslow")
  cat(sprintf("  Train C-index: %.4f\n", summary(coxFit_m7)$concordance[1]))
  
  lf_m7 <- if (length(latent_in_model)>0)
    as.formula(paste("MMSE ~ Years_bl +", paste(latent_in_model, collapse=" + ")))
  else as.formula("MMSE ~ Years_bl")
  
  lmeFit_m7 <- tryCatch(
    lme(lf_m7, random=~Years_bl|PTID, data=train_long_m7,
        control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE)),
    error=function(e) {
      cat("  ⚠ Random slope failed — intercept only\n")
      lme(lf_m7, random=~1|PTID, data=train_long_m7,
          control=lmeControl(opt="optim",maxIter=500,returnObject=TRUE))
    }
  )
  cat("  ✓ LME fitted\n")
  
  # For the joint model use single best univariate feature — same pattern as
  # Methods 4-6. Using all 3 Cox covariates dilutes the MMSE trajectory signal
  # and produces near-random dynamic AUC. The full coxFit_m7 (all features)
  # is still used for C-index, RMST, calibration, and Brier score.
  uni_c_jm_m7 <- sapply(final_features_m7, function(f) tryCatch(
    summary(coxph(as.formula(paste("Surv(time,event) ~", f)),
                  data = train_surv_m7))$concordance[1],
    error = function(e) 0.5
  ))
  top1_m7 <- final_features_m7[which.max(uni_c_jm_m7)]
  cat(sprintf("  JM Cox feature (single): %s (C=%.4f)\n",
              top1_m7, max(uni_c_jm_m7)))
  
  coxFit_m7_jm <- coxph(
    as.formula(paste("Surv(time,event) ~", top1_m7)),
    data = train_surv_m7, x = TRUE, method = "breslow"
  )
  lf_m7_jm <- as.formula(paste("MMSE ~ Years_bl +", top1_m7))
  lmeFit_m7_jm <- tryCatch(
    lme(lf_m7_jm, random = ~Years_bl | PTID, data = train_long_m7,
        control = lmeControl(opt = "optim", maxIter = 500, returnObject = TRUE)),
    error = function(e) lme(lf_m7_jm, random = ~1 | PTID, data = train_long_m7,
                            control = lmeControl(opt = "optim", maxIter = 500,
                                                 returnObject = TRUE))
  )
  
  cat("  Fitting joint model (single Cox feature, 20k iter)...\n")
  jointFit_m7 <- tryCatch(
    jm(coxFit_m7_jm, lmeFit_m7_jm, time_var = "Years_bl",
       functional_forms = ~value(MMSE) + slope(MMSE),
       n_iter = 20000, n_burnin = 5000, n_thin = 5, n_chains = 3, seed = 42),
    error = function(e) {
      cat(sprintf("  ⚠ Full JM failed (%s) — value only\n", e$message))
      tryCatch(
        jm(coxFit_m7_jm, lmeFit_m7_jm, time_var = "Years_bl",
           functional_forms = ~value(MMSE),
           n_iter = 12000, n_burnin = 3000, n_thin = 5, n_chains = 2, seed = 42),
        error = function(e2) { cat("  ⚠ JM failed entirely\n"); NULL }
      )
    }
  )
  if (!is.null(jointFit_m7)) cat("  ✓ Joint model fitted\n")
  
  metrics_thesis_latent <- compute_all_metrics_with_figures(
    coxFit_m7, test_surv_m7, test_long_m7, jointFit_m7, "MAE-JM",
    patient_level_data = test_pl_m7,
    cox_vars_hint      = top1_m7   # single feature matches JM Cox submodel
  )
  cat("\n✓ METHOD 7 COMPLETED\n")
}, error=function(e) {
  cat(sprintf("\n⚠ ERROR in MAE-JM: %s\n", e$message)); traceback()
})

# ============================================================
# FINAL COMPARISON TABLE
# ============================================================

cat("\n\n============================================================\n")
cat("FINAL COMPARISON\n")
cat("Each method uses its own CSV split column for train/test.\n")
cat("============================================================\n\n")

all_metrics <- Filter(Negate(is.null), list(
  metrics_baseline, metrics_image, metrics_tabular,
  metrics_concat_latent, metrics_noae_latent,
  metrics_no_mtl_latent, metrics_thesis_latent
))

if (length(all_metrics) > 0) {
  gv <- function(x, f) { v <- x[[f]]; if (is.null(v)) NA else v }
  
  comparison <- data.frame(
    Method         = sapply(all_metrics, gv, "name"),
    C_index        = sapply(all_metrics, gv, "cindex"),
    CI_Lower       = sapply(all_metrics, gv, "cindex_ci_lower"),
    CI_Upper       = sapply(all_metrics, gv, "cindex_ci_upper"),
    Brier_3yr      = sapply(all_metrics, gv, "brier_3yr"),
    Cal_Slope      = sapply(all_metrics, gv, "calibration_slope"),
    RMST_Diff      = sapply(all_metrics, gv, "rmst_diff"),
    RMST_pval      = sapply(all_metrics, gv, "rmst_pval"),
    MAE_yr         = sapply(all_metrics, gv, "mae_time"),
    JM_Dynamic_AUC = sapply(all_metrics, gv, "jm_dynamic_auc"),
    JM_alpha_val   = sapply(all_metrics, gv, "jm_association_val"),
    JM_alpha_slope = sapply(all_metrics, gv, "jm_association_slope"),
    N_Patients     = sapply(all_metrics, gv, "n_patients"),
    N_Events       = sapply(all_metrics, gv, "n_events"),
    stringsAsFactors = FALSE
  )
  
  bc <- comparison$C_index[comparison$Method == "Clinical Cox"]
  comparison$Improvement_Pct <- if (length(bc) > 0 && !is.na(bc))
    ((comparison$C_index - bc) / bc) * 100 else NA
  
  comparison <- comparison[order(-comparison$C_index), ]
  
  cat("DISCRIMINATION:\n")
  print(comparison[, c("Method","C_index","CI_Lower","CI_Upper","Improvement_Pct")])
  cat("\nCLINICAL UTILITY:\n")
  print(comparison[, c("Method","Brier_3yr","Cal_Slope","RMST_Diff","RMST_pval")])
  cat("\nJOINT MODEL:\n")
  print(comparison[, c("Method","JM_Dynamic_AUC","JM_alpha_val","JM_alpha_slope")])
  
  write.csv(comparison, "thesis_figures/FINAL_COMPARISON.csv", row.names=FALSE)
  
  pub_table <- comparison %>% mutate(
    `C-index (95% CI)` = ifelse(!is.na(CI_Lower),
                                sprintf("%.3f (%.3f-%.3f)", C_index, CI_Lower, CI_Upper),
                                sprintf("%.3f", C_index)),
    Improvement    = sprintf("%+.1f%%", Improvement_Pct),
    `Brier (3yr)`  = sprintf("%.3f", Brier_3yr),
    `Cal. Slope`   = sprintf("%.3f", Cal_Slope),
    `Dynamic AUC`  = ifelse(!is.na(JM_Dynamic_AUC),
                            sprintf("%.3f", JM_Dynamic_AUC), "—"),
    `alpha(value)` = ifelse(!is.na(JM_alpha_val),
                            sprintf("%.3f", JM_alpha_val), "—"),
    `alpha(slope)` = ifelse(!is.na(JM_alpha_slope),
                            sprintf("%.3f", JM_alpha_slope), "—")
  ) %>% select(Method, `C-index (95% CI)`, Improvement,
               `Brier (3yr)`, `Cal. Slope`, `Dynamic AUC`,
               `alpha(value)`, `alpha(slope)`)
  
  write.csv(pub_table, "thesis_figures/PUBLICATION_TABLE.csv", row.names=FALSE)
  
  best <- comparison[1, ]
  cat(sprintf("\n BEST: %s | C-index=%.4f\n", best$Method, best$C_index))
  cat("\n✓ ANALYSIS COMPLETE\n")
  cat("  thesis_figures/FINAL_COMPARISON.csv\n")
  cat("  thesis_figures/PUBLICATION_TABLE.csv\n")
  cat("  thesis_figures/ALL_METHODS_RESULTS.csv\n")
}
