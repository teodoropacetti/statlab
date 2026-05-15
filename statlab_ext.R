# =============================================================================
# statlab_ext.R - Extension to statlab.R
# Covers exercises from weeks 5-8 not yet present in the original library.
#
# Usage:
#   source("statlab.R")    # load original first
#   source("statlab_ext.R")
#   sl_setup_ext()
#
# NOTE: the last ~132 lines of statlab.R were not fetchable.
# The following functions may already be defined there:
#   sl_boot_plot, sl_lm_diagnostics, sl_vif, sl_step_aic, sl_predict_lm
# If R complains about duplicate definitions just remove the duplicates here.
#
# NOTE ON SEEDS: never set inside functions; always call set.seed() in your
# script before any random-number-generating call.
#
# NEW PACKAGES: pROC, mclust, car
# (mvtnorm, corrplot, skimr, dplyr already loaded by sl_setup)
# =============================================================================

# =============================================================================
# 0. SETUP EXTENSION
# =============================================================================

#' Install and load packages needed by statlab_ext (pROC, mclust, car)
#'
#' Call this once after sl_setup().
#'
#' @examples
#' sl_setup()
#' sl_setup_ext()
sl_setup_ext <- function() {
  pkgs <- c("pROC", "mclust", "car")
  new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(new)) install.packages(new, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
  cat("Extension packages loaded:", paste(pkgs, collapse = ", "), "\n")
}


# =============================================================================
# 1. BOOTSTRAP PLOT  (may already be in statlab.R - skip if so)
# =============================================================================

#' Histogram of a bootstrap distribution with estimate and CI bounds (Ex 12.7, 13.3)
#'
#' Plots the bootstrap distribution, marks the original estimate with a
#' vertical line, and adds the percentile CI bounds. A legend explains each
#' line. Matches the style used in the official solutions.
#'
#' @param boot_result  Output of sl_bootstrap(), sl_bootstrap2(), or
#'                     sl_bootstrap_gen().
#' @param level        Confidence level for the percentile CI (default 0.95).
#' @param title        Plot title.
#' @param breaks       Number of histogram bins.
#' @param col          Bar fill colour.
#'
#' @examples
#' set.seed(42)
#' boot12 <- sl_bootstrap2(serverA, serverB, fun = mean)
#' sl_boot_plot(boot12, level = 0.95, title = "Bootstrap: diff of means")
#'
#' set.seed(42)
#' boot13 <- sl_bootstrap(calorie, fun = sd, B = 2000)
#' sl_boot_plot(boot13, level = 0.99, title = "Bootstrap: sd of calorie intake")
sl_boot_plot <- function(boot_result, level = 0.95,
                         title = "Bootstrap distribution",
                         breaks = 30, col = "dodgerblue") {
  d     <- boot_result$distribution
  est   <- boot_result$estimate
  alpha <- 1 - level
  ci    <- quantile(d, probs = c(alpha / 2, 1 - alpha / 2))

  hist(d, breaks = breaks, freq = FALSE,
       main = title, xlab = "Bootstrap replications",
       col = col, border = "white")
  lines(density(d), col = "black", lwd = 1.5)
  abline(v = est,   col = "black",  lwd = 2,   lty = 1)
  abline(v = ci[1], col = "darkorange", lwd = 2, lty = 2)
  abline(v = ci[2], col = "darkorange", lwd = 2, lty = 2)
  legend("topright",
         legend = c("Original estimate",
                    sprintf("%.0f%% percentile CI", level * 100)),
         col    = c("black", "darkorange"),
         lty    = c(1, 2), lwd = 2, bty = "n")
  invisible(ci)
}


# =============================================================================
# 2. LINEAR REGRESSION EXTENSIONS
#    (sl_lm_diagnostics / sl_vif / sl_step_aic / sl_predict_lm may already
#    exist in the last section of statlab.R — skip duplicates as needed)
# =============================================================================

#' 4-panel residual diagnostic plots + Shapiro-Wilk test (Ex 14.8, 16.5, 19.3)
#'
#' Produces the four standard R diagnostic plots (Residuals vs Fitted,
#' Normal QQ, Scale-Location, Residuals vs Leverage) in a 2x2 grid,
#' then runs the Shapiro-Wilk normality test on the residuals.
#'
#' @param mod lm object
#'
#' @examples
#' sl_lm_diagnostics(mod14)
sl_lm_diagnostics <- function(mod) {
  old_par <- par(mfrow = c(2, 2))
  on.exit(par(old_par))
  plot(mod)
  sw <- shapiro.test(residuals(mod))
  cat("\n=== Shapiro-Wilk test on residuals ===\n")
  cat(sprintf("W = %.4f,  p-value = %.4f\n", sw$statistic, sw$p.value))
  if (sw$p.value < 0.05)
    cat("-> Normality assumption is REJECTED at 5% level.\n")
  else
    cat("-> No evidence against normality at 5% level.\n")
  invisible(sw)
}


#' Variance Inflation Factors for all predictors (Ex 17.5, 19.3, 20.5)
#'
#' A VIF > 5 is often considered a warning; > 10 indicates severe
#' multicollinearity.
#'
#' @param mod lm object
#'
#' @examples
#' sl_vif(mod17)
sl_vif <- function(mod) {
  if (!requireNamespace("car", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  v <- car::vif(mod)
  cat("=== Variance Inflation Factors ===\n")
  print(round(v, 4))
  high <- names(v)[v > 5]
  if (length(high))
    cat(sprintf("\nVIF > 5: %s\n", paste(high, collapse = ", ")))
  else
    cat("\nNo VIF > 5 detected.\n")
  invisible(v)
}


#' AIC-based stepwise variable selection for a linear model (Ex 19.4, 20.6)
#'
#' Wraps step() with direction = "backward" by default and prints a clean
#' summary of the selected model including R^2 and AIC.
#'
#' @param mod       lm object (full model to start from)
#' @param direction "backward" (default), "forward", or "both"
#' @param trace     Print each AIC step? (default TRUE)
#'
#' @examples
#' mod_full <- sl_lm(rating ~ profitability + capital + flexibility, ratings)
#' mod_sel  <- sl_step_aic(mod_full)
sl_step_aic <- function(mod, direction = "backward", trace = TRUE) {
  mod_sel <- step(mod, direction = direction, trace = trace)
  sm      <- summary(mod_sel)
  cat("\n=== Selected model ===\n")
  cat("Formula: "); print(formula(mod_sel))
  cat(sprintf("AIC  : %.4f\n", AIC(mod_sel)))
  cat(sprintf("R^2  : %.4f\n", sm$r.squared))
  cat(sprintf("R^2adj: %.4f\n", sm$adj.r.squared))
  invisible(mod_sel)
}


#' Confidence intervals for regression coefficients (Ex 17.4, 19.5, 22.1)
#'
#' @param mod   lm object
#' @param level confidence level (default 0.95)
#'
#' @examples
#' sl_lm_confint(mod14, level = 0.95)
#' sl_lm_confint(mod_sel, level = 0.99)
sl_lm_confint <- function(mod, level = 0.95) {
  ci <- confint(mod, level = level)
  cat(sprintf("=== %.0f%% Confidence intervals for coefficients ===\n",
              level * 100))
  print(round(ci, 4))
  invisible(ci)
}


#' Predict for new observations with confidence or prediction intervals (Ex 19.5, 21, 22)
#'
#' @param mod      lm object
#' @param new_data data.frame with the same predictor names as in mod
#' @param interval "prediction" (default) or "confidence"
#' @param level    confidence level (default 0.95)
#'
#' @examples
#' new_firm <- data.frame(profitability = 3.9, capital = 6.02, flexibility = 1.43)
#' sl_predict_lm(mod_sel, new_firm, interval = "prediction", level = 0.95)
#' sl_predict_lm(mod_sel, new_firm, interval = "confidence", level = 0.99)
sl_predict_lm <- function(mod, new_data, interval = "prediction", level = 0.95) {
  pred <- predict(mod, newdata = new_data, interval = interval, level = level)
  cat(sprintf("=== %.0f%% %s interval ===\n",
              level * 100, interval))
  print(round(pred, 4))
  invisible(pred)
}


#' Bootstrap confidence intervals for all regression coefficients (Ex 21.2)
#'
#' Resamples rows of the model data with replacement and refits the model
#' B times. Returns percentile CIs and optionally plots the bootstrap
#' distribution of each coefficient.
#'
#' NOTE: call set.seed() in your script before this function.
#'
#' @param mod   lm object
#' @param B     number of bootstrap replications (default 1000)
#' @param level confidence level (default 0.95)
#' @param plot  produce one histogram per coefficient? (default TRUE)
#'
#' @examples
#' set.seed(42)
#' sl_bootstrap_lm_ci(mod_sel, B = 1000, level = 0.95)
sl_bootstrap_lm_ci <- function(mod, B = 1000, level = 0.95, plot = TRUE) {
  df_orig <- model.frame(mod)
  form    <- formula(mod)
  n       <- nrow(df_orig)
  alpha   <- 1 - level

  boot_mat <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    coef(lm(form, data = df_orig[idx, , drop = FALSE]))
  })

  ci_mat <- t(apply(boot_mat, 1, quantile,
                    probs = c(alpha / 2, 1 - alpha / 2),
                    na.rm = TRUE))
  colnames(ci_mat) <- c("lwr", "upr")

  ci_ref <- confint(mod, level = level)

  cat(sprintf("=== Bootstrap %.0f%% CI vs confint() ===\n", level * 100))
  comparison <- data.frame(
    boot_lwr   = round(ci_mat[, "lwr"], 4),
    boot_upr   = round(ci_mat[, "upr"], 4),
    confint_lwr = round(ci_ref[, 1], 4),
    confint_upr = round(ci_ref[, 2], 4)
  )
  print(comparison)

  if (plot) {
    p <- nrow(boot_mat)
    old_par <- par(mfrow = c(ceiling(p / 2), 2))
    on.exit(par(old_par))
    for (i in seq_len(p)) {
      nm <- rownames(boot_mat)[i]
      hist(boot_mat[i, ], breaks = 30, freq = FALSE,
           main = paste("Bootstrap dist:", nm),
           xlab = nm, col = "dodgerblue", border = "white")
      abline(v = coef(mod)[i], col = "black", lwd = 2)
      abline(v = ci_mat[i, "lwr"], col = "darkorange", lwd = 2, lty = 2)
      abline(v = ci_mat[i, "upr"], col = "darkorange", lwd = 2, lty = 2)
    }
  }

  invisible(list(ci = as.data.frame(ci_mat), boot_matrix = t(boot_mat)))
}


#' Residuals plotted against each covariate + against fitted values (Ex 16.5)
#'
#' Useful for checking linearity and homoscedasticity assumption-by-assumption.
#'
#' @param mod lm object
#'
#' @examples
#' sl_residuals_vs_covariates(mod16)
sl_residuals_vs_covariates <- function(mod) {
  r    <- residuals(mod)
  mf   <- model.frame(mod)
  pred <- model.matrix(mod)[, -1, drop = FALSE]   # drop intercept column
  p    <- ncol(pred)

  old_par <- par(mfrow = c(ceiling((p + 1) / 3), 3))
  on.exit(par(old_par))

  plot(fitted(mod), r,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted", col = "dodgerblue", pch = 16)
  abline(h = 0, col = "red", lty = 2)

  for (j in seq_len(p)) {
    plot(pred[, j], r,
         xlab = colnames(pred)[j], ylab = "Residuals",
         main = paste("Residuals vs", colnames(pred)[j]),
         col = "dodgerblue", pch = 16)
    abline(h = 0, col = "red", lty = 2)
    lines(lowess(pred[, j], r), col = "darkorange", lwd = 2)
  }
  invisible(r)
}


#' Scatterplot with parallel regression lines for a factor covariate (Ex 23.6)
#'
#' Plots the response against a numeric predictor, coloring points by the
#' binary factor. Adds the two estimated parallel regression lines (same
#' slope, different intercepts) and prints their equations.
#'
#' @param mod         lm object fitted with a numeric predictor + a factor
#' @param x_var       name of the numeric predictor (string)
#' @param group_var   name of the factor/group variable (string)
#' @param y_var       name of the response variable (string, optional;
#'                    extracted from formula if NULL)
#' @param xlab, ylab  axis labels
#'
#' @examples
#' mod23 <- lm(Gas ~ Temp + Insul, data = whiteside)
#' sl_lm_parallel_lines(mod23, x_var = "Temp", group_var = "Insul")
sl_lm_parallel_lines <- function(mod, x_var, group_var,
                                  y_var   = NULL,
                                  xlab    = x_var,
                                  ylab    = NULL,
                                  title   = "Parallel regression lines",
                                  cols    = c("dodgerblue", "darkorange")) {
  df     <- model.frame(mod)
  if (is.null(y_var))
    y_var <- as.character(formula(mod)[[2L]])
  if (is.null(ylab)) ylab <- y_var

  y     <- df[[y_var]]
  x     <- df[[x_var]]
  grp   <- as.factor(df[[group_var]])
  levs  <- levels(grp)
  coefs <- coef(mod)

  plot(x, y, col = cols[as.integer(grp)], pch = 16,
       xlab = xlab, ylab = ylab, main = title)

  # intercept for group 1 = beta_0 ; group 2 = beta_0 + beta_factor
  # slope = beta for x_var (assumed constant)
  slope <- coefs[x_var]
  int1  <- coefs["(Intercept)"]

  # find the factor dummy coefficient (could be named "InsulAfter", etc.)
  factor_coef_name <- grep(group_var, names(coefs), value = TRUE)
  int2 <- if (length(factor_coef_name))
    int1 + coefs[factor_coef_name[1]]
  else
    int1   # fallback if factor is the reference

  abline(a = int1, b = slope, col = cols[1], lwd = 2)
  abline(a = int2, b = slope, col = cols[2], lwd = 2)

  legend("topright", legend = levs, col = cols,
         pch = 16, lty = 1, lwd = 2, bty = "n")

  cat("=== Parallel line equations ===\n")
  cat(sprintf("%s = %.4f %+.4f * %s\n", y_var, int1, slope, x_var))
  cat(sprintf("    [%s = %s]\n", group_var, levs[1]))
  cat(sprintf("%s = %.4f %+.4f * %s\n", y_var, int2, slope, x_var))
  cat(sprintf("    [%s = %s]\n", group_var, levs[2]))
  invisible(list(intercept1 = int1, intercept2 = int2, slope = slope))
}


# =============================================================================
# 3. LOGISTIC REGRESSION  (Ex 29, 30)
# =============================================================================

#' Fit a logistic regression model and print a full summary (Ex 29.3, 30)
#'
#' Mirrors the style of sl_lm(): fits the model, prints summary(),
#' and returns the glm object invisibly. For coefficient interpretation
#' use sl_logistic_or().
#'
#' @param formula R formula; response must be 0/1 integer or factor
#' @param data    data.frame
#'
#' @examples
#' mod29 <- sl_logistic(cancer ~ yrsmoke + bird + highstatus, data = bird)
sl_logistic <- function(formula, data) {
  mod <- glm(formula, data = data, family = binomial(link = "logit"))
  print(summary(mod))
  cat(sprintf("\nAIC: %.4f\n", AIC(mod)))
  cat(sprintf("Null deviance    : %.4f (df = %d)\n",
              mod$null.deviance, mod$df.null))
  cat(sprintf("Residual deviance: %.4f (df = %d)\n",
              mod$deviance, mod$df.residual))
  invisible(mod)
}


#' Print the logistic regression equation on the log-odds scale (Ex 29.4)
#'
#' @param mod glm (logistic) object
#'
#' @examples
#' sl_logistic_equation(mod29)
sl_logistic_equation <- function(mod) {
  coefs  <- coef(mod)
  nomi   <- names(coefs)
  terms  <- mapply(function(nm, val) {
    if (nm == "(Intercept)") sprintf("%.4f", val)
    else sprintf("%+.4f * %s", val, nm)
  }, nomi, coefs)
  y_name <- as.character(formula(mod)[[2L]])
  eq     <- paste0("logit[P(", y_name, "=1)] = ",
                   paste(terms, collapse = " "))
  cat("Fitted logistic equation (log-odds):\n  ", eq, "\n")
  invisible(eq)
}


#' AIC-based stepwise selection for a logistic regression model (Ex 29.3, 30)
#'
#' @param mod       glm (logistic) object (full model)
#' @param direction "backward" (default), "forward", or "both"
#' @param trace     Print each step? (default TRUE)
#'
#' @examples
#' mod29_sel <- sl_step_aic_logistic(mod29_full)
sl_step_aic_logistic <- function(mod, direction = "backward", trace = TRUE) {
  mod_sel <- step(mod, direction = direction, trace = trace)
  cat("\n=== Selected logistic model ===\n")
  cat("Formula: "); print(formula(mod_sel))
  cat(sprintf("AIC: %.4f\n", AIC(mod_sel)))
  invisible(mod_sel)
}


#' Odds ratios and confidence intervals for a logistic model (Ex 29.4, 29.5, 30)
#'
#' Exponentiates the coefficients and their profile likelihood CIs.
#' Coefficients with OR > 1 increase the odds of the event; OR < 1 decrease it.
#'
#' @param mod   glm (logistic) object
#' @param level confidence level for the CIs (default 0.95)
#'
#' @examples
#' sl_logistic_or(mod29_sel, level = 0.95)
sl_logistic_or <- function(mod, level = 0.95) {
  if (!requireNamespace("stats", quietly = TRUE)) stop("stats not available")
  ci     <- suppressMessages(confint(mod, level = level))
  coefs  <- coef(mod)
  out    <- data.frame(
    log_OR = round(coefs, 4),
    OR     = round(exp(coefs), 4),
    ci_lwr = round(exp(ci[, 1]), 4),
    ci_upr = round(exp(ci[, 2]), 4)
  )
  cat(sprintf("=== Odds Ratios (%.0f%% CI) ===\n", level * 100))
  print(out)
  invisible(out)
}


#' Predicted probabilities for all units in the training set (Ex 30)
#'
#' Returns and prints a summary of the fitted probability P(Y=1|X).
#'
#' @param mod glm (logistic) object
#'
#' @examples
#' probs29 <- sl_logistic_predict_prob(mod29_sel)
#' hist(probs29, breaks = 20, col = "dodgerblue", main = "Predicted P(cancer=1)")
sl_logistic_predict_prob <- function(mod) {
  p <- fitted(mod)
  cat("=== Summary of predicted probabilities ===\n")
  print(summary(p))
  invisible(p)
}


# =============================================================================
# 4. CLASSIFICATION METRICS  (Ex 31)
# =============================================================================

#' Confusion matrix from a fitted logistic model (Ex 31.1)
#'
#' @param mod       glm (logistic) object
#' @param threshold probability cut-off (default 0.5)
#'
#' @examples
#' sl_confusion_matrix(mod29_sel)
#' sl_confusion_matrix(mod29_sel, threshold = 0.4)
sl_confusion_matrix <- function(mod, threshold = 0.5) {
  y_true <- model.response(model.frame(mod))
  y_pred <- as.integer(fitted(mod) >= threshold)

  # ensure both factor levels present even if all predictions are one class
  y_true <- factor(y_true, levels = c(0, 1))
  y_pred <- factor(y_pred, levels = c(0, 1))

  cm <- table(Predicted = y_pred, Actual = y_true)
  cat(sprintf("=== Confusion matrix (threshold = %.2f) ===\n", threshold))
  print(cm)
  invisible(cm)
}


#' Sensitivity, specificity, and accuracy from a confusion matrix (Ex 31.2)
#'
#' @param cm        2x2 table from sl_confusion_matrix() or table()
#'                  Rows = Predicted (0/1), Columns = Actual (0/1)
#'
#' @examples
#' cm29 <- sl_confusion_matrix(mod29_sel)
#' sl_classification_metrics(cm29)
sl_classification_metrics <- function(cm) {
  TP <- cm["1", "1"]; TN <- cm["0", "0"]
  FP <- cm["1", "0"]; FN <- cm["0", "1"]

  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  acc  <- (TP + TN) / sum(cm)
  ppv  <- TP / (TP + FP)    # positive predictive value / precision
  npv  <- TN / (TN + FN)    # negative predictive value

  cat("=== Classification metrics ===\n")
  cat(sprintf("Sensitivity (recall)     : %.4f  (TP / (TP+FN))\n", sens))
  cat(sprintf("Specificity              : %.4f  (TN / (TN+FP))\n", spec))
  cat(sprintf("Accuracy                 : %.4f  ((TP+TN) / N)\n",  acc))
  cat(sprintf("Precision (PPV)          : %.4f  (TP / (TP+FP))\n", ppv))
  cat(sprintf("NPV                      : %.4f  (TN / (TN+FN))\n", npv))

  invisible(list(sensitivity = sens, specificity = spec,
                 accuracy = acc, precision = ppv, npv = npv))
}


#' ROC curve and AUC for a fitted logistic model (Ex 31.3, 31.4)
#'
#' Requires the pROC package (loaded by sl_setup_ext()).
#'
#' @param mod   glm (logistic) object
#' @param title plot title
#'
#' @examples
#' sl_roc_plot(mod29_sel, title = "ROC - bird keeping study")
#' sl_roc_plot(mod30_sel, title = "ROC - urine calcium")
sl_roc_plot <- function(mod, title = "ROC curve") {
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  y_true  <- as.numeric(model.response(model.frame(mod)))
  y_prob  <- fitted(mod)
  roc_obj <- pROC::roc(y_true, y_prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))

  pROC::plot.roc(roc_obj, col = "dodgerblue", lwd = 2,
                 main = title, print.auc = TRUE,
                 auc.polygon = TRUE,
                 auc.polygon.col = adjustcolor("dodgerblue", 0.15))

  cat(sprintf("AUC = %.4f\n", auc_val))
  invisible(roc_obj)
}


# =============================================================================
# 5. MIXTURE MODELS  (Ex 32, 33)
# =============================================================================

#' Select number of components + variance structure via BIC (Ex 32.2, 33.2)
#'
#' Prints the top BIC values and plots the BIC surface. The model with the
#' highest BIC is selected (mclust convention: higher is better).
#'
#' @param x          numeric vector (univariate) or matrix (multivariate)
#' @param G          integer vector of component numbers to try (default 1:9)
#' @param model_names character vector of mclust model names to try, or NULL
#'                   for all univariate models ("E", "V")
#'
#' @examples
#' sl_mclust_select(haemoglobin)
#' sl_mclust_select(snapper, G = 1:6)
sl_mclust_select <- function(x, G = 1:9, model_names = NULL) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  args <- list(data = x, G = G)
  if (!is.null(model_names)) args$modelNames <- model_names
  bic <- do.call(mclust::mclustBIC, args)
  cat("=== BIC values (higher is better) ===\n")
  print(bic)
  plot(bic, main = "BIC for Gaussian mixture models")
  best <- summary(bic)
  cat(sprintf("\nBest model: G = %d, model = %s\n",
              best$G, best$modelName))
  invisible(bic)
}


#' Fit a Gaussian mixture model with Mclust (Ex 32.3, 33.2)
#'
#' @param x          numeric vector or matrix
#' @param G          number of components
#' @param model_name mclust model name (e.g. "V" for varying variance,
#'                   "E" for equal variance). NULL = auto.
#'
#' @examples
#' mod32 <- sl_mclust_fit(haemoglobin, G = 2)
#' mod33 <- sl_mclust_fit(snapper, G = 4, model_name = "V")
sl_mclust_fit <- function(x, G, model_name = NULL) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  args <- list(data = x, G = G)
  if (!is.null(model_name)) args$modelName <- model_name
  mod <- do.call(mclust::Mclust, args)
  print(summary(mod, parameters = TRUE))
  invisible(mod)
}


#' Print estimated parameters of each mixture component (Ex 32.4, 33.3)
#'
#' @param mod Mclust object from sl_mclust_fit()
#'
#' @examples
#' sl_mclust_params(mod32)
sl_mclust_params <- function(mod) {
  par <- mod$parameters
  G   <- mod$G
  cat("=== Mixture model parameters ===\n")
  for (k in seq_len(G)) {
    cat(sprintf("\nComponent %d (pi = %.4f):\n", k, par$pro[k]))
    cat(sprintf("  Mean : %.4f\n", par$mean[k]))
    if (!is.null(par$variance$sigmasq))
      cat(sprintf("  Var  : %.4f\n",
                  if (length(par$variance$sigmasq) == 1)
                    par$variance$sigmasq
                  else
                    par$variance$sigmasq[k]))
  }
  invisible(par)
}


#' Classification of each observation (Ex 32.5, 33.4)
#'
#' Prints a frequency table of the MAP classification and returns a
#' data.frame with observation index, value, class, and posterior
#' probability of the assigned class.
#'
#' @param mod Mclust object
#'
#' @examples
#' sl_mclust_classify(mod32)
sl_mclust_classify <- function(mod) {
  cls  <- mod$classification
  z    <- mod$z
  prob <- apply(z, 1, max)

  cat("=== Classification (MAP rule) ===\n")
  print(table(Class = cls))

  out <- data.frame(
    id          = seq_along(cls),
    observation = as.vector(mod$data),
    class       = cls,
    post_prob   = round(prob, 4)
  )
  invisible(out)
}


#' Plot the clustering partition (density rug coloured by class) (Ex 32.6, 33.5)
#'
#' @param mod Mclust object
#'
#' @examples
#' sl_mclust_plot_partition(mod32)
sl_mclust_plot_partition <- function(mod) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  mclust::mclust2Dplot(data   = mod$data,
                        what   = "classification",
                        classification = mod$classification,
                        parameters = mod$parameters)
  # For univariate data, use mclustDensity plot
  if (is.null(dim(mod$data))) {
    plot(mod, what = "classification",
         main = "Mixture clustering partition")
  }
  invisible(mod)
}


#' Histogram overlaid with the estimated mixture density (Ex 32.7, 33.6)
#'
#' Works for univariate Mclust models. Computes the weighted sum of
#' component densities on a fine grid and overlays it on the histogram.
#'
#' @param mod   Mclust object (univariate)
#' @param label x-axis label
#' @param breaks histogram bins
#'
#' @examples
#' sl_mclust_density_plot(mod32, label = "Haemoglobin")
#' sl_mclust_density_plot(mod33, label = "Fish length (inches)")
sl_mclust_density_plot <- function(mod, label = "x", breaks = 30) {
  x     <- as.vector(mod$data)
  pars  <- mod$parameters
  props <- pars$pro
  means <- pars$mean
  vars  <- pars$variance$sigmasq
  if (is.null(vars)) vars <- rep(pars$variance$sigma2, length(props))

  xgrid <- seq(min(x), max(x), length.out = 500)
  dens  <- Reduce(`+`, lapply(seq_along(props), function(k) {
    v_k <- if (length(vars) == 1) vars else vars[k]
    props[k] * dnorm(xgrid, mean = means[k], sd = sqrt(v_k))
  }))

  hist(x, freq = FALSE, breaks = breaks,
       col = "dodgerblue", border = "white",
       main = paste("Estimated mixture density -", label),
       xlab = label)
  lines(xgrid, dens, col = "red", lwd = 2)

  # add component curves
  cols <- c("darkorange", "darkgreen", "purple", "brown", "darkblue")
  for (k in seq_along(props)) {
    v_k <- if (length(vars) == 1) vars else vars[k]
    lines(xgrid,
          props[k] * dnorm(xgrid, mean = means[k], sd = sqrt(v_k)),
          col = cols[(k - 1) %% length(cols) + 1],
          lwd = 1.5, lty = 2)
  }
  legend("topright",
         legend = c("Mixture", paste0("Component ", seq_along(props))),
         col    = c("red", cols[seq_along(props)]),
         lty    = c(1, rep(2, length(props))), lwd = 2, bty = "n")

  invisible(dens)
}


#' Predict class for new observations (Ex 33.7)
#'
#' @param mod      Mclust object
#' @param new_data numeric vector or matrix of new values
#'
#' @examples
#' sl_mclust_predict(mod33, new_data = 7)
sl_mclust_predict <- function(mod, new_data) {
  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Run sl_setup_ext() first.")
  pred <- predict(mod, newdata = new_data)
  out  <- data.frame(
    observation = new_data,
    class       = pred$classification,
    post_prob   = round(apply(pred$z, 1, max), 4)
  )
  cat("=== Out-of-sample predictions ===\n")
  print(out)
  invisible(out)
}


# =============================================================================
# 6. PROBABILITY UTILITIES  (Ex 25, 26, 28)
# =============================================================================

#' Binomial probability: P(X = k) or P(X <= k) for X ~ Bin(n, p) (Ex 25, 26)
#'
#' @param k     number of successes (scalar or vector)
#' @param n     number of trials
#' @param p     probability of success per trial
#' @param cumul if TRUE computes P(X <= k); default FALSE (exact probability)
#'
#' @examples
#' sl_binom_prob(k = 2, n = 10, p = 0.25)          # Ex 26
#' sl_binom_prob(k = 195, n = 200, p = 195/200, cumul = TRUE)  # Ex 25
sl_binom_prob <- function(k, n, p, cumul = FALSE) {
  fn <- if (cumul) pbinom else dbinom
  pr <- fn(k, size = n, prob = p)
  lbl <- if (cumul) sprintf("P(X <= %s)", paste(k, collapse = ","))
         else       sprintf("P(X = %s)",  paste(k, collapse = ","))
  cat(sprintf("X ~ Bin(n = %d, p = %.4g)\n", n, p))
  cat(sprintf("%s = %.7f\n", lbl, pr))
  invisible(pr)
}


#' Build a 2x2 contingency table from four cell counts (Ex 28)
#'
#' Cell layout (epidemiological convention):
#'   rows    = Exposure (yes = row 1, no = row 2)
#'   columns = Outcome  (yes = col 1, no = col 2)
#'
#' @param a exposed & outcome present
#' @param b exposed & outcome absent
#' @param c unexposed & outcome present
#' @param d unexposed & outcome absent
#' @param row_labels labels for the two rows
#' @param col_labels labels for the two columns
#'
#' @examples
#' tbl28 <- sl_contingency_table(a=115, b=4677, c=53, d=4739,
#'                                row_labels = c("IVF yes", "IVF no"),
#'                                col_labels  = c("Defect yes", "Defect no"))
sl_contingency_table <- function(a, b, c, d,
                                  row_labels = c("Exposed", "Unexposed"),
                                  col_labels  = c("Outcome yes", "Outcome no")) {
  m <- matrix(c(a, c, b, d), nrow = 2,
              dimnames = list(row_labels, col_labels))
  cat("=== 2x2 Contingency table ===\n")
  print(m)
  cat(sprintf("Row totals : %d, %d\n", a + b, c + d))
  cat(sprintf("Col totals : %d, %d\n", a + c, b + d))
  cat(sprintf("Grand total: %d\n",     a + b + c + d))
  invisible(m)
}


#' Odds, odds ratio, and relative risk from a 2x2 table (Ex 28)
#'
#' @param tbl 2x2 matrix (rows = exposure, cols = outcome) from
#'            sl_contingency_table() or table()
#'
#' @examples
#' tbl28 <- sl_contingency_table(115, 4677, 53, 4739)
#' sl_odds_ratio(tbl28)
sl_odds_ratio <- function(tbl) {
  stopifnot(all(dim(tbl) == c(2, 2)))
  a <- tbl[1, 1]; b <- tbl[1, 2]
  c <- tbl[2, 1]; d <- tbl[2, 2]

  odds_exp   <- a / b
  odds_unexp <- c / d
  OR         <- (a * d) / (b * c)
  RR         <- (a / (a + b)) / (c / (c + d))

  cat("=== Odds and Odds Ratio ===\n")
  cat(sprintf("Odds (exposed)  : %.4f  [P(event|exp) / P(no event|exp)]\n",
              odds_exp))
  cat(sprintf("Odds (unexposed): %.4f\n", odds_unexp))
  cat(sprintf("Odds Ratio (OR) : %.4f  (>1 = higher risk for exposed)\n", OR))
  cat(sprintf("Relative Risk   : %.4f\n", RR))

  invisible(list(odds_exposed   = odds_exp,
                 odds_unexposed = odds_unexp,
                 OR             = OR,
                 RR             = RR))
}


# =============================================================================
# SUMMARY — printed when this file is sourced
# =============================================================================
cat("\n")
cat("statlab_ext.R loaded. New functions available:\n")
cat("  Setup .............. sl_setup_ext\n")
cat("  Bootstrap plot ..... sl_boot_plot\n")
cat("  LM extras .......... sl_lm_diagnostics, sl_vif, sl_step_aic\n")
cat("                       sl_lm_confint, sl_predict_lm\n")
cat("                       sl_bootstrap_lm_ci\n")
cat("                       sl_residuals_vs_covariates\n")
cat("                       sl_lm_parallel_lines\n")
cat("  Logistic regression  sl_logistic, sl_logistic_equation\n")
cat("                       sl_step_aic_logistic, sl_logistic_or\n")
cat("                       sl_logistic_predict_prob\n")
cat("  Classification ..... sl_confusion_matrix, sl_classification_metrics\n")
cat("                       sl_roc_plot\n")
cat("  Mixture models ..... sl_mclust_select, sl_mclust_fit, sl_mclust_params\n")
cat("                       sl_mclust_classify, sl_mclust_plot_partition\n")
cat("                       sl_mclust_density_plot, sl_mclust_predict\n")
cat("  Probability ........ sl_binom_prob, sl_contingency_table, sl_odds_ratio\n")
cat("\n")
cat("NOTE: sl_boot_plot, sl_lm_diagnostics, sl_vif, sl_step_aic, sl_predict_lm\n")
cat("may already be defined in the last section of statlab.R.\n")
cat("Remove duplicates here if R raises a warning.\n")
