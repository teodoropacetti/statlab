# =============================================================================
#  statlab.R  -  Utility library for Statistical Modelling
#
#  Usage:   source("statlab.R")
#           sl_setup()           # install / load all dependencies
#
#  NOTE ON SEEDS: seeds are NOT managed inside these functions.
#                 Always call set.seed() in your script before generating
#                 random numbers, so your code is explicit and reproducible.
#                 Example:
#                   set.seed(625)
#                   X <- sl_bvnorm(n = 3000, mu1 = 0, mu2 = 0, s1 = 3,
#                                  s2 = 3, rho = 0)
#
#  DEPENDENCIES: mvtnorm, scatterplot3d, corrplot, skimr, dplyr
# =============================================================================


# =============================================================================
#  0.  SETUP
# =============================================================================

#' Install (if needed) and load all required packages
sl_setup <- function() {
  pkgs  <- c("mvtnorm", "scatterplot3d", "corrplot", "skimr", "dplyr")
  new   <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(new)) install.packages(new, quiet = TRUE)
  invisible(lapply(pkgs, library, character.only = TRUE))
  cat("Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}


# =============================================================================
#  1.  DESCRIPTIVE STATISTICS
# =============================================================================

#' Full descriptive statistics using skimr (as in official solutions)
#'
#' @param df  data.frame
#'
#' @examples
#'   sl_describe(school)
sl_describe <- function(df) {
  if (!requireNamespace("skimr", quietly = TRUE))
    stop("Run sl_setup() first.")
  out <- skimr::skim_without_charts(df) |>
    dplyr::select(-dplyr::any_of("complete_rate"))
  print(out)
  invisible(out)
}

#' Coefficient of variation for each numeric column (as in Ex. 5 solutions)
#'
#' @param df  data.frame
#'
#' @examples
#'   sl_cv(school)
sl_cv <- function(df) {
  cv_fn  <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  out    <- round(apply(num_df, 2, cv_fn), 4)
  print(out)
  invisible(out)
}

#' Manual descriptive statistics (fallback when skimr is unavailable)
#'
#' @param x  numeric vector or data.frame
sl_summary <- function(x) {
  if (is.numeric(x)) x <- data.frame(value = x)
  x <- x[, sapply(x, is.numeric), drop = FALSE]
  do.call(rbind, lapply(names(x), function(col) {
    v <- x[[col]][!is.na(x[[col]])]
    data.frame(
      variable = col, n = length(v),
      mean = mean(v), sd = sd(v),
      cv   = sd(v) / abs(mean(v)),
      min  = min(v),  Q1 = quantile(v, .25),
      median = median(v), Q3 = quantile(v, .75),
      max  = max(v),  IQR = IQR(v),
      row.names = NULL
    )
  }))
}

#' Descriptive statistics of model residuals (Ex. 14.4)
#'
#' Prints summary(), variance, and a brief normality indicator.
#'
#' @param mod  lm object
#'
#' @examples
#'   sl_residuals_describe(mod14)
sl_residuals_describe <- function(mod) {
  r <- residuals(mod)
  cat("=== Descriptive statistics of residuals ===\n")
  print(summary(r))
  cat(sprintf("\nVariance : %.6f\n", var(r)))
  cat(sprintf("Std. dev.: %.6f\n", sd(r)))
  cat(sprintf("Skewness : %.4f\n",  mean((r - mean(r))^3) / sd(r)^3))
  cat(sprintf("Kurtosis : %.4f  (excess)\n", mean((r - mean(r))^4) / sd(r)^4 - 3))
  invisible(r)
}

#' mu +/- k*sigma interval and count of observations outside it (Ex. 3.1)
#'
#' @param x  numeric vector
#' @param k  number of standard deviations (default 3)
#'
#' @examples
#'   sl_sigma_rule(x)       # checks 99.7% rule
#'   sl_sigma_rule(x, k=1)  # checks 68% rule
sl_sigma_rule <- function(x, k = 3) {
  x_min  <- mean(x) - k * sd(x)
  x_max  <- mean(x) + k * sd(x)
  n_low  <- length(x[x < x_min])
  n_high <- length(x[x > x_max])
  cat(sprintf("mu +/- %d*sigma = [%.4f, %.4f]\n", k, x_min, x_max))
  cat(sprintf("Obs. below lower bound: %d\n",  n_low))
  cat(sprintf("Obs. above upper bound: %d\n",  n_high))
  invisible(list(lower = x_min, upper = x_max,
                 n_low = n_low,  n_high = n_high))
}


# =============================================================================
#  2.  PROBABILITY: Normal and Student-t
# =============================================================================

#' P(a < X < b) for X ~ N(mu, sigma^2)
#'
#' IMPORTANT: R's pnorm() takes the standard deviation (sd), not the variance.
#'   Notation N(mu, sigma^2) -> pnorm(q, mean = mu, sd = sqrt(sigma^2))
#'   The function takes sigma (sd) directly, so N(5, 4) -> sigma = sqrt(4) = 2.
#'
#' @param a      lower bound (default -Inf)
#' @param b      upper bound (default +Inf)
#' @param mu     mean
#' @param sigma  standard deviation (NOT variance)
#'
#' @examples
#'   sl_prob_normal(b = 7,  mu = 5, sigma = 2)      # P(X < 7),    X ~ N(5, 4)
#'   sl_prob_normal(a = 3, b = 6, mu = 5, sigma = 2) # P(3 < X < 6)
#'   sl_prob_normal(a = -1, b = 1)                   # P(-1 < X < 1), X ~ N(0,1)
sl_prob_normal <- function(a = -Inf, b = Inf, mu = 0, sigma = 1) {
  p <- pnorm(b, mean = mu, sd = sigma) - pnorm(a, mean = mu, sd = sigma)
  cat(sprintf("X ~ N(mu = %.4g, sigma^2 = %.4g)\n", mu, sigma^2))
  cat(sprintf("P(%.4g < X < %.4g) = %.7f\n", a, b, p))
  invisible(p)
}

#' P(a < X < b) for X ~ T_nu (Student's t with nu degrees of freedom)
#'
#' @param a   lower bound
#' @param b   upper bound
#' @param nu  degrees of freedom
#'
#' @examples
#'   sl_prob_t(a = -1, b = 1, nu = 1)
#'   sl_prob_t(a = -1, b = 1, nu = 10)
sl_prob_t <- function(a = -Inf, b = Inf, nu = 1) {
  p <- pt(b, df = nu) - pt(a, df = nu)
  cat(sprintf("X ~ T(nu = %d)\n", as.integer(nu)))
  cat(sprintf("P(%.4g < X < %.4g) = %.7f\n", a, b, p))
  invisible(p)
}

#' Print the 68-95-99.7 rules for N(0,1) (Ex. 6)
sl_normal_rules <- function() {
  cat("Empirical rules for X ~ N(0, 1):\n")
  for (k in 1:3) {
    p <- pnorm(k) - pnorm(-k)
    cat(sprintf("  P(-%d < X < %d) = %.7f  (%.4f%%)\n", k, k, p, p * 100))
  }
}


# =============================================================================
#  3.  UNIVARIATE PLOTS
# =============================================================================

#' Histogram with overlaid density curve
#'
#' @param x       numeric vector
#' @param title   plot title
#' @param breaks  number of bins
#' @param col     bar colour (default "dodgerblue")
#'
#' @examples
#'   set.seed(130); x <- rnorm(1000, 7, sqrt(2))
#'   sl_hist_density(x, title = "Histogram of the realizations from N(7, 2)")
sl_hist_density <- function(x, title = "Histogram", breaks = 16,
                            col = "dodgerblue") {
  hist(x, breaks = breaks, freq = FALSE,
       main = title, xlab = deparse(substitute(x)), col = col)
  lines(density(x), col = "darkorange", lwd = 2)
}

#' Side-by-side histograms for two samples (Ex. 3.2)
#'
#' Mirrors the par(mfrow = c(1, 2)) + hist() pattern from the solutions.
#' Use xlim to restrict the axis when one sample has extreme outliers.
#'
#' @param x1, x2   numeric vectors
#' @param title1, title2  plot titles
#' @param xlim1, xlim2    optional x-axis limits for each panel
#'
#' @examples
#'   set.seed(130)
#'   t1 <- rt(1000, df = 1); t2 <- rt(1000, df = 10)
#'   sl_hist_compare(t1, t2, title1 = "v = 1", title2 = "v = 10",
#'                   xlim1 = c(-100, 100))
sl_hist_compare <- function(x1, x2, title1 = "Series 1", title2 = "Series 2",
                            breaks1 = 500, breaks2 = 30,
                            xlim1 = NULL, xlim2 = NULL) {
  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))

  args1 <- list(x1, main = title1, breaks = breaks1, freq = FALSE)
  if (!is.null(xlim1)) args1$xlim <- xlim1
  do.call(hist, args1)

  args2 <- list(x2, main = title2, breaks = breaks2, freq = FALSE)
  if (!is.null(xlim2)) args2$xlim <- xlim2
  do.call(hist, args2)
}

#' ECDF with optional theoretical curve overlay (Ex. 3.1 / 5.2 style)
#'
#' @param x      numeric vector
#' @param title  plot title
#' @param dist   "norm" to overlay N(mean(x), sd(x)) CDF;
#'               "t" to overlay T(df_t) CDF; NULL for none
#' @param df_t   degrees of freedom (required when dist = "t")
#' @param xlim   optional x-axis limits (useful for T(1) with extreme tails)
#' @param col    ECDF line colour
#'
#' @examples
#'   sl_ecdf(x, dist = "norm", title = "Empirical vs Theoretical CDFs")
#'   sl_ecdf(t1, dist = "t", df_t = 1,  xlim = c(-100, 100), title = "v = 1")
#'   sl_ecdf(t2, dist = "t", df_t = 10, title = "v = 10")
sl_ecdf <- function(x, title = "ECDF", dist = NULL, df_t = NULL,
                    xlim = NULL, col = "dodgerblue") {
  args <- list(ecdf(x), col = col, lwd = 3, do.points = FALSE,
               main = title, xlab = deparse(substitute(x)))
  if (!is.null(xlim)) args$xlim <- xlim
  do.call(plot, args)
  
  if (identical(dist, "norm")) {
    m <- mean(x); s <- sd(x)
    curve(pnorm(x, mean = m, sd = s), col = "red", lwd = 1, add = TRUE)
  } else if (identical(dist, "t")) {
    if (is.null(df_t)) stop("Specify df_t when using dist = 't'")
    curve(pt(x, df = df_t), col = "red", lwd = 1, add = TRUE)
  }
}

#' Compare ECDFs of two samples + N(0,1) reference line (Ex. 3.2)
#'
#' @param x1, x2         numeric vectors
#' @param label1, label2 legend labels
#' @param xlim           optional x-axis limits
#'
#' @examples
#'   sl_ecdf_compare(t1, t2, xlim = c(-100, 100),
#'                   label1 = "Tv, v = 1", label2 = "Tv, v = 10")
sl_ecdf_compare <- function(x1, x2, xlim = NULL,
                            title  = "Comparison of the ECDFs",
                            label1 = "Series 1", label2 = "Series 2") {
  args <- list(ecdf(x1), do.points = FALSE, col = "blue", lwd = 2, main = title)
  if (!is.null(xlim)) args$xlim <- xlim
  do.call(plot, args)
  plot(ecdf(x2), do.points = FALSE, col = "orange", lwd = 2, add = TRUE)
  curve(pnorm(x), lty = 2, lwd = 1, add = TRUE)
  legend("topleft",
    legend = c(label1, label2, "N(0, 1)"),
    lty    = c(1, 1, 2), col = c("blue", "orange", "black"), bty = "n")
}

#' Q-Q plot against the normal distribution (style of official solutions)
#'
#' Uses qqnorm() + qqline() exactly as in the solutions.
#'
#' @param x      numeric vector
#' @param title  plot title
#' @param col    point colour
#'
#' @examples
#'   sl_qqplot(x, title = "QQ-plot")
sl_qqplot <- function(x, title = "QQ-plot", col = "dodgerblue") {
  qqnorm(x, col = col, main = title)
  qqline(x, col = "black", lwd = 1.5)
}

#' ECDF + QQ side by side for a single sample (Ex. 3.1 panel)
#'
#' @param x      numeric vector
#' @param mu     theoretical mean for the normal overlay (default = mean(x))
#' @param sigma  theoretical sd  for the normal overlay (default = sd(x))
#'
#' @examples
#'   set.seed(130); x <- rnorm(1000, 7, sqrt(2))
#'   sl_panel_ecdf_qq(x, mu = 7, sigma = sqrt(2))
sl_panel_ecdf_qq <- function(x, mu = NULL, sigma = NULL,
                             title    = "Empirical vs Theoretical CDFs",
                             title_qq = "QQ-plot",
                             col      = "dodgerblue") {
  if (is.null(mu))    mu    <- mean(x)
  if (is.null(sigma)) sigma <- sd(x)

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))

  plot(ecdf(x), col = col, lwd = 3, do.points = FALSE, main = title)
  curve(pnorm(x, mean = mu, sd = sigma), col = "red", lwd = 1, add = TRUE)

  sl_qqplot(x, title = title_qq, col = col)
}

#' Side-by-side QQ plots for two t-distributed samples (Ex. 3.2)
#'
#' @examples
#'   sl_qqplot_t_compare(t1, t2)
sl_qqplot_t_compare <- function(t1, t2,
                                label1 = "QQ-plot for Tv, v = 1",
                                label2 = "QQ-plot for Tv, v = 10") {
  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))
  qqnorm(t1, col = "blue",   main = label1); qqline(t1)
  qqnorm(t2, col = "orange", main = label2); qqline(t2)
}

#' ECDF panel for all numeric columns of a data.frame (Ex. 5.2)
#'
#' Loops over columns and adds the fitted normal CDF in blue,
#' exactly as in the official solutions.
#'
#' @param df    data.frame
#' @param nrow  number of rows in the plot grid
#' @param ncol  number of columns in the plot grid (auto if NULL)
#'
#' @examples
#'   sl_ecdf_panel(school, nrow = 2, ncol = 3)
sl_ecdf_panel <- function(df, nrow = 2, ncol = NULL) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  n      <- ncol(num_df)
  if (is.null(ncol)) ncol <- ceiling(n / nrow)

  old_par <- par(mfrow = c(nrow, ncol))
  on.exit(par(old_par))

  for (i in seq_len(n)) {
    v <- num_df[, i]
    plot(ecdf(v), do.points = FALSE, main = names(num_df)[i],
         xlab = "x", ylab = "Fn(x)")
    curve(pnorm(x, mean = mean(v), sd = sd(v)),
          col = "blue", lwd = 1, add = TRUE)
  }
}

#' QQ-plot panel for all numeric columns of a data.frame (Ex. 5.3)
#'
#' @param df    data.frame
#' @param nrow  rows in the grid
#' @param ncol  cols in the grid (auto if NULL)
#'
#' @examples
#'   sl_qqplot_panel(school, nrow = 2, ncol = 3)
sl_qqplot_panel <- function(df, nrow = 2, ncol = NULL) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  n      <- ncol(num_df)
  if (is.null(ncol)) ncol <- ceiling(n / nrow)

  old_par <- par(mfrow = c(nrow, ncol))
  on.exit(par(old_par))

  for (i in seq_len(n)) {
    qqnorm(num_df[, i], main = names(num_df)[i], col = "black")
    qqline(num_df[, i], col = "darkorange", lwd = 1.5)
  }
}

#' Univariate plot panel for all numeric columns: boxplot + hist + density (Ex. 4.2)
#'
#' For each numeric column produces a row of three plots:
#'   boxplot  |  histogram with density curve  |  density plot
#'
#' @param df   data.frame
#'
#' @examples
#'   sl_univariate_panel_df(profits)
sl_univariate_panel_df <- function(df) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  n      <- ncol(num_df)

  old_par <- par(mfrow = c(n, 3), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  for (i in seq_len(n)) {
    v   <- num_df[, i]
    nm  <- names(num_df)[i]

    # boxplot
    boxplot(v, main = paste(nm, "- Boxplot"),
            col = "dodgerblue", horizontal = TRUE)

    # histogram + density curve
    hist(v, freq = FALSE, main = paste(nm, "- Histogram"),
         xlab = nm, col = "dodgerblue", border = "white")
    lines(density(v), col = "darkorange", lwd = 2)

    # standalone density plot
    plot(density(v), main = paste(nm, "- Density"),
         xlab = nm, col = "darkorange", lwd = 2)
    rug(v, col = adjustcolor("black", 0.3))
  }
}

#' Single boxplot, optionally grouped (for Ex. 12.1 server comparison)
#'
#' @param x      numeric vector
#' @param group  optional grouping vector (factor or character)
#'
#' @examples
#'   sl_boxplot(c(A, B), group = rep(c("A","B"), each=8),
#'              title = "Processing times by server")
sl_boxplot <- function(x, group = NULL, title = "Boxplot", ylab = "") {
  if (is.null(group)) {
    boxplot(x, main = title, ylab = ylab,
            col = "dodgerblue", notch = TRUE)
  } else {
    boxplot(x ~ group, main = title, ylab = ylab,
            col = c("dodgerblue", "darkorange"), notch = FALSE)
  }
}


# =============================================================================
#  4.  BIVARIATE NORMAL DISTRIBUTION
# =============================================================================

#' Build a 2x2 covariance matrix from standard deviations and correlation
#'
#' Avoids writing matrix(c(...), nrow = 2) by hand every time.
#'
#' @param s1   standard deviation of X1
#' @param s2   standard deviation of X2
#' @param rho  correlation coefficient (-1 < rho < 1)
#'
#' @examples
#'   Sigma <- sl_sigma2(s1 = 3, s2 = 3, rho = 0)
#'   Sigma <- sl_sigma2(s1 = 10, s2 = 10, rho = -0.2)
sl_sigma2 <- function(s1, s2, rho) {
  matrix(c(s1^2,         rho * s1 * s2,
           rho * s1 * s2, s2^2), nrow = 2)
}

#' Generate n samples from a bivariate normal (uses mvtnorm::rmvnorm)
#'
#' Returns a matrix with columns X1 and X2 - access with X[, 1] and X[, 2],
#' exactly as in the official solutions.
#'
#' NOTE: call set.seed() in your script BEFORE this function.
#'
#' @param n    number of observations
#' @param mu1  mean of X1
#' @param mu2  mean of X2
#' @param s1   standard deviation of X1
#' @param s2   standard deviation of X2
#' @param rho  correlation
#'
#' @examples
#'   set.seed(625)
#'   X <- sl_bvnorm(n = 3000, mu1 = 0, mu2 = 0, s1 = 3, s2 = 3, rho = 0)
#'   summary(X)
#'   sd(X[, 1]); cov(X[, 1], X[, 2])
sl_bvnorm <- function(n = 1000, mu1 = 0, mu2 = 0, s1 = 1, s2 = 1, rho = 0) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Run sl_setup() first.")
  X <- mvtnorm::rmvnorm(n = n, mean = c(mu1, mu2),
                        sigma = sl_sigma2(s1, s2, rho))
  colnames(X) <- c("X1", "X2")
  X
}

#' Bivariate scatter plot (Ex. 7.2)
#'
#' @param X      matrix or data.frame with two columns
#' @param title  plot title
#' @param col    point colour
#' @param asp    aspect ratio (use asp = 1 for undistorted ellipses)
#'
#' @examples
#'   sl_scatter_bv(X, title = "Scatter plot")
sl_scatter_bv <- function(X, title = "Scatter plot",
                          col = "blue", asp = NULL) {
  args <- list(X[, 1], X[, 2], col = col, main = title,
               xlab = "First component (X)",
               ylab = "Second component (Y)")
  if (!is.null(asp)) args$asp <- asp
  do.call(plot, args)
}

#' 3D wireframe density plot for a bivariate normal (Ex. 7.4)
#'
#' Reproduces the scatterplot3d + double for-loop approach from the solutions.
#'
#' @param mu1, mu2    means
#' @param s1, s2      standard deviations
#' @param rho         correlation
#' @param n_grid      grid points per axis (default 50, as in solutions)
#' @param angle       viewing angle (default 70)
#' @param col         wireframe colour (default "darkorange")
#' @param xrange, yrange  axis ranges (auto = mean +/- 4*sd)
#'
#' @examples
#'   sl_bvnorm_3d(mu1 = 0, mu2 = 0, s1 = 3, s2 = 3, rho = 0)
sl_bvnorm_3d <- function(mu1 = 0, mu2 = 0, s1 = 1, s2 = 1, rho = 0,
                         n_grid = 50, angle = 70, col = "darkorange",
                         title  = "Scatter Plot 3D",
                         xrange = NULL, yrange = NULL) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE) ||
      !requireNamespace("mvtnorm",       quietly = TRUE))
    stop("Run sl_setup() first.")

  mu    <- c(mu1, mu2)
  sigma <- sl_sigma2(s1, s2, rho)
  r1    <- if (is.null(xrange)) c(mu1 - 4*s1, mu1 + 4*s1) else xrange
  r2    <- if (is.null(yrange)) c(mu2 - 4*s2, mu2 + 4*s2) else yrange

  z1      <- seq(r1[1], r1[2], length.out = n_grid)
  z2      <- seq(r2[1], r2[2], length.out = n_grid)
  griglia <- expand.grid(z1, z2)
  dens    <- matrix(mvtnorm::dmvnorm(griglia, mean = mu, sigma = sigma),
                    ncol = length(z1), byrow = TRUE)

  p3d <- scatterplot3d::scatterplot3d(
    z1, z2, seq(min(dens), max(dens), length.out = n_grid),
    xlim = r1, ylim = r2, type = "n", angle = angle, grid = FALSE,
    main = title,
    xlab = expression(X[1]), ylab = expression(X[2]), zlab = "f(X1, X2)")

  for (i in length(z1):1)
    p3d$points3d(rep(z1[i], length(z2)), z2, dens[i, ],
                 type = "l", col = col)
  for (i in length(z2):1)
    p3d$points3d(z1, rep(z2[i], length(z1)), dens[, i],
                 type = "l", col = col)

  invisible(dens)
}

#' Contour plot of the theoretical bivariate normal density (Ex. 7.5)
#'
#' @param asp  aspect ratio: use asp = 1 for undistorted ellipses (as in solutions)
#'
#' @examples
#'   sl_bvnorm_contour(mu1 = 0, mu2 = 0, s1 = 3, s2 = 3, rho = 0, asp = 1)
sl_bvnorm_contour <- function(mu1 = 0, mu2 = 0, s1 = 1, s2 = 1, rho = 0,
                              n_grid  = 300, nlevels = 10,
                              title   = "Contour plot",
                              xlab    = "First component (X)",
                              ylab    = "Second component (Y)",
                              col     = "blue", asp = 1,
                              xrange  = NULL, yrange = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Run sl_setup() first.")

  mu    <- c(mu1, mu2)
  sigma <- sl_sigma2(s1, s2, rho)
  r1    <- if (is.null(xrange)) c(mu1 - 4*s1, mu1 + 4*s1) else xrange
  r2    <- if (is.null(yrange)) c(mu2 - 4*s2, mu2 + 4*s2) else yrange

  z1      <- seq(r1[1], r1[2], length.out = n_grid)
  z2      <- seq(r2[1], r2[2], length.out = n_grid)
  griglia <- expand.grid(z1, z2)
  dens    <- matrix(mvtnorm::dmvnorm(griglia, mean = mu, sigma = sigma),
                    ncol = length(z1))

  contour(z1, z2, dens, col = col, main = title,
          xlab = xlab, ylab = ylab, xlim = r1, ylim = r2,
          asp = asp, nlevels = nlevels)
  invisible(dens)
}

#' Contour plot overlaid on observed data (Ex. 10.3)
#'
#' Computes the sample mean vector and covariance matrix, draws the
#' theoretical normal contour lines, then adds the data points.
#' Standard method to visually assess bivariate normality.
#'
#' @param X               matrix or 2-column data.frame
#' @param xrange, yrange  axis ranges (auto from data range if NULL)
#' @param col_contour     contour line colour
#' @param col_points      data point colour
#'
#' @examples
#'   sl_contour_data(f1sim,
#'                   xlab = "Speed (km/h)",
#'                   ylab = "Lateral acceleration (g)")
sl_contour_data <- function(X, title = "Contour plot",
                            xlab = "Variable 1", ylab = "Variable 2",
                            n_grid = 300, nlevels = 10,
                            col_contour = "blue", col_points = "orange",
                            xrange = NULL, yrange = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("Run sl_setup() first.")

  X     <- as.matrix(X[, 1:2])
  mu    <- colMeans(X)
  sigma <- cov(X)
  r1    <- if (is.null(xrange)) range(X[, 1]) else xrange
  r2    <- if (is.null(yrange)) range(X[, 2]) else yrange

  z1      <- seq(r1[1], r1[2], length.out = n_grid)
  z2      <- seq(r2[1], r2[2], length.out = n_grid)
  griglia <- expand.grid(z1, z2)
  dens    <- matrix(mvtnorm::dmvnorm(griglia, mean = mu, sigma = sigma),
                    ncol = length(z1))

  contour(z1, z2, dens, col = col_contour, main = title,
          xlab = xlab, ylab = ylab, nlevels = nlevels)
  points(X[, 1], X[, 2], col = col_points, pch = 16, cex = 0.8)
  invisible(list(mu = mu, sigma = sigma))
}

#' Spectral decomposition of a covariance matrix (Ex. 7.6)
#'
#' Prints eigenvalues, eigenvectors, and verifies the reconstruction
#' V %*% diag(lambda) %*% t(V) == Sigma.
#'
#' @param sigma  2x2 covariance matrix
#'
#' @examples
#'   sl_spectral(sl_sigma2(3, 3, 0))
#'   sl_spectral(sl_sigma2(10, 10, -0.2))
sl_spectral <- function(sigma) {
  eig <- eigen(sigma)
  cat("=== Eigenvalues ===\n");  print(eig$values)
  cat("=== Eigenvectors (columns) ===\n"); print(eig$vectors)
  cat("=== Reconstruction: V %*% diag(lambda) %*% t(V) ===\n")
  print(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
  invisible(eig)
}


# =============================================================================
#  5.  CORRELATIONS
# =============================================================================

#' Sample correlation matrix, rounded (Ex. 4.1 / 10.2)
#'
#' @param df      data.frame (only numeric columns used)
#' @param digits  decimal places
#'
#' @examples
#'   sl_corr(profits)
sl_corr <- function(df, digits = 2) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  r      <- cor(num_df, use = "complete.obs")
  print(round(r, digits))
  invisible(r)
}

#' Partial correlation matrix via precision matrix inversion (Ex. 4.3)
#'
#' @param df      data.frame
#' @param digits  decimal places
#'
#' @examples
#'   sl_partial_corr(profits)
sl_partial_corr <- function(df, digits = 4) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  R      <- cor(num_df, use = "complete.obs")
  P      <- solve(R)                       # precision matrix
  D      <- diag(1 / sqrt(diag(P)))
  pc     <- -D %*% P %*% D
  diag(pc) <- 1
  rownames(pc) <- colnames(pc) <- colnames(R)
  cat("Partial correlation matrix:\n")
  print(round(pc, digits))
  invisible(pc)
}

#' Raw and partial corrplots side by side (Ex. 4.4 / 15.5)
#'
#' @examples
#'   sl_corrplot_dual(profits)
sl_corrplot_dual <- function(df) {
  if (!requireNamespace("corrplot", quietly = TRUE))
    stop("Run sl_setup() first.")
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  R      <- cor(num_df, use = "complete.obs")
  pc     <- sl_partial_corr(invisible(num_df))

  old_par <- par(mfrow = c(1, 2))
  on.exit(par(old_par))
  corrplot::corrplot(R,  type = "upper",
                     main = "Raw correlations",    mar = c(0,0,2,0))
  corrplot::corrplot(pc, type = "upper",
                     main = "Partial correlations", mar = c(0,0,2,0))
}

#' Scatterplot matrix with correlation coefficients (Ex. 4.5 / 15.3)
#'
#' Upper panel: scatter + r value.
#' Diagonal:    histogram.
#' Lower panel: scatter.
#'
#' @examples
#'   sl_pairs(profits)
sl_pairs <- function(df, ...) {
  num_df <- df[, sapply(df, is.numeric), drop = FALSE]
  pairs(num_df,
    pch = 16, cex = 0.4, col = adjustcolor("dodgerblue", 0.5),
    upper.panel = function(x, y, ...) {
      points(x, y, pch = 16, cex = 0.4,
             col = adjustcolor("dodgerblue", 0.5))
      r <- cor(x, y, use = "complete.obs")
      legend("topleft", legend = paste("r =", round(r, 2)),
             bty = "n", cex = 0.85)
    },
    diag.panel = function(x, ...) {
      usr <- par("usr"); on.exit(par(usr = usr))
      par(usr = c(usr[1:2], 0, 1.5))
      h <- hist(x, plot = FALSE)
      y <- h$counts / max(h$counts)
      rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], y,
           col = "dodgerblue", border = "white")
    },
    ...)
}


# =============================================================================
#  6.  BOOTSTRAP
# =============================================================================

#' Non-parametric bootstrap for a single-sample estimator
#'
#' NOTE: call set.seed() in your script before this function.
#'
#' @param x    numeric vector of data
#' @param fun  estimator function: mean | median | sd | var | custom
#' @param B    number of bootstrap replications (default 1000)
#'
#' @return invisible list: $estimate, $se, $bias, $distribution
#'
#' @examples
#'   set.seed(42)
#'   boot11 <- sl_bootstrap(times, fun = mean, B = 1000)
#'   boot13 <- sl_bootstrap(calorie, fun = sd, B = 2000)
sl_bootstrap <- function(x, fun = mean, B = 1000) {
  theta_hat <- fun(x)
  boot_dist <- replicate(B, fun(sample(x, size = length(x), replace = TRUE)))
  se   <- sd(boot_dist)
  bias <- mean(boot_dist) - theta_hat
  cat(sprintf("Original estimate : %.6f\n", theta_hat))
  cat(sprintf("Bootstrap SE      : %.6f\n", se))
  cat(sprintf("Bootstrap bias    : %.6f\n", bias))
  invisible(list(estimate = theta_hat, se = se, bias = bias,
                 distribution = boot_dist))
}

#' Non-parametric bootstrap for the difference between two groups
#'
#' NOTE: call set.seed() in your script before this function.
#'
#' @param x1, x2  numeric vectors for the two groups
#' @param fun      estimator applied to each group (default mean)
#' @param B        number of replications
#'
#' @examples
#'   A <- c(176, 125, 152, 180, 159, 168, 160, 151)
#'   B <- c(164, 121, 137, 169, 144, 145, 156, 139)
#'   set.seed(42)
#'   boot12_mean <- sl_bootstrap2(A, B, fun = mean)
#'   set.seed(42)
#'   boot12_med  <- sl_bootstrap2(A, B, fun = median)
sl_bootstrap2 <- function(x1, x2, fun = mean, B = 1000) {
  theta_hat <- fun(x1) - fun(x2)
  boot_dist <- replicate(B, {
    fun(sample(x1, length(x1), replace = TRUE)) -
    fun(sample(x2, length(x2), replace = TRUE))
  })
  se   <- sd(boot_dist)
  bias <- mean(boot_dist) - theta_hat
  cat(sprintf("Original difference: %.6f\n", theta_hat))
  cat(sprintf("Bootstrap SE       : %.6f\n", se))
  cat(sprintf("Bootstrap bias     : %.6f\n", bias))
  invisible(list(estimate = theta_hat, se = se, bias = bias,
                 distribution = boot_dist))
}

#' Percentile confidence interval from a bootstrap result (Ex. 12.4 / 12.6 / 13.3)
#'
#' @param boot_result  output of sl_bootstrap() or sl_bootstrap2()
#' @param level        confidence level (default 0.95)
#'
#' @examples
#'   sl_boot_ci(boot12_mean, level = 0.95)
#'   sl_boot_ci(boot12_med,  level = 0.90)
sl_boot_ci <- function(boot_result, level = 0.95) {
  alpha <- 1 - level
  ci    <- quantile(boot_result$distribution,
                    probs = c(alpha / 2, 1 - alpha / 2))
  cat(sprintf("%.0f%% percentile CI: [%.4f, %.4f]\n",
              level * 100, ci[1], ci[2]))
  cat(sprintf("CI length: %.4f\n", diff(ci)))
  invisible(ci)
}

#' Bootstrap standard error for the sample correlation coefficient (Ex. 10.2)
#'
#' Wrapper around sl_bootstrap_gen() with paired = TRUE, since (x, y) pairs
#' must be resampled jointly by row to preserve the bivariate structure.
#'
#' @param x  numeric vector (first variable)
#' @param y  numeric vector (second variable, same length as x)
#' @param B  number of bootstrap replications (default 1000)
#'
#' @return invisible list: $estimate, $se, $bias, $distribution
#'         Compatible with sl_boot_ci() and sl_boot_plot().
#'
#' @examples
#'   set.seed(42)
#'   boot_cor <- sl_bootstrap_cor(f1sim$speed, f1sim$lateralAcc, B = 1000)
#'   sl_boot_ci(boot_cor, level = 0.95)
#'   sl_boot_plot(boot_cor, title = "Bootstrap - Correlation coefficient")
sl_bootstrap_cor <- function(x, y, B = 1000) {
  if (length(x) != length(y))
    stop("x and y must have the same length.")
  
  sl_bootstrap_gen(
    data      = list(x = x, y = y),
    statistic = function(d) cor(d$x, d$y),
    B         = B,
    paired    = TRUE
  )
}

#' Generalised non-parametric bootstrap
#'
#' A single function that subsumes sl_bootstrap(), sl_bootstrap2(), and
#' sl_bootstrap_cor(). Data are passed as a named list; the statistic is
#' any function that accepts that list and returns a single number.
#'
#' Resampling strategies (controlled by 'paired'):
#'   paired = FALSE  each element of 'data' is resampled independently.
#'                   Use for one-sample estimators and two-group differences.
#'   paired = TRUE   one shared index vector is drawn and applied to every
#'                   element of 'data'. Use whenever the pairing between
#'                   vectors must be preserved (e.g. correlation, regression
#'                   coefficients, any row-wise statistic).
#'
#' @param data       named list of numeric vectors, e.g.
#'                     list(x = times)
#'                     list(a = serverA, b = serverB)
#'                     list(x = speed,   y = acc)
#' @param statistic  function(d) -> scalar, where d has the same structure
#'                   as 'data'. Use d$name or d[[1]], d[[2]], etc.
#' @param B          number of replications (default 1000)
#' @param paired     logical: resample jointly by row? (default FALSE)
#'
#' @return invisible list: $estimate, $se, $bias, $distribution
#'         Compatible with sl_boot_ci() and sl_boot_plot().
#'
#' @examples
#'   # --- one sample: mean ---
#'   set.seed(42)
#'   sl_bootstrap_gen(list(x = times),
#'                    statistic = function(d) mean(d$x))
#'
#'   # --- one sample: sd ---
#'   set.seed(42)
#'   sl_bootstrap_gen(list(x = calorie),
#'                    statistic = function(d) sd(d$x), B = 2000)
#'
#'   # --- two independent groups: difference of means ---
#'   set.seed(42)
#'   sl_bootstrap_gen(list(a = serverA, b = serverB),
#'                    statistic = function(d) mean(d$a) - mean(d$b))
#'
#'   # --- two independent groups: difference of medians ---
#'   set.seed(42)
#'   sl_bootstrap_gen(list(a = serverA, b = serverB),
#'                    statistic = function(d) median(d$a) - median(d$b))
#'
#'   # --- paired: correlation (must use paired = TRUE) ---
#'   set.seed(42)
#'   boot_cor <- sl_bootstrap_gen(
#'     list(x = f1sim$speed, y = f1sim$lateralAcc),
#'     statistic = function(d) cor(d$x, d$y),
#'     paired = TRUE
#'   )
#'   sl_boot_ci(boot_cor, level = 0.95)
#'   sl_boot_plot(boot_cor, title = "Bootstrap - Correlation")
#'
#'   # --- paired: custom statistic (ratio of means) ---
#'   set.seed(42)
#'   sl_bootstrap_gen(list(x = speed, y = acc),
#'                    statistic = function(d) mean(d$x) / mean(d$y),
#'                    paired = TRUE)
sl_bootstrap_gen <- function(data, statistic, B = 1000, paired = FALSE) {
  if (!is.list(data))
    stop("'data' must be a list of numeric vectors, e.g. list(x = myvar).")
  if (!is.function(statistic))
    stop("'statistic' must be a function(d) -> scalar.")
  
  # lengths of all elements
  lens <- vapply(data, length, integer(1))
  
  if (paired && length(unique(lens)) > 1)
    stop("When paired = TRUE all elements of 'data' must have the same length.")
  
  # point estimate on the original data
  theta_hat <- statistic(data)
  
  # bootstrap loop
  boot_dist <- replicate(B, {
    if (paired) {
      # draw one shared index, apply to every element
      idx      <- sample(lens[1], replace = TRUE)
      resampled <- lapply(data, function(v) v[idx])
    } else {
      # resample each element independently
      resampled <- lapply(data, function(v) sample(v, replace = TRUE))
    }
    statistic(resampled)
  })
  
  se   <- sd(boot_dist)
  bias <- mean(boot_dist) - theta_hat
  
  cat(sprintf("Original estimate : %.6f\n", theta_hat))
  cat(sprintf("Bootstrap SE      : %.6f\n", se))
  cat(sprintf("Bootstrap bias    : %.6f\n", bias))
  
  invisible(list(estimate     = theta_hat,
                 se           = se,
                 bias         = bias,
                 distribution = boot_dist))
}


# =============================================================================
#  7.  LINEAR REGRESSION
# =============================================================================

#' Fit a linear model and print a full summary (Ex. 14.2 / 15.6)
#'
#' @param formula  R formula, e.g. y ~ x1 + x2
#' @param data     data.frame
#' @return         invisible lm object
#'
#' @examples
#'   mod14 <- sl_lm(rating ~ profitability + capital + flexibility,
#'                  data = ratings)
sl_lm <- function(formula, data) {
  mod <- lm(formula, data = data)
  print(summary(mod))
  cat(sprintf("\nRSE   : %.4f  (df = %d)\n", sigma(mod), df.residual(mod)))
  cat(sprintf("R^2   : %.4f\n", summary(mod)$r.squared))
  cat(sprintf("R^2adj: %.4f\n", summary(mod)$adj.r.squared))
  invisible(mod)
}

#' Print the fitted model equation as a readable string (Ex. 14.3)
#'
#' @examples
#'   sl_lm_equation(mod14)
sl_lm_equation <- function(mod) {
  coefs   <- coef(mod)
  nomi    <- names(coefs)
  termini <- mapply(function(nm, val) {
    if (nm == "(Intercept)") sprintf("%.4f", val)
    else                     sprintf("%+.4f * %s", val, nm)
  }, nomi, coefs)
  eq <- paste(deparse(formula(mod)[[2L]]), "=",
              paste(termini, collapse = " "))
  cat("Fitted equation:\n  ", eq, "\n")
  invisible(eq)
}

#' Residual for a specific observation (Ex. 14.7)
#'
#' @param mod  lm object
#' @param id   row index (R integer index, not a variable ID column)
#'
#' @examples
#'   sl_residual(mod14, id = 29)
sl_residual <- function(mod, id) {
  r <- residuals(mod)[id]
  cat(sprintf("Residual for observation %d: %.6f\n", id, r))
  invisible(r)
}

#' Residual diagnostic plots + Shapiro-Wilk test (Ex. 14.8)
#'
#' Produces the four standard diagnostic plots (Residuals vs Fitted,
#' Normal QQ, Scale-Location, Cook's distance) and reports the
#' Shapiro-Wilk test for normality.
#'
#' @examples
#'   sl_residual_diag(mod14)
sl_residual_diag <- function(mod) {
  old_par <- par(mfrow = c(2, 2))
  on.exit(par(old_par))
  plot(mod, which = 1:4)

  sw <- shapiro.test(residuals(mod))
  cat("\n=== Shapiro-Wilk test on residuals ===\n")
  cat(sprintf("W = %.4f,  p-value = %.4f\n", sw$statistic, sw$p.value))
  if (sw$p.value > 0.05)
    cat("Do not reject H0: residuals are consistent with normality (alpha = 5%)\n")
  else
    cat("Reject H0: possible violation of normality in residuals\n")
  invisible(sw)
}

#' Prediction with confidence or prediction interval
#'
#' @param mod      lm object
#' @param newdata  data.frame with the new observation(s)
#' @param type     "prediction" (default) or "confidence"
#'
#' @examples
#'   sl_predict(mod14, newdata = ratings[ratings$id == 29, ])
sl_predict <- function(mod, newdata, type = "prediction") {
  pred <- predict(mod, newdata = newdata, interval = type)
  cat(sprintf("Prediction (%s interval):\n", type))
  print(round(pred, 4))
  invisible(pred)
}


# =============================================================================
#  8.  UTILITIES
# =============================================================================

#' Load an .Rdata file into the global environment
#'
#' Equivalent of pandas read_csv() but for R's native format.
#' Prints the names of all loaded objects.
#'
#' @examples
#'   sl_load("school.Rdata")   # then use head(school)
#'   sl_load("profits.Rdata")
sl_load <- function(path) {
  env  <- new.env()
  nomi <- load(path, envir = env)
  for (nm in nomi)
    assign(nm, get(nm, envir = env), envir = .GlobalEnv)
  cat("Objects loaded into .GlobalEnv:", paste(nomi, collapse = ", "), "\n")
  invisible(nomi)
}

#' Frequency table with counts and proportions
#'
#' @examples
#'   sl_freq(school$calworks > 20)
sl_freq <- function(x) {
  tab <- table(x)
  out <- data.frame(value = names(tab),
                    n     = as.integer(tab),
                    prop  = as.numeric(prop.table(tab)))
  print(out)
  invisible(out)
}


# =============================================================================
#  Banner
# =============================================================================
cat("\n")
cat("===========================================================================\n")
cat("  statlab.R  *  Statistical Modelling  *  A.Y. 2025-2026\n")
cat("===========================================================================\n")
cat("  SETUP         sl_setup()\n")
cat("\n")
cat("  DESCRIPTIVES  sl_describe(df)............................................# skimr output\n")
cat("                sl_cv(df)..................................................# coeff. of variation\n")
cat("                sl_summary(df).............................................# manual fallback\n")
cat("                sl_residuals_describe(mod).................................# Ex 14.4\n")
cat("                sl_sigma_rule(x, k)........................................# 68/95/99.7 check\n")
cat("\n")
cat("  PROBABILITY   sl_prob_normal(a, b, mu, sigma)\n")
cat("                sl_prob_t(a, b, nu)\n")
cat("                sl_normal_rules()\n")
cat("\n")
cat("  UNIVARIATE    sl_hist_density(x).........................................# single histogram\n")
cat("  PLOTS         sl_hist_compare(x1, x2)....................................# side-by-side  Ex 3.2\n")
cat("                sl_ecdf(x, dist)...........................................# ECDF + theory\n")
cat("                sl_ecdf_compare(x1, x2)....................................# two ECDFs     Ex 3.2\n")
cat("                sl_qqplot(x)...............................................# QQ normal\n")
cat("                sl_panel_ecdf_qq(x)........................................# ECDF+QQ panel Ex 3.1\n")
cat("                sl_qqplot_t_compare(t1,t2).................................# QQ compare    Ex 3.2\n")
cat("                sl_ecdf_panel(df)..........................................# all cols ECDF  Ex 5.2\n")
cat("                sl_qqplot_panel(df)........................................# all cols QQ    Ex 5.3\n")
cat("                sl_univariate_panel_df(df).................................# box+hist+dens  Ex 4.2\n")
cat("                sl_boxplot(x, group).......................................# single boxplot\n")
cat("\n")
cat("  BIVARIATE     sl_sigma2(s1, s2, rho).....................................# build Sigma\n")
cat("  NORMAL        sl_bvnorm(n,mu1,mu2,s1,s2,rho).............................# generate\n")
cat("                sl_scatter_bv(X)...........................................# scatter       Ex 7.2\n")
cat("                sl_bvnorm_3d(...)..........................................# wireframe 3D  Ex 7.4\n")
cat("                sl_bvnorm_contour(...).....................................# contour       Ex 7.5\n")
cat("                sl_contour_data(X).........................................# contour+data  Ex 10.3\n")
cat("                sl_spectral(sigma).........................................# eigen decomp  Ex 7.6\n")
cat("\n")
cat("  CORRELATIONS  sl_corr(df)................................................# raw\n")
cat("                sl_partial_corr(df)........................................# partial\n")
cat("                sl_corrplot_dual(df).......................................# both corrplots\n")
cat("                sl_pairs(df)...............................................# scatterplot matrix\n")
cat("\n")
cat("  BOOTSTRAP     sl_bootstrap(x, fun, B)....................................# single group\n")
cat("                sl_bootstrap2(x1,x2,fun,B).................................# two groups\n")
cat("                sl_boot_ci(res, level).....................................# percentile CI\n")
cat("                sl_boot_plot(res, level)...................................# histogram + CI\n")
cat("                sl_bootstrap_cor(x, y, B)..................................# bootstrap for standard error\n")
cat("                sl_bootstrap_gen(data, statistic, B, paired)...............# generalized bootstrap\n")
cat("\n")
cat("  REGRESSION    sl_lm(formula, data).......................................# fit + summary\n")
cat("                sl_lm_equation(mod)........................................# readable eq.\n")
cat("                sl_residual(mod, id).......................................# single residual\n")
cat("                sl_residuals_describe(mod).................................# residual stats\n")
cat("                sl_residual_diag(mod)......................................# 4 plots + SW\n")
cat("                sl_predict(mod, newdata)...................................# with interval\n")
cat("\n")
cat("  UTILITIES     sl_load(path)    sl_freq(x)\n")
cat("===========================================================================\n")
cat("  -> Start with:  sl_setup()\n\n")
