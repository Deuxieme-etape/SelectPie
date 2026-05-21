#' XbetaFunc
#'
#' A helper function producing the predicted values of Y from X*hat_beta.
#'
#' @param data A numeric matrix or data frame of regressors.
#' @param beta A numeric vector of estimated coefficients.
#'
#' @return A numeric vector of predicted values (X \%*\% beta).
#'
#' @export
XbetaFunc <- function(data, beta) {
  as.numeric(data) %*% beta
}


#' cdfevd
#'
#' Marginal CDF helper function (Gumbel / extreme value distribution).
#' Kept for completeness; the main objective uses the log form directly.
#'
#' @param x Numeric vector. Values to be transformed via the standardized
#'   Gumbel CDF to the (0, 1) interval.
#'
#' @return Numeric vector of CDF values in (0, 1).
#'
#' @export
cdfevd <- function(x) {
  x <- x - digamma(1)
  exp(-exp(-x))
}


#' int_evd_closed
#'
#' Closed-form marginal integral for the bivariate Gumbel-logistic copula.
#'
#' Computes \eqn{\int_{a}^{\infty} f_{X,Y}(x, y)\, dx} in closed form for
#' the bivariate extreme value distribution with dependence parameter
#' \eqn{m \ge 1}.
#'
#' @param x Numeric vector. Lower limit of integration (the selection index).
#' @param y Numeric vector. Value of the outcome variable.
#' @param m Numeric scalar. Dependence parameter (\eqn{m \ge 1}); larger
#'   values imply stronger dependence.
#'
#' @return A strictly positive numeric vector (floored at
#'   \code{.Machine$double.xmin}) so that \code{log()} is safe to apply.
#'
#' @export
int_evd_closed <- function(x, y, m) {
  x <- x - digamma(1)
  y <- y - digamma(1)

  ex <- exp(-m * x)
  ey <- exp(-m * y)
  A  <- ex + ey
  V  <- A^(1 / m)

  # Conditional survival contribution at lower limit x = a
  dFdy_a <- exp(-V) * A^(1 / m - 1) * ey

  # Gumbel marginal pdf at x = Inf
  fy <- exp(-y) * exp(-exp(-y))

  out <- fy - dFdy_a

  # Numerical guards so log(out) is always finite
  out[!is.finite(out)] <- 0
  out[out < .Machine$double.xmin] <- .Machine$double.xmin
  out
}


#' SelectPie
#'
#' Maximum likelihood estimation for compositional outcome models with
#' selection into the sample, with optional bootstrapped standard errors.
#'
#' @param data A data.frame containing all variables.
#' @param y1 Character string. Name of the binary selection indicator.
#' @param x1 Character vector of regressor names for the selection equation.
#' @param y2 Character string. Name of the continuous compositional outcome.
#' @param x2 Character vector of regressor names for the outcome equation.
#' @param method Optimization method passed to optim. Default is "BFGS".
#' @param maxit Maximum number of iterations. Default is 1000.
#' @param multistart_corr Logical. Use multiple starting values for dependence parameter.
#' @param B Integer. Number of bootstrap replications. Set to 0 (default) to skip (only obtain point estimates).
#'
#' @return A matrix with point estimates and (if B > 0) bootstrapped standard errors.
#' @export
SelectPie <- function(data, y1, x1, y2, x2,
                      method          = "BFGS",
                      maxit           = 1000,
                      multistart_corr = TRUE,
                      B               = 0) {
  
  # ---- Internal estimation workhorse (no bootstrap) ----
  .fit <- function(data) {
    
    y1v <- data[[y1]]
    y2v <- data[[y2]]
    
    X1 <- stats::model.matrix(stats::reformulate(x1, response = NULL), data = data)
    X2 <- stats::model.matrix(stats::reformulate(x2, response = NULL), data = data)
    
    p1 <- ncol(X1)
    p2 <- ncol(X2)
    
    idx_sel <- (y1v == 1)
    
    OLSStage2 <- stats::lm.fit(X2[idx_sel, , drop = FALSE], y2v[idx_sel])
    b2 <- stats::coef(OLSStage2)
    b2[is.na(b2)] <- 0
    
    if (length(unique(y1v)) == 1L) {
      corr_initial <- -2
      para0 <- c(rep(0, p1), b2, corr_initial)
    } else {
      OLSStage1 <- stats::lm.fit(X1, y1v)
      b1 <- stats::coef(OLSStage1)
      b1[is.na(b1)] <- 0
      
      r <- suppressWarnings(
        stats::cor(OLSStage2$residuals, OLSStage1$residuals[idx_sel])
      )
      if (!is.finite(r)) r <- 0
      
      corr_initial <- log(1 / (1 - abs(r)) - 1)
      corr_initial <- max(min(corr_initial, 2), -2)
      
      para0 <- c(b1, b2, corr_initial)
    }
    
    objFunc <- function(para) {
      beta1 <- para[seq_len(p1)]
      beta2 <- para[(p1 + 1L):(p1 + p2)]
      m     <- exp(para[length(para)]) + 1
      
      xb1 <- as.vector(-(X1 %*% beta1))
      xb2 <- as.vector( (X2 %*% beta2))
      u   <- y2v - xb2
      
      x1s   <- xb1 - digamma(1)
      logF1 <- -exp(-x1s)
      obj1  <- (1 - y1v) * logF1
      obj1[obj1 < -600] <- -600
      
      integ <- int_evd_closed(xb1, u, m)
      logI  <- log(integ)
      logI[logI < -600] <- -600
      obj2  <- y1v * logI
      
      -sum(obj1 + obj2)
    }
    
    res <- stats::optim(para0, objFunc,
                        method  = method,
                        control = list(maxit = maxit))
    
    if (isTRUE(multistart_corr)) {
      for (c0 in c(-2, 2)) {
        ptry <- para0
        ptry[length(ptry)] <- c0
        r2 <- stats::optim(ptry, objFunc,
                           method  = method,
                           control = list(maxit = maxit))
        if (is.finite(r2$value) && r2$value < res$value) res <- r2
      }
    }
    
    out <- res$par
    out[length(out)] <- 1 - 1 / (1 + exp(out[length(out)]))
    out
  }
  
  # ---- Point estimates on full data ----
  point_est <- .fit(data)
  
  # ---- Internal fit statistics (EVD copula likelihood) ----
  .fit_stats <- function(point_est) {
    
    k1 <- length(x1) + 1
    k2 <- length(x2) + 1
    
    beta1 <- point_est[seq_len(k1)]
    beta2 <- point_est[(k1 + 1L):(k1 + k2)]
    rho   <- point_est[k1 + k2 + 1L]
    
    # Recover m from back-transformed rho
    m <- 1 / (1 - rho)
    
    X1 <- stats::model.matrix(stats::reformulate(x1, response = NULL), data = data)
    X2 <- stats::model.matrix(stats::reformulate(x2, response = NULL), data = data)
    
    y1v <- data[[y1]]
    y2v <- data[[y2]]
    
    xb1 <- as.vector(-(X1 %*% beta1))
    xb2 <- as.vector( (X2 %*% beta2))
    u   <- y2v - xb2
    
    ll <- rep(NA_real_, nrow(data))
    
    observed     <- y1v == 1 & !is.na(y2v)
    not_observed <- y1v == 0
    
    # Not selected: log F_1(xb1)
    x1s <- xb1 - digamma(1)
    ll[not_observed] <- -exp(-x1s[not_observed])
    
    # Selected: log integral of bivariate EVD
    integ <- int_evd_closed(xb1[observed], u[observed], m)
    ll[observed] <- log(pmax(integ, .Machine$double.xmin))
    
    logLik_val <- sum(ll, na.rm = TRUE)
    
    n_params <- length(point_est)
    n_obs    <- sum(!is.na(y1v))
    
    data.frame(
      logLik = round(logLik_val,          3),
      AIC    = round(-2 * logLik_val + 2 * n_params,        3),
      BIC    = round(-2 * logLik_val + log(n_obs) * n_params, 3),
      n      = n_obs,
      k      = n_params
    )
  }
  
  # ---- Fit Statistics ----
  fit_stats <- .fit_stats(point_est)
  
  # ---- Row names for output matrix ----
  x1_full <- as.vector(rbind(c("(Intercept1)", x1), ""))
  x2_full <- as.vector(rbind(c("(Intercept2)", x2, "Corr."), ""))
  row_names <- c(x1_full, x2_full)
  
  k <- length(x1) + length(x2) + 3 # plus 3 because of the two intercepts and Corr.
  
  # ---- Bootstrap (only if B > 0) ----
  if (B > 0) {
    n     <- nrow(data)
    index <- seq_len(n)
    
    boot_ests <- matrix(NA, nrow = B, ncol = k)
    
    for (b in seq_len(B)) {
      sample_index   <- sample(index, n, replace = TRUE)
      bootstrap_data <- data[sample_index, ]
      
      boot_ests[b, ] <- tryCatch(
        .fit(bootstrap_data),
        error = function(e) rep(NA_real_, k)
      )
    }
    
    boot_se <- apply(boot_ests, 2, stats::sd, na.rm = TRUE)
    
    # ---- Format estimates with significance stars and SEs in parentheses ----
    results <- matrix(NA_character_, nrow = 2 * k, ncol = 1,
                      dimnames = list(row_names, y1))
    
    for (r in seq_len(k)) {
      est  <- point_est[r]
      se   <- boot_se[r]
      tstat <- abs(est / se)
      
      # Stars based on t-statistic
      stars <- if (tstat > 2.575) {
        "$^{***}$"
      } else if (tstat > 1.96) {
        "$^{**}$"
      } else if (tstat > 1.645) {
        "$^{*}$"
      } else {
        ""
      }
      
      # Estimate row (odd): rounded value + stars
      results[2 * r - 1, 1] <- paste0(round(est, 3), stars)
      
      # SE row (even): rounded value in parentheses
      results[2 * r,     1] <- paste0("(", round(se, 3), ")")
    }
    
  } else {
    
    # No bootstrap — just return rounded point estimates (odd rows only)
    results <- matrix(round(point_est, 3), nrow = k, ncol = 1,
                      dimnames = list(row_names[c(TRUE, FALSE)], y1))
  }
  
  list(
    results = results, 
    fit_stats = fit_stats
  )
  
}


#' latex_table
#'
#' Format a SelectPie results list or a matrix as a LaTeX \code{table} environment.
#' When passed a SelectPie results list (with \code{results} and \code{fit_stats}),
#' the function automatically splits the table at the outcome equation and appends
#' fit statistics at the bottom.
#'
#' @param x A SelectPie results list (from \code{\link{SelectPie}}) or a matrix
#'   coercible via \code{as.matrix}.
#' @param caption Character string or \code{NULL}. Table caption.
#' @param label Character string or \code{NULL}. LaTeX label for
#'   \code{\\ref\{\}} cross-referencing.
#' @param align Character string or \code{NULL}. Column alignment
#'   specification, e.g. \code{"lcccc"}. Defaults to left-aligning the
#'   row-name column and centring all data columns.
#' @param booktabs Logical. If \code{TRUE} (the default), uses
#'   \code{\\toprule}, \code{\\midrule}, and \code{\\bottomrule} from
#'   the \pkg{booktabs} LaTeX package instead of \code{\\hline}.
#'
#' @return Invisibly returns the character vector of LaTeX lines; called
#'   primarily for the side effect of printing via \code{cat}.
#'
#' @export
latex_table <- function(x,
                        caption  = NULL,
                        label    = NULL,
                        align    = NULL,
                        booktabs = TRUE) {
  
  # ---- Unpack SelectPie list if supplied ----
  if (is.list(x) && all(c("results", "fit_stats") %in% names(x))) {
    fit_stats <- x$fit_stats
    x         <- x$results
  } else {
    fit_stats <- NULL
  }
  
  x  <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  
  # ---- Default alignment ----
  if (is.null(align)) {
    align <- paste0("l", paste(rep("c", nc), collapse = ""))
  }
  
  # ---- Detect split point from row names ----
  intercept2_row <- which(rownames(x) == "(Intercept2)")
  split_at <- if (length(intercept2_row) == 1L) intercept2_row else NULL
  
  # ---- Rules ----
  top_rule <- if (booktabs) "\\toprule"    else "\\hline"
  mid_rule <- if (booktabs) "\\midrule"    else "\\hline"
  bot_rule <- if (booktabs) "\\bottomrule" else "\\hline"
  
  # ---- Section header helper ----
  n_cols <- nc + 1L
  make_section_header <- function(label) {
    paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", label, "}} \\\\")
  }
  
  # ---- Build LaTeX lines ----
  out <- c("\\begin{table}[!htbp]", "\\centering")
  
  if (!is.null(caption)) {
    out <- c(out, paste0("\\caption{", caption, "}"))
  }
  if (!is.null(label)) {
    out <- c(out, paste0("\\label{", label, "}"))
  }
  
  out <- c(out, paste0("\\begin{tabular}{", align, "}"))
  
  # Header row
  header <- paste(c("", colnames(x)), collapse = " & ")
  out <- c(out, top_rule, paste0(header, " \\\\"), mid_rule)
  
  # ---- Data rows ----
  for (i in seq_len(nr)) {
    
    # Selection equation header before row 1
    if (!is.null(split_at) && i == 1L) {
      out <- c(out, make_section_header("Selection Equation"))
    }
    
    # Outcome equation header at split point
    if (!is.null(split_at) && i == split_at) {
      out <- c(out, mid_rule, make_section_header("Outcome Equation"))
    }
    
    row_i <- paste(c(rownames(x)[i], x[i, ]), collapse = " & ")
    out   <- c(out, paste0(row_i, " \\\\"))
  }
  
  # ---- Fit statistics ----
  if (!is.null(fit_stats)) {
    fit_row <- paste(
      c("",
        paste0("logLik: ", fit_stats$logLik,
               " $|$ AIC: ",    fit_stats$AIC,
               " $|$ BIC: ",    fit_stats$BIC,
               " $|$ $n$: ",    fit_stats$n)),
      collapse = " & "
    )
    out <- c(out,
             mid_rule,
             paste0("\\multicolumn{", n_cols, "}{l}{\\textit{Fit Statistics}} \\\\"),
             paste0(fit_row, " \\\\"))
  }
  
  out <- c(out, bot_rule, "\\end{tabular}", "\\end{table}")
  
  cat(paste(out, collapse = "\n"))
  invisible(out)
}
