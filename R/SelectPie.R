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
  
  dFdy_a <- exp(-V) * A^(1 / m - 1) * ey
  fy     <- exp(-y) * exp(-exp(-y))
  out    <- fy - dFdy_a
  
  out[!is.finite(out)] <- 0
  out[out < .Machine$double.xmin] <- .Machine$double.xmin
  out
}


#' SelectPie
#'
#' Maximum likelihood estimation for compositional outcome models with
#' selection into the sample, with optional bootstrapped standard errors
#' and predicted log-ratios.
#'
#' Supports three estimators: full MLE under the bivariate extreme value
#' (Gumbel-logistic) copula (\code{distr = "ev"}), full MLE under bivariate
#' normality (\code{distr = "normal"}), and the Heckman two-step procedure
#' (\code{estimator = "heckman"}).
#'
#' @param data A \code{data.frame} containing all variables.
#' @param y1 Character string. Name of the binary selection indicator
#'   (1 = selected, 0 = not selected).
#' @param x1 Character vector of regressor names for the selection equation.
#'   An intercept is added automatically via \code{model.matrix}.
#' @param y2 Character string. Name of the continuous compositional outcome
#'   (observed only when \code{y1 == 1}).
#' @param x2 Character vector of regressor names for the outcome equation.
#'   An intercept is added automatically via \code{model.matrix}.
#' @param distr Character string. Distribution assumption for the MLE.
#'   One of \code{"ev"} (default, Gumbel-logistic copula) or
#'   \code{"normal"} (bivariate normal). Ignored when
#'   \code{estimator = "heckman"}.
#' @param estimator Character string. One of \code{"mle"} (default) for full
#'   maximum likelihood, or \code{"heckman"} for the two-step Heckman
#'   procedure.
#' @param method Optimization method passed to \code{\link[stats]{optim}}.
#'   Default is \code{"BFGS"}.
#' @param maxit Integer. Maximum number of iterations. Default is \code{1000}.
#' @param multistart_corr Logical. If \code{TRUE} (the default), the
#'   optimiser is restarted from two additional starting values.
#' @param B Integer. Number of bootstrap replications. Set to \code{0}
#'   (default) to return point estimates only.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{results}{A character matrix of point estimates (with significance
#'       stars when \code{B > 0}) and bootstrapped standard errors in
#'       parentheses (when \code{B > 0}).}
#'     \item{fit_stats}{A \code{data.frame} with \code{logLik}, \code{AIC},
#'       \code{BIC}, \code{n}, and \code{k}.}
#'     \item{predictions}{A \code{data.frame} with \code{log_ratio} (point
#'       estimates of \eqn{X_2 \hat\beta_2}) and, when \code{B > 0},
#'       \code{log_ratio_se} (bootstrap standard errors).}
#'   }
#'
#' @seealso \code{\link{int_evd_closed}}, \code{\link{cdfevd}},
#'   \code{\link{compose_shares}}, \code{\link{latex_table}}
#'
#' @export
SelectPie <- function(data, y1, x1, y2, x2,
                      distr           = c("ev", "normal"),
                      estimator       = c("mle", "heckman"),
                      method          = "BFGS",
                      maxit           = 1000,
                      multistart_corr = TRUE,
                      B               = 0) {
  
  distr     <- match.arg(distr)
  estimator <- match.arg(estimator)
  
  # ---- Internal estimation workhorse ----
  .fit <- function(data) {
    
    y1v <- data[[y1]]
    y2v <- data[[y2]]
    
    X1 <- stats::model.matrix(stats::reformulate(x1, response = NULL), data = data)
    X2 <- stats::model.matrix(stats::reformulate(x2, response = NULL), data = data)
    
    p1 <- ncol(X1)
    p2 <- ncol(X2)
    
    idx_sel <- (y1v == 1)
    
    # ---- Heckman two-step ----
    if (estimator == "heckman") {
      probit_fit <- stats::glm(
        stats::reformulate(x1, response = y1),
        family = stats::binomial(link = "probit"),
        data   = data
      )
      xb_probit     <- as.vector(X1 %*% stats::coef(probit_fit))
      imr           <- stats::dnorm(xb_probit) / stats::pnorm(xb_probit)
      data$.__imr__ <- y1v * imr
      
      outcome_fit <- stats::lm(
        stats::reformulate(c(x2, ".__imr__"), response = y2),
        data = data[idx_sel, ]
      )
      return(c(stats::coef(probit_fit), stats::coef(outcome_fit)))
    }
    
    # ---- MLE starting values ----
    OLSStage2 <- stats::lm.fit(X2[idx_sel, , drop = FALSE], y2v[idx_sel])
    b2 <- stats::coef(OLSStage2)
    b2[is.na(b2)] <- 0
    
    if (distr == "ev") {
      
      if (length(unique(y1v)) == 1L) {
        para0 <- c(rep(0, p1), b2, -2)
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
      
    } else if (distr == "normal") {
      
      if (length(unique(y1v)) == 1L) {
        para0 <- c(rep(0, p1), b2, rho0 = 0, logsigma0 = 0)
      } else {
        OLSStage1 <- stats::lm.fit(X1, y1v)
        b1 <- stats::coef(OLSStage1)
        b1[is.na(b1)] <- 0
        r <- suppressWarnings(
          stats::cor(OLSStage2$residuals, OLSStage1$residuals[idx_sel])
        )
        if (!is.finite(r)) r <- 0
        r         <- max(min(r, 0.95), -0.95)
        rho0      <- log((r + 1) / (1 - r))
        sigma0    <- max(stats::sd(OLSStage2$residuals, na.rm = TRUE), 1e-4)
        logsigma0 <- if (is.finite(log(sigma0))) log(sigma0) else 0
        para0 <- c(b1, b2, rho0, logsigma0)
      }
    }
    
    # ---- Objective function ----
    objFunc <- function(para) {
      beta1 <- para[seq_len(p1)]
      beta2 <- para[(p1 + 1L):(p1 + p2)]
      
      xb1 <- as.vector(-(X1 %*% beta1))
      xb2 <- as.vector( (X2 %*% beta2))
      u   <- y2v - xb2
      
      if (distr == "ev") {
        m     <- exp(para[length(para)]) + 1
        x1s   <- xb1 - digamma(1)
        logF1 <- -exp(-x1s)
        obj1  <- (1 - y1v) * logF1
        obj1[obj1 < -600] <- -600
        integ <- int_evd_closed(xb1, u, m)
        logI  <- log(integ)
        logI[logI < -600] <- -600
        obj2  <- y1v * logI
        
      } else if (distr == "normal") {
        rho   <- 2 / (1 + exp(-para[length(para) - 1L])) - 1
        sigma <- exp(para[length(para)])
        logF1 <- stats::pnorm(xb1, log.p = TRUE)
        obj1  <- (1 - y1v) * logF1
        obj1[obj1 < -600] <- -600
        z                     <- u / sigma
        log_density_u         <- stats::dnorm(z, log = TRUE) - log(sigma)
        log_selection_given_u <- stats::pnorm(
          (-xb1 + rho * z) / sqrt(1 - rho^2), log.p = TRUE
        )
        logI <- log_density_u + log_selection_given_u
        logI[logI < -600] <- -600
        obj2 <- y1v * logI
      }
      
      -sum(obj1 + obj2)
    }
    
    # ---- Optimisation ----
    res <- stats::optim(para0, objFunc,
                        method  = method,
                        control = list(maxit = maxit))
    
    if (isTRUE(multistart_corr)) {
      for (c0 in c(-2, 2)) {
        ptry <- para0
        if (distr == "ev") {
          ptry[length(ptry)] <- c0
        } else if (distr == "normal") {
          ptry[length(ptry) - 1L] <- c0
        }
        r2 <- stats::optim(ptry, objFunc,
                           method  = method,
                           control = list(maxit = maxit))
        if (is.finite(r2$value) && r2$value < res$value) res <- r2
      }
    }
    
    # ---- Back-transform parameters ----
    out <- res$par
    if (distr == "ev") {
      out[length(out)] <- 1 - 1 / (1 + exp(out[length(out)]))
    } else if (distr == "normal") {
      out[length(out) - 1L] <- 2 / (1 + exp(-out[length(out) - 1L])) - 1
      out <- out[-length(out)]
    }
    out
  }
  
  # ---- Internal: extract beta2 from parameter vector ----
  .get_beta2 <- function(par) {
    par[(length(par) - length(x2) - 1L):(length(par) - 1L)]
  }
  
  # ---- Internal: predicted log-ratios from original data ----
  .predict_lr <- function(par) {
    X2_full <- stats::model.matrix(
      stats::reformulate(x2, response = NULL), data = data
    )
    beta2 <- .get_beta2(par)
    as.vector(X2_full %*% beta2)
  }
  
  # ---- Internal fit statistics ----
  .fit_stats <- function(point_est) {
    
    n_obs    <- sum(!is.na(data[[y1]]))
    n_params <- length(point_est)
    
    if (estimator == "heckman") {
      return(data.frame(logLik = NA_real_, AIC = NA_real_,
                        BIC = NA_real_, n = n_obs, k = n_params))
    }
    
    k1 <- length(x1) + 1L
    k2 <- length(x2) + 1L
    
    beta1 <- point_est[seq_len(k1)]
    beta2 <- point_est[(k1 + 1L):(k1 + k2)]
    rho   <- point_est[k1 + k2 + 1L]
    
    X1 <- stats::model.matrix(stats::reformulate(x1, response = NULL), data = data)
    X2 <- stats::model.matrix(stats::reformulate(x2, response = NULL), data = data)
    
    y1v <- data[[y1]]
    y2v <- data[[y2]]
    
    xb1 <- as.vector(-(X1 %*% beta1))
    xb2 <- as.vector( (X2 %*% beta2))
    u   <- y2v - xb2
    
    ll           <- rep(NA_real_, nrow(data))
    observed     <- y1v == 1 & !is.na(y2v)
    not_observed <- y1v == 0
    
    if (distr == "ev") {
      m <- 1 / (1 - rho)
      x1s <- xb1 - digamma(1)
      ll[not_observed] <- -exp(-x1s[not_observed])
      integ <- int_evd_closed(xb1[observed], u[observed], m)
      ll[observed] <- log(pmax(integ, .Machine$double.xmin))
      
    } else if (distr == "normal") {
      ll[not_observed] <- stats::pnorm(xb1[not_observed], log.p = TRUE)
      z                     <- u[observed]
      log_density_u         <- stats::dnorm(z, log = TRUE)
      log_selection_given_u <- stats::pnorm(
        (-xb1[observed] + rho * z) / sqrt(1 - rho^2), log.p = TRUE
      )
      ll[observed] <- log_density_u + log_selection_given_u
    }
    
    logLik_val <- sum(ll, na.rm = TRUE)
    data.frame(
      logLik = round(logLik_val,                               3),
      AIC    = round(-2 * logLik_val + 2 * n_params,          3),
      BIC    = round(-2 * logLik_val + log(n_obs) * n_params, 3),
      n      = n_obs,
      k      = n_params
    )
  }
  
  # ---- Point estimates and fit statistics on full data ----
  point_est <- .fit(data)
  fit_stats <- .fit_stats(point_est)
  
  # ---- Point estimate log-ratios ----
  lr_point <- .predict_lr(point_est)
  
  # ---- Row names for results matrix ----
  last_term <- if (estimator == "heckman") "IMR" else "Corr."
  x1_full   <- as.vector(rbind(c("(Intercept1)", x1), ""))
  x2_full   <- as.vector(rbind(c("(Intercept2)", x2, last_term), ""))
  row_names <- c(x1_full, x2_full)
  k         <- length(x1) + length(x2) + 3
  
  # ---- Bootstrap ----
  if (B > 0) {
    n     <- nrow(data)
    index <- seq_len(n)
    
    boot_ests <- matrix(NA_real_, nrow = B, ncol = k)
    boot_lr   <- matrix(NA_real_, nrow = B, ncol = n)  # predictions per draw
    
    for (b in seq_len(B)) {
      boot_result <- tryCatch({
        sample_idx     <- sample(index, n, replace = TRUE)
        bootstrap_data <- data[sample_idx, ]
        par_b          <- .fit(bootstrap_data)
        
        list(
          coefs = par_b,
          lr    = .predict_lr(par_b)  # predicted on ORIGINAL data
        )
      }, error = function(e) list(coefs = rep(NA_real_, k),
                                  lr    = rep(NA_real_, n)))
      
      boot_ests[b, ] <- boot_result$coefs
      boot_lr[b, ]   <- boot_result$lr
    }
    
    boot_se    <- apply(boot_ests, 2, stats::sd, na.rm = TRUE)
    boot_lr_se <- apply(boot_lr,   2, stats::sd, na.rm = TRUE)
    
    # ---- Format results matrix ----
    results <- matrix(NA_character_, nrow = 2 * k, ncol = 1,
                      dimnames = list(row_names, y1))
    
    for (r in seq_len(k)) {
      est   <- point_est[r]
      se    <- boot_se[r]
      tstat <- abs(est / se)
      
      stars <- if (tstat > 2.575) "$^{***}$" else
        if (tstat > 1.96)  "$^{**}$"  else
          if (tstat > 1.645) "$^{*}$"   else ""
      
      results[2 * r - 1, 1] <- paste0(round(est, 3), stars)
      results[2 * r,     1] <- paste0("(", round(se, 3), ")")
    }
    
    predictions <- data.frame(
      log_ratio    = lr_point,
      log_ratio_se = boot_lr_se
    )
    
  } else {
    
    # No bootstrap — point estimates only
    results <- matrix(round(point_est, 3), nrow = k, ncol = 1,
                      dimnames = list(row_names[c(TRUE, FALSE)], y1))
    
    predictions <- data.frame(
      log_ratio = lr_point
    )
  }
  
  list(
    results     = results,
    fit_stats   = fit_stats,
    predictions = predictions
  )
}


#' latex_table
#'
#' Format a SelectPie results list or a matrix as a LaTeX \code{table}
#' environment. When passed a SelectPie results list, the function
#' automatically splits the table at the outcome equation and appends
#' fit statistics at the bottom.
#'
#' @param x A SelectPie results list (from \code{\link{SelectPie}}) or a
#'   matrix coercible via \code{as.matrix}.
#' @param caption Character string or \code{NULL}. Table caption.
#' @param label Character string or \code{NULL}. LaTeX label.
#' @param align Character string or \code{NULL}. Column alignment
#'   specification. Defaults to left-aligning the row-name column and
#'   centring all data columns.
#' @param booktabs Logical. If \code{TRUE} (the default), uses
#'   \code{\\toprule}, \code{\\midrule}, and \code{\\bottomrule}.
#'
#' @return Invisibly returns the character vector of LaTeX lines.
#'
#' @export
latex_table <- function(x,
                        caption  = NULL,
                        label    = NULL,
                        align    = NULL,
                        booktabs = TRUE) {
  
  # ---- Unpack SelectPie list ----
  if (is.list(x) && all(c("results", "fit_stats") %in% names(x))) {
    fit_stats <- x$fit_stats
    x         <- x$results
  } else {
    fit_stats <- NULL
  }
  
  x  <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  
  if (is.null(align)) {
    align <- paste0("l", paste(rep("c", nc), collapse = ""))
  }
  
  intercept2_row <- which(rownames(x) == "(Intercept2)")
  split_at       <- if (length(intercept2_row) == 1L) intercept2_row else NULL
  
  top_rule <- if (booktabs) "\\toprule"    else "\\hline"
  mid_rule <- if (booktabs) "\\midrule"    else "\\hline"
  bot_rule <- if (booktabs) "\\bottomrule" else "\\hline"
  
  n_cols <- nc + 1L
  make_section_header <- function(label) {
    paste0("\\multicolumn{", n_cols, "}{l}{\\textit{", label, "}} \\\\")
  }
  
  out <- c("\\begin{table}[!htbp]", "\\centering")
  if (!is.null(caption)) out <- c(out, paste0("\\caption{", caption, "}"))
  if (!is.null(label))   out <- c(out, paste0("\\label{",   label,   "}"))
  out <- c(out, paste0("\\begin{tabular}{", align, "}"))
  
  header <- paste(c("", colnames(x)), collapse = " & ")
  out <- c(out, top_rule, paste0(header, " \\\\"), mid_rule)
  
  for (i in seq_len(nr)) {
    if (!is.null(split_at) && i == 1L)
      out <- c(out, make_section_header("Selection Equation"))
    if (!is.null(split_at) && i == split_at)
      out <- c(out, mid_rule, make_section_header("Outcome Equation"))
    row_i <- paste(c(rownames(x)[i], x[i, ]), collapse = " & ")
    out   <- c(out, paste0(row_i, " \\\\"))
  }
  
  if (!is.null(fit_stats)) {
    fmt <- function(val) if (is.na(val)) "---" else as.character(val)
    fit_row <- paste(
      c("", paste0("logLik: ", fmt(fit_stats$logLik),
                   " $|$ AIC: ", fmt(fit_stats$AIC),
                   " $|$ BIC: ", fmt(fit_stats$BIC),
                   " $|$ $n$: ", fit_stats$n)),
      collapse = " & "
    )
    out <- c(out, mid_rule, make_section_header("Fit Statistics"),
             paste0(fit_row, " \\\\"))
  }
  
  out <- c(out, bot_rule, "\\end{tabular}", "\\end{table}")
  out <- gsub("_", " ", out)
  
  cat(paste(out, collapse = "\n"))
  invisible(out)
}


#' compose_shares
#'
#' Convert predicted log-ratios from multiple \code{\link{SelectPie}} and/or
#' \code{lm} objects into compositional vote shares using the additive
#' log-ratio (ALR) inverse transformation. Propagates uncertainty via the
#' delta method when standard errors are available.
#'
#' The reference category receives the residual share
#' \eqn{1 / (1 + \sum_j \exp(\ell_j))}, and each modelled category receives
#' \eqn{\exp(\ell_j) / (1 + \sum_j \exp(\ell_j))}.
#'
#' For \code{\link{SelectPie}} objects, prediction SEs come from the
#' bootstrap distribution stored in \code{$predictions$log_ratio_se}.
#' For \code{lm} objects, prediction SEs are computed analytically from the
#' variance-covariance matrix using the efficient formula
#' \eqn{\sqrt{\mathrm{rowSums}((X \Sigma) \circ X)}}, where \eqn{\Sigma} is
#' \code{vcov(model)}.
#'
#' Standard errors for each share are then propagated via the delta method:
#' \deqn{
#'   \mathrm{SE}(\hat{s}_j) \approx \sqrt{
#'     \hat{s}_j^2 (1 - \hat{s}_j)^2 \, \mathrm{SE}(\ell_j)^2 +
#'     \sum_{k \neq j} \hat{s}_j^2 \hat{s}_k^2 \, \mathrm{SE}(\ell_k)^2
#'   }
#' }
#'
#' @param models A named list where each element is either a
#'   \code{\link{SelectPie}} results list or a fitted \code{lm} object.
#'   Names identify the categories and become column names in the output.
#' @param newdata A \code{data.frame} of observations for which to compute
#'   predictions. Required when any element of \code{models} is an \code{lm}
#'   object; ignored for \code{\link{SelectPie}} objects which already carry
#'   their predictions internally.
#' @param reference Character string. Name for the reference category column.
#'   Default is \code{"reference"}.
#'
#' @return A \code{data.frame} with one share column per modelled category,
#'   one share column for the reference category, and (when SEs are available
#'   for all models) one \code{_se} column per category containing
#'   delta-method standard errors.
#'
#' @seealso \code{\link{SelectPie}}, \code{\link{simulate_shock}}
#'
#' @export
compose_shares <- function(models, newdata = NULL, reference = "reference") {
  
  if (!is.list(models) || is.null(names(models))) {
    stop("models must be a named list of SelectPie result objects or lm objects.")
  }
  
  # ---- Helper: extract log-ratio and SE from a single model ----
  .extract_lr <- function(model, nm) {
    
    # --- lm object ---
    if (inherits(model, "lm")) {
      if (is.null(newdata)) {
        stop("newdata must be supplied when models contains lm objects.")
      }
      X2  <- stats::model.matrix(model$terms, data = newdata)
      lr  <- as.vector(X2 %*% stats::coef(model))
      # Efficient diagonal of X2 %*% vcov %*% t(X2)
      V   <- stats::vcov(model)
      se  <- sqrt(rowSums((X2 %*% V) * X2))
      return(list(log_ratio = lr, log_ratio_se = se))
    }
    
    # --- SelectPie object ---
    if (all(c("results", "fit_stats", "predictions") %in% names(model))) {
      lr <- model$predictions$log_ratio
      se <- if ("log_ratio_se" %in% names(model$predictions))
        model$predictions$log_ratio_se
      else
        NULL
      return(list(log_ratio = lr, log_ratio_se = se))
    }
    
    stop("'", nm, "' must be either a SelectPie results list or an lm object.")
  }
  
  # ---- Extract from all models ----
  extracted <- mapply(.extract_lr, models, names(models),
                      SIMPLIFY = FALSE)
  
  cat_names <- names(models)
  n_obs     <- length(extracted[[1]]$log_ratio)
  
  # Check consistent number of observations
  lens <- vapply(extracted, function(e) length(e$log_ratio), integer(1))
  if (any(lens != n_obs)) {
    stop("All models must produce the same number of predicted values. ",
         "Check that newdata has the correct number of rows.")
  }
  
  # ---- Build log-ratio and SE matrices ----
  lr_matrix <- matrix(NA_real_, nrow = n_obs, ncol = length(cat_names),
                      dimnames = list(NULL, cat_names))
  for (nm in cat_names) lr_matrix[, nm] <- extracted[[nm]]$log_ratio
  
  has_se <- all(vapply(extracted,
                       function(e) !is.null(e$log_ratio_se),
                       logical(1)))
  
  se_matrix <- if (has_se) {
    m <- matrix(NA_real_, nrow = n_obs, ncol = length(cat_names),
                dimnames = list(NULL, cat_names))
    for (nm in cat_names) m[, nm] <- extracted[[nm]]$log_ratio_se
    m
  } else NULL
  
  # ---- ALR inverse transformation ----
  exp_lr       <- exp(lr_matrix)
  denom        <- 1 + rowSums(exp_lr)
  share_matrix <- exp_lr / denom
  ref_share    <- 1 / denom
  
  # ---- Assemble output ----
  out <- as.data.frame(share_matrix)
  colnames(out) <- cat_names
  out[[reference]] <- ref_share
  
  # ---- Delta method SEs ----
  if (has_se) {
    
    all_cats   <- c(cat_names, reference)
    all_shares <- cbind(share_matrix, ref_share)
    colnames(all_shares) <- all_cats
    
    for (j in all_cats) {
      sj    <- all_shares[, j]
      var_j <- rep(0, n_obs)
      
      for (nm in cat_names) {
        sk   <- all_shares[, nm]
        se_k <- se_matrix[, nm]
        
        if (j == nm) {
          # Own derivative: d(s_j)/d(lr_j) = s_j * (1 - s_j)
          var_j <- var_j + (sj * (1 - sj))^2 * se_k^2
        } else {
          # Cross derivative: d(s_j)/d(lr_k) = -s_j * s_k
          var_j <- var_j + (sj * sk)^2 * se_k^2
        }
      }
      
      out[[paste0(j, "_se")]] <- sqrt(var_j)
    }
  }
  
  out
}


#' simulate_shock
#'
#' Estimate the average marginal effect of a shock to a continuous variable
#' on compositional vote shares, with uncertainty quantified via a
#' nonparametric bootstrap that re-estimates all models on each draw.
#'
#' @param data A \code{data.frame} containing all variables.
#' @param models A named list of model specifications, one per modelled
#'   category. Each element must be a list with named elements:
#'   \code{y1}, \code{x1}, \code{y2}, and \code{x2}.
#' @param shock_var Character string. Name of the variable to shock.
#' @param shock_sd Numeric. Number of standard deviations to add. Default
#'   is \code{1}.
#' @param reference Character string. Name for the reference category.
#'   Default is \code{"reference"}.
#' @param B Integer. Number of bootstrap draws. Default is \code{1000}.
#' @param distr Character string. Passed to \code{\link{SelectPie}}.
#' @param estimator Character string. Passed to \code{\link{SelectPie}}.
#' @param method Character string. Passed to \code{\link{SelectPie}}.
#' @param maxit Integer. Passed to \code{\link{SelectPie}}.
#' @param multistart_corr Logical. Passed to \code{\link{SelectPie}}.
#' @param seed Integer or \code{NULL}. Random seed. Default is \code{NULL}.
#'
#' @return A \code{data.frame} with columns \code{category}, \code{mean},
#'   \code{sd}, \code{lb}, and \code{ub} summarising the bootstrap
#'   distribution of average differences (shocked minus baseline) for each
#'   category.
#'
#' @seealso \code{\link{SelectPie}}, \code{\link{compose_shares}}
#'
#' @export
simulate_shock <- function(data,
                           models,
                           shock_var,
                           shock_sd        = 1,
                           reference       = "reference",
                           B               = 1000,
                           distr           = "ev",
                           estimator       = "mle",
                           method          = "BFGS",
                           maxit           = 1000,
                           multistart_corr = TRUE,
                           seed            = NULL) {
  
  if (!is.list(models) || is.null(names(models)))
    stop("models must be a named list of model specifications.")
  
  for (nm in names(models)) {
    missing_fields <- setdiff(c("y1", "x1", "y2", "x2"), names(models[[nm]]))
    if (length(missing_fields) > 0)
      stop("Model '", nm, "' is missing: ", paste(missing_fields, collapse = ", "))
  }
  
  if (!shock_var %in% names(data))
    stop("shock_var '", shock_var, "' not found in data.")
  
  if (!is.null(seed)) set.seed(seed)
  
  cat_names  <- names(models)
  all_cats   <- c(cat_names, reference)
  n_obs      <- nrow(data)
  index      <- seq_len(n_obs)
  
  shock_size        <- shock_sd * stats::sd(data[[shock_var]], na.rm = TRUE)
  data_shocked      <- data
  data_shocked[[shock_var]] <- data_shocked[[shock_var]] + shock_size
  
  boot_diffs <- matrix(NA_real_, nrow = B, ncol = length(all_cats),
                       dimnames = list(NULL, all_cats))
  
  for (b in seq_len(B)) {
    
    if (b %% 100 == 0) message("Bootstrap draw ", b, " of ", B)
    
    result <- tryCatch({
      
      sample_idx     <- sample(index, n_obs, replace = TRUE)
      bootstrap_data <- data[sample_idx, ]
      
      # Re-estimate all models on bootstrap data
      fitted <- lapply(cat_names, function(nm) {
        spec <- models[[nm]]
        SelectPie(data            = bootstrap_data,
                  y1              = spec$y1,
                  x1              = spec$x1,
                  y2              = spec$y2,
                  x2              = spec$x2,
                  distr           = distr,
                  estimator       = estimator,
                  method          = method,
                  maxit           = maxit,
                  multistart_corr = multistart_corr,
                  B               = 0)
      })
      names(fitted) <- cat_names
      
      # Baseline shares on original data
      shares_base    <- compose_shares(fitted, reference = reference)
      
      # Shocked shares — temporarily replace predictions
      fitted_shocked <- lapply(cat_names, function(nm) {
        spec   <- models[[nm]]
        fit_nm <- fitted[[nm]]
        X2_s   <- stats::model.matrix(
          stats::reformulate(spec$x2, response = NULL), data = data_shocked
        )
        par      <- fit_nm$results[c(TRUE, FALSE), , drop = FALSE]
        par_num  <- as.numeric(gsub("\\$\\^\\{\\*+\\}\\$", "", par[, 1]))
        beta2    <- par_num[(length(par_num) - length(spec$x2) - 1L):
                              (length(par_num) - 1L)]
        fit_nm$predictions$log_ratio <- as.vector(X2_s %*% beta2)
        fit_nm
      })
      names(fitted_shocked) <- cat_names
      
      shares_shock <- compose_shares(fitted_shocked, reference = reference)
      
      # Mean difference across observations
      vapply(all_cats, function(cat) {
        mean(shares_shock[[cat]] - shares_base[[cat]], na.rm = TRUE)
      }, numeric(1))
      
    }, error = function(e) rep(NA_real_, length(all_cats)))
    
    boot_diffs[b, ] <- result
  }
  
  n_failed <- sum(rowSums(is.na(boot_diffs)) > 0)
  if (n_failed > 0) {
    message(n_failed, " bootstrap draw(s) failed and were excluded.")
    boot_diffs <- boot_diffs[rowSums(is.na(boot_diffs)) == 0, , drop = FALSE]
  }
  
  means <- apply(boot_diffs, 2, mean, na.rm = TRUE)
  sds   <- apply(boot_diffs, 2, stats::sd, na.rm = TRUE)
  
  data.frame(
    category = all_cats,
    mean     = round(means, 6),
    sd       = round(sds,   6),
    lb       = round(means - 1.96 * sds, 6),
    ub       = round(means + 1.96 * sds, 6),
    row.names = NULL
  )
}