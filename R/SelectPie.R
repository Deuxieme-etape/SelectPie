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
#' selection into the sample.
#'
#' Fits a two-equation system in which a binary selection equation
#' (equation 1) determines whether the continuous compositional outcome
#' (equation 2) is observed.  Joint dependence between the two equations
#' is modelled via a bivariate Gumbel-logistic (extreme value) copula with
#' dependence parameter \eqn{m \ge 1}, reparameterised as
#' \eqn{m = \exp(\gamma) + 1} so that optimisation is unconstrained.
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
#' @param method Character string passed to \code{\link[stats]{optim}}.
#'   Default is \code{"BFGS"}.
#' @param maxit Integer. Maximum number of iterations passed to
#'   \code{\link[stats]{optim}}. Default is \code{1000}.
#' @param multistart_corr Logical. If \code{TRUE} (the default), the
#'   optimiser is restarted from two additional starting values for the
#'   dependence parameter (\eqn{\gamma \in \{-2, 2\}}) and the best
#'   solution is retained.
#'
#' @return A named numeric vector containing:
#'   \describe{
#'     \item{Selection coefficients}{Entries \code{1:p1} correspond to
#'       \code{x1} (including intercept).}
#'     \item{Outcome coefficients}{Entries \code{(p1+1):(p1+p2)} correspond
#'       to \code{x2} (including intercept).}
#'     \item{rho}{The final entry is the back-transformed dependence
#'       parameter \eqn{\rho = 1 - 1/(1 + e^{\hat\gamma})} in (0, 1).}
#'   }
#'
#' @seealso \code{\link{int_evd_closed}}, \code{\link{cdfevd}}
#'
#' @export
SelectPie <- function(data, y1, x1, y2, x2,
                      method          = "BFGS",
                      maxit           = 1000,
                      multistart_corr = TRUE) {

  # ---- Extract outcome vectors ----
  y1v <- data[[y1]]
  y2v <- data[[y2]]

  # ---- Design matrices (intercept included automatically) ----
  X1 <- stats::model.matrix(stats::reformulate(x1, response = NULL), data = data)
  X2 <- stats::model.matrix(stats::reformulate(x2, response = NULL), data = data)

  p1 <- ncol(X1)
  p2 <- ncol(X2)

  idx_sel <- (y1v == 1)

  # ---- Starting values via fast OLS ----
  OLSStage2 <- stats::lm.fit(X2[idx_sel, , drop = FALSE], y2v[idx_sel])
  b2 <- stats::coef(OLSStage2)
  b2[is.na(b2)] <- 0

  if (length(unique(y1v)) == 1L) {
    # Degenerate selection: no variation in y1
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

  # ---- Log-likelihood (negated for minimisation) ----
  objFunc <- function(para) {
    beta1 <- para[seq_len(p1)]
    beta2 <- para[(p1 + 1L):(p1 + p2)]
    m     <- exp(para[length(para)]) + 1  # constrained to (1, Inf)

    xb1 <- as.vector(-(X1 %*% beta1))
    xb2 <- as.vector( (X2 %*% beta2))
    u   <- y2v - xb2

    # Non-selected log-contribution: log F_1(xb1)
    x1s   <- xb1 - digamma(1)
    logF1 <- -exp(-x1s)
    obj1  <- (1 - y1v) * logF1
    obj1[obj1 < -600] <- -600

    # Selected log-contribution: log integral
    integ <- int_evd_closed(xb1, u, m)
    logI  <- log(integ)
    logI[logI < -600] <- -600
    obj2  <- y1v * logI

    -sum(obj1 + obj2)
  }

  # ---- Primary optimisation ----
  res <- stats::optim(para0, objFunc,
                      method  = method,
                      control = list(maxit = maxit))

  # ---- Optional multi-start over dependence parameter ----
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

  # ---- Back-transform dependence parameter to rho in (0, 1) ----
  out <- res$par
  out[length(out)] <- 1 - 1 / (1 + exp(out[length(out)]))

  # ---- Attach names ----
  names(out) <- c(
    paste0("sel_",     colnames(X1)),
    paste0("outcome_", colnames(X2)),
    "rho"
  )

  out
}


#' latex_table
#'
#' Format a matrix or data frame as a LaTeX \code{table} environment.
#'
#' @param x A matrix or object coercible to one via \code{as.matrix}.
#' @param caption Character string or \code{NULL}. Table caption.
#' @param label Character string or \code{NULL}. LaTeX label for
#'   \code{\\ref\{\}} cross-referencing.
#' @param align Character string or \code{NULL}. Column alignment
#'   specification, e.g. \code{"lcccc"}.  Defaults to left-aligning the
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

  x  <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)

  # Default alignment: row-name column left, data columns centred
  if (is.null(align)) {
    align <- paste0("l", paste(rep("c", nc), collapse = ""))
  }

  # ---- Build LaTeX lines ----
  out <- c("\\begin{table}[!htbp]", "\\centering")

  if (!is.null(caption)) {
    out <- c(out, paste0("\\caption{", caption, "}"))
  }
  if (!is.null(label)) {
    out <- c(out, paste0("\\label{",   label,   "}"))
  }

  out <- c(out, paste0("\\begin{tabular}{", align, "}"))

  top_rule <- if (booktabs) "\\toprule" else "\\hline"
  mid_rule <- if (booktabs) "\\midrule" else "\\hline"
  bot_rule <- if (booktabs) "\\bottomrule" else "\\hline"

  # Header row
  header <- paste(c("", colnames(x)), collapse = " & ")
  out <- c(out, top_rule, paste0(header, " \\\\"), mid_rule)

  # Data rows
  for (i in seq_len(nr)) {
    row_i <- paste(c(rownames(x)[i], x[i, ]), collapse = " & ")
    out   <- c(out, paste0(row_i, " \\\\"))
  }

  out <- c(out, bot_rule, "\\end{tabular}", "\\end{table}")

  cat(paste(out, collapse = "\n"))
  invisible(out)
}
