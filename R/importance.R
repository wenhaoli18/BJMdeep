#' Permutation importance of the covariates of a DeepSurv survival sub-model
#'
#' @description A neural survival sub-model has no hazard ratios to report.
#' This measures how much each covariate contributes instead: the drop in
#' Harrell's concordance (C-index) of the fitted risk score when that
#' covariate's column is randomly permuted across subjects, averaged over
#' \code{n_perm} permutations. A larger drop means the network relies more
#' on that covariate; a drop near 0 means it barely uses it. Pass held-out
#' data to measure what generalizes rather than what was fit.
#'
#' @param object A \code{\link{survivalSubDeep}} or \code{\link{deepsurv}} fit.
#' @param data One row per subject, with the outcome and covariates of
#'   the network's formula.
#' @param n_perm Number of permutations per covariate.
#' @param seed Optional random seed.
#' @return A data frame with columns \code{variable}, \code{c_drop} (mean
#'   drop in C-index) and \code{sd}, sorted by \code{c_drop}, with the
#'   unpermuted C-index as attribute \code{"c_index"}.
#' @examples
#' \donttest{
#' data(pbc3, package = "BJM")
#' base <- pbc3[!duplicated(pbc3$id), ]
#' fit <- deepsurv(survival::Surv(years, status3) ~ age + sex + edema + ascites, base,
#'                 deepsurvControl(hidden = 8, n_ensemble = 2, epochs = 300, seed = 1))
#' survivalImportance(fit, base, n_perm = 5, seed = 1)
#' }
#' @export
survivalImportance <- function(object, data, n_perm = 20L, seed = NULL) {
  nn <- if (inherits(object, "survivalSubDeep")) object$deepsurv else object
  if (!inherits(nn, "deepsurv")) {
    stop("`object` must come from survivalSubDeep() or deepsurv().", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  assert_positive_integer(n_perm, "n_perm")
  form <- nn$formula
  missing_vars <- setdiff(all.vars(form), names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf("`data` lacks %s.", paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  y <- eval(form[[2]], data, environment(form))
  c_of <- function(d) {
    lp <- deepsurv_lp(nn, d)
    ok <- !is.na(lp) & !is.na(y[, "time"])
    unname(survival::concordance(y[ok] ~ lp[ok], reverse = TRUE)$concordance)
  }
  c0 <- c_of(data)
  vars <- all.vars(stats::delete.response(stats::terms(form)))
  res <- with_seed(seed, lapply(vars, function(v) {
    drops <- vapply(seq_len(n_perm), function(b) {
      d <- data
      d[[v]] <- d[[v]][sample.int(nrow(d))]
      c0 - c_of(d)
    }, numeric(1))
    data.frame(variable = v, c_drop = mean(drops), sd = stats::sd(drops))
  }))
  out <- do.call(rbind, res)
  out <- out[order(-out$c_drop), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "c_index") <- c0
  out
}
