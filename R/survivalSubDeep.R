#' DeepSurv survival sub-model for the backward joint model
#'
#' @description A drop-in replacement for \code{BJM::survivalSub()} whose
#' marginal survival model is DeepSurv (\code{\link{deepsurv}}) instead of a
#' linear Cox model. The backward joint model only asks this sub-model for
#' the marginal distribution of the event time given the baseline
#' covariates, and a DeepSurv model is still proportional in the hazard,
#' \eqn{S(t \mid x) = \exp\{-H_0(t) e^{g(x)}\}}. So the network's risk score
#' \eqn{g(x)} is stored as a data column (\code{score_name}) and entered
#' into an ordinary \code{coxph} model as its only covariate, with the
#' coefficient fixed at 1: its Breslow baseline hazard is then exactly
#' DeepSurv's, and the result is a regular \code{survivalSub.BJM} object
#' that \code{BJM::predictRisk()}, \code{BJM::predictLongitudinal()},
#' \code{BJM::simulateTrajectory()} and \code{BJM::performancePlot()} use
#' as is. BJM itself is not modified.
#'
#' \strong{The one extra step}: every data frame given to BJM's prediction
#' functions needs the score column; add it with \code{\link{addDeepScore}}.
#'
#' @param data_survival_fitting One row per subject, with the outcome and
#'   covariates of \code{form_marginal_surv} (and of
#'   \code{form_conditional_cr}).
#' @param form_marginal_surv A right-censored \code{Surv(time, status) ~
#'   covariates} formula: the network's inputs.
#' @param form_conditional_cr Optional competing-risks (event-type) formula,
#'   fit by logistic regression exactly as in \code{BJM::survivalSub()}.
#' @param control Output of \code{\link{deepsurvControl}}.
#' @param score_name Name of the score column added to the data.
#'
#' @return A \code{survivalSub.BJM} object (with the extra class
#'   \code{"survivalSubDeep"}) with, in addition, \code{deepsurv} (the
#'   \code{\link{deepsurv}} fit), \code{score_name} and
#'   \code{form_deepsurv} (\code{form_marginal_surv} as given).
#'   Its \code{coxph_fit} has the single coefficient 1 on the score; its
#'   standard error is not meaningful.
#' @details With a few hundred events a network rarely beats a
#' well-specified Cox model; compare the two on held-out data with
#' \code{BJM::performancePlot()} before relying on it.
#' @examples
#' \donttest{
#' library(BJM)
#' data(pbc3)
#' base <- pbc3[!duplicated(pbc3$id), ]
#' surv_fit <- survivalSubDeep(base,
#'   Surv(years, status3) ~ age + sex + drug + edema + ascites + hepatomegaly + spiders,
#'   control = deepsurvControl(hidden = c(16, 8), n_ensemble = 3, seed = 1))
#' surv_fit
#'
#' long_fit <- longitudinalSub(list(pbc3[pbc3$status3 == 1, ]),
#'   list("long1" = serBilir ~ year + age + sex + years), list("long1" = ~ year | id))
#' trans <- survivalTrans(c(1, 3, 5, 7))
#' d <- addDeepScore(pbc3[pbc3$id == 2 & pbc3$year <= 5, ], surv_fit)
#' predictRisk(list(d), long_fit, surv_fit, prediction_time = 5, horizon = 2,
#'   time_variable = "year", trans$survival_variable_all, trans$survival_trans_function)
#' }
#' @export
survivalSubDeep <- function(data_survival_fitting, form_marginal_surv, form_conditional_cr = NULL,
                            control = deepsurvControl(), score_name = "deep_score") {
  if (!is.data.frame(data_survival_fitting)) {
    stop("`data_survival_fitting` must be a data.frame.", call. = FALSE)
  }
  if (!inherits(form_marginal_surv, "formula")) {
    stop("`form_marginal_surv` must be a formula, e.g. Surv(time, status) ~ covariates.", call. = FALSE)
  }
  if (!is.character(score_name) || length(score_name) != 1 || !nzchar(score_name) ||
      make.names(score_name) != score_name) {
    stop("`score_name` must be a single syntactic column name.", call. = FALSE)
  }
  if (score_name %in% all.vars(form_marginal_surv)) {
    stop(sprintf("`score_name` (\"%s\") is already a variable of `form_marginal_surv`.", score_name),
         call. = FALSE)
  }
  nn <- deepsurv(form_marginal_surv, data_survival_fitting, control)

  d <- data_survival_fitting
  d[[score_name]] <- stats::predict(nn, d)
  surv_form <- stats::reformulate(score_name, response = form_marginal_surv[[2]])
  environment(surv_form) <- environment(form_marginal_surv)

  ### BJM builds the object (and the competing-risks model) the usual way;
  ### then the Cox coefficient on the score is fixed at 1, so the Breslow
  ### baseline hazard is the one given the network's score
  fit <- BJM::survivalSub(d, surv_form, form_conditional_cr)
  fit$coxph_fit <- suppressWarnings(survival::coxph(
    surv_form, data = d, ties = "breslow", init = 1,
    control = survival::coxph.control(iter.max = 0), x = TRUE, y = TRUE))
  fit$deepsurv <- nn
  fit$score_name <- score_name
  fit$form_deepsurv <- form_marginal_surv
  class(fit) <- c("survivalSubDeep", class(fit))
  fit
}

#' Add the DeepSurv score column for BJM's prediction functions
#'
#' @description Adds the network's risk score of each subject, as column
#' \code{survival_fit_all$score_name}, to a data frame or to every data frame
#' of a list (the \code{data_predict_all} of \code{BJM::predictRisk()} and
#' friends). The score depends only on the baseline covariates, so every row
#' of a subject gets the same value.
#'
#' @param data A data frame, or a list of data frames.
#' @param survival_fit_all A \code{\link{survivalSubDeep}} fit.
#' @return \code{data} with the score column added.
#' @export
addDeepScore <- function(data, survival_fit_all) {
  if (!inherits(survival_fit_all, "survivalSubDeep")) {
    stop("`survival_fit_all` must come from survivalSubDeep().", call. = FALSE)
  }
  add <- function(d) {
    if (!is.data.frame(d)) stop("`data` must be a data.frame or a list of data.frames.", call. = FALSE)
    d[[survival_fit_all$score_name]] <- stats::predict(survival_fit_all$deepsurv, d)
    d
  }
  if (is.data.frame(data)) add(data) else lapply(data, add)
}

#' @export
print.survivalSubDeep <- function(x, ...) {
  cat("BJM survival sub-model with a DeepSurv marginal survival model\n")
  cat(sprintf("(score column \"%s\"; add it to prediction data with addDeepScore())\n\n", x$score_name))
  print(x$deepsurv)
  if (!is.null(x$glm_fit)) {
    cat("\nCompeting-risks (event-type) model, logistic regression:\n")
    cat("Formula: "); print(x$form_conditional_cr)
    stats::printCoefmat(summary(x$glm_fit)$coefficients)
  }
  invisible(x)
}
