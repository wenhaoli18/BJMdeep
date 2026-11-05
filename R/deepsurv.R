#' Fit a DeepSurv model
#'
#' @description A proportional hazards model whose log relative hazard is a
#' neural network of the covariates (DeepSurv; Katzman et al., 2018):
#' \eqn{h(t \mid x) = h_0(t) \exp\{g(x)\}}, with \eqn{g} a multilayer
#' perceptron trained on the Breslow partial likelihood (with weight decay,
#' early stopping and an ensemble of networks, see
#' \code{\link{deepsurvControl}}) and \eqn{h_0} the Breslow estimate given
#' \eqn{g}. Implemented in base R. Most users want
#' \code{\link{survivalSubDeep}}, which wraps this as a BJM survival
#' sub-model.
#'
#' @param formula A right-censored \code{Surv(time, status) ~ x} formula,
#'   without \code{strata()} terms.
#' @param data One row per subject.
#' @param control Output of \code{deepsurvControl()}.
#' @return An object of class \code{"deepsurv"}; \code{predict()} gives the
#'   risk score \eqn{g(x)} of new subjects.
#' @references Katzman JL, Shaham U, Cloninger A, Bates J, Jiang T, Kluger Y
#'   (2018). DeepSurv: personalized treatment recommender system using a Cox
#'   proportional hazards deep neural network. \emph{BMC Medical Research
#'   Methodology}, 18, 24.
#' @examples
#' library(survival)
#' data(pbc3, package = "BJM")
#' base <- pbc3[!duplicated(pbc3$id), ]
#' fit <- deepsurv(Surv(years, status3) ~ age + sex + edema, base,
#'                 deepsurvControl(hidden = 8, n_ensemble = 2, epochs = 300, seed = 1))
#' fit
#' head(predict(fit, base))
#' @export
deepsurv <- function(formula, data, control = deepsurvControl()) {
  form_marginal_surv <- formula
  if (!inherits(control, "deepsurvControl")) {
    stop("`control` must come from deepsurvControl().", call. = FALSE)
  }
  tt <- stats::terms(form_marginal_surv, specials = "strata", data = data)
  if (length(attr(tt, "specials")$strata) > 0) {
    stop("deepsurv() does not support strata() terms; ",
         "add the stratification variable as a covariate instead.", call. = FALSE)
  }
  mf <- stats::model.frame(form_marginal_surv, data = data, na.action = stats::na.omit)
  n_dropped <- nrow(data) - nrow(mf)
  if (n_dropped > 0) {
    warning(sprintf("deepsurv(): %d subject(s) with missing covariates or outcome were left out of the fit.",
                    n_dropped), call. = FALSE)
  }
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv") || attr(y, "type") != "right") {
    stop("deepsurv() needs a right-censored Surv(time, status) outcome.",
         call. = FALSE)
  }
  time <- unname(y[, "time"])
  status <- unname(y[, "status"])
  if (sum(status) < 2) stop("deepsurv() needs at least 2 events.", call. = FALSE)

  rhs <- stats::delete.response(stats::terms(mf))
  x_raw <- stats::model.matrix(rhs, mf)
  x_raw <- x_raw[, colnames(x_raw) != "(Intercept)", drop = FALSE]
  if (ncol(x_raw) == 0) {
    stop("deepsurv() needs at least one covariate in `formula`.",
         call. = FALSE)
  }
  center <- colMeans(x_raw)
  scale <- apply(x_raw, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  x <- sweep(sweep(x_raw, 2, center, "-"), 2, scale, "/")

  fits <- with_seed(control$seed, lapply(seq_len(control$n_ensemble), function(k)
    deepsurv_train_one(x, time, status, control, k)))
  nets <- lapply(fits, `[[`, "params")
  history <- lapply(fits, `[[`, "history")

  eta <- deepsurv_ensemble_eta(nets, x, control$activation)
  eta_shift <- mean(eta)
  eta <- eta - eta_shift

  c_index <- unname(survival::concordance(survival::Surv(time, status) ~ eta, reverse = TRUE)$concordance)

  out <- list(formula = form_marginal_surv, nets = nets, activation = control$activation, control = control,
       terms = rhs, xlevels = stats::.getXlevels(rhs, mf),
       contrasts = attr(x_raw, "contrasts"), x_colnames = colnames(x_raw),
       center = center, scale = scale, eta_shift = eta_shift,
       cum_basehaz = breslow_cum_basehaz(time, status, eta),
       time_variable = all.vars(form_marginal_surv[[2]])[1],
       max_followup = max(time), n = length(time), nevent = as.integer(sum(status)),
       n_dropped = n_dropped, c_index = c_index,
       loglik = cox_partial_loglik(eta, time, status),
       history = history, x = x_raw)
  class(out) <- "deepsurv"
  out
}

#' @describeIn deepsurv Risk score \eqn{g(x)} (log relative hazard, on the
#'   scale of the tabulated baseline hazard) of each row of \code{newdata};
#'   \code{NA} for a row with a missing covariate. The outcome columns are
#'   not needed.
#' @param object A \code{"deepsurv"} fit.
#' @param newdata Data frame of subjects to score.
#' @param ... Unused.
#' @export
predict.deepsurv <- function(object, newdata, ...) {
  deepsurv_lp(object, newdata)
}

#' @export
print.deepsurv <- function(x, ...) {
  ctl <- x$control
  cat("DeepSurv model (proportional hazards, neural log relative hazard)\n")
  cat("Formula: "); print(x$formula)
  cat(sprintf("Network : %s (%s)%s\n", paste(c(length(x$x_colnames), ctl$hidden, 1), collapse = " -> "),
              ctl$activation, if (length(ctl$hidden) == 0) ", linear" else ""))
  cat(sprintf("Ensemble: %d network(s), weight decay %g, learning rate %g\n",
              length(x$nets), ctl$weight_decay, ctl$learning_rate))
  steps <- vapply(x$history, function(h)
    as.numeric(if (all(is.na(h$valid_loss))) nrow(h) else which.min(h$valid_loss)), numeric(1))
  cat(sprintf("Steps kept: %s (early stopping on %g%% held out)\n",
              paste(steps, collapse = ", "), 100 * ctl$valid_fraction))
  cat(sprintf("n = %d, events = %d%s\n", x$n, x$nevent,
              if (x$n_dropped > 0) sprintf(", %d left out for missing values", x$n_dropped) else ""))
  cat(sprintf("Concordance = %.3f (apparent), log partial likelihood = %.2f\n", x$c_index, x$loglik))
  invisible(x)
}

#' Train one DeepSurv network
#'
#' @description Full-batch Adam on the mean negative Breslow log partial
#' likelihood per event plus the L2 penalty, with early stopping on a
#' held-out, event-stratified subset of subjects.
#' @return A list with \code{params} (the weights of the best step) and
#'   \code{history} (a data frame of training and held-out loss per step).
#' @keywords internal
deepsurv_train_one <- function(x, time, status, control, k = 1) {
  n <- nrow(x)
  valid <- logical(n)
  if (control$valid_fraction > 0) {
    for (s in c(0, 1)) {
      idx <- which(status == s)
      n_take <- floor(control$valid_fraction * length(idx))
      if (n_take > 0) valid[idx[sample.int(length(idx), n_take)]] <- TRUE
    }
    ### early stopping needs events on both sides of the split
    if (sum(status[valid]) < 2 || sum(status[!valid]) < 2) valid[] <- FALSE
  }
  use_valid <- any(valid)
  tr <- !valid

  params <- mlp_init(ncol(x), control$hidden, control$activation)
  adam <- lapply(params, function(p) list(m = 0 * p, v = 0 * p))
  b1 <- 0.9; b2 <- 0.999; eps <- 1e-8
  best <- params
  best_loss <- Inf
  since_best <- 0L
  hist_train <- hist_valid <- rep(NA_real_, control$epochs)

  for (epoch in seq_len(control$epochs)) {
    fw <- mlp_forward(params, x[tr, , drop = FALSE], control$activation)
    lg <- cox_loss_grad(fw$eta, time[tr], status[tr])
    grads <- mlp_backward(params, fw, lg$grad, control$activation)
    penalty <- 0
    for (j in seq_along(params)) {
      if (mlp_is_weight(params, j)) {
        penalty <- penalty + 0.5 * control$weight_decay * sum(params[[j]]^2)
        grads[[j]] <- grads[[j]] + control$weight_decay * params[[j]]
      }
      adam[[j]]$m <- b1 * adam[[j]]$m + (1 - b1) * grads[[j]]
      adam[[j]]$v <- b2 * adam[[j]]$v + (1 - b2) * grads[[j]]^2
      m_hat <- adam[[j]]$m / (1 - b1^epoch)
      v_hat <- adam[[j]]$v / (1 - b2^epoch)
      params[[j]] <- params[[j]] - control$learning_rate * m_hat / (sqrt(v_hat) + eps)
    }
    hist_train[epoch] <- lg$loss + penalty

    monitor <- if (use_valid) {
      eta_v <- mlp_forward(params, x[valid, , drop = FALSE], control$activation)$eta
      cox_loss_grad(eta_v, time[valid], status[valid])$loss
    } else hist_train[epoch]
    hist_valid[epoch] <- if (use_valid) monitor else NA_real_
    if (is.finite(monitor) && monitor < best_loss - 1e-7) {
      best_loss <- monitor
      best <- params
      since_best <- 0L
    } else {
      since_best <- since_best + 1L
      if (use_valid && since_best >= control$patience) break
    }
    if (control$verbose && (epoch %% 100 == 0 || epoch == 1)) {
      message(sprintf("network %d, step %d: train loss %.4f%s", k, epoch, hist_train[epoch],
                      if (use_valid) sprintf(", held-out loss %.4f", monitor) else ""))
    }
  }
  keep <- seq_len(epoch)
  list(params = if (use_valid) best else params,
       history = data.frame(epoch = keep, train_loss = hist_train[keep], valid_loss = hist_valid[keep]))
}

#' Mean negative Breslow log partial likelihood and its gradient
#'
#' @description For risk scores \code{eta}, the loss
#' \eqn{-\frac{1}{D}\sum_{i: d_i = 1} \{\eta_i - \log \sum_{j: t_j \ge t_i}
#' e^{\eta_j}\}} (\eqn{D} the number of events; tied event times handled as
#' by \code{coxph(ties = "breslow")}) and its gradient with respect to
#' \code{eta}, \eqn{-\frac{1}{D}\{d_j - e^{\eta_j} \sum_{i: d_i = 1,
#' t_i \le t_j} 1 / R_i\}}, both in \eqn{O(n \log n)}.
#' @return A list with \code{loss} and \code{grad}.
#' @keywords internal
cox_loss_grad <- function(eta, time, status) {
  n_events <- sum(status)
  if (n_events == 0) return(list(loss = NA_real_, grad = 0 * eta))
  m <- max(eta)
  r <- exp(eta - m)
  ### risk-set sums: decreasing time, cumulative sum, read at the last tie
  o <- order(time, decreasing = TRUE)
  ts <- time[o]
  cs <- cumsum(r[o])
  risk <- numeric(length(eta))
  risk[o] <- cs[findInterval(-ts, -ts)]
  loglik <- sum(status * (eta - m - log(risk)))
  ### sum over events with t_i <= t_j of 1 / R_i: increasing time, cumulative
  a <- status / risk
  o2 <- order(time)
  ta <- time[o2]
  ca <- cumsum(a[o2])
  A <- numeric(length(eta))
  A[o2] <- ca[findInterval(ta, ta)]
  list(loss = -loglik / n_events, grad = -(status - r * A) / n_events)
}

#' Breslow log partial likelihood
#' @keywords internal
cox_partial_loglik <- function(eta, time, status) {
  -sum(status) * cox_loss_grad(eta, time, status)$loss
}

#' Breslow baseline cumulative hazard given risk scores
#'
#' @description Tabulated, like \code{survival::basehaz(centered = FALSE)},
#' at every distinct observed time: \eqn{H_0(t) = \sum_{t_k \le t} d_k /
#' \sum_{j: t_j \ge t_k} e^{\eta_j}}.
#' @return A data frame with columns \code{hazard} and \code{time}.
#' @keywords internal
breslow_cum_basehaz <- function(time, status, eta) {
  ut <- sort(unique(time))
  d <- as.numeric(tapply(status, factor(time, levels = ut), sum))
  r <- exp(eta)
  o <- order(time, decreasing = TRUE)
  cs <- cumsum(r[o])
  ts <- time[o]
  ### risk set at each distinct time: everyone with t_j >= t
  risk_at <- cs[findInterval(-ut, -ts)]
  data.frame(hazard = cumsum(d / risk_at), time = ut)
}

# -- A small multilayer perceptron with a scalar, bias-free output -------------
# Parameters are a flat list: W1, b1, W2, b2, ..., w_out. The output has no
# bias because the partial likelihood does not identify one.

#' @keywords internal
mlp_init <- function(p, hidden, activation) {
  params <- list()
  fan_in <- p
  for (h in hidden) {
    sd <- switch(activation, tanh = sqrt(2 / (fan_in + h)), selu = sqrt(1 / fan_in), sqrt(2 / fan_in))
    params[[length(params) + 1]] <- matrix(stats::rnorm(fan_in * h, sd = sd), fan_in, h)
    params[[length(params) + 1]] <- matrix(0, 1, h)
    fan_in <- h
  }
  ### a linear model starts at the null model; a network's output layer
  ### starts small so the first steps do not saturate the partial likelihood
  params[[length(params) + 1]] <- if (length(hidden) == 0) matrix(0, fan_in, 1) else
    matrix(stats::rnorm(fan_in, sd = 0.1 / sqrt(fan_in)), fan_in, 1)
  params
}

#' @keywords internal
mlp_is_weight <- function(params, j) {
  ### W1, b1, ..., w_out: the biases sit at even positions before w_out
  j == length(params) || j %% 2 == 1
}

#' @keywords internal
mlp_activation <- function(z, activation) {
  switch(activation,
         relu = pmax(z, 0),
         tanh = tanh(z),
         elu = ifelse(z > 0, z, exp(pmin(z, 0)) - 1),
         selu = 1.0507009873554805 * ifelse(z > 0, z, 1.6732632423543772 * (exp(pmin(z, 0)) - 1)))
}

#' @keywords internal
mlp_activation_grad <- function(z, activation) {
  switch(activation,
         relu = (z > 0) * 1,
         tanh = 1 - tanh(z)^2,
         elu = ifelse(z > 0, 1, exp(pmin(z, 0))),
         selu = 1.0507009873554805 * ifelse(z > 0, 1, 1.6732632423543772 * exp(pmin(z, 0))))
}

#' @keywords internal
mlp_forward <- function(params, x, activation) {
  n_hidden <- (length(params) - 1) / 2
  h <- x
  hs <- list(x)
  zs <- list()
  for (l in seq_len(n_hidden)) {
    z <- h %*% params[[2 * l - 1]] + matrix(params[[2 * l]], nrow(h), ncol(params[[2 * l]]), byrow = TRUE)
    h <- mlp_activation(z, activation)
    zs[[l]] <- z
    hs[[l + 1]] <- h
  }
  list(eta = drop(h %*% params[[length(params)]]), hs = hs, zs = zs)
}

#' @keywords internal
mlp_backward <- function(params, fw, grad_eta, activation) {
  n_hidden <- (length(params) - 1) / 2
  grads <- vector("list", length(params))
  g <- matrix(grad_eta, ncol = 1)
  grads[[length(params)]] <- crossprod(fw$hs[[n_hidden + 1]], g)
  g_h <- g %*% t(params[[length(params)]])
  for (l in rev(seq_len(n_hidden))) {
    g_z <- g_h * mlp_activation_grad(fw$zs[[l]], activation)
    grads[[2 * l - 1]] <- crossprod(fw$hs[[l]], g_z)
    grads[[2 * l]] <- matrix(colSums(g_z), nrow = 1)
    g_h <- g_z %*% t(params[[2 * l - 1]])
  }
  grads
}

#' @keywords internal
deepsurv_ensemble_eta <- function(nets, x, activation) {
  ### matrix(): vapply() returns a plain vector for a single patient
  etas <- vapply(nets, function(p) mlp_forward(p, x, activation)$eta, numeric(nrow(x)))
  rowMeans(matrix(etas, nrow = nrow(x)))
}

#' Risk score of a DeepSurv fit for new subjects
#'
#' @description \eqn{g(x)}, averaged over the ensemble and on the scale of
#' the tabulated baseline hazard, for each row of \code{newdata}; \code{NA}
#' for a row with a missing covariate, so rows stay aligned with patients
#' (as \code{survival_lp()} requires).
#' @keywords internal
deepsurv_lp <- function(nn_fit, newdata) {
  mf <- stats::model.frame(nn_fit$terms, newdata, na.action = stats::na.pass, xlev = nn_fit$xlevels)
  x <- stats::model.matrix(nn_fit$terms, mf, contrasts.arg = nn_fit$contrasts)
  x <- x[, nn_fit$x_colnames, drop = FALSE]
  x <- sweep(sweep(x, 2, nn_fit$center, "-"), 2, nn_fit$scale, "/")
  out <- rep(NA_real_, nrow(newdata))
  ok <- stats::complete.cases(x)
  if (any(ok)) {
    out[ok] <- deepsurv_ensemble_eta(nn_fit$nets, x[ok, , drop = FALSE], nn_fit$activation) -
      nn_fit$eta_shift
  }
  out
}

