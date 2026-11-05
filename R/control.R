#' Control parameters for a DeepSurv survival sub-model
#'
#' @description Settings for \code{\link{deepsurv}} and
#' \code{\link{survivalSubDeep}}, which replace the linear predictor of the
#' marginal Cox model by a neural network (DeepSurv; Katzman et al., 2018): the hazard is
#' \eqn{h(t \mid x) = h_0(t) \exp\{g(x)\}}, with \eqn{g} a multilayer
#' perceptron of the baseline covariates trained on the Cox partial
#' likelihood, and \eqn{h_0} the Breslow estimate given \eqn{g}.
#'
#' @param hidden Integer vector of hidden-layer widths, e.g. \code{c(32, 16)}
#'   for two hidden layers. \code{integer(0)} gives no hidden layer, i.e. an
#'   ordinary linear Cox model fit by gradient descent (useful as a check).
#' @param activation Hidden-layer activation: \code{"relu"}, \code{"elu"},
#'   \code{"selu"} or \code{"tanh"}.
#' @param weight_decay L2 penalty on the network weights (not the biases),
#'   added to the mean negative log partial likelihood per event.
#' @param learning_rate Step size of the Adam optimizer.
#' @param epochs Maximum number of full-batch Adam steps.
#' @param valid_fraction Share of subjects (stratified by event status) held
#'   out of each network's training for early stopping. \code{0} trains on
#'   everyone for \code{epochs} steps.
#' @param patience Early stopping: training stops when the held-out partial
#'   likelihood has not improved for this many steps, and the weights of
#'   the best step are kept.
#' @param n_ensemble Number of networks, each with its own random start and
#'   held-out split; their \eqn{g(x)} are averaged. Averaging reduces the
#'   variance of a network fit on the few hundred events typical of
#'   joint-model data.
#' @param seed Optional random seed, for reproducible fits. R's random
#'   number generator state is restored afterwards.
#' @param verbose Print the training progress of each network.
#'
#' @return A list of class \code{"deepsurvControl"}.
#' @references Katzman JL, Shaham U, Cloninger A, Bates J, Jiang T, Kluger Y
#'   (2018). DeepSurv: personalized treatment recommender system using a Cox
#'   proportional hazards deep neural network. \emph{BMC Medical Research
#'   Methodology}, 18, 24.
#' @seealso \code{\link{survivalSubDeep}}, \code{\link{deepsurv}}
#' @examples
#' deepsurvControl(hidden = c(16, 8), n_ensemble = 3, seed = 1)
#' @export
deepsurvControl <- function(hidden = c(32L, 16L), activation = c("relu", "elu", "selu", "tanh"),
                            weight_decay = 1e-3, learning_rate = 1e-3, epochs = 2000L,
                            valid_fraction = 0.2, patience = 100L, n_ensemble = 5L,
                            seed = NULL, verbose = FALSE) {
  activation <- match.arg(activation)
  if (!is.numeric(hidden) || anyNA(hidden) || any(hidden < 1) || any(hidden != round(hidden))) {
    stop("`hidden` must be a vector of positive whole numbers (or integer(0) for no hidden layer).",
         call. = FALSE)
  }
  if (!is.numeric(weight_decay) || length(weight_decay) != 1 || is.na(weight_decay) || weight_decay < 0) {
    stop("`weight_decay` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(learning_rate) || length(learning_rate) != 1 || is.na(learning_rate) ||
      learning_rate <= 0) {
    stop("`learning_rate` must be a single positive number.", call. = FALSE)
  }
  assert_positive_integer(epochs, "epochs")
  assert_positive_integer(patience, "patience")
  assert_positive_integer(n_ensemble, "n_ensemble")
  if (!is.numeric(valid_fraction) || length(valid_fraction) != 1 || is.na(valid_fraction) ||
      valid_fraction < 0 || valid_fraction >= 1) {
    stop("`valid_fraction` must be a single number in [0, 1).", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || is.na(seed))) {
    stop("`seed` must be NULL or a single number.", call. = FALSE)
  }
  out <- list(hidden = as.integer(hidden), activation = activation, weight_decay = weight_decay,
              learning_rate = learning_rate, epochs = as.integer(epochs),
              valid_fraction = valid_fraction, patience = as.integer(patience),
              n_ensemble = as.integer(n_ensemble), seed = seed, verbose = isTRUE(verbose))
  class(out) <- "deepsurvControl"
  out
}


#' @keywords internal
assert_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x != round(x) || x < 1) {
    stop(sprintf("`%s` must be a single positive integer.", arg_name), call. = FALSE)
  }
  invisible(TRUE)
}

#' Run code with a seed, restoring the caller's random number generator state
#' @keywords internal
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  code
}
