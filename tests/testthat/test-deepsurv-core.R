test_that("cox_loss_grad() matches coxph's Breslow partial likelihood, with ties", {
  time <- round(base$years)
  status <- base$status3
  eta <- 0.03 * (base$age - 50) + 0.4 * (base$sex == "female")
  fit <- coxph(Surv(time, status) ~ offset(eta), ties = "breslow")
  expect_equal(cox_partial_loglik(eta, time, status), fit$loglik, tolerance = 1e-8)
  expect_equal(cox_partial_loglik(eta + 7, time, status), fit$loglik, tolerance = 1e-8)
})

test_that("cox_loss_grad() gradient agrees with finite differences", {
  set.seed(1)
  n <- 40
  time <- sample(1:12, n, replace = TRUE)
  status <- rbinom(n, 1, 0.6)
  eta <- rnorm(n)
  expect_equal(cox_loss_grad(eta, time, status)$grad,
               numeric_grad(function(e) cox_loss_grad(e, time, status)$loss, eta), tolerance = 1e-6)
})

test_that("mlp_backward() agrees with finite differences for every activation", {
  set.seed(2)
  x <- matrix(rnorm(30 * 4), 30, 4)
  time <- rexp(30)
  status <- rbinom(30, 1, 0.7)
  for (act in c("relu", "elu", "selu", "tanh")) {
    params <- mlp_init(4, c(5, 3), act)
    params[[length(params)]] <- matrix(rnorm(3), 3, 1)
    fw <- mlp_forward(params, x, act)
    grads <- mlp_backward(params, fw, cox_loss_grad(fw$eta, time, status)$grad, act)
    for (j in seq_along(params)) {
      f <- function(v) {
        p <- params
        p[[j]][] <- v
        cox_loss_grad(mlp_forward(p, x, act)$eta, time, status)$loss
      }
      expect_equal(c(grads[[j]]), numeric_grad(f, c(params[[j]])), tolerance = 1e-5,
                   info = sprintf("%s, parameter block %d", act, j))
    }
  }
})

test_that("breslow_cum_basehaz() matches survival::basehaz(centered = FALSE)", {
  b <- c(0.035, -0.3)
  fit <- suppressWarnings(coxph(Surv(years, status3) ~ age + sex, data = base, ties = "breslow",
                                init = b, control = coxph.control(iter.max = 0)))
  eta <- drop(model.matrix(~ age + sex, base)[, -1] %*% b)
  ours <- breslow_cum_basehaz(base$years, base$status3, eta)
  theirs <- basehaz(fit, centered = FALSE)
  expect_equal(ours$time, theirs$time)
  expect_equal(ours$hazard, theirs$hazard, tolerance = 1e-8)
})

test_that("deepsurv() with no hidden layer is the Cox model", {
  form <- Surv(years, status3) ~ age + sex + edema
  cox <- coxph(form, base, ties = "breslow")
  lin <- deepsurv(form, base, linear_control())
  lp_cox <- unname(predict(cox, base, type = "lp"))
  lp_lin <- predict(lin, base)
  expect_equal(lp_lin - mean(lp_lin), lp_cox - mean(lp_cox), tolerance = 1e-3)
  expect_equal(lin$loglik, cox$loglik[2], tolerance = 1e-5)
})

test_that("deepsurv() predicts NA only for rows with a missing covariate", {
  fit <- deepsurv(Surv(years, status3) ~ age + sex + edema, base,
                  deepsurvControl(hidden = 4, epochs = 50, n_ensemble = 2, seed = 1))
  nd <- base[1:3, ]
  nd$edema[2] <- NA
  lp <- predict(fit, nd)
  expect_true(is.na(lp[2]))
  expect_equal(lp[c(1, 3)], predict(fit, base[c(1, 3), ]))
  expect_length(predict(fit, base[1, ]), 1)
  expect_output(print(fit), "DeepSurv")
})

test_that("fits are reproducible with a seed and leave the RNG alone", {
  ctl <- deepsurvControl(hidden = 4, epochs = 100, n_ensemble = 2, seed = 5)
  set.seed(42)
  before <- .Random.seed
  a <- deepsurv(Surv(years, status3) ~ age + sex + edema, base, ctl)
  expect_identical(.Random.seed, before)
  b <- deepsurv(Surv(years, status3) ~ age + sex + edema, base, ctl)
  expect_identical(predict(a, base), predict(b, base))
})

test_that("input errors are informative", {
  expect_error(deepsurvControl(hidden = 0), "positive whole numbers")
  expect_error(deepsurvControl(valid_fraction = 1), "valid_fraction")
  expect_error(deepsurv(Surv(years, status3) ~ age + strata(sex), base), "strata")
  expect_error(deepsurv(Surv(years, status3) ~ 1, base), "at least one covariate")
  expect_error(deepsurv(Surv(years, status3) ~ age, base, list(hidden = 4)), "deepsurvControl")
})
