test_that("survivalSubDeep() with a linear network reproduces BJM's Cox predictions", {
  fx <- setup_long_fixture()
  cox <- BJM::survivalSub(base, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex)
  lin <- survivalSubDeep(base, Surv(years, status3) ~ age + sex, status4 ~ years + age + sex,
                         control = linear_control())
  expect_s3_class(lin, c("survivalSubDeep", "survivalSub.BJM"))
  expect_equal(unname(coef(lin$coxph_fit)), 1)

  args <- list(prediction_time = 5, horizon = 1, time_variable = "year",
               survival_variable_all = fx$trans$survival_variable_all,
               survival_trans_function = fx$trans$survival_trans_function,
               bandcount1 = 20, bandcount2 = 40)
  r_cox <- do.call(BJM::predictRisk, c(list(fx$data_predict_all, fx$long_fit_all, cox), args))
  r_lin <- do.call(BJM::predictRisk, c(list(addDeepScore(fx$data_predict_all, lin),
                                            fx$long_fit_all, lin), args))
  expect_equal(unname(r_lin$risk_prob_1), unname(r_cox$risk_prob_1), tolerance = 1e-3)
  expect_equal(unname(r_lin$risk_prob_2), unname(r_cox$risk_prob_2), tolerance = 1e-3)
})

test_that("its baseline hazard is DeepSurv's Breslow estimate", {
  fit <- survivalSubDeep(base, Surv(years, status3) ~ age + sex + edema,
                         control = deepsurvControl(hidden = 4, epochs = 100, n_ensemble = 2, seed = 3))
  bh <- basehaz(fit$coxph_fit, centered = FALSE)
  expect_equal(bh$hazard, fit$deepsurv$cum_basehaz$hazard, tolerance = 1e-8)
  expect_equal(unname(predict(fit$coxph_fit, addDeepScore(base, fit), type = "lp", reference = "zero")),
               predict(fit$deepsurv, base))
})

test_that("a DeepSurv sub-model runs through BJM's prediction and evaluation", {
  fx <- setup_long_fixture()
  form <- Surv(years, status3) ~ age + sex + drug + edema + ascites + hepatomegaly + spiders
  fit <- survivalSubDeep(base, form, status4 ~ years + age + sex,
                         control = deepsurvControl(hidden = c(8, 4), epochs = 300, n_ensemble = 2,
                                                   seed = 11))
  expect_true(fit$deepsurv$c_index > 0.6)
  expect_output(print(fit), "addDeepScore")

  dp <- addDeepScore(fx$data_predict_all, fit)
  expect_true(all(vapply(dp, function(d) "deep_score" %in% names(d), logical(1))))
  risk <- BJM::predictRisk(dp, fx$long_fit_all, fit, prediction_time = 5, horizon = 1,
                           time_variable = "year", fx$trans$survival_variable_all,
                           fx$trans$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
  p <- unname(c(risk$risk_prob_1, risk$risk_prob_2))
  expect_true(all(p >= 0 & p <= 1))
  expect_true(sum(p) <= 1)

  imp <- survivalImportance(fit, base, n_perm = 3, seed = 1)
  expect_named(imp, c("variable", "c_drop", "sd"))
  expect_setequal(imp$variable, all.vars(form[[3]]))
})

test_that("survivalSubDeep() input errors are informative", {
  expect_error(survivalSubDeep(base, Surv(years, status3) ~ age, score_name = "age"), "already")
  expect_error(survivalSubDeep(base, Surv(years, status3) ~ age, score_name = "a b"), "syntactic")
  cox <- BJM::survivalSub(base, Surv(years, status3) ~ age, NULL)
  expect_error(addDeepScore(base, cox), "survivalSubDeep")
  expect_error(survivalImportance(cox, base), "survivalSubDeep")
})
