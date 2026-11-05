## Cox vs DeepSurv survival sub-models inside the backward joint model, pbc3
##
## 5-fold cross-validation by subject. In every fold the longitudinal
## sub-model (serBilir, albumin) is fit once on the training subjects, and
## three marginal survival sub-models are swapped in under it:
##   cox_basic : Cox,      age + sex                      (the vignette model)
##   cox_rich  : Cox,      age + sex + 7 baseline clinical variables
##   deepsurv  : DeepSurv, same covariates as cox_rich
## Held-out subjects are scored with performancePlot() (landmark AUC and
## Brier score, 2-year window), and the baseline survival model alone with
## Harrell's C on the held-out fold.
##
## BJM itself is unmodified: the DeepSurv sub-model comes from BJMdeep and
## the prediction data get its score column with addDeepScore().
## Run from the BJMdeep root after installing BJM and BJMdeep:
##   Rscript dev/deepsurv_pbc3_benchmark.R
## Writes dev/results/deepsurv_pbc3_*.csv and a figure.

suppressPackageStartupMessages({
  library(BJM)
  library(BJMdeep)
  library(ggplot2)
})
data(pbc3)

out_dir <- file.path("dev", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_folds <- 5
landmarks <- c(1, 3, 5, 7)
horizon <- 2
set.seed(2026)

ids <- unique(pbc3$id)
fold_of <- setNames(sample(rep(seq_len(n_folds), length.out = length(ids))), ids)

rich <- ~ age + sex + drug + edema + ascites + hepatomegaly + spiders + histologic
forms <- list(
  cox_basic = Surv(years, status3) ~ age + sex,
  cox_rich  = update(rich, Surv(years, status3) ~ .),
  deepsurv  = update(rich, Surv(years, status3) ~ .)
)
ds_control <- deepsurvControl(hidden = c(32, 16), activation = "elu", weight_decay = 1e-3,
                              learning_rate = 1e-3, epochs = 3000, valid_fraction = 0.2,
                              patience = 150, n_ensemble = 5, seed = 1)

long_sub_fixed <- list("serBilir" = serBilir ~ year + age + sex + years,
                       "albumin"  = albumin  ~ year + age + sex + years)
long_sub_random <- list("serBilir" = ~ year | id, "albumin" = ~ year | id)
trans <- survivalTrans(c(1, 3, 5, 7))

perf_all <- list()
cidx_all <- list()
for (k in seq_len(n_folds)) {
  cat(sprintf("fold %d/%d\n", k, n_folds))
  train <- pbc3[fold_of[as.character(pbc3$id)] != k, ]
  test  <- pbc3[fold_of[as.character(pbc3$id)] == k, ]
  train_base <- train[!duplicated(train$id), ]
  ### the complete-case rule of the covariates must not differ between models
  keep <- complete.cases(train_base[, all.vars(rich)])
  train_base <- train_base[keep, ]

  long_fit_all <- longitudinalSub(train[train$status3 == 1, ], long_sub_fixed, long_sub_random)

  for (m in names(forms)) {
    t0 <- Sys.time()
    fit <- if (m == "deepsurv") {
      survivalSubDeep(train_base, forms[[m]], NULL, control = ds_control)
    } else {
      survivalSub(train_base, forms[[m]], NULL)
    }
    test_m <- if (m == "deepsurv") addDeepScore(test, fit) else test
    test_base_m <- test_m[!duplicated(test_m$id), ]
    secs_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    lp <- predict(fit$coxph_fit, newdata = test_base_m, type = "lp", reference = "zero",
                  na.action = na.pass)
    ok <- !is.na(lp)
    cidx <- survival::concordance(Surv(years, status3) ~ lp, data = test_base_m[ok, ],
                                  reverse = TRUE)$concordance
    cidx_all[[length(cidx_all) + 1]] <- data.frame(fold = k, model = m, c_index = cidx,
                                                   fit_seconds = secs_fit)

    p <- performancePlot(test_m, long_fit_all, fit, prediction_time = landmarks, horizon = horizon,
                         time_variable = "year", trans$survival_variable_all,
                         trans$survival_trans_function, bandcount1 = 10, bandcount2 = 20)
    perf_all[[length(perf_all) + 1]] <- cbind(fold = k, model = m, p$data)
  }
}

perf <- do.call(rbind, perf_all)
cidx <- do.call(rbind, cidx_all)
write.csv(perf, file.path(out_dir, "deepsurv_pbc3_by_fold.csv"), row.names = FALSE)

summ <- aggregate(value ~ model + measure + landmark, perf, function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
summ <- do.call(data.frame, summ)
names(summ)[names(summ) == "value.mean"] <- "mean"
names(summ)[names(summ) == "value.se"] <- "se"
write.csv(summ, file.path(out_dir, "deepsurv_pbc3_summary.csv"), row.names = FALSE)

c_summ <- aggregate(cbind(c_index, fit_seconds) ~ model, cidx, mean)
write.csv(c_summ, file.path(out_dir, "deepsurv_pbc3_cindex.csv"), row.names = FALSE)

cat("\nHeld-out C-index of the baseline survival sub-model (mean over folds):\n")
print(c_summ, row.names = FALSE, digits = 3)
cat("\nHeld-out landmark AUC / Brier of the joint-model predictions (mean over folds):\n")
print(reshape(summ[, c("model", "measure", "landmark", "mean")], idvar = c("measure", "landmark"),
              timevar = "model", direction = "wide"), row.names = FALSE, digits = 3)

g <- ggplot(summ, aes(landmark, mean, colour = model)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  facet_wrap(~ measure, scales = "free_y") +
  labs(x = "Landmark time (years)", y = "Held-out value (mean ± SE over 5 folds)",
       colour = "Survival sub-model",
       title = "pbc3: backward joint model with a Cox vs DeepSurv survival sub-model",
       subtitle = sprintf("%d-year prediction window, biomarkers serBilir + albumin", horizon)) +
  theme_bw()
ggsave(file.path(out_dir, "deepsurv_pbc3.png"), g, width = 8, height = 4, dpi = 150)
