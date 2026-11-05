# BJMdeep

**Experimental.** A neural (DeepSurv) survival sub-model for the
[BJM](https://github.com/wenhaoli18/BJM) backward joint model, packaged
separately so BJM itself stays unchanged.

The backward joint model factorizes as `f(T | X) · f(Y | T, X)`. The
survival sub-model only has to supply the marginal distribution of the event
time `T` given baseline covariates `X`. BJMdeep replaces BJM's linear Cox
model there with DeepSurv (Katzman et al., 2018): `h(t | x) = h0(t) exp{g(x)}`,
where `g` is a multilayer perceptron trained on the Cox partial likelihood
(weight decay, early stopping, and an ensemble of networks), and `h0` is the
Breslow estimate given `g`.

The model stays proportional in the hazard, so the network score `g(x)` is
stored as a data column and entered into an ordinary `coxph` model as its only
covariate, with the coefficient fixed at 1. The result is a regular
`survivalSub.BJM` object. `predictRisk()`, `predictLongitudinal()`,
`simulateTrajectory()` and `performancePlot()` use it unchanged. Everything is
base R; no deep-learning library is needed.

## Installation

```r
# install.packages("BJM")
devtools::install_github("wenhaoli18/BJMdeep")
```

## Usage

```r
library(BJM)
library(BJMdeep)
data(pbc3)
base <- pbc3[!duplicated(pbc3$id), ]

## 1. DeepSurv survival sub-model (instead of survivalSub())
surv_fit <- survivalSubDeep(
  base,
  Surv(years, status3) ~ age + sex + drug + edema + ascites + hepatomegaly + spiders + histologic,
  control = deepsurvControl(hidden = c(32, 16), activation = "elu", n_ensemble = 5, seed = 1)
)
surv_fit

## 2. Longitudinal sub-model: as usual
long_fit <- longitudinalSub(
  list(pbc3[pbc3$status3 == 1, ]),
  list("serBilir" = serBilir ~ year + age + sex + years),
  list("serBilir" = ~ year | id)
)

## 3. Prediction: the one extra step is addDeepScore() on the prediction data
trans <- survivalTrans(c(1, 3, 5, 7))
d <- addDeepScore(pbc3[pbc3$id == 2 & pbc3$year <= 5, ], surv_fit)
predictRisk(list(d), long_fit, surv_fit, prediction_time = 5, horizon = 2,
            time_variable = "year", trans$survival_variable_all,
            trans$survival_trans_function)

## Which covariates does the network use? (a network has no hazard ratios)
survivalImportance(surv_fit, base, n_perm = 20, seed = 1)
```

A network with no hidden layer, `deepsurvControl(hidden = integer(0))`, is a
linear Cox model. The tests check that it reproduces `coxph()` and BJM's
dynamic predictions.

## Does it help? pbc3, 5-fold cross-validation

`dev/deepsurv_pbc3_benchmark.R` keeps the longitudinal sub-model fixed
(serBilir and albumin) and swaps the survival sub-model. Results are held-out
means over 5 folds, with a 2-year prediction window.

| | Cox, age + sex | Cox, + 6 baseline variables | DeepSurv, same variables |
|---|---|---|---|
| Survival sub-model C-index | 0.574 | 0.761 | 0.763 |
| Joint-model AUC at landmarks 1 / 3 / 5 / 7 | .887 / .872 / .807 / .878 | .889 / .881 / .807 / .849 | .890 / .879 / .798 / .850 |
| Joint-model Brier score at landmarks 1 / 3 / 5 / 7 | .092 / .091 / .117 / .137 | .098 / .090 / .118 / .143 | .096 / .089 / .120 / .142 |

On pbc3, DeepSurv is no better than a linear Cox model with the same
covariates: every paired fold difference is within its standard error. Richer
baseline covariates improve the survival sub-model a lot, but they barely
change the dynamic predictions, which come mostly from the biomarkers. Each
fold has only 3–9 events per landmark window. A network is more likely to pay
off with large cohorts, high-dimensional baseline data (omics, imaging or text
embeddings), or strongly non-linear covariate effects.

## Limitations

- Right-censored outcomes only. No `strata()` and no interval censoring.
- The competing-risks (event-type) model is still BJM's logistic regression.
- The `coxph_fit` inside the returned object has a fixed coefficient, so its
  standard error and the forest `plot()` are not meaningful.

## Reference

Katzman JL, Shaham U, Cloninger A, Bates J, Jiang T, Kluger Y (2018).
DeepSurv: personalized treatment recommender system using a Cox proportional
hazards deep neural network. *BMC Medical Research Methodology*, 18, 24.
