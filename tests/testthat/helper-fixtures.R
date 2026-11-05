library(survival)
data(pbc3, package = "BJM")
base <- pbc3[!duplicated(pbc3$id), ]

numeric_grad <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(i) {
    e <- replace(numeric(length(x)), i, h)
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))
}

setup_long_fixture <- function() {
  long_fit_all <- BJM::longitudinalSub(
    list(pbc3[pbc3$status3 == 1, ], pbc3[pbc3$status3 == 1, ]),
    list("long1" = serBilir ~ year + age + sex + (years) + (years) * year,
         "long2" = albumin ~ year + age + sex + (years) + (years) * year),
    list("long1" = ~ year | id, "long2" = ~ year | id))
  trans <- BJM::survivalTrans(c(1, 3, 5, 7))
  d <- pbc3[pbc3$id == 2 & pbc3$year <= 5, ]
  list(long_fit_all = long_fit_all, trans = trans, data_predict_all = list(d, d))
}

linear_control <- function() {
  deepsurvControl(hidden = integer(0), weight_decay = 0, learning_rate = 0.02, epochs = 4000,
                  valid_fraction = 0, n_ensemble = 1, seed = 1)
}
