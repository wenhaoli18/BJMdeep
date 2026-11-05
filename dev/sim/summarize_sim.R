## Summarizes dev/sim/results/*.rds into tables and a figure.
## Run from the BJMdeep root:  Rscript dev/sim/summarize_sim.R
suppressPackageStartupMessages(library(ggplot2))
files <- list.files(file.path("dev", "sim", "results"), "^cell_.*\\.rds$", full.names = TRUE)
res <- do.call(rbind, lapply(files, readRDS))
res$measure[res$measure == "C-index (survival sub-model)"] <- "C-index"
res$scenario <- factor(res$scenario, c("linear", "nonlinear", "highdim"))

## AUC and Brier score averaged over the landmarks within each replicate
per_rep <- aggregate(value ~ rep + n_train + scenario + alpha + model + measure, res, mean)
wide <- reshape(per_rep, idvar = c("rep", "n_train", "scenario", "alpha", "measure"),
                timevar = "model", direction = "wide")
wide$ds_minus_cox <- wide$value.deepsurv - wide$value.cox
wide$oracle_minus_cox <- wide$value.oracle - wide$value.cox

mse <- function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
tab <- aggregate(cbind(value.cox, value.deepsurv, value.oracle, ds_minus_cox, oracle_minus_cox) ~
                   measure + scenario + n_train + alpha, wide, mean)
se <- aggregate(ds_minus_cox ~ measure + scenario + n_train + alpha, wide,
                function(v) sd(v) / sqrt(length(v)))
tab$ds_minus_cox_se <- se$ds_minus_cox
## share of the Cox-to-oracle gap that DeepSurv closes (only meaningful when the gap is clear)
tab$gap_closed <- ifelse(abs(tab$oracle_minus_cox) > 0.005, tab$ds_minus_cox / tab$oracle_minus_cox, NA)
tab$n_reps <- aggregate(ds_minus_cox ~ measure + scenario + n_train + alpha, wide, length)$ds_minus_cox
tab <- tab[order(tab$measure, tab$scenario, tab$alpha, tab$n_train), ]
write.csv(tab, file.path("dev", "sim", "sim_summary.csv"), row.names = FALSE)

fmt <- function(m) {
  t <- tab[tab$measure == m, ]
  data.frame(scenario = t$scenario, alpha = t$alpha, n = t$n_train,
             cox = sprintf("%.3f", t$value.cox), deepsurv = sprintf("%.3f", t$value.deepsurv),
             oracle = sprintf("%.3f", t$value.oracle),
             `deepsurv - cox` = sprintf("%+.3f (%.3f)", t$ds_minus_cox, t$ds_minus_cox_se),
             `gap closed` = ifelse(is.na(t$gap_closed), "-", sprintf("%.0f%%", 100 * t$gap_closed)),
             check.names = FALSE)
}
for (m in c("C-index", "AUC", "Brier score")) {
  cat(sprintf("\n== %s%s ==\n", m, if (m == "C-index") " of the survival sub-model (held out)"
              else " of the joint-model dynamic predictions, mean over landmarks 1/3/5 (held out)"))
  print(fmt(m), row.names = FALSE)
}

plot_dat <- aggregate(value ~ scenario + n_train + alpha + model + measure, per_rep,
                      function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
plot_dat <- do.call(data.frame, plot_dat)
plot_dat$biomarker <- factor(ifelse(plot_dat$alpha < 0.3, "weak biomarker (alpha = 0.1)",
                                    "strong biomarker (alpha = 0.5)"),
                             c("weak biomarker (alpha = 0.1)", "strong biomarker (alpha = 0.5)"))
plot_dat$model <- factor(plot_dat$model, c("cox", "deepsurv", "oracle"))
for (m in c("C-index", "AUC")) {
  d <- plot_dat[plot_dat$measure == m, ]
  g <- ggplot(d, aes(factor(n_train), value.mean, colour = model, group = model)) +
    geom_point(position = position_dodge(0.4), size = 2) +
    geom_errorbar(aes(ymin = value.mean - 2 * value.se, ymax = value.mean + 2 * value.se),
                  width = 0.2, position = position_dodge(0.4)) +
    facet_grid(biomarker ~ scenario) +
    scale_colour_manual(values = c(cox = "#1f77b4", deepsurv = "#d62728", oracle = "grey40")) +
    labs(x = "Training sample size", y = sprintf("Held-out %s (mean ± 2 SE)", m), colour = NULL,
         title = if (m == "C-index") "Survival sub-model alone" else
           "Joint-model dynamic prediction (mean over landmarks 1, 3, 5; 2-year window)") +
    theme_bw()
  ggsave(file.path("dev", "sim", sprintf("sim_%s.png", tolower(gsub("[^A-Za-z]", "", m)))),
         g, width = 9, height = 5.5, dpi = 150)
}
