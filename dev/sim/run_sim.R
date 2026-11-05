## Runs every cell of the simulation design, one .rds per cell in
## dev/sim/results/ (cells already done are skipped, so it can be resumed).
## Run from the BJMdeep root:  Rscript dev/sim/run_sim.R [n_reps] [n_cores]
source("dev/sim/sim_functions.R")
args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 20
n_cores <- if (length(args) >= 2) as.integer(args[2]) else parallel::detectCores()

design <- expand.grid(scenario = c("linear", "nonlinear", "highdim"), n_train = c(500, 2000),
                      alpha = c(0.1, 0.5), stringsAsFactors = FALSE)
out_dir <- file.path("dev", "sim", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(design))) {
  cell <- design[i, ]
  file <- file.path(out_dir, sprintf("cell_%s_n%d_a%g.rds", cell$scenario, cell$n_train, cell$alpha))
  if (file.exists(file)) next
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_len(n_reps), function(r)
    tryCatch(one_replicate(r, cell$n_train, cell$scenario, cell$alpha, n_test = 500),
             error = function(e) { message("rep ", r, ": ", conditionMessage(e)); NULL }),
    mc.cores = n_cores)
  saveRDS(do.call(rbind, res), file)
  cat(sprintf("%s  cell %d/%d  %s n=%d alpha=%g  (%d reps ok, %.1f min)\n", format(Sys.time(), "%H:%M"),
              i, nrow(design), cell$scenario, cell$n_train, cell$alpha,
              sum(!vapply(res, is.null, logical(1))),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
