###############################################################################
##  sim_part1a_q20.R — main simulation for q=20 only
###############################################################################
suppressPackageStartupMessages({
  library(MASS); library(Matrix); library(mvtnorm)
})
source("BSL_core.R")
set.seed(20260522)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

N_REP <- 3; N_ITER <- 4000; BURNIN <- 1000; EDGE_PRIOR <- 0.15
ns_main <- c(100, 200, 300, 500); q <- 20

run_one_rep <- function(q, n, seed_offset = 0) {
  set.seed(seed_offset)
  A_true <- random_dag(q, prob = 3 / (2 * q - 2))
  dat <- simulate_dag_data(n, A_true)
  alpha_p <- n + q - 5
  fit_n <- mcmc_dag(dat$X, n_iter = N_ITER, burn_in = BURNIN,
                    prior = "nCPNG", g = 2, alpha = alpha_p,
                    edge_prior = EDGE_PRIOR)
  A_n <- (fit_n$post_prob > 0.5) * 1
  m_n <- perf_metrics(A_n, A_true)
  roc_n <- roc_from_probs(fit_n$post_prob, A_true)
  auc_n <- auc_from_roc(roc_n)
  fit_c <- mcmc_dag(dat$X, n_iter = N_ITER, burn_in = BURNIN,
                    prior = "CPNIG", g = 2, edge_prior = EDGE_PRIOR)
  A_c <- (fit_c$post_prob > 0.5) * 1
  m_c <- perf_metrics(A_c, A_true)
  roc_c <- roc_from_probs(fit_c$post_prob, A_true)
  auc_c <- auc_from_roc(roc_c)
  A_g <- ges_like(dat$X, max_iter = 30)
  m_g <- perf_metrics(A_g, A_true)
  A_p <- pc_like(dat$X, alpha = 0.05)
  m_p <- perf_metrics(A_p, A_true)
  A_w <- tryCatch(notears_like(dat$X, lambda = 0.05),
                  error = function(e) matrix(0, q, q))
  m_w <- perf_metrics(A_w, A_true)
  list(metrics = list(nCPNG = m_n, CPNIG = m_c, GES = m_g,
                      PC = m_p, NOTEARS = m_w),
       roc = list(nCPNG = roc_n, CPNIG = roc_c),
       auc = c(nCPNG = auc_n, CPNIG = auc_c),
       q = q, n = n)
}

cat("\n=== q=20 simulation ===\n")
t0 <- Sys.time()
results_q20 <- list()
for (n in ns_main) {
  cat(sprintf("  (q=20, n=%d) ...", n))
  reps <- vector("list", N_REP)
  for (r in seq_len(N_REP)) {
    reps[[r]] <- run_one_rep(q, n, seed_offset = 1000 * q + 100 * n + r)
  }
  results_q20[[length(results_q20) + 1]] <- list(q = q, n = n, reps = reps)
  cat(sprintf(" done (%.1fs total)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
saveRDS(results_q20, "results_q20.rds")
cat("Saved results_q20.rds\n")
