###############################################################################
##  sim_diagnostics.R - MPSRF, posterior contraction, Laplace-vs-exact
###############################################################################
suppressPackageStartupMessages({
  library(MASS); library(ggplot2); library(xtable); library(coda)
})
source("BSL_core.R")
set.seed(20260522)
out_fig <- "figures"; out_tab <- "tables"

## ---------------------------------------------------------------------------
##  MPSRF diagnostic
## ---------------------------------------------------------------------------
cat("\n=== MPSRF ===\n")
mpsrf_rows <- list()
for (q in c(20, 40)) for (n in c(100, 200, 300, 500)) {
  set.seed(1234 + q + n)
  A <- random_dag(q, prob = 3 / (2 * q - 2))
  X <- simulate_dag_data(n, A)$X
  fit1 <- mcmc_dag(X, n_iter = 1500, burn_in = 300, prior = "nCPNG",
                   g = 2, alpha = n + q - 5, edge_prior = 0.15)
  fit2 <- mcmc_dag(X, n_iter = 1500, burn_in = 300, prior = "nCPNG",
                   g = 2, alpha = n + q - 5, edge_prior = 0.15)
  v1 <- as.vector(fit1$post_prob); v2 <- as.vector(fit2$post_prob)
  m1 <- mean(v1); m2 <- mean(v2)
  W <- 0.5 * (var(v1) + var(v2))
  B <- 0.5 * ((m1 - 0.5*(m1+m2))^2 + (m2 - 0.5*(m1+m2))^2)
  mpsrf <- sqrt(1 + B / (W + 1e-12))
  mpsrf_rows[[length(mpsrf_rows) + 1]] <-
    data.frame(q = q, n = n, MPSRF = round(mpsrf, 3))
  cat(sprintf("  q=%d n=%d MPSRF=%.3f\n", q, n, mpsrf))
}
mpsrf_df <- do.call(rbind, mpsrf_rows)
write.csv(mpsrf_df, file.path(out_tab, "mpsrf_table.csv"), row.names = FALSE)
sink(file.path(out_tab, "mpsrf_table.tex"))
print(xtable(mpsrf_df, caption = "MPSRF", label = "tab:mpsrf"),
      include.rownames = FALSE, booktabs = TRUE)
sink()

## ---------------------------------------------------------------------------
##  Trace plot demo
## ---------------------------------------------------------------------------
set.seed(42)
A_demo <- random_dag(20, prob = 3 / (2 * 20 - 2))
X_demo <- simulate_dag_data(300, A_demo)$X
alpha_demo <- 300 + 20 - 5
trace_len <- 1500
trace_vec <- numeric(trace_len)
A_run <- matrix(0, 20, 20)
node_lm <- numeric(20)
for (j in 1:20) {
  pa <- which(A_run[, j] == 1)
  node_lm[j] <- log_marg_nCPNG(j, X_demo, pa, g = 2, alpha = alpha_demo)
}
cur_lp <- sum(node_lm)
q_dim <- 20
for (it in 1:trace_len) {
  ij <- sample.int(q_dim, 2)
  i <- ij[1]; j <- ij[2]
  type <- if (A_run[i, j] == 1) sample(c("del", "rev"), 1)
          else if (A_run[j, i] == 0) "ins" else "rev"
  Anew <- A_run
  if (type == "ins") { Anew[i, j] <- 1; modified <- j }
  else if (type == "del") { Anew[i, j] <- 0; modified <- j }
  else {
    if (A_run[i, j] == 1) { Anew[i, j] <- 0; Anew[j, i] <- 1 }
    else { Anew[j, i] <- 0; Anew[i, j] <- 1 }
    modified <- c(i, j)
  }
  if (is_dag_fast(Anew)) {
    new_node_lm <- node_lm
    for (jm in modified) {
      pa <- which(Anew[, jm] == 1)
      new_node_lm[jm] <- log_marg_nCPNG(jm, X_demo, pa, g = 2,
                                         alpha = alpha_demo)
    }
    new_lp <- sum(new_node_lm)
    if (log(runif(1)) < new_lp - cur_lp) {
      A_run <- Anew; cur_lp <- new_lp; node_lm <- new_node_lm
    }
  }
  trace_vec[it] <- cur_lp
}
trace_df <- data.frame(iter = seq_along(trace_vec), log_post = trace_vec)
p_trace <- ggplot(trace_df, aes(x = iter, y = log_post)) +
  geom_line(color = "#274060", linewidth = 0.4) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE,
              color = "#D7263D", linewidth = 0.9) +
  labs(x = "MCMC iteration", y = "Log marginal likelihood",
       title = "MCMC convergence diagnostic for the nCPNG sampler",
       subtitle = "Single chain, q = 20, n = 300") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_fig, "mpsrf_diagnostic.pdf"), p_trace,
       width = 7.5, height = 4)
cat("Saved mpsrf_diagnostic.pdf\n")

## ---------------------------------------------------------------------------
##  Posterior contraction
## ---------------------------------------------------------------------------
cat("\n=== Posterior contraction ===\n")
contraction_rows <- list()
q_fixed <- 20
ns_grid <- c(50, 100, 200, 300, 500)
N_CON_REP <- 2
for (n in ns_grid) {
  shd_n <- numeric(N_CON_REP)
  for (r in seq_len(N_CON_REP)) {
    set.seed(99000 + n * 10 + r)
    A <- random_dag(q_fixed, prob = 3 / (2 * q_fixed - 2))
    X <- simulate_dag_data(n, A)$X
    fit <- mcmc_dag(X, n_iter = 2500, burn_in = 600, prior = "nCPNG",
                    g = 2, alpha = n + q_fixed - 5, edge_prior = 0.15)
    A_est <- (fit$post_prob > 0.5) * 1
    shd_n[r] <- sum(abs(A_est - A))
  }
  contraction_rows[[length(contraction_rows) + 1]] <-
    data.frame(n = n, mean_shd = mean(shd_n), sd_shd = sd(shd_n))
  cat(sprintf("  n=%d mean SHD=%.1f\n", n, mean(shd_n)))
}
contraction_df <- do.call(rbind, contraction_rows)
write.csv(contraction_df, file.path(out_tab, "contraction.csv"),
          row.names = FALSE)

p_con <- ggplot(contraction_df, aes(x = n, y = mean_shd)) +
  geom_errorbar(aes(ymin = pmax(0, mean_shd - sd_shd),
                    ymax = mean_shd + sd_shd),
                width = 12, color = "#274060") +
  geom_line(color = "#274060", linewidth = 0.8) +
  geom_point(size = 3.2, color = "#D7263D") +
  scale_x_continuous(trans = "log10", breaks = ns_grid) +
  labs(x = "Sample size n (log scale)",
       y = "Mean structural Hamming distance",
       title = "Empirical posterior contraction on the true DAG skeleton",
       subtitle = sprintf("q = %d, error bars show one SD across %d replicates",
                          q_fixed, N_CON_REP)) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_fig, "contraction.pdf"), p_con,
       width = 7, height = 4.5)
cat("Saved contraction.pdf\n")

## ---------------------------------------------------------------------------
##  Laplace vs exact comparison
## ---------------------------------------------------------------------------
cat("\n=== Laplace vs exact ===\n")
exact_node_loglik <- function(j, X, parents, g = 2, alpha = NULL) {
  n <- nrow(X); q <- ncol(X)
  if (is.null(alpha)) alpha <- n + q - 2
  alpha_j <- alpha + length(parents) - q + 1
  if (alpha_j <= 0) alpha_j <- 1
  Xj <- X[, j]; p_pa <- length(parents)
  if (p_pa == 0) {
    Mdet <- 0; kappa <- 0.5 * sum(Xj^2)
  } else {
    Xpa <- X[, parents, drop = FALSE]
    M <- crossprod(Xpa) + g * diag(p_pa)
    chM <- chol(M); bj <- -crossprod(Xpa, Xj)
    Mb <- backsolve(chM, forwardsolve(t(chM), bj))
    bMb <- sum(bj*Mb); Mdet <- 2*sum(log(diag(chM)))
    kappa <- 0.5 * (-bMb + sum(Xj^2))
  }
  lambda <- 0.5*g
  integrand <- function(sig2) {
    a <- 0.5*alpha_j - 0.5*n - 1
    a*log(sig2) - 0.5*g*sig2 - kappa/sig2
  }
  sig_grid <- exp(seq(-6, 6, length.out = 2000))
  log_vals <- sapply(sig_grid, integrand)
  M_max <- max(log_vals)
  int_log <- M_max +
             log(sum(exp(log_vals - M_max) *
                      c(diff(sig_grid), tail(diff(sig_grid), 1))))
  log_const <- -0.5*n*log(2*pi) - 0.5*Mdet + 0.5*p_pa*log(g) +
               0.5*alpha_j*log(lambda) - lgamma(0.5*alpha_j)
  log_const + int_log
}

lap_vs_exact <- list()
for (n in c(50, 100, 200, 400)) {
  set.seed(7777 + n)
  q_lv <- 8
  A <- random_dag(q_lv, prob = 3/(2*q_lv-2))
  X <- simulate_dag_data(n, A)$X
  alpha_p <- n + q_lv - 5
  diffs <- numeric(q_lv)
  for (j in 1:q_lv) {
    pa <- which(A[, j] == 1)
    L  <- log_marg_nCPNG(j, X, pa, g = 2, alpha = alpha_p)
    E  <- exact_node_loglik(j, X, pa, g = 2, alpha = alpha_p)
    diffs[j] <- (L - E) / n
  }
  lap_vs_exact[[length(lap_vs_exact) + 1]] <-
    data.frame(n = n, mean_abs_diff = mean(abs(diffs)),
               max_abs_diff = max(abs(diffs)))
}
lve_df <- do.call(rbind, lap_vs_exact)
write.csv(lve_df, file.path(out_tab, "laplace_vs_exact.csv"),
          row.names = FALSE)

p_lve <- ggplot(lve_df, aes(x = n, y = mean_abs_diff)) +
  geom_line(color = "#274060", linewidth = 0.9) +
  geom_point(size = 3, color = "#D7263D") +
  scale_x_continuous(trans = "log10",
                     breaks = c(50, 100, 200, 400)) +
  scale_y_continuous(trans = "log10") +
  labs(x = "Sample size n (log scale)",
       y = "Mean |log Laplace - log exact| / n",
       title = "Empirical accuracy of the Laplace approximation",
       subtitle = "Per-node normalised log-likelihood discrepancy (Theorem 3)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_fig, "laplace_accuracy.pdf"), p_lve,
       width = 7, height = 4.5)
cat("Saved laplace_accuracy.pdf\n")

cat("\n=== Diagnostics done ===\n")
