###############################################################################
##  sim_runtime.R - runtime scaling experiment
###############################################################################
suppressPackageStartupMessages({
  library(MASS); library(ggplot2); library(reshape2); library(xtable)
})
source("BSL_core.R")
set.seed(20260522)
out_fig <- "figures"; out_tab <- "tables"

cat("\n=== Runtime scaling ===\n")
runtime_rows <- list()
RT_REPS <- 1
for (q in c(10, 20, 40, 60)) {
  for (n in c(100, 500)) {
    cat(sprintf("  runtime q=%d n=%d\n", q, n))
    t_n <- t_c <- t_g <- t_w <- numeric(RT_REPS)
    for (r in seq_len(RT_REPS)) {
      set.seed(50000 + 1000 * q + 100 * n + r)
      A <- random_dag(q, prob = 3 / (2 * q - 2))
      X <- simulate_dag_data(n, A)$X
      n_it <- if (q <= 20) 400 else 200
      alpha_p <- n + q - 5
      t_n[r] <- system.time(mcmc_dag(X, n_iter = n_it, burn_in = 30,
                                      prior = "nCPNG", g = 2,
                                      alpha = alpha_p,
                                      edge_prior = 0.15))["elapsed"] /
                  n_it * 1000
      t_c[r] <- system.time(mcmc_dag(X, n_iter = n_it, burn_in = 30,
                                      prior = "CPNIG", g = 2,
                                      edge_prior = 0.15))["elapsed"] /
                  n_it * 1000
      t_g[r] <- system.time(ges_like(X,
                                     max_iter = if (q <= 20) 12 else 5))["elapsed"]
      t_w[r] <- system.time(tryCatch(notears_like(X, lambda = 0.05),
                              error = function(e) NULL))["elapsed"]
    }
    runtime_rows[[length(runtime_rows) + 1]] <- data.frame(
      q = q, n = n,
      nCPNG_ms = round(mean(t_n), 3),
      CPNIG_ms = round(mean(t_c), 3),
      GES_s = round(mean(t_g), 3),
      NOTEARS_s = round(mean(t_w), 3),
      nCPNG_sd = round(sd(t_n), 4),
      CPNIG_sd = round(sd(t_c), 4))
  }
}
runtime_df <- do.call(rbind, runtime_rows)
write.csv(runtime_df, file.path(out_tab, "runtime_table.csv"), row.names = FALSE)
sink(file.path(out_tab, "runtime_table.tex"))
print(xtable(runtime_df, caption = "Mean runtimes",
             label = "tab:runtime"),
      include.rownames = FALSE, booktabs = TRUE)
sink()
cat("\n--- Runtime ---\n"); print(runtime_df)

rt_long <- melt(runtime_df[, c("q", "n", "nCPNG_ms", "CPNIG_ms")],
                id.vars = c("q", "n"),
                variable.name = "method", value.name = "ms_per_iter")
rt_long$method <- gsub("_ms", "", rt_long$method)
rt_long$n <- factor(rt_long$n, levels = c(100, 500))

p_rt <- ggplot(rt_long, aes(x = q, y = ms_per_iter, color = method,
                            shape = n, linetype = n)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.6) +
  scale_color_manual(values = c("nCPNG" = "#D7263D", "CPNIG" = "#1B98E0")) +
  scale_x_continuous(breaks = c(10, 20, 40, 60)) +
  labs(x = "Number of nodes (q)",
       y = "Runtime per MCMC iteration (ms)",
       color = "Method", shape = "Sample size", linetype = "Sample size",
       title = "Computational scaling of nCPNG and CPNIG samplers") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_fig, "runtime_scaling.pdf"), p_rt,
       width = 7.5, height = 4.5)
cat("Saved runtime_scaling.pdf\n")
