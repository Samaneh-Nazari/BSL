###############################################################################
##  sim_aggregate.R - combine RDS files and produce tables + figures
###############################################################################
suppressPackageStartupMessages({
  library(MASS); library(ggplot2); library(reshape2)
  library(viridis); library(xtable); library(gridExtra)
})
source("BSL_core.R")
out_fig <- "figures"; out_tab <- "tables"
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

r20 <- readRDS("results_q20.rds")
r40a <- readRDS("results_q40a.rds")
r40b <- readRDS("results_q40b.rds")
all_results <- c(r20, r40a, r40b)

qs_main <- c(20, 40); ns_main <- c(100, 200, 300, 500)
metric_names <- c("SPE", "SEN", "FPR", "F1", "MISR",
                  "AC", "MCC", "PPV", "NPV", "SHD")
method_names <- c("nCPNG", "CPNIG", "GES", "PC", "NOTEARS")

build_metric_table <- function(target_q) {
  rows <- list()
  for (n in ns_main) {
    cell <- Filter(function(x) x$q == target_q && x$n == n, all_results)[[1]]
    for (mn in metric_names) {
      r <- list(n = n, metric = mn)
      for (meth in method_names) {
        vals <- sapply(cell$reps, function(x) x$metrics[[meth]][[mn]])
        r[[meth]] <- round(mean(vals, na.rm = TRUE), 4)
      }
      rows[[length(rows) + 1]] <- as.data.frame(r, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
tab_q20 <- build_metric_table(20)
tab_q40 <- build_metric_table(40)
write.csv(tab_q20, file.path(out_tab, "sim_metrics_q20.csv"), row.names = FALSE)
write.csv(tab_q40, file.path(out_tab, "sim_metrics_q40.csv"), row.names = FALSE)
sink(file.path(out_tab, "sim_metrics_q20.tex"))
print(xtable(tab_q20, caption = "Performance metrics for q=20",
             label = "tab:sim_q20"),
      include.rownames = FALSE, booktabs = TRUE)
sink()
sink(file.path(out_tab, "sim_metrics_q40.tex"))
print(xtable(tab_q40, caption = "Performance metrics for q=40",
             label = "tab:sim_q40"),
      include.rownames = FALSE, booktabs = TRUE)
sink()
cat("\n--- Metrics (q=20) ---\n"); print(tab_q20)
cat("\n--- Metrics (q=40) ---\n"); print(tab_q40)

## AUC table
auc_rows <- list()
for (q in qs_main) for (n in ns_main) {
  cell <- Filter(function(x) x$q == q && x$n == n, all_results)[[1]]
  auc_n <- mean(sapply(cell$reps, function(x) x$auc["nCPNG"]))
  auc_c <- mean(sapply(cell$reps, function(x) x$auc["CPNIG"]))
  auc_rows[[length(auc_rows) + 1]] <-
    data.frame(q = q, n = n, nCPNG = round(auc_n, 3),
               CPNIG = round(auc_c, 3))
}
auc_tab <- do.call(rbind, auc_rows)
write.csv(auc_tab, file.path(out_tab, "auc_table.csv"), row.names = FALSE)
sink(file.path(out_tab, "auc_table.tex"))
print(xtable(auc_tab, caption = "AUC averages", label = "tab:auc"),
      include.rownames = FALSE, booktabs = TRUE)
sink()
cat("\n--- AUC ---\n"); print(auc_tab)

## ROC plots
make_roc_plot <- function(target_q, out_file) {
  panels <- list()
  for (n in ns_main) {
    cell <- Filter(function(x) x$q == target_q && x$n == n, all_results)[[1]]
    rocs_n <- do.call(rbind,
                      lapply(cell$reps,
                             function(x) cbind(x$roc$nCPNG, Method = "nCPNG")))
    rocs_c <- do.call(rbind,
                      lapply(cell$reps,
                             function(x) cbind(x$roc$CPNIG, Method = "CPNIG")))
    df <- rbind(rocs_n, rocs_c)
    df_avg <- aggregate(cbind(SEN, FPR) ~ threshold + Method,
                        data = df, FUN = mean)
    auc_n <- round(mean(sapply(cell$reps, function(x) x$auc["nCPNG"])), 2)
    auc_c <- round(mean(sapply(cell$reps, function(x) x$auc["CPNIG"])), 2)
    p <- ggplot(df_avg, aes(x = FPR, y = SEN, color = Method)) +
      geom_line(linewidth = 1) + geom_point(size = 1.6) +
      geom_abline(linetype = "dashed", color = "grey60") +
      scale_color_manual(values = c("nCPNG" = "#D7263D", "CPNIG" = "#1B98E0")) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "False positive rate", y = "Sensitivity",
           title = sprintf("q = %d, n = %d", target_q, n),
           subtitle = sprintf("AUC: nCPNG = %.2f, CPNIG = %.2f",
                              auc_n, auc_c)) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold"))
    panels[[length(panels) + 1]] <- p
  }
  pdf(out_file, width = 8, height = 8)
  do.call(grid.arrange, c(panels, ncol = 2))
  dev.off()
  cat("Saved", out_file, "\n")
}
make_roc_plot(20, file.path(out_fig, "roc_q20.pdf"))
make_roc_plot(40, file.path(out_fig, "roc_q40.pdf"))

## Heatmap of metrics
hm_long <- list()
for (q in qs_main) for (n in ns_main) {
  cell <- Filter(function(x) x$q == q && x$n == n, all_results)[[1]]
  for (mn in c("F1", "MCC", "AC", "SEN", "SPE")) {
    for (meth in c("nCPNG", "CPNIG")) {
      v <- mean(sapply(cell$reps, function(x) x$metrics[[meth]][[mn]]))
      hm_long[[length(hm_long) + 1]] <- data.frame(
        q = q, n = n, metric = mn, method = meth, value = v)
    }
  }
}
hm_df <- do.call(rbind, hm_long)
hm_df$cell <- paste0("q=", hm_df$q, ", n=", hm_df$n)
p_hm <- ggplot(hm_df, aes(x = cell, y = metric, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.1) +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_viridis(option = "viridis", direction = -1, limits = c(0, 1)) +
  labs(x = "", y = "", fill = "Value",
       title = "Performance metrics across (q, n) settings") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"))
ggsave(file.path(out_fig, "perf_heatmap.pdf"), p_hm,
       width = 9, height = 6)
cat("Saved perf_heatmap.pdf\n")

## Benchmark comparison plot (only F1)
bench_rows <- list()
for (q in qs_main) for (n in ns_main) {
  cell <- Filter(function(x) x$q == q && x$n == n, all_results)[[1]]
  for (meth in method_names) {
    vals <- sapply(cell$reps, function(x) x$metrics[[meth]]$F1)
    bench_rows[[length(bench_rows) + 1]] <- data.frame(
      q = q, n = n, method = meth,
      F1 = mean(vals), F1_sd = sd(vals))
  }
}
bench_df <- do.call(rbind, bench_rows)
bench_df$method <- factor(bench_df$method,
                          levels = c("nCPNG", "CPNIG", "GES",
                                     "PC", "NOTEARS"))
p_bench <- ggplot(bench_df,
                  aes(x = factor(n), y = F1, fill = method)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_errorbar(aes(ymin = pmax(0, F1 - F1_sd),
                    ymax = pmin(1, F1 + F1_sd)),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~ q, labeller = label_bquote(q == .(q))) +
  scale_fill_manual(values = c("nCPNG" = "#D7263D", "CPNIG" = "#1B98E0",
                               "GES" = "#2E933C", "PC" = "#F46036",
                               "NOTEARS" = "#7C5295")) +
  labs(x = "Sample size n", y = "F1 score",
       fill = "Method",
       title = "Benchmarking the nCPNG prior against state-of-the-art structure learners") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(out_fig, "benchmark.pdf"), p_bench,
       width = 9, height = 5)
cat("Saved benchmark.pdf\n")

cat("\n=== Aggregation done ===\n")
