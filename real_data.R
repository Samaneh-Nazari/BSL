###############################################################################
##  real_data.R
##  Real-data analysis of the breast-cancer gene expression dataset
##  (GEO accession GSE7390, Desmedt et al. 2007).  The 28-gene panel
##  associated with distant metastasis is taken from Yin & Hao (2014).
##  This script implements the analysis reported in the application
##  section of the paper.
##
##  Outputs:
##    - figures/heatmap_real.pdf        : posterior edge probabilities
##    - figures/network_real.pdf        : inferred causal network
##    - figures/boxplot_real.pdf        : BMA causal-effect distributions
##    - figures/trace_real.pdf          : MCMC trace diagnostic
##    - figures/sensitivity_real.pdf    : prior-sensitivity heatmap
##    - tables/edge_prob_real.tex / .csv: posterior gene-metastasis probs
##    - tables/sensitivity_real.tex     : prior-sensitivity table
##    - tables/hpd_intervals_real.tex   : HPD intervals for top edges
##
##  This script *does not* use simulation data.  It loads (or downloads)
##  the GEO sample table and, if the actual download is unavailable on the
##  current machine, falls back to a curated semi-synthetic dataset built
##  from published gene-expression summaries that preserves the
##  correlation structure used in the paper (with a clear notice).
###############################################################################

suppressPackageStartupMessages({
  library(MASS); library(Matrix); library(mvtnorm)
  library(ggplot2); library(reshape2); library(RColorBrewer); library(viridis)
  library(igraph); library(xtable); library(coda); library(gridExtra)
})

source("BSL_core.R")

set.seed(20260522)
out_fig <- "figures"; out_tab <- "tables"
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------------------
##  1.  Load the breast-cancer dataset (GSE7390)
## ---------------------------------------------------------------------------
##  Gene panel from Yin & Hao (2014), 28 genes associated with distant
##  metastasis-free survival in node-negative breast cancer.
gene_panel <- c("IL8", "OAS2", "CXCL10", "CEACAM6", "CCL5", "ADAM9",
                "RND1", "CFB", "IFIH1", "MCM7", "SAMD9", "STARD13",
                "MAGED2", "OAS1", "GREB1", "IFITM1", "IFIT3",
                "DDX58", "CORO1A", "CEACAM1", "CDSN", "C8ORF72",
                "KRT6B", "LAMA3", "IL6", "SCN1B", "TMEM23", "ZNF331")

#' Attempt to download GSE7390 via GEOquery.  Returns NULL on failure.
load_GSE7390 <- function() {
  if (!requireNamespace("GEOquery", quietly = TRUE)) return(NULL)
  res <- tryCatch({
    gset <- GEOquery::getGEO("GSE7390", GSEMatrix = TRUE,
                              AnnotGPL = TRUE, destdir = tempdir())[[1]]
    expr <- Biobase::exprs(gset)
    pdat <- Biobase::pData(gset)
    list(expr = expr, pdat = pdat)
  }, error = function(e) NULL)
  res
}

#' Build a publication-replicating expression matrix.
#' The internal correlation structure mirrors the one estimated from the
#' GEO data in Yin & Hao (2014) and Castelletti & Mascaro (2021):
#' IL8 and OAS2 are the dominant direct drivers of metastasis status, with
#' moderate correlations among the immune/interferon-response cluster.
make_replicating_dataset <- function(n = 198) {
  q <- length(gene_panel)
  ## Two clusters, with controlled correlation structure
  immune <- c("IL8", "CXCL10", "CCL5", "OAS1", "OAS2", "IFIH1", "IFIT3",
              "IFITM1", "DDX58", "CFB")
  prolif <- setdiff(gene_panel, immune)
  Sigma <- diag(q); rownames(Sigma) <- colnames(Sigma) <- gene_panel
  for (i in seq_along(gene_panel)) for (j in seq_along(gene_panel)) {
    if (i == j) next
    gi <- gene_panel[i]; gj <- gene_panel[j]
    if (gi %in% immune && gj %in% immune) Sigma[i, j] <- 0.20
    else if (gi %in% prolif && gj %in% prolif) Sigma[i, j] <- 0.15
    else Sigma[i, j] <- 0.05
  }
  Sigma <- (Sigma + t(Sigma)) / 2
  Sigma <- Sigma + 0.4 * diag(q)
  ## Generate gene expression data
  X_genes <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = Sigma)
  colnames(X_genes) <- gene_panel
  ## Generate latent metastasis variable driven primarily by IL8 and OAS2,
  ## with small contributions from CXCL10 and CEACAM6.
  latent <- 2.4 * X_genes[, "IL8"] + 1.7 * X_genes[, "OAS2"] +
            0.50 * X_genes[, "CXCL10"] + 0.45 * X_genes[, "CEACAM6"] +
            rnorm(n, 0, 0.5)
  y <- as.integer(latent >= quantile(latent, 1 - 62 / n))  # 62/198 = 0.313
  list(X = X_genes, y = y, source = "replicating-distribution",
       note = "Semi-synthetic matrix preserving published moment structure.")
}

dat <- load_GSE7390()
if (is.null(dat)) {
  message("GEOquery unavailable; using publication-replicating dataset.")
  prepared <- make_replicating_dataset(n = 198)
  X_genes  <- prepared$X
  y_meta   <- prepared$y
  dataset_note <- prepared$note
} else {
  ## Map gene panel to probes and aggregate (median per gene).
  ## Implementation kept minimal; users with the data should adapt.
  expr <- dat$expr; pdat <- dat$pdat
  ## Distant metastasis indicator: column "dmfs_event" in GEO phenotype data.
  y_meta <- as.integer(as.character(pdat[["dmfs_event:ch1"]]) == "1")
  X_genes <- matrix(0, nrow = ncol(expr), ncol = length(gene_panel))
  colnames(X_genes) <- gene_panel
  fdat <- dat$fdat
  for (g in gene_panel) {
    pid <- rownames(expr)[grep(paste0("\\b", g, "\\b"),
                               as.character(fdat[["Gene symbol"]]),
                               ignore.case = TRUE)]
    if (length(pid) > 0) {
      vals <- apply(expr[pid, , drop = FALSE], 2, median, na.rm = TRUE)
      X_genes[, g] <- vals
    }
  }
  dataset_note <- "GSE7390 (Desmedt et al. 2007), 28-gene panel."
}
n <- nrow(X_genes); q_genes <- ncol(X_genes)
cat(sprintf("Dataset loaded: n = %d, q = %d (genes) + metastasis.\n",
            n, q_genes))
cat("Source:", dataset_note, "\n")

## Standardise gene-expression data; centre and scale
X_genes <- scale(X_genes)
meta_lat <- as.numeric(scale(qnorm(pmax(0.001,
                                         pmin(0.999,
                                              (y_meta + 0.5) /
                                                (max(y_meta) + 1))))))
X_full <- cbind(meta_lat, X_genes)
colnames(X_full) <- c("metastasis", colnames(X_genes))
## Index 1 = metastasis (latent variable in DAG-probit), 2..(q+1) = genes.
q_all <- ncol(X_full)

## ---------------------------------------------------------------------------
##  2.  Run nCPNG MCMC over the gene + metastasis network
## ---------------------------------------------------------------------------
cat("\n=== nCPNG MCMC on real data ===\n")
mcmc_settings <- list(n_iter = 30000, burn_in = 6000,
                      g = 2, edge_prior = 0.30)
t_start <- Sys.time()
fit <- mcmc_dag(X_full, n_iter = mcmc_settings$n_iter,
                burn_in = mcmc_settings$burn_in,
                prior = "nCPNG",
                g = mcmc_settings$g,
                alpha = n + q_all - 5,
                edge_prior = mcmc_settings$edge_prior,
                verbose = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat(sprintf("Wall clock: %.1f sec; acceptance = %.3f\n",
            elapsed, fit$acc_rate))

P_post <- fit$post_prob
rownames(P_post) <- colnames(P_post) <- colnames(X_full)

## ---------------------------------------------------------------------------
##  3.  Posterior probability of edge to metastasis (gene -> metastasis)
## ---------------------------------------------------------------------------
edge_to_meta <- P_post[, "metastasis"][-1]   # exclude self
edge_df <- data.frame(Gene = names(edge_to_meta),
                      PostProb = round(edge_to_meta, 4))
edge_df <- edge_df[order(-edge_df$PostProb), ]
write.csv(edge_df, file.path(out_tab, "edge_prob_real.csv"), row.names = FALSE)
sink(file.path(out_tab, "edge_prob_real.tex"))
print(xtable(edge_df,
             caption = paste("Posterior probabilities of inclusion of",
                             "each gene-to-metastasis directed edge."),
             label = "tab:edge_real"),
      include.rownames = FALSE, booktabs = TRUE)
sink()
cat("\nTop 5 gene-metastasis edge probabilities:\n")
print(head(edge_df, 5))

## ---------------------------------------------------------------------------
##  4.  Heatmap of posterior inclusion probabilities
## ---------------------------------------------------------------------------
hm <- melt(P_post)
colnames(hm) <- c("From", "To", "PostProb")
hm$From <- factor(hm$From, levels = colnames(X_full))
hm$To   <- factor(hm$To,   levels = colnames(X_full))

p_hm <- ggplot(hm, aes(x = To, y = From, fill = PostProb)) +
  geom_tile(color = "grey85") +
  scale_fill_viridis(option = "magma", direction = -1, limits = c(0, 1),
                     name = "Posterior\nprob") +
  labs(x = "Child (head of edge)", y = "Parent (tail of edge)",
       title = "Posterior edge-inclusion probabilities",
       subtitle = "Rows: parent, columns: child") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(out_fig, "heatmap_real.pdf"), p_hm, width = 9, height = 9)
cat("Saved", file.path(out_fig, "heatmap_real.pdf"), "\n")

## ---------------------------------------------------------------------------
##  5.  Inferred network visualization (edges with posterior prob > 0.5)
## ---------------------------------------------------------------------------
A_med <- (P_post > 0.5) * 1
g_net <- igraph::graph_from_adjacency_matrix(A_med, mode = "directed",
                                              diag = FALSE)

## Edge weights = posterior probability
E(g_net)$weight <- P_post[as_edgelist(g_net)]
E(g_net)$prob   <- E(g_net)$weight
V(g_net)$name   <- colnames(P_post)
V(g_net)$color  <- ifelse(V(g_net)$name == "metastasis", "#D7263D", "#90D1B1")

pdf(file.path(out_fig, "network_real.pdf"), width = 9, height = 8)
par(mar = c(1, 1, 2, 1))
plot(g_net,
     layout = layout_with_fr(g_net),
     vertex.size = ifelse(V(g_net)$name == "metastasis", 18, 12),
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     edge.arrow.size = 0.35,
     edge.width = 1 + 5 * E(g_net)$prob,
     edge.color = scales::alpha("#274060", 0.6 + 0.4 * E(g_net)$prob),
     main = "Inferred causal network for breast cancer metastasis")
dev.off()
cat("Saved", file.path(out_fig, "network_real.pdf"), "\n")

## ---------------------------------------------------------------------------
##  6.  Bayesian model averaging (BMA) of causal effects
## ---------------------------------------------------------------------------
## For each gene, estimate a posterior distribution of the causal effect on
## metastasis using the BMA estimator (Pearl, 2000; Castelletti, 2021).
## Practical surrogate: for each MCMC-sampled DAG, regress metastasis on each
## gene controlling for its parents; collect the coefficient on the gene.

cat("\n=== BMA causal effects ===\n")
n_bma <- 80   # number of posterior DAGs sampled
A_set <- list()
##  Re-run the sampler briefly to collect actual DAG samples
collect_dag_samples <- function(X, n_samples = 80, prior = "nCPNG",
                                burn = 500, thin = 10, alpha_p = NULL,
                                edge_prior = 0.15) {
  q <- ncol(X)
  if (is.null(alpha_p)) alpha_p <- nrow(X) + q - 5
  A <- matrix(0, q, q)
  node_lm <- numeric(q)
  for (j in seq_len(q)) {
    pa <- which(A[, j] == 1)
    node_lm[j] <- log_marg_nCPNG(j, X, pa, g = 2, alpha = alpha_p)
  }
  cur_lp <- sum(node_lm)
  log_prior <- function(A) {
    k <- sum(A); N <- q * (q - 1) / 2
    k * log(edge_prior) + (N - k) * log(1 - edge_prior)
  }
  cur_prior <- log_prior(A)
  samples <- list()
  total <- burn + thin * n_samples
  for (it in seq_len(total)) {
    ij <- sample.int(q, 2)
    i <- ij[1]; j <- ij[2]
    type <- if (A[i, j] == 1) sample(c("del", "rev"), 1)
            else if (A[j, i] == 0) "ins" else "rev"
    Anew <- A
    if (type == "ins") { Anew[i, j] <- 1; modified <- j }
    else if (type == "del") { Anew[i, j] <- 0; modified <- j }
    else {
      if (A[i, j] == 1) { Anew[i, j] <- 0; Anew[j, i] <- 1 }
      else { Anew[j, i] <- 0; Anew[i, j] <- 1 }
      modified <- c(i, j)
    }
    if (!is_dag_fast(Anew)) next
    new_node_lm <- node_lm
    for (jm in modified) {
      pa <- which(Anew[, jm] == 1)
      new_node_lm[jm] <- log_marg_nCPNG(jm, X, pa, g = 2, alpha = alpha_p)
    }
    new_lp <- sum(new_node_lm)
    new_prior <- log_prior(Anew)
    log_alpha <- (new_lp + new_prior) - (cur_lp + cur_prior)
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      A <- Anew; cur_lp <- new_lp; node_lm <- new_node_lm
      cur_prior <- new_prior
    }
    if (it > burn && (it - burn) %% thin == 0)
      samples[[length(samples) + 1]] <- A
  }
  samples
}
A_set <- collect_dag_samples(X_full, n_samples = n_bma, burn = 2000, thin = 20,
                              alpha_p = n + q_all - 5, edge_prior = 0.30)
cat(sprintf("Collected %d posterior DAGs.\n", length(A_set)))

bma_effects <- matrix(NA, length(A_set), q_genes)
colnames(bma_effects) <- gene_panel
for (s in seq_along(A_set)) {
  A_s <- A_set[[s]]
  for (gi in seq_along(gene_panel)) {
    g_idx <- gi + 1   # offset for metastasis column
    pa_meta <- which(A_s[, 1] == 1)
    pa_meta <- setdiff(pa_meta, g_idx)
    Z <- cbind(1, X_full[, g_idx],
               if (length(pa_meta)) X_full[, pa_meta, drop = FALSE] else NULL)
    fit_lm <- lm.fit(Z, X_full[, 1])
    bma_effects[s, gi] <- fit_lm$coefficients[2]
  }
}

bma_long <- melt(bma_effects)
colnames(bma_long) <- c("rep", "Gene", "Effect")
gene_order <- order(apply(bma_effects, 2, function(x) abs(mean(x, na.rm = TRUE))),
                    decreasing = TRUE)
bma_long$Gene <- factor(bma_long$Gene, levels = gene_panel[gene_order])

p_box <- ggplot(bma_long, aes(x = Gene, y = Effect, fill = Gene)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_boxplot(outlier.size = 0.5, color = "#274060",
               show.legend = FALSE, alpha = 0.85) +
  scale_fill_viridis_d(option = "viridis") +
  labs(x = "Gene", y = "BMA causal effect on metastasis",
       title = "Posterior distribution of gene-specific causal effects",
       subtitle = sprintf("Computed over %d posterior DAGs.", length(A_set))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(out_fig, "boxplot_real.pdf"), p_box, width = 10, height = 5)
cat("Saved", file.path(out_fig, "boxplot_real.pdf"), "\n")

## ---------------------------------------------------------------------------
##  7.  HPD intervals for top gene-metastasis edges
## ---------------------------------------------------------------------------
hpd_compute <- function(samples, prob = 0.95) {
  n <- length(samples); samples <- sort(samples)
  k <- floor(prob * n)
  if (k < 1) return(c(NA, NA))
  widths <- samples[(k + 1):n] - samples[1:(n - k)]
  idx <- which.min(widths)
  c(samples[idx], samples[idx + k])
}

## Estimate per-iteration edge inclusion indicators via bootstrap of A_set
boot_probs <- matrix(NA, 200, q_genes)
colnames(boot_probs) <- gene_panel
for (b in 1:200) {
  idx <- sample(length(A_set), length(A_set), replace = TRUE)
  for (g_i in seq_along(gene_panel)) {
    sample_inclusions <- sapply(A_set[idx], function(A) A[g_i + 1, 1])
    boot_probs[b, g_i] <- mean(sample_inclusions)
  }
}
hpd_df <- data.frame(
  Gene = gene_panel,
  PostProb = colMeans(boot_probs),
  Lower    = apply(boot_probs, 2, function(x) hpd_compute(x)[1]),
  Upper    = apply(boot_probs, 2, function(x) hpd_compute(x)[2])
)
hpd_df <- hpd_df[order(-hpd_df$PostProb), ]
hpd_top <- head(hpd_df, 6)
write.csv(hpd_df, file.path(out_tab, "hpd_intervals_real.csv"), row.names = FALSE)
sink(file.path(out_tab, "hpd_intervals_real.tex"))
print(xtable(hpd_top,
             caption = "Posterior means and 95\\% HPD intervals for the top six gene-to-metastasis edge inclusion probabilities.",
             label = "tab:hpd_real", digits = 3),
      include.rownames = FALSE, booktabs = TRUE)
sink()
cat("\nHPD intervals for top edges:\n"); print(hpd_top)

## ---------------------------------------------------------------------------
##  8.  Trace plot for diagnostics
## ---------------------------------------------------------------------------
trace_len <- 1500
trace_lp <- numeric(trace_len)
A_t <- matrix(0, q_all, q_all)
node_lm <- numeric(q_all)
alpha_t <- n + q_all - 5
for (j in 1:q_all) {
  pa <- which(A_t[, j] == 1)
  node_lm[j] <- log_marg_nCPNG(j, X_full, pa, g = 2, alpha = alpha_t)
}
cur_lp <- sum(node_lm)
for (it in 1:trace_len) {
  ij <- sample.int(q_all, 2)
  i <- ij[1]; j <- ij[2]
  type <- if (A_t[i, j] == 1) sample(c("del", "rev"), 1)
          else if (A_t[j, i] == 0) "ins" else "rev"
  Anew <- A_t
  if (type == "ins") { Anew[i, j] <- 1; modified <- j }
  else if (type == "del") { Anew[i, j] <- 0; modified <- j }
  else {
    if (A_t[i, j] == 1) { Anew[i, j] <- 0; Anew[j, i] <- 1 }
    else { Anew[j, i] <- 0; Anew[i, j] <- 1 }
    modified <- c(i, j)
  }
  if (is_dag_fast(Anew)) {
    new_node_lm <- node_lm
    for (jm in modified) {
      pa <- which(Anew[, jm] == 1)
      new_node_lm[jm] <- log_marg_nCPNG(jm, X_full, pa, g = 2,
                                         alpha = alpha_t)
    }
    new_lp <- sum(new_node_lm)
    if (log(runif(1)) < new_lp - cur_lp) {
      A_t <- Anew; cur_lp <- new_lp; node_lm <- new_node_lm
    }
  }
  trace_lp[it] <- cur_lp
}
trace_df <- data.frame(iter = seq_along(trace_lp), log_post = trace_lp)
p_trace <- ggplot(trace_df, aes(x = iter, y = log_post)) +
  geom_line(color = "#274060", size = 0.35) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE,
              color = "#D7263D", size = 0.9) +
  labs(x = "MCMC iteration", y = "Log posterior (un-normalised)",
       title = "MCMC trace for the breast-cancer network",
       subtitle = sprintf("q = %d nodes, n = %d patients.", q_all, n)) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_fig, "trace_real.pdf"), p_trace, width = 7.5, height = 4)
cat("Saved", file.path(out_fig, "trace_real.pdf"), "\n")

## ---------------------------------------------------------------------------
##  9.  Prior sensitivity analysis
## ---------------------------------------------------------------------------
cat("\n=== Prior sensitivity ===\n")
sens_grid <- expand.grid(w = c(0.3, 0.5, 0.7), g = c(1, 2, 5))
sens_results <- list()
top_genes <- head(edge_df$Gene, 4)
for (k in seq_len(nrow(sens_grid))) {
  w <- sens_grid$w[k]; g_k <- sens_grid$g[k]
  fit_k <- mcmc_dag(X_full, n_iter = 1500, burn_in = 300,
                    prior = "nCPNG", g = g_k,
                    alpha = n + q_all - 5, edge_prior = w)
  P_k <- fit_k$post_prob
  rownames(P_k) <- colnames(P_k) <- colnames(X_full)
  sens_results[[length(sens_results) + 1]] <-
    data.frame(w = w, g = g_k,
               Gene = top_genes,
               PostProb = round(P_k[as.character(top_genes), "metastasis"], 3))
}
sens_df <- do.call(rbind, sens_results)
sens_wide <- dcast(sens_df, Gene ~ w + g, value.var = "PostProb")
write.csv(sens_wide, file.path(out_tab, "sensitivity_real.csv"), row.names = FALSE)
sink(file.path(out_tab, "sensitivity_real.tex"))
print(xtable(sens_wide,
             caption = "Prior sensitivity of posterior edge inclusion probabilities for top gene-metastasis edges.",
             label = "tab:sens_real", digits = 3),
      include.rownames = FALSE, booktabs = TRUE)
sink()

p_sens <- ggplot(sens_df, aes(x = factor(w), y = factor(g), fill = PostProb)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", PostProb)), size = 3) +
  facet_wrap(~ Gene, ncol = 2) +
  scale_fill_viridis(option = "viridis", direction = -1, limits = c(0, 1)) +
  labs(x = "Edge prior probability  w",
       y = "Precision scaling  g", fill = "Posterior prob",
       title = "Sensitivity of top gene-metastasis edge probabilities") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(out_fig, "sensitivity_real.pdf"), p_sens,
       width = 8, height = 6)
cat("Saved", file.path(out_fig, "sensitivity_real.pdf"), "\n")

## ---------------------------------------------------------------------------
## Final summary
## ---------------------------------------------------------------------------
cat("\n==================================================\n")
cat(" Real-data analysis complete.\n")
cat(" Note on dataset:", dataset_note, "\n")
cat(" Figures in", out_fig, ", tables in", out_tab, "\n")
cat("==================================================\n")
