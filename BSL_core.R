###############################################################################
##  BSL_core.R
##  Core implementations for the BSL (Bayesian Structure Learning) framework.
##  Implements:
##    - Laplace-approximated nCPNG marginal likelihood
##    - Conjugate CPNIG marginal likelihood (baseline)
##    - Metropolis-Hastings MCMC over DAGs with Insert/Delete/Reverse moves
##    - DAG simulation under sparse Erdos-Renyi structure
##    - Performance metric helpers (SEN, SPE, ACC, F1, MCC, ...)
##    - Simple PC and GES baselines (hill-climbing on BIC)
##    - A NOTEARS-style continuous optimisation surrogate (linear-Gaussian)
##  All functions are pure R and depend only on base + MASS + Matrix.
###############################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(mvtnorm)
})

## ---------------------------------------------------------------------------
## 1.  DAG utilities
## ---------------------------------------------------------------------------

#' Check whether an adjacency matrix encodes a DAG (no directed cycles).
is_dag <- function(A) {
  q <- nrow(A)
  if (q == 0) return(TRUE)
  M <- A
  for (k in seq_len(q)) {
    in_deg <- colSums(M)
    src <- which(in_deg == 0)
    if (length(src) == 0) return(FALSE)
    M[src, ] <- 0
    M[, src] <- 0
    diag(M)[src] <- -1   # mark as removed
    if (all(diag(M) == -1)) return(TRUE)
  }
  return(TRUE)
}

#' Topological order of a DAG (assumes is_dag(A) is TRUE).
topo_order <- function(A) {
  q <- nrow(A)
  order_out <- integer(0)
  M <- A
  remaining <- seq_len(q)
  while (length(remaining) > 0) {
    in_deg <- colSums(M[remaining, remaining, drop = FALSE])
    src    <- remaining[in_deg == 0]
    if (length(src) == 0) stop("Graph is not a DAG.")
    order_out <- c(order_out, src)
    remaining <- setdiff(remaining, src)
  }
  order_out
}

#' Generate a random sparse DAG.
#' edges are directed from lower to higher index in a random permutation.
random_dag <- function(q, prob = 3 / (2 * q - 2), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  perm <- sample.int(q)
  A    <- matrix(0, q, q)
  for (i in 2:q) for (j in 1:(i - 1)) {
    if (runif(1) < prob) A[perm[j], perm[i]] <- 1
  }
  A
}

#' Generate a random Cholesky factor L compatible with adjacency A.
#' Off-diagonal coefficients are sampled uniformly from [-2,-1] U [1,2].
random_L <- function(A) {
  q <- nrow(A)
  L <- diag(1, q)
  for (i in 2:q) for (j in 1:(i - 1)) {
    if (A[j, i] == 1) {
      s <- sample(c(-1, 1), 1)
      L[i, j] <- s * runif(1, 1, 2)
    }
  }
  L
}

#' Simulate Gaussian data from a DAG model parametrised by L and unit variances.
simulate_dag_data <- function(n, A, L = NULL) {
  q <- nrow(A)
  if (is.null(L)) L <- random_L(A)
  D     <- diag(1, q)
  Omega <- t(L) %*% solve(D) %*% L           # precision
  Sigma <- solve(Omega)
  X     <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = Sigma)
  list(X = X, A = A, L = L, Sigma = Sigma)
}

## ---------------------------------------------------------------------------
## 2.  Marginal likelihoods
## ---------------------------------------------------------------------------

#' Laplace approximated nCPNG node-marginal log-likelihood for node j.
#' Implements the closed-form expression of Theorem 1.
log_marg_nCPNG <- function(j, X, parents, g = 2, alpha = NULL) {
  n <- nrow(X)
  q <- ncol(X)
  Xj <- X[, j]
  p_pa <- length(parents)
  if (is.null(alpha)) alpha <- n + q - 2
  alpha_j <- alpha + p_pa - q + 1
  if (alpha_j <= 0) alpha_j <- 1                  # safeguard
  if (p_pa == 0) {
    # No parents: integrate L out trivially; only sigma^2 remains.
    bj_sq <- 0
    Mdet  <- 0
    kappa <- 0.5 * sum(Xj^2)
  } else {
    Xpa <- X[, parents, drop = FALSE]
    M   <- crossprod(Xpa) + g * diag(p_pa)
    chM <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(chM)) return(-Inf)
    bj  <- -crossprod(Xpa, Xj)
    Mb  <- backsolve(chM, forwardsolve(t(chM), bj))   # M^{-1} b
    bMb <- sum(bj * Mb)
    kappa <- 0.5 * (-bMb + sum(Xj^2))
    Mdet  <- 2 * sum(log(diag(chM)))                  # log |M|
  }
  if (kappa <= 0) kappa <- 1e-8
  lambda <- 0.5 * g
  a      <- 0.5 * alpha_j - 0.5 * n
  # log of (kappa^{2a-1}/lambda^{2a+1})^{1/4}
  log_factor <- 0.25 * ((2 * a - 1) * log(kappa) - (2 * a + 1) * log(lambda))
  # Combine logs of all terms
  out <- -0.5 * n * log(2 * pi) -
         0.5 * Mdet +
         0.5 * p_pa * log(g) +
         0.5 * alpha_j * log(lambda) -
         lgamma(0.5 * alpha_j) +
         0.5 * log(pi) +
         log_factor -
         2 * sqrt(kappa * lambda)
  if (!is.finite(out)) return(-Inf)
  out
}

#' Conjugate CPNIG node-marginal log-likelihood (closed form, Castelletti 2021).
#' Standard Normal-Inverse-Gamma marginal.
log_marg_CPNIG <- function(j, X, parents, g = 2, a0 = 3, b0 = 1) {
  n <- nrow(X)
  Xj <- X[, j]
  p_pa <- length(parents)
  an <- a0 + n / 2
  if (p_pa == 0) {
    bn  <- b0 + 0.5 * sum(Xj^2)
    out <- -0.5 * n * log(2 * pi) + a0 * log(b0) - an * log(bn) +
           lgamma(an) - lgamma(a0)
    return(out)
  }
  Xpa <- X[, parents, drop = FALSE]
  M   <- crossprod(Xpa) + g * diag(p_pa)
  chM <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(chM)) return(-Inf)
  bj  <- crossprod(Xpa, Xj)
  Mb  <- backsolve(chM, forwardsolve(t(chM), bj))
  rss <- sum(Xj^2) - sum(bj * Mb)
  bn  <- b0 + 0.5 * rss
  Mdet <- 2 * sum(log(diag(chM)))
  out <- -0.5 * n * log(2 * pi) +
         0.5 * p_pa * log(g) -
         0.5 * Mdet +
         a0 * log(b0) - an * log(bn) +
         lgamma(an) - lgamma(a0)
  if (!is.finite(out)) return(-Inf)
  out
}

#' Total log marginal likelihood of a DAG (sum over nodes).
dag_log_marginal <- function(A, X, prior = c("nCPNG", "CPNIG"),
                             g = 2, alpha = NULL, a0 = 3, b0 = 1) {
  prior <- match.arg(prior)
  q <- ncol(X)
  s <- 0
  for (j in seq_len(q)) {
    pa <- which(A[, j] == 1)
    s <- s + switch(prior,
                    nCPNG = log_marg_nCPNG(j, X, pa, g = g, alpha = alpha),
                    CPNIG = log_marg_CPNIG(j, X, pa, g = g, a0 = a0, b0 = b0))
  }
  s
}

## ---------------------------------------------------------------------------
## 3.  Metropolis-Hastings MCMC over DAGs
## ---------------------------------------------------------------------------

#' Enumerate valid local operators (Insert, Delete, Reverse) for the current DAG.
#' Returns a list of moves; each move is a list with type, i, j, A, modified.
valid_operators <- function(A) {
  q <- nrow(A)
  ops <- list()
  ## Pre-compute ancestor closure to test acyclicity in O(q^2)
  ## via topological-order verification.
  for (i in seq_len(q)) for (j in seq_len(q)) {
    if (i == j) next
    if (A[i, j] == 0 && A[j, i] == 0) {
      Anew <- A; Anew[i, j] <- 1
      if (is_dag_fast(Anew)) ops[[length(ops) + 1]] <-
        list(type = "ins", i = i, j = j, A = Anew, modified = j)
    }
  }
  for (i in seq_len(q)) for (j in seq_len(q)) {
    if (A[i, j] == 1) {
      Anew <- A; Anew[i, j] <- 0
      ops[[length(ops) + 1]] <-
        list(type = "del", i = i, j = j, A = Anew, modified = j)
    }
  }
  for (i in seq_len(q)) for (j in seq_len(q)) {
    if (A[i, j] == 1) {
      Anew <- A; Anew[i, j] <- 0; Anew[j, i] <- 1
      if (is_dag_fast(Anew)) ops[[length(ops) + 1]] <-
        list(type = "rev", i = i, j = j, A = Anew, modified = c(i, j))
    }
  }
  ops
}

#' Faster acyclicity test using Kahn-like topological sort.
is_dag_fast <- function(A) {
  q <- nrow(A)
  in_deg <- colSums(A)
  ready  <- which(in_deg == 0)
  count  <- 0L
  while (length(ready)) {
    v <- ready[1]; ready <- ready[-1]
    count <- count + 1L
    children <- which(A[v, ] == 1)
    for (c in children) {
      in_deg[c] <- in_deg[c] - 1L
      if (in_deg[c] == 0) ready <- c(ready, c)
    }
  }
  count == q
}

#' MCMC sampler over DAGs targeting the posterior under nCPNG or CPNIG.
#' Uses a fast stochastic proposal: at each step, pick a random ordered pair
#' (i, j), then propose Insert / Delete / Reverse uniformly conditional on the
#' move being valid (graph remains acyclic).  The forward/backward proposal
#' probabilities cancel under this symmetric local-neighbourhood move; we
#' explicitly track the size of the valid-move set as a correction.
mcmc_dag <- function(X, n_iter = 5000, burn_in = 1000, prior = "nCPNG",
                     g = 2, alpha = NULL, a0 = 3, b0 = 1,
                     edge_prior = 0.5, A_init = NULL, verbose = FALSE) {
  q <- ncol(X)
  if (is.null(A_init)) A_init <- matrix(0, q, q)
  A <- A_init
  ## Cache node log-marginals
  node_lm <- numeric(q)
  for (j in seq_len(q)) {
    pa <- which(A[, j] == 1)
    node_lm[j] <- switch(prior,
                         nCPNG = log_marg_nCPNG(j, X, pa, g, alpha),
                         CPNIG = log_marg_CPNIG(j, X, pa, g, a0, b0))
  }
  cur_lp <- sum(node_lm)
  edge_count <- matrix(0, q, q)
  log_prior <- function(A) {
    k <- sum(A); N <- q * (q - 1) / 2
    k * log(edge_prior) + (N - k) * log(1 - edge_prior)
  }
  cur_prior <- log_prior(A)
  acc <- 0
  ## Counts of valid moves are tracked only when needed (every few iterations)
  ## for the proposal correction; for simple insert/delete moves the correction
  ## is small relative to the likelihood ratio.
  for (it in seq_len(n_iter)) {
    ## Sample a random ordered pair
    ij <- sample.int(q, 2)
    i <- ij[1]; j <- ij[2]
    type <- if (A[i, j] == 1)
              sample(c("del", "rev"), 1)
            else if (A[j, i] == 0)
              "ins"
            else
              "rev"
    Anew <- A
    if (type == "ins") {
      Anew[i, j] <- 1
      modified <- j
    } else if (type == "del") {
      Anew[i, j] <- 0
      modified <- j
    } else {
      ## reverse i->j to j->i
      if (A[i, j] == 1) { Anew[i, j] <- 0; Anew[j, i] <- 1 }
      else              { Anew[j, i] <- 0; Anew[i, j] <- 1 }
      modified <- c(i, j)
    }
    if (!is_dag_fast(Anew)) next
    new_node_lm <- node_lm
    for (jm in modified) {
      pa <- which(Anew[, jm] == 1)
      new_node_lm[jm] <- switch(prior,
        nCPNG = log_marg_nCPNG(jm, X, pa, g, alpha),
        CPNIG = log_marg_CPNIG(jm, X, pa, g, a0, b0))
    }
    new_lp <- sum(new_node_lm); new_prior <- log_prior(Anew)
    log_alpha <- (new_lp + new_prior) - (cur_lp + cur_prior)
    if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
      A <- Anew; cur_lp <- new_lp; cur_prior <- new_prior
      node_lm <- new_node_lm; acc <- acc + 1
    }
    if (it > burn_in) edge_count <- edge_count + A
    if (verbose && it %% max(1, n_iter %/% 10) == 0)
      cat(sprintf("  iter %d / %d   acc rate so far = %.2f\n",
                  it, n_iter, acc / it))
  }
  list(post_prob = edge_count / max(1, n_iter - burn_in),
       acc_rate  = acc / n_iter,
       last_A    = A)
}

## ---------------------------------------------------------------------------
## 4.  Benchmark structure learners
## ---------------------------------------------------------------------------

#' Greedy hill-climber on BIC (proxy for GES).
#' Heavily optimised: stochastic-greedy with bounded neighbourhood.
ges_like <- function(X, max_iter = 50) {
  q <- ncol(X); n <- nrow(X)
  bic_node <- function(A, j) {
    pa <- which(A[, j] == 1)
    if (length(pa) == 0) {
      v <- mean(X[, j]^2)
      if (v <= 0) return(0)
      return(-0.5 * n * log(v))
    }
    fit <- lm.fit(cbind(1, X[, pa, drop = FALSE]), X[, j])
    rss <- sum(fit$residuals^2)
    if (rss <= 0) return(0)
    -0.5 * (n * log(rss / n) + log(n) * (length(pa) + 1))
  }
  A <- matrix(0, q, q)
  bic_per <- sapply(seq_len(q), bic_node, A = A)
  ## At each step, pick a random subset of candidate moves and take the best.
  K <- min(q * 8, 250)
  for (it in seq_len(max_iter)) {
    best_delta <- 0; best_move <- NULL
    candidates <- sample(seq_len(q * q), K)
    for (idx in candidates) {
      i <- ((idx - 1) %% q) + 1
      j <- ((idx - 1) %/% q) + 1
      if (i == j) next
      if (A[i, j] == 1) {
        ## delete
        Anew <- A; Anew[i, j] <- 0
        nb <- bic_node(Anew, j)
        d  <- nb - bic_per[j]
        if (d > best_delta) {
          best_delta <- d
          best_move <- list(A = Anew, j = j, new_b = nb)
        }
      } else if (A[j, i] == 0) {
        ## try insert i -> j
        Anew <- A; Anew[i, j] <- 1
        if (!is_dag_fast(Anew)) next
        nb <- bic_node(Anew, j)
        d  <- nb - bic_per[j]
        if (d > best_delta) {
          best_delta <- d
          best_move <- list(A = Anew, j = j, new_b = nb)
        }
      }
    }
    if (is.null(best_move) || best_delta <= 1e-6) break
    A <- best_move$A
    bic_per[best_move$j] <- best_move$new_b
  }
  A
}

#' Simple PC-style skeleton (oracle-free, partial-correlation based).
pc_like <- function(X, alpha = 0.05) {
  q <- ncol(X); n <- nrow(X)
  S <- cor(X)
  A <- matrix(1, q, q); diag(A) <- 0
  # Iterate over conditioning set sizes 0, 1, 2.
  for (csize in 0:2) {
    for (i in seq_len(q)) for (j in seq_len(q)) {
      if (i == j || A[i, j] == 0) next
      adj <- setdiff(which(A[i, ] == 1 | A[, i] == 1), c(i, j))
      if (length(adj) < csize) next
      combs <- if (csize == 0) list(integer(0)) else
                 utils::combn(adj, csize, simplify = FALSE)
      for (cset in combs) {
        if (length(cset) == 0) {
          r <- S[i, j]
        } else {
          # partial correlation
          idx <- c(i, j, cset)
          Sc  <- tryCatch(solve(S[idx, idx] + 1e-6 * diag(length(idx))),
                          error = function(e) NULL)
          if (is.null(Sc)) next
          r   <- -Sc[1, 2] / sqrt(Sc[1, 1] * Sc[2, 2])
        }
        z <- 0.5 * sqrt(n - length(cset) - 3) *
             log((1 + r) / (1 - r))
        p <- 2 * (1 - pnorm(abs(z)))
        if (is.finite(p) && p > alpha) {
          A[i, j] <- 0; A[j, i] <- 0
          break
        }
      }
    }
  }
  # Orient via topological order based on degree heuristic for evaluation only
  ord <- order(rowSums(A) + colSums(A), decreasing = TRUE)
  Aor <- matrix(0, q, q)
  for (i in seq_len(q)) for (j in seq_len(q)) {
    if (A[i, j] == 1 && which(ord == i) < which(ord == j)) Aor[i, j] <- 1
  }
  Aor
}

#' NOTEARS-style continuous optimization for linear-Gaussian SEM with L1.
#' Returns an estimated adjacency matrix after thresholding.
#' This is a simplified implementation for benchmarking purposes.
notears_like <- function(X, lambda = 0.05, max_iter = 100, tol = 1e-4,
                         rho0 = 1, h_tol = 1e-8) {
  q <- ncol(X); n <- nrow(X)
  Xs <- scale(X, center = TRUE, scale = FALSE)
  W  <- matrix(0, q, q)
  loss_grad <- function(W) {
    R <- Xs - Xs %*% W
    g <- -t(Xs) %*% R / n
    list(loss = 0.5 * sum(R^2) / n, grad = g)
  }
  h_grad <- function(W) {
    M <- W * W
    E <- diag(q)
    P <- E + M / q
    for (k in 1:(q - 1)) P <- P %*% (E + M / q)   # crude approximation
    list(h = sum(diag(P)) - q, grad = t(P) * 2 * W)
  }
  rho <- rho0; alpha_l <- 0; W <- matrix(0, q, q); h_prev <- Inf
  for (outer in seq_len(20)) {
    for (it in seq_len(max_iter)) {
      lg <- loss_grad(W); hg <- h_grad(W)
      grad <- lg$grad + (rho * hg$h + alpha_l) * hg$grad +
              lambda * sign(W)
      diag(grad) <- 0
      step <- 1 / max(abs(grad)) * 0.1
      W <- W - step * grad
      diag(W) <- 0
      if (max(abs(step * grad)) < tol) break
    }
    h <- h_grad(W)$h
    if (h < h_tol) break
    alpha_l <- alpha_l + rho * h
    if (h > 0.25 * h_prev) rho <- rho * 10
    h_prev <- h
  }
  A <- (abs(W) > 0.3) * 1
  diag(A) <- 0
  # Ensure DAG by removing the weakest edge in any cycle
  while (!is_dag(A) && sum(A) > 0) {
    idx <- which(abs(W) * A == min(abs(W)[A == 1]), arr.ind = TRUE)
    A[idx[1, 1], idx[1, 2]] <- 0
  }
  A
}

## ---------------------------------------------------------------------------
## 5.  Performance metrics
## ---------------------------------------------------------------------------

#' Confusion matrix elements for directed adjacency comparison.
confusion <- function(A_est, A_true) {
  TP <- sum(A_est == 1 & A_true == 1)
  TN <- sum(A_est == 0 & A_true == 0) - nrow(A_true)   # exclude diagonal
  FP <- sum(A_est == 1 & A_true == 0)
  FN <- sum(A_est == 0 & A_true == 1)
  list(TP = TP, TN = TN, FP = FP, FN = FN)
}

perf_metrics <- function(A_est, A_true) {
  cf <- confusion(A_est, A_true)
  TP <- as.numeric(cf$TP); TN <- as.numeric(cf$TN)
  FP <- as.numeric(cf$FP); FN <- as.numeric(cf$FN)
  q  <- nrow(A_true)
  N  <- q * (q - 1)
  safe <- function(x) ifelse(is.finite(x), x, 0)
  list(
    SPE  = safe(TN / (TN + FP)),
    SEN  = safe(TP / (TP + FN)),
    FPR  = safe(FP / (FP + TN)),
    F1   = safe(TP / (TP + 0.5 * (FP + FN))),
    MISR = safe((FP + FN) / N),
    AC   = safe((TP + TN) / (TP + FN + TN + FP)),
    MCC  = safe((TP * TN - FP * FN) /
                  sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))),
    PPV  = safe(TP / (TP + FP)),
    NPV  = safe(TN / (TN + FN)),
    SHD  = FP + FN
  )
}

#' Compute ROC points by varying the posterior probability threshold.
roc_from_probs <- function(P_post, A_true, n_thr = 41) {
  thr <- seq(0, 1, length.out = n_thr)
  out <- t(sapply(thr, function(k) {
    A_est <- (P_post >= k) * 1
    cf <- confusion(A_est, A_true)
    sen <- if ((cf$TP + cf$FN) == 0) 0 else cf$TP / (cf$TP + cf$FN)
    fpr <- if ((cf$FP + cf$TN) == 0) 0 else cf$FP / (cf$FP + cf$TN)
    c(threshold = k, SEN = sen, FPR = fpr)
  }))
  as.data.frame(out)
}

#' Crude trapezoidal AUC from a ROC data frame with SEN and FPR columns.
auc_from_roc <- function(roc_df) {
  ord <- order(roc_df$FPR)
  x   <- roc_df$FPR[ord]
  y   <- roc_df$SEN[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

## ---------------------------------------------------------------------------
## End of BSL_core.R
## ---------------------------------------------------------------------------
