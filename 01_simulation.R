# ======================================================================
# Monte Carlo Simulation for Evaluating the Shortfall Probability Estimator
# ======================================================================

library(frontier)
library(dplyr)
library(boot)
library(parallel)

set.seed(12345)

# ----------------------------------------------------------------------
# Scenario Definitions
# ----------------------------------------------------------------------

scenarios <- list(
  Low = list(sigma_u = 0.08, sigma_v = 0.040, delta = 0.15),
  Medium = list(sigma_u = 0.15, sigma_v = 0.133, delta = 0.30),
  High = list(sigma_u = 0.25, sigma_v = 0.332, delta = 0.40)
)

n_values <- c(50, 100, 200, 500)
R <- 1000
B <- 999

# ----------------------------------------------------------------------
# Helper Function: Shortfall Probability
# ----------------------------------------------------------------------

compute_sp <- function(delta, mu_u, sigma_v) {
  1 - pnorm((delta - mu_u) / sigma_v)
}

# ----------------------------------------------------------------------
# Core Simulation Function
# ----------------------------------------------------------------------

run_simulation <- function(sigma_u, sigma_v, delta, n, R, B) {
  
  mu_u_true <- sigma_u * sqrt(2 / pi)
  sp_true <- compute_sp(delta, mu_u_true, sigma_v)
  
  sp_est <- numeric(R)
  cover_delta <- logical(R)
  cover_boot <- logical(R)
  
  for (r in 1:R) {
    
    # Generate data
    x <- runif(n, 0, 10)
    v <- rnorm(n, 0, sigma_v)
    u <- abs(rnorm(n, 0, sigma_u))
    y <- 1 + 0.5 * x + v - u
    d <- data.frame(y = y, x = x)
    
    # Fit SFA model
    m <- tryCatch(
      sfa(y ~ x, data = d, ineffDecrease = TRUE, truncNorm = FALSE),
      error = function(e) NULL
    )
    if (is.null(m)) next
    
    coef_m <- summary(m)$Coef
    sv <- exp(coef_m["sigmaSqV", "Estimate"] / 2)
    su <- exp(coef_m["sigmaSqU", "Estimate"] / 2)
    mu_hat <- su * sqrt(2 / pi)
    
    sp_hat <- compute_sp(delta, mu_hat, sv)
    sp_est[r] <- sp_hat
    
    # Delta method (approximate)
    se_delta <- sqrt(
      (1 / sv^2) * (su^2 / (2 * pi)) +
        ((delta - mu_hat)^2 / (2 * sv^4)) * sv^2
    )
    ci_l <- sp_hat - 1.96 * se_delta
    ci_u <- sp_hat + 1.96 * se_delta
    cover_delta[r] <- (ci_l <= sp_true & ci_u >= sp_true)
    
    # Bootstrap
    boot_sp <- function(data, idx, delta) {
      bd <- data[idx, ]
      bm <- tryCatch(
        sfa(y ~ x, data = bd, ineffDecrease = TRUE, truncNorm = FALSE),
        error = function(e) NULL
      )
      if (is.null(bm)) return(NA)
      bc <- summary(bm)$Coef
      bsv <- exp(bc["sigmaSqV", "Estimate"] / 2)
      bsu <- exp(bc["sigmaSqU", "Estimate"] / 2)
      bmu <- bsu * sqrt(2 / pi)
      compute_sp(delta, bmu, bsv)
    }
    
    boot_res <- tryCatch(
      boot(d, boot_sp, R = B, delta = delta),
      error = function(e) NULL
    )
    if (!is.null(boot_res)) {
      bvals <- boot_res$t[!is.na(boot_res$t)]
      if (length(bvals) > 0) {
        ci_lb <- quantile(bvals, 0.025)
        ci_ub <- quantile(bvals, 0.975)
        cover_boot[r] <- (ci_lb <= sp_true & ci_ub >= sp_true)
      }
    }
  }
  
  # Remove NA values
  sp_est <- sp_est[!is.na(sp_est)]
  cover_delta <- cover_delta[!is.na(cover_delta)]
  cover_boot <- cover_boot[!is.na(cover_boot)]
  
  data.frame(
    true_sp = sp_true,
    bias = mean(sp_est) - sp_true,
    rmse = sqrt(mean((sp_est - sp_true)^2)),
    cover_delta = mean(cover_delta),
    cover_boot = mean(cover_boot),
    n_eff = length(sp_est)
  )
}

# ----------------------------------------------------------------------
# Run Simulations
# ----------------------------------------------------------------------

cat("\n========================================\n")
cat("MONTE CARLO SIMULATION\n")
cat("========================================\n\n")

results_list <- list()

for (sc_name in names(scenarios)) {
  sc <- scenarios[[sc_name]]
  cat(sc_name, "scenario:\n")
  for (n in n_values) {
    cat("  n =", n, "... ")
    res <- run_simulation(sc$sigma_u, sc$sigma_v, sc$delta, n, R, B)
    results_list[[paste(sc_name, n, sep = "_")]] <- res
    cat("done\n")
  }
}

# ----------------------------------------------------------------------
# Build Results Table
# ----------------------------------------------------------------------

sim_table <- do.call(rbind, lapply(names(results_list), function(id) {
  parts <- strsplit(id, "_")[[1]]
  sc <- parts[1]
  n <- as.integer(parts[2])
  r <- results_list[[id]]
  data.frame(
    Scenario = sc,
    n = n,
    Delta_Method = round(r$cover_delta, 3),
    Bootstrap = round(r$cover_boot, 3),
    Bias = round(r$bias, 4),
    RMSE = round(r$rmse, 4)
  )
}))

print(sim_table, row.names = FALSE)

# ----------------------------------------------------------------------
# Save Results
# ----------------------------------------------------------------------

write.csv(sim_table, "output/simulation_results.csv", row.names = FALSE)

cat("\nSimulation results saved to output/simulation_results.csv\n")
cat("All simulations completed successfully.\n")
