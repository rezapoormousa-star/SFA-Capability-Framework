# ======================================================================
# Steel Industry Data Analysis
# Data source: Sathishkumar et al. (2021), Building Research & Information
# ======================================================================

library(frontier)
library(dplyr)
library(boot)
library(parallel)

set.seed(12345)

# ----------------------------------------------------------------------
# Load Data Directly from Public Repository
# ----------------------------------------------------------------------

cat("\n========================================\n")
cat("STEEL INDUSTRY DATA ANALYSIS\n")
cat("========================================\n")
steel_url <- "https://raw.githubusercontent.com/Hokfu/Energy-Consumption-Model/main/Steel_industry_data.csv"
steel_data <- read.csv(steel_url, header = TRUE)

cat("Observations:", nrow(steel_data), "\n")

# Prepare log variables
steel_data$log_energy <- log(steel_data$energy)
steel_data$log_reactive <- log(steel_data$reactive_power)

# ----------------------------------------------------------------------
# Fit SFA Model
# ----------------------------------------------------------------------

model_steel <- sfa(
  log_energy ~ log_reactive + power_factor + load_type,
  data = steel_data,
  ineffDecrease = TRUE,
  truncNorm = FALSE
)

cat("\n--- Model Summary ---\n")
print(summary(model_steel))

# ----------------------------------------------------------------------
# Extract Parameters
# ----------------------------------------------------------------------

coef_steel <- summary(model_steel)$Coef

beta0 <- coef_steel["(Intercept)", "Estimate"]
beta1 <- coef_steel["log_reactive", "Estimate"]
beta2 <- coef_steel["power_factor", "Estimate"]
beta3 <- coef_steel["load_type", "Estimate"]

sigma_v <- exp(coef_steel["sigmaSqV", "Estimate"] / 2)
sigma_u <- exp(coef_steel["sigmaSqU", "Estimate"] / 2)

se_beta0 <- coef_steel["(Intercept)", "Std. Error"]
se_beta1 <- coef_steel["log_reactive", "Std. Error"]
se_beta2 <- coef_steel["power_factor", "Std. Error"]
se_beta3 <- coef_steel["load_type", "Std. Error"]

se_sigma_v <- sigma_v * coef_steel["sigmaSqV", "Std. Error"] / 2
se_sigma_u <- sigma_u * coef_steel["sigmaSqU", "Std. Error"] / 2

mu_u <- sigma_u * sqrt(2 / pi)

cat("\n--- Parameter Estimates ---\n")
cat(sprintf("beta_0 = %8.3f (SE = %6.4f)\n", beta0, se_beta0))
cat(sprintf("beta_1 = %8.3f (SE = %6.4f)\n", beta1, se_beta1))
cat(sprintf("beta_2 = %8.3f (SE = %6.4f)\n", beta2, se_beta2))
cat(sprintf("beta_3 = %8.3f (SE = %6.4f)\n", beta3, se_beta3))
cat(sprintf("sigma_v = %8.3f (SE = %6.4f)\n", sigma_v, se_sigma_v))
cat(sprintf("sigma_u = %8.3f (SE = %6.4f)\n", sigma_u, se_sigma_u))
cat(sprintf("mu_u = %8.3f\n", mu_u))

# ----------------------------------------------------------------------
# Compute SP(0.10) and OCI(0.10)
# ----------------------------------------------------------------------

compute_sp <- function(delta, mu_u, sigma_v) {
  1 - pnorm((delta - mu_u) / sigma_v)
}

delta_steel <- 0.10
SP_steel <- compute_sp(delta_steel, mu_u, sigma_v)
OCI_steel <- (delta_steel - mu_u) / sigma_v

cat("\n--- Results for delta = 0.10 ---\n")
cat(sprintf("SP(0.10) = %.3f\n", SP_steel))
cat(sprintf("OCI(0.10) = %.3f\n", OCI_steel))

# ----------------------------------------------------------------------
# Bootstrap for SP(0.10)
# ----------------------------------------------------------------------

cat("\n--- Bootstrap (B = 999) ---\n")

boot_sp_steel <- function(data, indices, delta) {
  boot_data <- data[indices, ]
  m <- tryCatch({
    sfa(log_energy ~ log_reactive + power_factor + load_type,
        data = boot_data, ineffDecrease = TRUE, truncNorm = FALSE)
  }, error = function(e) NULL)
  if (is.null(m)) return(NA)
  cm <- summary(m)$Coef
  sv <- exp(cm["sigmaSqV", "Estimate"] / 2)
  su <- exp(cm["sigmaSqU", "Estimate"] / 2)
  mu <- su * sqrt(2 / pi)
  return(compute_sp(delta, mu, sv))
}

set.seed(54321)
boot_steel <- boot(
  data = steel_data,
  statistic = boot_sp_steel,
  R = 999,
  delta = 0.10,
  parallel = "multicore",
  ncpus = max(1, detectCores() - 1)
)

boot_vals <- boot_steel$t[!is.na(boot_steel$t)]
ci_boot <- quantile(boot_vals, probs = c(0.025, 0.975))

cat(sprintf("Bootstrap mean = %.3f\n", mean(boot_vals)))
cat(sprintf("95%% CI = [%.3f, %.3f]\n", ci_boot[1], ci_boot[2]))

# ----------------------------------------------------------------------
# Sensitivity Analysis
# ----------------------------------------------------------------------

cat("\n--- Sensitivity Analysis ---\n")

delta_values <- c(0.05, 0.08, 0.10, 0.15, 0.20)

sensitivity_table <- data.frame(
  delta = delta_values,
  SP = sapply(delta_values, function(d) compute_sp(d, mu_u, sigma_v)),
  OCI = sapply(delta_values, function(d) (d - mu_u) / sigma_v)
)

print(sensitivity_table, row.names = FALSE)

# ----------------------------------------------------------------------
# Wald Test for Constant Mean Inefficiency
# ----------------------------------------------------------------------

cat("\n--- Wald Test for Stability of E[u_i] ---\n")

u_hat <- efficiencies(model_steel, asInEff = TRUE)
u_hat <- as.numeric(u_hat)
u_bar <- mean(u_hat, na.rm = TRUE)
u_var <- var(u_hat, na.rm = TRUE)

W_stat <- sum((u_hat - u_bar)^2 / u_var, na.rm = TRUE)
df <- length(u_hat) - 1
p_value <- pchisq(W_stat, df = df, lower.tail = FALSE)

cat(sprintf("W statistic = %.2f\n", W_stat))
cat(sprintf("df = %d\n", df))
cat(sprintf("p-value = %.4f\n", p_value))

if (p_value > 0.05) {
  cat("Conclusion: Fail to reject H0 (E[u_i] is stable)\n")
} else {
  cat("Conclusion: Reject H0 (E[u_i] is not stable)\n")
}

# ----------------------------------------------------------------------
# Save Results
# ----------------------------------------------------------------------

write.csv(sensitivity_table, "output/steel_results.csv", row.names = FALSE)

cat("\nSteel industry results saved to output/steel_results.csv\n")
cat("Steel industry analysis completed successfully.\n")
