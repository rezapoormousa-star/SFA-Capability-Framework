# ======================================================================
# AutoMPG Data Analysis
# Data source: UCI Machine Learning Repository
# ======================================================================

library(frontier)
library(dplyr)
library(boot)

set.seed(12345)

# ----------------------------------------------------------------------
# Load Data Directly from UCI Repository
# ----------------------------------------------------------------------

cat("\n========================================\n")
cat("AUTO MPG DATA ANALYSIS\n")
cat("========================================\n")

auto_mpg_url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"

column_names <- c("mpg", "cylinders", "displacement", "horsepower", 
                  "weight", "acceleration", "model_year", "origin", "car_name")

automp_data <- read.table(auto_mpg_url, header = FALSE, 
                          col.names = column_names, 
                          na.strings = "?", 
                          stringsAsFactors = FALSE)

automp_data <- na.omit(automp_data)

cat("Observations:", nrow(automp_data), "\n")

# Prepare variables
automp_data <- automp_data %>%
  mutate(
    log_mpg = log(mpg),
    log_weight = log(weight),
    log_disp = log(displacement),
    log_hp = log(horsepower),
    age = 82 - model_year,
    origin_Japan = ifelse(origin == 3, 1, 0),
    origin_Europe = ifelse(origin == 2, 1, 0)
  )

# ----------------------------------------------------------------------
# Fit SFA Model
# ----------------------------------------------------------------------

model_automp <- sfa(
  log_mpg ~ log_weight + log_disp + log_hp,
  data = automp_data,
  ineffDecrease = TRUE,
  truncNorm = FALSE,
  logDepVar = TRUE,
  z = ~ age + origin_Japan + origin_Europe
)

cat("\n--- Model Summary ---\n")
print(summary(model_automp))

# ----------------------------------------------------------------------
# Extract Parameters
# ----------------------------------------------------------------------

coef_automp <- summary(model_automp)$Coef

sigma_v <- exp(coef_automp["sigmaSqV", "Estimate"] / 2)
sigma_u_base <- exp(coef_automp["sigmaSqU", "Estimate"] / 2)
gamma_age <- coef_automp["z_age", "Estimate"]
gamma_japan <- coef_automp["z_origin_Japan", "Estimate"]
gamma_europe <- coef_automp["z_origin_Europe", "Estimate"]

mu_u_base <- sigma_u_base * sqrt(2 / pi)

# ----------------------------------------------------------------------
# Compute SP(0.15) and OCI(0.15)
# ----------------------------------------------------------------------

compute_sp <- function(delta, mu_u, sigma_v) {
  1 - pnorm((delta - mu_u) / sigma_v)
}

delta_automp <- 0.15
SP_automp <- compute_sp(delta_automp, mu_u_base, sigma_v)
OCI_automp <- (delta_automp - mu_u_base) / sigma_v

cat("\n--- Results for delta = 0.15 ---\n")
cat(sprintf("SP(0.15) = %.3f\n", SP_automp))
cat(sprintf("OCI(0.15) = %.3f\n", OCI_automp))

# ----------------------------------------------------------------------
# SP by Origin
# ----------------------------------------------------------------------

sigma_u_japan <- sigma_u_base * exp(gamma_japan)
sigma_u_europe <- sigma_u_base * exp(gamma_europe)

mu_u_japan <- sigma_u_japan * sqrt(2 / pi)
mu_u_europe <- sigma_u_europe * sqrt(2 / pi)

SP_japan <- compute_sp(delta_automp, mu_u_japan, sigma_v)
SP_europe <- compute_sp(delta_automp, mu_u_europe, sigma_v)
SP_usa <- SP_automp

OCI_japan <- (delta_automp - mu_u_japan) / sigma_v
OCI_europe <- (delta_automp - mu_u_europe) / sigma_v
OCI_usa <- OCI_automp

cat("\n--- SP by Origin ---\n")
cat(sprintf("American:  SP = %.3f, OCI = %.3f, N = %d\n",
            SP_usa, OCI_usa, sum(automp_data$origin == 1)))
cat(sprintf("Japanese:  SP = %.3f, OCI = %.3f, N = %d\n",
            SP_japan, OCI_japan, sum(automp_data$origin == 3)))
cat(sprintf("European:  SP = %.3f, OCI = %.3f, N = %d\n",
            SP_europe, OCI_europe, sum(automp_data$origin == 2)))

# ----------------------------------------------------------------------
# Sensitivity Analysis
# ----------------------------------------------------------------------

cat("\n--- Sensitivity Analysis ---\n")

delta_values <- c(0.05, 0.10, 0.15, 0.20, 0.25)

sensitivity_table <- data.frame(
  delta = delta_values,
  SP = sapply(delta_values, function(d) compute_sp(d, mu_u_base, sigma_v)),
  OCI = sapply(delta_values, function(d) (d - mu_u_base) / sigma_v)
)

print(sensitivity_table, row.names = FALSE)

# ----------------------------------------------------------------------
# Wald Test for Constant Mean Inefficiency
# ----------------------------------------------------------------------

cat("\n--- Wald Test for Stability of E[u_i] ---\n")

u_hat <- efficiencies(model_automp, asInEff = TRUE)
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

write.csv(sensitivity_table, "output/automp_results.csv", row.names = FALSE)

cat("\nAutoMPG results saved to output/automp_results.csv\n")
cat("AutoMPG analysis completed successfully.\n")
