# ======================================================================
# Install Required Packages (Run Once)
# ======================================================================

packages <- c("frontier", "dplyr", "boot", "parallel")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

cat("All required packages are installed.\n")
