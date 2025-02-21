.onLoad <- function(libname, pkgname) {
  # Load your libraries or do other initialization here
  required_packages <- c(
    "GGally", "GIGrvg", "car", "class", "dplyr", "expm", "fMultivar",
    "ggnetwork", "ggplot2", "grpreg", "gtools", "hetGP", "infotheo",
    "invgamma", "laGP", "mvtnorm", "network", "quantreg", "robustbase",
    "rpart", "sna", "splines", "tidyr", "tidyverse"
  )

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("The package '%s' is not installed. Please install it to use this library.", pkg))
    }
  }
}
