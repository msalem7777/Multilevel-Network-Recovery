.onAttach <- function(libname, pkgname) {
  # Load your libraries or do other initialization here
  required_packages <- c(
    "GGally", "GIGrvg", "car", "class", "doParallel", "dplyr", "expm", "fMultivar", "foreach",
    "ggnetwork", "ggplot2", "gtools", "infotheo",
    "invgamma", "laGP", "mvtnorm", "network", "parallel", "plgp", "quantreg", "robustbase",
    "rpart", "sna", "splines", "tidyr", "tidyverse"
  )

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      warning(sprintf("The package '%s' is not installed. Please install it to use this library.", pkg))
    }
  }
}
