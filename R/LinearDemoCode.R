#' Linear Model Demo Code
#'
#' This function demonstrates how to generate simulated data, fit an MLNR model,
#' and visualize the results.
#'
#' @param num_sets Number of top-level sets.
#' @param n_samples Number of observations.
#' @param distn_mean Mean for the distribution (default 0).
#' @param distn_sd Standard deviation for the distribution (default 1).
#' @param skew Skewness for the distribution (default 0.5).
#' @param sigmasq_y Variance for the output variable (default 1).
#' @param a, b Hyperparameters for the model.
#' @param penalty Type of penalty for the model (default "function").
#' @param dist Distribution type (default "mvn").
#' @param mthd Method for fitting (default "VB").
#' @param Restarts Number of restarts for the fitting algorithm (default 1).
#' @param rel_method Method for relevance (default "gp").
#' @param w_set Weighting scheme (default "avg.corr").
#' @return A list containing model results, intermediate outputs, and visualizations.
#'
#' @import GGally
#' @import GIGrvg
#' @import class
#' @import doParallel
#' @import dplyr
#' @import expm
#' @import fMultivar
#' @import foreach
#' @import ggnetwork
#' @import ggplot2
#' @import gtools
#' @import infotheo
#' @import invgamma
#' @import laGP
#' @import mvtnorm
#' @import network
#' @import parallel
#' @import rpart
#' @import sna
#' @import tidyr
#' @import tidyverse
#' @importFrom Matrix forceSymmetric
#' @importFrom doParallel registerDoParallel
#' @importFrom ClusterR GMM
#' @importFrom ClusterR predict_GMM
#' @importFrom ald rALD dALD
#' @importFrom plgp covar
#' @importFrom plgp covar.sep
#' @importFrom utils combn
#' @export
demo_code_Linear <- function(num_sets = 5, n_samples = 100, mlnr_rho = 1, distn_mean = 0, distn_sd = 1, skew = 0.5,
                                sigmasq_y = 1, a = 1, b = 1, penalty = "function", dist = "mvn", mthd = "VB",
                                Restarts = 1, rel_method = "gp", w_set = "avg.corr") {

  # Generate the simulated data
  dat = GenSimDatLinear("mvn", n_samples, num_sets, distn_mean = distn_mean, distn_sd = distn_sd, skew = skew)
  dat_pred = GenSimDatLinear("mvn", 3, num_sets, distn_mean = distn_mean, distn_sd = distn_sd, skew = skew)
  y_test = dat_pred[1:(dim(dat_pred)[1]-1), 1]
  dat_pred = dat_pred[, 2:dim(dat_pred)[2]]  # Removes response from testing data

  # Fit model
  model = MLNR(dat, num_sets, mlnr_rho = 1, skipper = 300, smpl.sz = 2, N_norm = 2000, level_1_connected = 0,
               sigmasq_y = sigmasq_y, a = a, b = b, ald_p = 0.5, n0 = 1, s0 = 1, pi_0 = 0.5,
               a_al = 0.5, b_al = 0.5, sigmasq_alpha = 1, penalty = penalty, dist = dist, mthd = mthd,
               Restarts = Restarts, rel_method = rel_method, w_set = w_set)

  # Create a list to store the results
  results <- list()

  # Plot Actual vs. Predicted
  plot_obj <- plot(model$yhat, model$y, main = "Actual vs. Predicted")
  results$plot_actual_vs_predicted <- plot_obj

  # Create visualization
  names_vector <- paste0("S_{", 1:num_sets, "}")
  network_viz_obj <- CreateNetworkViz(dat, model[["gamma"]], model[["all_xi"]], num_sets, set_names = names_vector)
  results$network_viz <- network_viz_obj

  # Apply prediction
  predictions <- MLNR.predict(dat_pred, model, cov_transform = "scale")
  results$predictions <- predictions

  # Store the model results and any other relevant output
  results$model <- model
  results$y_test <- y_test
  results$dat_pred <- dat_pred

  # Return the list with all results
  return(results)
}
