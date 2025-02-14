#' Prediction Wrapper for MLNR
#'
#' This function predicts a continuous response using the MLNR model for the input dataset
#'
#' @param dat Data frame containing the covariates for which the response is to be predicted. Should be formatted similarly to input data to the model, with the first column dropped.
#' @param model list output of function `MLNR'containing all the saved model parameters.
#'
#' @return A vector containing predicted values of the response.
#'
#' @export
#'
#' @importFrom stats cor rexp scale
#' @importFrom mvtnorm rmvnorm
#' @importFrom infotheo mutinformation discretize
#' @importFrom Matrix forceSymmetric
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom invgamma rinvgamma
#' @importFrom expm sqrtm
#' @importFrom robustbase covMcd
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr filter select mutate arrange
#'
MLNR.predict = function(dat_pred, model, cov_transform = "scale"){

  # Creating a y
  alpha_mats = model[["alpha.mats"]]
  kmat_dfs = model[["kmats"]]
  num_pwy = model[["num_pwy"]]
  gam_mod = model[["gamma"]]
  selected_indcs = gam_mod*seq(1,num_pwy)

  f_xi = list()

  # applying the above function
  pwy_dfs = pathway_creator(dat, num_pwy, transform = "scale")
  pwy_dfs_pred = pathway_creator_pred(dat_pred, dat, num_pwy, transform = "scale")

  y_hat = rep(0, (dim(dat_pred)[1]-1))

  for(i in 1:num_pwy){
    if(i %in% selected_indcs){
      stringr_xi = paste0("xi.",1)
      f_xi[[i]] = model[[stringr_xi]]
      X = as.matrix(pwy_dfs[[i]][, f_xi[[i]]])
      XX = as.matrix(pwy_dfs_pred[[i]][, f_xi[[i]]])
      Cn = plgp::covar(XX, X,d=(4/(3*nrow(X)))^(0.2)*sqrt(1), g = 0.001)
      y_hat = y_hat + Cn%*%as.matrix(alpha_mats[[i]])
    }
  }

  y_hat = y_hat*sd(y) + mean(y)

  return(y_hat)

}
