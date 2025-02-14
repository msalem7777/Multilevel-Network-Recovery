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
MLNR.predict = function(dat, model){

  dat_pred = dat[-c((length(y)))]

  # Creating a y
  y = model[["y"]]
  X = model[["X"]]
  alpha_mats = model[["alpha.mats"]]
  kmat_dfs = model[["kmats"]]
  num_pwy = model[["num_pwy"]]

  # applying the above function
  pwy_dfs = pathway_creator(dat[1:(dim(dat)[1]-1),], num_pwy)

  y_hat = rep(0, (dim(dat)[1]-1))

  for(i in 1:num_pwy){
    XX = (pwy_dfs[[i]]-mean(X[[i]]))/sd(X[[i]])
    Cn = plgp::covar(XX, X[i],d=(4/(3*nrow(X[[i]])))^(0.2)*sqrt(1))
    y_hat =+ Cn%*%alpha_mats[[i]]
  }

  return(y_hat)

}
