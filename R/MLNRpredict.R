#' Prediction Wrapper for MLNR
#'
#' This function predicts a continuous response using the MLNR model for the input dataset.
#'
#' @param dat_pred Data frame containing the covariates for which the response is to be predicted. Should be formatted similarly to the input data to the model, with the first column dropped.
#' @param model List output of function `MLNR` containing all the saved model parameters.
#' @param cov_transform Character string specifying the transformation to apply to covariates. Default is "scale".
#'
#' @return A vector containing predicted values of the response.
#'
#' @import GGally
#' @import GIGrvg
#' @import class
#' @import dplyr
#' @import expm
#' @import fMultivar
#' @import ggnetwork
#' @import ggplot2
#' @import gtools
#' @import infotheo
#' @import invgamma
#' @import laGP
#' @import mvtnorm
#' @import network
#' @import rpart
#' @import sna
#' @import tidyr
#' @import tidyverse
#' @importFrom Matrix forceSymmetric
#' @importFrom ClusterR GMM
#' @importFrom ClusterR predict_GMM
#' @importFrom ald rALD dALD
#' @importFrom plgp covar
#' @importFrom plgp covar.sep
#' @importFrom utils combn
#' @export
MLNR.predict = function(dat_pred, model, cov_transform = "scale"){

  # Creating a y
  y = model[["y"]]
  dat = model[["data"]]
  alpha_mats = model[["alpha.mats"]]
  kmat_dfs = model[["kmats"]]
  num_pwy = model[["num_sets"]]
  gam_mod = model[["gamma"]]
  selected_indcs = gam_mod*seq(1,num_pwy)

  f_xi = list()

  # applying the above function
  pwy_dfs = pathway_creator(dat, num_pwy, transform = cov_transform)
  pwy_dfs_pred = pathway_creator_pred(dat_pred, dat, num_pwy, transform = cov_transform)

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
