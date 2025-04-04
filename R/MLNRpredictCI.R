#' MLNR.predictCI
#'
#' Prediction Wrapper with Credible Interval for MLNR
#'
#' This function predicts a continuous response using the MLNR model for the input dataset.
#'
#' @param dat_pred Data frame containing the covariates for which the response is to be predicted. Should be formatted similarly to the input data to the model, with the first column dropped.
#' @param model List output of function `MLNR` containing all the saved model parameters.
#' @param cov_transform Character string specifying the transformation to apply to covariates. Default is "scale".
#' @param scale_up Boolean. Uses mean and variance of training data to Upscale response (default = FALSE)
#' @param CI_level The level of the credible interval. Between 0 and 1 (default = 0.95)
#' @param sampler Number of samples for CI estimation (default = 1000)
#'
#' @return A dataframe containing containing predicted values of the response, lower bound and upper bound of CI
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
MLNR.predictCI = function(dat_pred, model, interval="credible", cov_transform = "none", scale_up=FALSE, CI_level=0.95, sampler = 1000){

  # Creating a y
  y = model[["y"]]
  dat = model[["data"]]
  alpha_mats = model[["alpha.mats"]]
  kmat_dfs = model[["kmats"]]
  num_pwy = model[["num_sets"]]
  gam_mod = model[["gamma"]]
  mlnr_rho = model[["mlnr_rho"]]
  mu_alpha_post = model[["mu_alpha"]]
  sigma_alpha_post = model[["sigma_alpha"]]
  sigmasq_eps = model[["sigmasq_eps"]]
  selected_indcs = gam_mod*seq(1,num_pwy)
  selected_indcs = selected_indcs[selected_indcs!=0]

  # applying the above function
  pwy_dfs = pathway_creator(dat[2:dim(dat)[2]], num_pwy, transform = cov_transform)
  pwy_dfs_pred = pathway_creator_pred(dat_pred, dat[2:dim(dat)[2]], num_pwy, transform = cov_transform)

  y_hat = rep(0, (dim(dat_pred)[1]-1))
  y_hat_samples = matrix(0, nrow = length(y_hat), ncol = sampler)

  cntr = 1
  for(i in selected_indcs){
    stringr_xi = paste0("xi.",i)
    xi_indcs = model[[stringr_xi]]*seq(1,ncol(pwy_dfs[[i]]))
    X = as.matrix(pwy_dfs[[i]][, xi_indcs])
    XX = as.matrix(pwy_dfs_pred[[i]][, xi_indcs])
    Cn = plgp::covar(XX, X,d = mlnr_rho, g = 0.001)
    if(sum(model[[stringr_xi]])==0){
      y_hat = y_hat + 0
      # eps_sampler = rmvnorm(sampler, rep(0,length(y_hat)), sigmasq_eps*diag(length(y_hat)))
      # y_hat_samples = y_hat_samples + t(eps_sampler) # Compute y_hat for all 1000 samples
    } else {
      y_hat = y_hat + Cn%*%as.matrix(alpha_mats[[cntr]])
      alpha_sampler = rmvnorm(sampler, mu_alpha_post[[cntr]], sigma_alpha_post[[cntr]])
      # eps_sampler = rmvnorm(sampler, rep(0,length(y_hat)), sigmasq_eps*diag(length(y_hat)))
      y_hat_samples = y_hat_samples + Cn %*% t(alpha_sampler) #+ t(eps_sampler) # Compute y_hat for all 1000 samples
    }

    cntr = cntr + 1
  }

  pL = (1 - CI_level)/2
  pU = 1 - pL

  if(interval == "predictive"){
    y_hat_L = apply(y_hat_samples, 1, quantile, probs = pL) + qnorm(pL,0, sd = sqrt(sigmasq_eps))
    y_hat_U = apply(y_hat_samples, 1, quantile, probs = pU) + qnorm(pU,0, sd = sqrt(sigmasq_eps))
  } else if(interval == "credible"){
    y_hat_L = apply(y_hat_samples, 1, quantile, probs = pL)
    y_hat_U = apply(y_hat_samples, 1, quantile, probs = pU)
  }


  if(scale_up == TRUE){
    y_hat = y_hat*sd(y) + mean(y)
    y_hat_L = y_hat_L*sd(y) + mean(y)
    y_hat_U = y_hat_U*sd(y) + mean(y)
  }

  df <- data.frame(lower = y_hat_L, pred = y_hat, upper = y_hat_U)

  return(df)



}
