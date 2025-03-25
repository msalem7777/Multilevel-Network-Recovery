#' MLNR
#'
#' Helper functions for Bayesian Multilevel Network Recovery
#'
#' @name MLNR
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
NULL


#' nlsep
#'
#' Computes seperable negative log-likelihood --see Gramacy 2020
#'
#' @param par Parameter vector for constructed GP.
#' @param X A dataframe or data matrix.
#' @param Y A continuous response variable as a numeric vector.
#' @return A numeric value for negative log-likelihood.
nlsep <- function(par, X, Y){
  theta <- par[1:ncol(X)]
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  return(-ll)
}

#' gradnlsep
#'
#' Computes the gradient for a seperable negative log-likelihood --see Gramacy 2020
#'
#' @param par Parameter vector for constructed GP.
#' @param X A dataframe or data matrix.
#' @param Y A continuous response variable as a numeric vector.
#' @return A numeric value for negative log-likelihood.
gradnlsep <- function(par, X, Y){
  theta <- par[1:ncol(X)]
  g <- par[ncol(X)+1]
  n <- length(Y)
  K <- covar.sep(X, d=theta, g=g)
  Ki <- solve(K)
  KiY <- Ki %*% Y
  ## loop over theta components
  dlltheta <- rep(NA, length(theta))
  for(k in 1:length(dlltheta)) {
    dotK <- K * distance(X[,k])/(theta[k]^2)
    dlltheta[k] <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
      (1/2)*sum(diag(Ki %*% dotK))
  }
  ## for g
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  return(-c(dlltheta, dllg))
}


#' log_sum_exp
#'
#' Computes the log-sum-exponential trrick.
#'
#' @param x A numeric vector.
#' @return A numeric value.
log_sum_exp <- function(x) {
  max_x <- max(x)
  sum_exp <- sum(exp(x - max_x))
  log_sum_exp <- log(sum_exp) + max_x
  return(log_sum_exp)
}


#' log_ratio_d_e
#'
#' Computes the log-sum-exp trick for two vectors simultaneously.
#'
#' @param log_d Numeric vector..
#' @param log_e Numeric vector.
#' @return A numeric value.
log_ratio_d_e <- function(log_d, log_e) {
  max_log <- max(log_d, log_e)
  log_denominator <- log_sum_exp(c(log_d - max_log, log_e - max_log)) + max_log
  log_numerator <- log_d
  log_ratio <- log_numerator - log_denominator
  return(log_ratio)
}

#' sigmoid
#'
#' Function to compute the sigmoid
#'
#' @param x A numeric vector containing the values to which we want to apply a sigmoid transformation.
#' @return a transformed version of the input.
sigmoid = function(x){
  return(1/(1+exp(-x)))
}

#' dcorr
#'
#' Function to compute density of the estimated squared-correlation statistic
#'
#' @param c2 A numeric vector containing squared correlation values.
#' @param n the number of observations.
#' @param ncp (default=0) A logical value to specify the non-centrality parameter if desired. Seeting a non-zero `ncp` switches to an F-distribution while zero is a beta distribution
#' @return a numeric value for density.
dcorr = function(c2, n, ncp=0){
  if(ncp==0){
    dens = dbeta(c2, 0.5, (n-2)/2)
  } else {
    dens = df(c2, 0.5, (n-2)/2, ncp = ((n*c2)/(1-c2)))
  }

  return(dens)
}

#' compute_correlation
#'
#' Function to compute correlation between two sets
#'
#' @param df1 A dataframe containing the first set
#' @param df2 A dataframe containing the second set
#' @return A scalar value for the average of correlations between the input dataframes
compute_correlation <- function(df1, df2) {
  cor_matrix <- cor(df1, df2)
  avg_correlation <- mean(cor_matrix)
  return(avg_correlation)
}


#' min_max_normalization_df
#'
#' Function to apply min-max normalization to all columns of a dataframe
#'
#' @param df A dataframe containing numeric values
min_max_normalization_df <- function(df) {
  # Apply min-max normalization to each column
  normalized_df <- as.data.frame(lapply(df, function(x) {
    min_val <- min(x)
    max_val <- max(x)
    normalized <- (x - min_val) / (max_val - min_val)
    return(normalized)
  }))

  return(normalized_df)
}

#' prv.alpha
#'
#' Function to extract last sampled value from a matrix of sampled values. Finds the last sampled non-NA entry for each column
#'
#' @param x A matrix containing sampled numeric values
#' @return A numeric vector of last sampled values

prv.alpha = function(x, i){

  if(is.matrix(x)==TRUE){

    lst.out = c()
    for(colmnr in 1:(dim(x)[2])){

      cnt = i

      if(is.na(x[cnt-1,colmnr])==FALSE){
        lst.out = append(lst.out, x[cnt-1,colmnr])
      } else {
        while(is.na(x[cnt-1,colmnr])==TRUE)
          cnt = cnt-1
        lst.out = append(lst.out, x[cnt-1,colmnr])
      }

    }

    return(lst.out)

  } else {

    cnt = i

    if(is.na(x[cnt-1,1])==FALSE){
      return(x[cnt-1,])
    } else {
      while(is.na(x[cnt-1,1])==TRUE)
        cnt = cnt-1
      return(x[cnt-1,])
    }

  }
}


#' pathway_creator
#'
#' Function to create a list of dataframes, one for each set.
#'
#' @param dat A dataframe containing the data in its perferred format. Response value in first column, set number in last row.
#' @param num_pw A scalar numeric for the number of sets in the data
#' @param transform (default="scale") A string value indicating the desired tranformation on the covariates. "sclae" for a scale transformation. "minmax" for min-max normalization
#' @return A list of set dataframes

pathway_creator = function(dat, num_pw, transform = "scale"){
  group_list = c(unique(unlist(dat[(dim(dat)[1]),2:(dim(dat)[2])])))
  listr = list()
  for(i in 1:num_pw){
    if(transform=="scale"){
      listr[[i]] = as.data.frame(scale(dat[1:(dim(dat)[1]-1),c((dat[dim(dat)[1],])==group_list[i])]))
    } else if(transform == "minmax"){
      listr[[i]] = min_max_normalization_df(dat[1:(dim(dat)[1]-1),c((dat[dim(dat)[1],])==group_list[i])])
    } else {
      listr[[i]] = dat[1:(dim(dat)[1]-1),c((dat[dim(dat)[1],])==group_list[i])]
    }
  }

  pwy_dfs = listr

  return(pwy_dfs)
}

#' corr_mat_creator
#'
#' Function to create a list of correlation dataframes each holding correlations between elements of a set, one for each set.
#'
#' @param df_list A list of dataframes of sets
#' @param num_pw A scalar numeric for the number of sets in the data
#' @return A list of set dataframes

corr_mat_creator = function(df_list, num_pw){
  listr = list()
  for(i in 1:num_pw){
    listr[[i]] = as.data.frame(cor(df_list[[i]]))
    diag(listr[[i]]) = 0
  }

  corr_mats = listr

  return(corr_mats)
}

#' kmat_creator
#'
#' Function to create a list of kernel matrices dataframes each holding the Gaussian kernel matrix from covariates of a set, one for each set.
#'
#' @param df_list A list of dataframes of sets
#' @param num_pw A scalar numeric for the number of sets in the data
#' @return A list of set dataframes

kmat_creator = function(df_list, num_pw, mlnr_rho){
  listr = list()
  for(i in 1:num_pw){
      # listr[[i]] = listr[[i]] + hetGP::cov_gen(as.matrix(df_list[[i]]), theta=(4/(3*nrow(df_list[[i]])))^(0.2)*sqrt(1), type = "Gaussian")
      # listr[[i]] =  plgp::covar(as.matrix(df_list[[i]]),d=(4/(3*nrow(df_list[[i]])))^(0.2)*sqrt(1), g = 0.00001)
      listr[[i]] =  plgp::covar(as.matrix(df_list[[i]]),d = mlnr_rho, g = 0.00001)

      if (all(listr[[i]] == 1)) {
        listr[[i]] = matrix(0, nrow = nrow(listr[[i]]), ncol = ncol(listr[[i]]))
      }
  }

  return(listr)
}

#' bigB_creator
#'
#' Function to extract diagonals from list of kernel matrices dataframes for asymmetric LaPlace B parameter
#'
#' @param kmat_dfs A list of dataframes of Gaussian kernel matrices for sets
#' @param num_pw A scalar numeric for the number of sets in the data
#' @return A list of lists, each containing the diagonal entries of a kernel matrix

bigB_creator = function(kmat_dfs, num_pw){
  listr = list()
  for(i in 1:num_pw){
    listr[[i]] = diag(dim(kmat_dfs[[i]])[1])
  }

  ald_bigB_inv = listr

  return(ald_bigB_inv)
}

#' alpha_creator
#'
#' Function to solve for initial kernel weights given kernel matrices and initialize a placeholder list of dataframes for kernel weights
#'
#' @param kmat_dfs A list of dataframes of Gaussian kernel matrices for sets
#' @param N A scalar numeric for the number of iterations in the Gibbs sampler
#' @param distn A string for the choice of distribution "mvn" for multivariate normal (default). "ald" for asymmetric LaPlace.
#' @return A list of matrices, each containing the initial kernel weights

alpha_creator <- function(y, kmat_dfs, alpha_prior_V, N, dist = "mvn",
                          sigmasq = NULL, ald_tau = NULL,
                          ald_theta = NULL, ald_z_vec = NULL,
                          ald_z_mat = NULL, ald_bigB_inv = NULL){
  num_pw = length(kmat_dfs)
  listr = list()
  for(i in 1:num_pw){
    listr[[i]] = matrix(NA, nrow = N, ncol = nrow(kmat_dfs[[i]]))
    # initializing alphas
    if(dist == "mvn"){

      alph_var = sigmasq[1,1]*solve(sigmasq[1,1]*solve(alpha_prior_V[[i]])+t(kmat_dfs[[i]])%*%(kmat_dfs[[i]]))
      alph_var = as.matrix(Matrix::forceSymmetric(alph_var))
      alph_mean = (1/sigmasq[1,1])*alph_var%*%t(kmat_dfs[[i]])%*%(y - mean(y))

      listr[[i]][1,] = rmvnorm(1, alph_mean, alph_var)

    } else if(dist == "ald"){
      # Check if all required arguments are provided
      if (is.null(ald_tau) || is.null(ald_theta) || is.null(ald_z_vec) ||
          is.null(ald_z_mat) || is.null(ald_bigB_inv)) {
        stop("Error: Missing required arguments for 'ald' distribution.")
      }
      ald_bigB_inv[[i]] = t(kmat_dfs[[i]])%*%solve(ald_z_mat)%*%kmat_dfs[[i]]/(sqrt(sigmasq[1,])*ald_tau^2) + solve(alpha_prior_V[[i]])
      ald_bigB = solve(ald_bigB_inv[[i]])
      alph_mean = ald_bigB%*%(kmat_dfs[[i]]%*%(solve(ald_z_mat))%*%as.matrix(y-mean(y)-ald_theta*ald_z_vec)/(sqrt(sigmasq[1,])*ald_tau^2))
      alph_var = as.matrix(Matrix::forceSymmetric(ald_bigB))
      listr[[i]][1,] = rmvnorm(1, alph_mean,alph_var)
    }
  }

  return(listr)
}

#' prob_notnull_creator
#'
#' Function to initialize set probabilities and create placeholders for computed probabilities
#'
#' @param pwy_dfs A list of dataframes of sets
#' @param N A scalar numeric for the number of iterations in the Gibbs sampler
#' @param init A numeric value for the initial probability of inclusion of sets. Default is 1
#' @return A list of lists containing the initial inclusion probability for each set

prob_notnull_creator = function(pwy_dfs, N, init = 1){
  num_pw = length(pwy_dfs)
  listr = list()
  for(i in 1:num_pw){
    listr[[i]] = rep(init,N)
  }

  prob_not_null_vecs = listr

  return(prob_not_null_vecs)
}


#' curr_xi_creator
#'
#' Function to initialize element inclusion variables and create placeholders for inclusion variables
#'
#' @param pwy_dfs A list of dataframes of sets
#' @param N A scalar numeric for the number of iterations in the Gibbs sampler
#' @return A list of matrices containing the initial inclusion variables for each set's elements

curr_xi_creator = function(pwy_dfs, N){
  listr = list()
  num_pw = length(pwy_dfs)
  for(i in 1:num_pw){
    listr[[i]] = as.matrix(rbind(sample(c(1), dim(pwy_dfs[[i]])[2], replace = TRUE),matrix(NA, nrow=N-1, ncol=dim(pwy_dfs[[i]])[2])))
  }

  curr_xi_dfs = listr

  return(curr_xi_dfs)
}

#' pathway_creator_pred
#'
#' Function to create a list of dataframes, one for each set, using the scaling of the training data
#'
#' @param dat A dataframe containing the data in its perferred format. Response value in first column, set number in last row.
#' @param num_pw A scalar numeric for the number of sets in the data
#' @param transform (default="scale") A string value indicating the desired tranformation on the covariates. "sclae" for a scale transformation. "minmax" for min-max normalization
#' @param mu mean of training covariates
#' @param sigma standard deviation of training covariates
#' @param xmin minimum of training covariates
#' @param xmax maximum of training covariates
#' @return A list of set dataframes

pathway_creator_pred = function(dat_pred, dat, num_pw, transform){
  group_list = c(unique(unlist(dat[(dim(dat)[1]),2:(dim(dat)[2])])))
  listr = list()
  for(i in 1:num_pw){
    if(transform == "scale"){
      # Extract columns based on the group
      plc_hld_pred = dat_pred[1:(dim(dat_pred)[1]-1), c((dat_pred[dim(dat_pred)[1], ]) == group_list[i])]
      plc_hld_dat = dat[1:(dim(dat)[1]-1), c((dat[dim(dat)[1], ]) == group_list[i])]

      mu = colMeans(plc_hld_dat, na.rm = TRUE)
      sigma = apply(plc_hld_dat, 2, sd, na.rm = TRUE)

      # Scale dat_pred using mu and sigma from dat
      listr[[i]] = sweep(plc_hld_pred, 2, mu, "-") / sigma
    } else if(transform == "minmax"){
      # Perform Min-Max normalization using dat's min and max
      plc_hld_pred = dat_pred[1:(dim(dat_pred)[1]-1), c((dat_pred[dim(dat_pred)[1], ]) == group_list[i])]
      plc_hld_dat = dat[1:(dim(dat)[1]-1), c((dat[dim(dat)[1], ]) == group_list[i])]

      xmin = apply(plc_hld_dat, 2, min, na.rm = TRUE)
      xmax = apply(plc_hld_dat, 2, max, na.rm = TRUE)

      # Apply min-max normalization
      listr[[i]] = sweep(plc_hld_pred, 2, xmin, "-") / (xmax - xmin)
    } else {
      # No transformation, directly use the data
      listr[[i]] = dat_pred[1:(dim(dat_pred)[1]-1), c((dat[dim(dat)[1], ]) == group_list[i])]
    }
  }
  return(listr)
}
