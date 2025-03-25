#' MLNR: Bayesian Multilevel Network Recovery
#'
#' @description Implements the MLNR model for Bayesian inference in multilevel network regression.
#'
#' @param dat Input dataset.
#' @param num_pwy Number of pathways.
#' @param skipper Number of Gibbs sampling iterations for the kernel weights between each element selection run.
#' @param smpl.sz Sample size of pathways tested at once.
#' @param reg_results If 1, the model samples final weights given selected elements and sets, and computes MSE.
#' @param N_norm Number of total iterations for normalization.
#' @param level_1_connected If 1, supernodes are connected; if 0, they are not.
#' @param sigmasq_y Initial overall variance.
#' @param a Shape parameter for the inverse gamma prior on the overall variance.
#' @param b Scale parameter for the inverse gamma prior on the overall variance.
#' @param ald_p Asymmetric Laplace Distribution (ALD) quantile parameter.
#' @param n0 ALD scale prior parameter.
#' @param s0 ALD scale prior parameter.
#' @param pi_0 Binomial probability for pathway selection.
#' @param a_al Shape parameter for the inverse gamma prior on the coefficient variance.
#' @param b_al Scale parameter for the inverse gamma prior on the coefficient variance.
#' @param sigmasq_alpha Initial variance for coefficients.
#' @param penalty Type of penalty, either "weights" (penalizes only kernel weights) or "function" (penalizes function complexity).
#' @param dist Distribution type, either "mvn" (multivariate normal) or "ald" (asymmetric Laplace).
#' @param mthd Method for inference, either "VB" (variational Bayes) or "MCMC" (Markov Chain Monte Carlo).
#' @param Restarts Number of restarts to help the Ising model escape local modes.
#' @param rel_method Relevance method, either "gp" (uses actual observations), "mi" (mutual information), or "cor" (correlation squared).
#' @param w_set Weight setting, such as "avg.cor" for using average correlation.
#' @return MLNR model results.
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
MLNR = function(dat, num_pwy, mlnr_rho = 1, skipper = 300, smpl.sz = 2, N_norm = 2000, level_1_connected = 1, sigmasq_y = 1,
                a = 1, b = 1, ald_p = 0.5, n0=1, s0=1, pi_0 = 0.5, a_al=0.5, b_al=0.5, sigmasq_alpha=1,
                penalty = "weights", dist = "mvn", mthd = "VB", Restarts = 1, rel_method = "mi", w_set="avg.corr"){

  library(foreach)
  library(doParallel)

  prcent = 1/(2*N_norm)*100

  eps = sqrt(.Machine$double.eps)

  out_list = list()

  # Num cores for parallelization
  num_cores <- min(num_pwy,(parallel::detectCores() - 1))  # Reserve one core for system processes

  # Placeholder matrices
  MLN_G_mat = as.data.frame(matrix(0, nrow = num_pwy, ncol=2))
  MLN_mat = as.data.frame(matrix(0, nrow = (dim(dat)[2]-1), ncol=1))

  MLN_vec_mse = c()

  # Creating a y
  y = as.data.frame(dat[,1])
  og_df = dat[-c((length(y))),]
  group_identifier_vector = c(dat[(length(y)),])
  group_list = paste(c(1:num_pwy))
  y = y[-c(nrow(y)),1]
  nx = length(y)

  out_list[["y"]] = y
  out_list[["X"]] = dat[1:(dim(dat)[1]-1),2:dim(dat)[2]]
  # Standardizing y
  y = scale(y)
  y =  as.numeric(y)

  # applying the above function
  pwy_dfs = pathway_creator(dat[2:dim(dat)[2]], num_pwy)

  if(rel_method == "cor"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = as.data.frame(as.numeric(cor(y,pwy_dfs[[i]])^2))
    }
    GP_status = 0
  } else if(rel_method == "mi"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = as.data.frame(matrix(NA, nrow = (dim(pwy_dfs[[i]])[2]), ncol = 1))
      for(j in 1:(dim(pwy_dfs[[i]])[2])){
        y_tilde[[i]][j,1] = infotheo::mutinformation(infotheo::discretize(y),infotheo::discretize(pwy_dfs[[i]][j]), method="emp")
      }
    }
    GP_status = 0
  } else if(rel_method == "gp"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = y
    }
    GP_status = 1
  }

  # applying the above function
  corr_mats = corr_mat_creator(pwy_dfs, num_pwy)

  # applying the above function
  kmat_dfs = kmat_creator(pwy_dfs, num_pwy, mlnr_rho)

  # applying the above function
  ald_bigB_inv = bigB_creator(kmat_dfs, num_pwy)

  # Setting up overall variance placeholder
  sigmasq = data.frame(rep(NA, N_norm))
  sigmasq[1,] = sigmasq_y

  # Setting Asymmetric LaPlace Tau
  ald_tau = sqrt(2/(ald_p*(1-ald_p)))

  # Setting Asymmetric LaPlace Theta
  ald_theta = (1-2*ald_p)/(ald_p*(1-ald_p))

  # Setting z initial
  ald_z_vec = rexp(length(y), rate=1)
  ald_z_mat = diag(ald_z_vec, nrow = length(y))

  alpha_prior_V = list()
  for(i in 1:num_pwy){
    if(penalty == "function"){
      alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*kmat_dfs[[i]]+sqrt(.Machine$double.eps)*diag(nrow = nx)))
    } else if(penalty == "weights"){
      alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*diag(nrow=nx)))
    }
  }
  # applying the above function
  if (dist == "mvn"){
    alpha_mats = alpha_creator(y, kmat_dfs, alpha_prior_V, N_norm, dist = dist, sigmasq = sigmasq)
  } else if(dist == "ald"){
    alpha_mats = alpha_creator(y, kmat_dfs, alpha_prior_V, N_norm, dist = dist, sigmasq = sigmasq, ald_tau=ald_tau, ald_theta=ald_theta, ald_z_vec=ald_z_vec, ald_z_mat=ald_z_mat, ald_bigB_inv=ald_bigB_inv)
  }

  # Setting up a placeholder for path selection variable (10 pathways)
  gamma = matrix(NA, nrow = N_norm, ncol = num_pwy)
  # Initializing path selection variable
  gamma[1,] = rep(1,num_pwy)

  # applying the above function
  prob_not_null_vecs = prob_notnull_creator(pwy_dfs, N_norm)

  # applying the above function
  curr_xi_dfs = curr_xi_creator(pwy_dfs, N_norm)

  if (w_set=="avg.corr"){
    # Compute average correlation between each pair of data frames
    W_mat <- matrix(NA, nrow = num_pwy, ncol = num_pwy)

    for (i in 1:num_pwy) {
      for (j in 1:num_pwy) {
        if (i != j) {
          correlation <- compute_correlation(pwy_dfs[[i]], pwy_dfs[[j]])
          W_mat[i, j] <- correlation
        } else {
          W_mat[i, j] <- 1  # Set diagonal elements to 1
        }
      }
    }

    W_mat = 100*W_mat^2
  }

  gam_mod = rep(1,num_pwy)
  omega = 1.0

  for (i in 2:N_norm){

    alpha_mats_k = lapply(alpha_mats, function(x) prv.alpha(x, i))

    selected_all = as.double(sample(seq(1, num_pwy), num_pwy))

    if(dist == "mvn"){

      # Gamma by K matrices
      gamma_lstr = as.list(matrix(prv.alpha(gamma, i),nrow = length(gamma[i-1,]), ncol=1))
      kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)

      for (j in 1:num_pwy){

        # sampling first path coefficients from spike slab prior
        alph_var = sigmasq[i-1,1]*solve(sigmasq[i-1,1]*solve(alpha_prior_V[[j]])+t(kmat_dfs[[j]])%*%(kmat_dfs[[j]]))
        alph_var = as.matrix(Matrix::forceSymmetric(alph_var))
        alph_mean = (1/sigmasq[i-1,1])*alph_var%*%t(kmat_dfs[[j]])%*%(y - mean(y)-omega*Reduce("+",(Map('%*%',kmat_dfs_gammad[-j],alpha_mats_k[-j]))))
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)

      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] = rinvgamma(1,length(y)/2+a/2,0.5*sum((y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)+b/2)

    } else if(dist == "ald"){

      # Gamma by K matrices
      gamma_lstr = as.list(matrix(prv.alpha(gamma,i),nrow = length(gamma[i-1,]), ncol=1))
      kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)


      for (j in 1:num_pwy){

        # sampling first path coefficients from spike slab prior
        ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(sqrt(sigmasq[i-1,])*ald_tau^2) + solve(alpha_prior_V[[j]])
        ald_bigB = solve(ald_bigB_inv[[j]])
        alph_mean = ald_bigB%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix((y- mean(y)-omega*Reduce("+",(Map('%*%',kmat_dfs_gammad[-j],alpha_mats_k[-j]))))-ald_theta*ald_z_vec)/(sqrt(sigmasq[i-1,])*ald_tau^2))
        alph_var = as.matrix(Matrix::forceSymmetric(ald_bigB))
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)
      }

      ald_delta_sq = ((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)/(ald_tau^2*sqrt(sigmasq[i-1,]))
      ald_gamma_sq = 2/sqrt(sigmasq[i-1,]) + (ald_theta^2)/(sqrt(sigmasq[i-1,])*ald_tau^2)
      ald_delta = sqrt(ald_delta_sq)
      ald_gamma = sqrt(ald_gamma_sq)

      plc_hld_z = numeric(1)
      for (k in 1:length(ald_z_vec)){
        plc_hld_z = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
        ald_z_vec[k] = mean(plc_hld_z, na.rm=TRUE)
        ald_z_mat[k,k] = ald_z_vec[k]
      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] =  (rinvgamma(1,n0/2+3*length(y)/2, s0/2+sum(ald_z_vec)+0.5*sum((solve(ald_z_mat))%*%(as.matrix(y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k)))-ald_theta*ald_z_vec)^2)/(ald_tau^2))))^2
    }

    selected = as.double(sample(seq(1, num_pwy), smpl.sz))

    if(dist == "mvn"){
      for (j in selected){

        gamma_flipr = gamma_lstr
        gamma_flipr[[j]] = 1*(gamma_lstr[[j]]==0)

        kmat_reducer = kmat_dfs[[j]]%*%alpha_mats[[j]][i,]
        if(gamma_lstr[[j]]==0){
          kmat_dfs_tstr = Map('*',kmat_dfs,gamma_flipr)
        } else if(gamma_lstr[[j]]==1){
          kmat_dfs_tstr = Map('*',kmat_dfs,gamma_lstr)
        }

        if(level_1_connected == 1){
          gamma_lstr_neg_coding = replace(unlist(gamma_lstr), unlist(gamma_lstr)==0, -1)
          delta_mat = (as.matrix(gamma_lstr_neg_coding)%*%t(as.matrix(gamma_lstr_neg_coding)))[selected, selected]
          upper_tri_delta_mat = delta_mat*upper.tri(delta_mat, diag = FALSE)
          if(gamma_flipr[[j]]==1){
            prior_ising_e = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          } else if(gamma_flipr[[j]]==0){
            prior_ising_d = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          }

          gamma_flipr_neg_coding = replace(unlist(gamma_flipr), unlist(gamma_flipr)==0, -1)
          delta_mat = (as.matrix(gamma_flipr_neg_coding)%*%t(as.matrix(gamma_flipr_neg_coding)))[selected, selected]
          upper_tri_delta_mat = delta_mat*upper.tri(delta_mat, diag = FALSE)
          if(gamma_flipr[[j]]==1){
            prior_ising_d = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          } else if(gamma_flipr[[j]]==0){
            prior_ising_e = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          }

        } else{
          prior_ising_d = 0.5
          prior_ising_e = 0.5
        }

        d = dmvnorm(y, (mean(y) +(Reduce("+",(Map('%*%',kmat_dfs_tstr[selected],alpha_mats_k[selected]))))),sigmasq[i,]*diag(length(y)))*prior_ising_d

        e = dmvnorm(y, (mean(y) +(Reduce("+",(Map('%*%',kmat_dfs_tstr[selected],alpha_mats_k[selected]))))-kmat_reducer),sigmasq[i,]*diag(length(y)))*prior_ising_e

        # Computing the probability that gamma j equals 1
        prob_not_null_vecs[[j]][i] = exp(log_ratio_d_e(log(d), log(e)))

        if(is.na(prob_not_null_vecs[[j]][i]) == TRUE){
          prob_not_null_vecs[[j]][i] = 0
        } else {}

        # Updating our value for gamma 1
        if(prob_not_null_vecs[[j]][i]>=0.5){
          gamma[i,j] = 1
        } else {gamma[i,j] = 0}
      }

    } else if(dist == "ald"){

      for (j in selected){

        gamma_flipr = gamma_lstr
        gamma_flipr[[j]] = 1*(gamma_lstr[[j]]==0)

        kmat_reducer = kmat_dfs[[j]]%*%alpha_mats[[j]][i,]
        if(gamma_lstr[[j]]==0){
          kmat_dfs_tstr = Map('*',kmat_dfs,gamma_flipr)
        } else if(gamma_lstr[[j]]==1){
          kmat_dfs_tstr = Map('*',kmat_dfs,gamma_lstr)
        }

        if(level_1_connected == 1){
          gamma_lstr_neg_coding = replace(unlist(gamma_lstr), unlist(gamma_lstr)==0, -1)
          delta_mat = (as.matrix(gamma_lstr_neg_coding)%*%t(as.matrix(gamma_lstr_neg_coding)))[selected, selected]
          upper_tri_delta_mat = delta_mat*upper.tri(delta_mat, diag = FALSE)
          if(gamma_flipr[[j]]==1){
            prior_ising_e = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          } else if(gamma_flipr[[j]]==0){
            prior_ising_d = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          }


          gamma_flipr_neg_coding = replace(unlist(gamma_flipr), unlist(gamma_flipr)==0, -1)
          delta_mat = (as.matrix(gamma_flipr_neg_coding)%*%t(as.matrix(gamma_flipr_neg_coding)))[selected, selected]
          upper_tri_delta_mat = delta_mat*upper.tri(delta_mat, diag = FALSE)
          if(gamma_flipr[[j]]==1){
            prior_ising_d = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          } else if(gamma_flipr[[j]]==0){
            prior_ising_e = exp(sum(W_mat[selected, selected]*upper_tri_delta_mat))
          }

        } else{
          prior_ising_d = 0.5
          prior_ising_e = 0.5
        }

        d = e = c()
        m1_ald = (mean(y)+Reduce("+",(Map('%*%',kmat_dfs_tstr[selected],alpha_mats_k[selected]))))
        m2_ald = (mean(y)+Reduce("+",(Map('%*%',kmat_dfs_tstr[selected],alpha_mats_k[selected])))-kmat_reducer)
        for(lp in 1:nx){
          d = c(d, dALD(y[lp], m1_ald[lp]))
          e = c(e, dALD(y[lp], m2_ald[lp]))
        }

        d = prod(d)*(prior_ising_d)

        e = prod(e)*(prior_ising_e)

        # Computing the probability that gamma j equals 1
        prob_not_null_vecs[[j]][i] = exp(log_ratio_d_e(log(d), log(e)))

        if(is.na(prob_not_null_vecs[[j]][i]) == TRUE){
          prob_not_null_vecs[[j]][i] = 0
        } else {}

        # Updating our value for gamma 1
        if(prob_not_null_vecs[[j]][i]>=0.5){
          gamma[i,j] = 1
        } else {gamma[i,j] = 0}
      }
    }


    # Gene Selection
    if(mthd == 'MCMC'){


      if(i%%skipper==0){

        # initialize parallelization clusters
        cl <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cl)

        results <- foreach(j = 1:num_pwy, .packages = c("plgp")) %dopar% {
          if (rel_method == "gp") {
            Glstr <- as.list(matrix(prv.alpha(gamma, i), nrow = length(gamma[i-1,]), ncol = 1))
            KDFG <- Map('*', kmat_dfs, Glstr)
            AMK <- Map('*', alpha_mats_k, Glstr)
            y_tilde_j <- y - mean(y) - Reduce("+", Map('%*%', KDFG, AMK)) + KDFG[[j]] %*% AMK[[j]]
          } else {
            y_tilde_j <- unlist(y_tilde[j])  # Default assignment if rel_method is not "gp"
          }

          p_gam <- mean(gamma[(max((i - skipper), 1)):i, j], na.rm = TRUE)

          elbo_vec <- numeric(Restarts)
          mu_list <- vector("list", Restarts)

          for (r in 1:Restarts) {
            outc <- BIGM(
              y = y_tilde_j,
              X = pwy_dfs[[j]],
              d_opt = mlnr_rho,
              corrmat = corr_mats[[j]],
              num_locs = ncol(pwy_dfs[[j]]),
              k_on = 1,
              N = 300,
              rand = TRUE,
              prior_scaler = 1.0,
              lkli = dist
            )
            elbo_vec[r] <- outc$elbo
            mu_list[[r]] <- outc$mu_r
          }

          elbo_vec[is.nan(elbo_vec)] <- 0
          curr_xi_dfs_j <- ((p_gam * mu_list[[which.max(elbo_vec)]] + (1 - p_gam) * -1) > 0) * 1

          idx <- replace(curr_xi_dfs_j, is.na(curr_xi_dfs_j), 0)
          if (sum(idx) == 0) {
            idx <- sample(1:ncol(pwy_dfs[[j]]), sample(1:ncol(pwy_dfs[[j]]), 1))
          }
          # kmat_dfs_j <- plgp::covar(as.matrix(pwy_dfs[[j]][, idx]), d = (4 / (3 * nrow(pwy_dfs[[j]])))^(0.2) * sqrt(1), g = 0.00001)
          kmat_dfs_j <- plgp::covar(as.matrix(pwy_dfs[[j]][, idx]), d = mlnr_rho, g = 0.00001)

          list(curr_xi_dfs_j = curr_xi_dfs_j, kmat_dfs_j = kmat_dfs_j)
        }
        for (j in 1:num_pwy) {
          curr_xi_dfs[[j]][i, ] <- results[[j]]$curr_xi_dfs_j
          kmat_dfs[[j]] <- results[[j]]$kmat_dfs_j
        }
        # Terminate parallelization clusters
        parallel::stopCluster(cl)

      } else {
        for(j in 1:num_pwy){
          curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
        }
      }

      # variational bayes implementation
    } else if(mthd == "VB"){
      if(i%%skipper==0){

        # initialize parallelization clusters
        cl <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cl)

        results <- foreach(j = 1:num_pwy, .packages = c("plgp")) %dopar% {
          if (rel_method == "gp") {
            Glstr <- as.list(matrix(prv.alpha(gamma, i), nrow = length(gamma[i-1,]), ncol = 1))
            KDFG <- Map('*', kmat_dfs, Glstr)
            AMK <- Map('*', alpha_mats_k, Glstr)
            y_tilde_j <- y - mean(y) - Reduce("+", Map('%*%', KDFG, AMK)) + KDFG[[j]] %*% AMK[[j]]
          } else {
            y_tilde_j <- unlist(y_tilde[j])  # Default assignment if rel_method is not "gp"
          }

          p_gam <- mean(gamma[(max((i - skipper), 1)):i, j], na.rm = TRUE)

          elbo_vec <- numeric(Restarts)
          mu_list <- vector("list", Restarts)

          for (r in 1:Restarts) {
            outc <- VB(
              y = y_tilde_j,
              X = pwy_dfs[[j]],
              d_opt = mlnr_rho,
              corrmat = corr_mats[[j]],
              num_locs = ncol(pwy_dfs[[j]]),
              Sigmat = sigmasq[i,],
              GP = GP_status,
              lambda = 1,
              iters = 30,
              rand = TRUE,
              prior_scaler = 1.0,
              lkli = dist
            )
            elbo_vec[r] <- outc$elbo
            mu_list[[r]] <- outc$mu_r
          }

          elbo_vec[is.nan(elbo_vec)] <- 0
          selected_mu <- mu_list[[which.max(elbo_vec)]]
          curr_xi_dfs_j <- ((p_gam * selected_mu + (1 - p_gam) * -1) > 0) * 1

          idx <- replace(curr_xi_dfs_j, is.na(curr_xi_dfs_j), 0)
          if (sum(idx) == 0) {
            idx <- sample(1:ncol(pwy_dfs[[j]]), sample(1:ncol(pwy_dfs[[j]]), 1))
          }
          # kmat_dfs_j <- plgp::covar(as.matrix(pwy_dfs[[j]][, idx]), d = (4 / (3 * nrow(pwy_dfs[[j]])))^(0.2) * sqrt(1), g = 0.00001)
          kmat_dfs_j <- plgp::covar(as.matrix(pwy_dfs[[j]][, idx]), d = mlnr_rho, g = 0.00001)

          list(curr_xi_dfs_j = curr_xi_dfs_j, kmat_dfs_j = kmat_dfs_j)
        }
        for (j in 1:num_pwy) {
          curr_xi_dfs[[j]][i, ] <- results[[j]]$curr_xi_dfs_j
          kmat_dfs[[j]] <- results[[j]]$kmat_dfs_j
        }
        # Terminate parallelization clusters
        parallel::stopCluster(cl)

      } else {
        for(j in 1:num_pwy){
          curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
        }
      }

    }

    prcent = prcent + 1/(2*N_norm)*100
    cat(sprintf("\rProgress: %.1f%%", prcent))
  }

  MLN_gamma_results = numeric(num_pwy)
  for(j in 1:num_pwy){
    MLN_gamma_results[j] = mean(gamma[((N_norm*0.7)):(N_norm),j], na.rm=TRUE)
  }

  #min-max normalisation
  gam_mod = (((MLN_gamma_results - min(MLN_gamma_results, na.rm=TRUE))/(max(MLN_gamma_results, na.rm=TRUE) - min(MLN_gamma_results, na.rm=TRUE)))>0.5)*1

  f_xi = list()

  # Gene Selection
  if(mthd == 'MCMC'){

    # initialize parallelization clusters
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    results <- foreach(j = 1:num_pwy, .packages = c("plgp")) %dopar% {
      if (rel_method == "gp") {
        Glstr <- gam_mod
        KDFG <- Map('*', kmat_dfs, Glstr)
        AMK <- Map('*', alpha_mats_k, Glstr)
        y_tilde_j <- y
      } else {
        y_tilde_j <- unlist(y_tilde[j])  # Default assignment if rel_method is not "gp"
      }

      p_gam <- gam_mod[j]

      elbo_vec <- numeric(Restarts)
      mu_list <- vector("list", Restarts)

      for (r in 1:Restarts) {
        outc <- BIGM(
          y = y_tilde_j,
          X = pwy_dfs[[j]],
          d_opt = mlnr_rho,
          corrmat = corr_mats[[j]],
          num_locs = ncol(pwy_dfs[[j]]),
          k_on = 1,
          N = 300,
          rand = TRUE,
          prior_scaler = 1.0,
          lkli = dist
        )
        elbo_vec[r] <- outc$elbo
        mu_list[[r]] <- outc$mu_r
      }

      elbo_vec[is.nan(elbo_vec)] <- 0
      curr_xi_dfs_j <- ((p_gam * mu_list[[which.max(elbo_vec)]] + (1 - p_gam) * -1) > 0) * 1

      f_xi_j <- ((p_gam * mu_list[[which.max(elbo_vec)]] + (1 - p_gam) * -1) > 0) * 1

      list(curr_xi_dfs_j = curr_xi_dfs_j, f_xi_j = f_xi_j)
    }
    for (j in 1:num_pwy) {
      curr_xi_dfs[[j]][i, ] <- results[[j]]$curr_xi_dfs_j
      f_xi[[j]] <- results[[j]]$f_xi_j
    }
    # Terminate parallelization clusters
    parallel::stopCluster(cl)
    # variational bayes implementation
  } else if(mthd == "VB"){

    # initialize parallelization clusters
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)

    results <- foreach(j = 1:num_pwy, .packages = c("plgp")) %dopar% {
      if (rel_method == "gp") {
        Glstr <- gam_mod
        KDFG <- Map('*', kmat_dfs, Glstr)
        AMK <- Map('*', alpha_mats_k, Glstr)
        y_tilde_j <- y
      } else {
        y_tilde_j <- unlist(y_tilde[j])  # Default assignment if rel_method is not "gp"
      }

      p_gam <- gam_mod[j]

      elbo_vec <- numeric(Restarts)
      mu_list <- vector("list", Restarts)

      for (r in 1:Restarts) {
        outc <- VB(
          y = y_tilde_j,
          X = pwy_dfs[[j]],
          d_opt = mlnr_rho,
          corrmat = corr_mats[[j]],
          num_locs = ncol(pwy_dfs[[j]]),
          Sigmat = sigmasq[i, ],
          GP = GP_status,
          lambda = 1,
          iters = 30,
          rand = TRUE,
          prior_scaler = 1.0,
          lkli = dist
        )
        elbo_vec[r] <- outc$elbo
        mu_list[[r]] <- outc$mu_r
      }

      elbo_vec[is.nan(elbo_vec)] <- 0
      f_xi_j <- ((p_gam * mu_list[[which.max(elbo_vec)]] + (1 - p_gam) * -1) > 0) * 1

      list(f_xi_j = f_xi_j)
    }
    f_xi <- lapply(results, function(res) res$f_xi_j)
    # Terminate parallelization clusters
    parallel::stopCluster(cl)
  }

  mln_indic = 1
  MLN_results = currxi_results = numeric(ncol(dat)-1)
  for(j in 1:num_pwy){
    for(k in 1:ncol(pwy_dfs[[j]])){
      MLN_results[mln_indic] = f_xi[[j]][k]
      currxi_results[mln_indic] = round(f_xi[[j]][k])
      mln_indic = mln_indic+1
    }
  }

  # Rebuilding Pathways using only selected genes
  for(j in 1:num_pwy){
    xi_vec = f_xi[[j]]
    kmat_dfs[[j]] = plgp::covar(as.matrix(pwy_dfs[[j]][, xi_vec]), d = mlnr_rho, g = 0.00001)
  }

  selected_indcs = gam_mod*seq(1,num_pwy)
  kmat_dfs_fin = kmat_dfs[selected_indcs]


  alpha_prior_V = list()
  for(i in 1:sum(gam_mod)){
    if(penalty == "function"){
        alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*kmat_dfs_fin[[i]]+sqrt(.Machine$double.eps)*diag(nrow = nx)))
    } else if(penalty == "weights"){
        alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*diag(nrow=nx)))
    }
  }

  # applying the above function
  if (dist == "mvn"){
    alpha_mats = alpha_creator(y, kmat_dfs_fin, alpha_prior_V, N_norm, dist = dist, sigmasq = sigmasq)
  } else if(dist == "ald"){
    alpha_mats = alpha_creator(y, kmat_dfs_fin, alpha_prior_V, N_norm, dist = dist, sigmasq = sigmasq, ald_tau=ald_tau, ald_theta=ald_theta, ald_z_vec=ald_z_vec, ald_z_mat=ald_z_mat, ald_bigB_inv=ald_bigB_inv)
  }

  ll_vec = c()
  # YOURE USING WRONG INDEXING
  for(i in 2:N_norm){

    alpha_mats_k = lapply(alpha_mats, function(x) prv.alpha(x, i))

    if(dist == "mvn"){

      for (j in 1:sum(gam_mod)){

        if(sum(gam_mod)==1){
          Rem_term = 0
        } else {
          Rem_term = Reduce("+",(Map('%*%',kmat_dfs_fin[-j],alpha_mats_k[-j])))
        }

        # sampling first path coefficients from spike slab prior
        alph_var = sigmasq[i-1,1]*solve(sigmasq[i-1,1]*solve(alpha_prior_V[[j]])+t(kmat_dfs_fin[[j]])%*%(kmat_dfs_fin[[j]]))
        alph_var = as.matrix(Matrix::forceSymmetric(alph_var))
        alph_mean = (1/sigmasq[i-1,1])*alph_var%*%t(kmat_dfs_fin[[j]])%*%(y - mean(y)-omega*Rem_term)
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)

      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] = rinvgamma(1,length(y)/2+a/2,0.5*sum((y - mean(y) - Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))^2)+b/2)

    } else if(dist == "ald"){

      for (j in 1:sum(gam_mod)){

        if(sum(gam_mod)==1){
          Rem_term = 0
        } else {
          Rem_term = Reduce("+",(Map('%*%',kmat_dfs_fin[-j],alpha_mats_k[-j])))
        }

        # sampling first path coefficients from spike slab prior
        ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(sqrt(sigmasq[i-1,])*ald_tau^2) + solve(alpha_prior_V[[j]])
        ald_bigB = solve(ald_bigB_inv[[j]])
        alph_mean = ald_bigB%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix((y- mean(y)-omega*Rem_term)-ald_theta*ald_z_vec)/(sqrt(sigmasq[i-1,])*ald_tau^2))
        alph_var = as.matrix(Matrix::forceSymmetric(ald_bigB))
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)
      }

      ald_delta_sq = ((y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))^2)/(ald_tau^2*sqrt(sigmasq[i-1,]))
      ald_gamma_sq = 2/sqrt(sigmasq[i-1,]) + (ald_theta^2)/(sqrt(sigmasq[i-1,])*ald_tau^2)
      ald_delta = sqrt(ald_delta_sq)
      ald_gamma = sqrt(ald_gamma_sq)

      plc_hld_z = numeric(1)
      for (k in 1:length(ald_z_vec)){
        plc_hld_z = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
        ald_z_vec[k] = mean(plc_hld_z, na.rm=TRUE)
        ald_z_mat[k,k] = ald_z_vec[k]
      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] =  (rinvgamma(1,n0/2+3*length(y)/2, s0/2+sum(ald_z_vec)+0.5*sum((solve(ald_z_mat))%*%(as.matrix(y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k)))-ald_theta*ald_z_vec)^2)/(ald_tau^2))))^2
    }

    yhat = as.numeric(Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))
    ll_vec = c(ll_vec, sum(dnorm(y, yhat, sigmasq[i,], log=TRUE)))

    prcent = prcent + 1/(2*N_norm)*100
    cat(sprintf("\rProgress: %.1f%%", prcent))
  }

  posterior_mean_alpha = list()
  posterior_var_alpha = list()
  # posterior_mean
  if(dist == "mvn"){

    for (j in 1:sum(gam_mod)){

      if(sum(gam_mod)==1){
        Rem_term = 0
      } else {
        Rem_term = Reduce("+",(Map('%*%',kmat_dfs_fin[-j],alpha_mats_k[-j])))
      }

      # sampling first path coefficients from spike slab prior
      alph_var = sigmasq[(N_norm-1),1]*solve(sigmasq[(N_norm-1),1]*solve(alpha_prior_V[[j]])+t(kmat_dfs_fin[[j]])%*%(kmat_dfs_fin[[j]]))
      alph_var = as.matrix(Matrix::forceSymmetric(alph_var))
      posterior_mean_alpha[[j]] = (1/sigmasq[(N_norm-1),1])*alph_var%*%t(kmat_dfs_fin[[j]])%*%(y - mean(y)-omega*Rem_term)
      posterior_var_alpha[[j]] = alph_var
    }
  } else if(dist == "ald"){

    for (j in 1:sum(gam_mod)){

      if(sum(gam_mod)==1){
        Rem_term = 0
      } else {
        Rem_term = Reduce("+",(Map('%*%',kmat_dfs_fin[-j],alpha_mats_k[-j])))
      }

      # sampling first path coefficients from spike slab prior
      ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(sqrt(sigmasq[(N_norm-1),])*ald_tau^2) + solve(alpha_prior_V[[j]])
      ald_bigB = solve(ald_bigB_inv[[j]])
      posterior_mean_alpha[[j]] = ald_bigB%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix((y- mean(y)-omega*Rem_term)-ald_theta*ald_z_vec)/(sqrt(sigmasq[(N_norm-1),])*ald_tau^2))
      posterior_var_alpha[[j]] = as.matrix(Matrix::forceSymmetric(ald_bigB))
    }
  }

  alpha_mats_k = posterior_mean_alpha

  yhat = as.numeric(Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))*sd(out_list[['y']])+mean(out_list[['y']])
  MLN_mse = mean((y-Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))^2)


  ll_met = sum(dnorm(y, yhat, sigmasq[N_norm,], log=TRUE))
  model_metric = ll_met

  out_list[["yhat"]] = yhat
  out_list[["gamma"]] = gam_mod
  out_list[["gamma_prob"]] = MLN_gamma_results
  cntr = 1
  for(j in 1:num_pwy){
    # stringr_alpha = paste0("alpha.",j)
    stringr_xi = paste0("xi.",j)
    out_list[[stringr_xi]] = f_xi[[j]]
    # if(gam_mod[j]==1){
    #   out_list[[stringr_alpha]] = alpha_mats_k[[cntr]]
    #   cntr = cntr+1
    # } else {
    #   out_list[[stringr_alpha]] = rep(0, nx)
    # }
  }
  out_list[["mse"]] = MLN_mse
  out_list[["model.metric"]] = model_metric
  out_list[["kmats"]] = kmat_dfs_fin
  out_list[["alpha.mats"]] = alpha_mats_k
  out_list[["num_sets"]] = num_pwy
  out_list[["data"]] = dat
  out_list[["all_xi"]] = MLN_results
  out_list[["mlnr_rho"]] = mlnr_rho
  out_list[["mu_alpha"]] = posterior_mean_alpha
  out_list[["sigma_alpha"]] = posterior_var_alpha

  print("Done!")
  cat("\n")

  return(out_list)

}

