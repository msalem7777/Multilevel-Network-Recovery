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
#'
#' @return MLNR model results.
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
MLNR = function(dat, num_pwy, skipper = 300, smpl.sz = 2, N_norm = 2000, level_1_connected = 1, sigmasq_y = 1,
                a = 1, b = 1, ald_p = 0.5, n0=1, s0=1, pi_0 = 0.5, a_al=0.5, b_al=0.5, sigmasq_alpha=1,
                penalty = "function", dist = "mvn", mthd = "MCMC", Restarts = 1, rel_method = "cor", w_set="avg.corr"){

  eps = sqrt(.Machine$double.eps)

  # Placeholder matrices
  MLN_G_mat = as.data.frame(matrix(0, nrow = num_pwy, ncol=2))
  MLN_mat = as.data.frame(matrix(0, nrow = (dim(dat)[2]-1), ncol=1))

  MLN_vec_mse = c()

  # Creating a y
  y = dat[,1]
  og_df = dat[-c((length(y))),]
  group_identifier_vector = c(dat[(length(y)),])
  group_list = paste(c(1:num_pwy))
  y = y[-c(length(y))]
  y =  as.numeric(y)

  # Standardizing y
  y = scale(y)

  # applying the above function
  pwy_dfs = pathway_creator(dat[2:dim(dat)[2]], num_pwy)

  if(rel_method == "cor"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = as.data.frame(as.numeric(cor(y,pwy_dfs[[i]])^2))
    }
  } else if(rel_method == "mi"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = as.data.frame(matrix(NA, nrow = (dim(pwy_dfs[[i]])[2]), ncol = 1))
      for(j in 1:(dim(pwy_dfs[[i]])[2])){
        y_tilde[[i]][j,1] = mutinformation(discretize(y),discretize(pwy_dfs[[i]][j]), method="emp")
      }
    }
  } else if(rel_method == "gp"){
    y_tilde = list()
    for(i in 1:num_pwy){
      y_tilde[[i]] = y
    }
  }

  # applying the above function
  corr_mats = corr_mat_creator(pwy_dfs, num_pwy)

  # applying the above function
  kmat_dfs = kmat_creator(pwy_dfs, num_pwy)

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
  alpha_mats = alpha_creator(kmat_dfs, N_norm, dist = dist)

  # Setting up a placeholder for path selection variable (10 pathways)
  gamma = matrix(NA, nrow = N_norm, ncol = num_pwy)
  # Initializing path selection variable
  gamma[1,] = rep(1,num_pwy)

  # applying the above function
  prob_not_null_vecs = prob_notnull_creator(pwy_dfs, N_norm)

  # applying the above function
  curr_xi_dfs = curr_xi_creator(pwy_dfs, N_norm)

  if (w_set=="avg.cor"){
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
    print(paste("i is ",i))

    alpha_mats_k = lapply(alpha_mats, prv.alpha)

    selected_all = as.double(sample(seq(1, num_pwy), num_pwy))

    if(distn == "mvn"){

      # Gamma by K matrices
      gamma_lstr = as.list(matrix(prv.alpha(gamma),nrow = length(gamma[i-1,]), ncol=1))
      kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)

      for (j in 1:num_pwy){

        # sampling first path coefficients from spike slab prior
        alph_var = sigmasq[i-1,1]*solve(sigmasq[i-1,1]*solve(alpha_prior_V[[j]])+t(kmat_dfs[[j]])%*%(kmat_dfs[[j]]))
        alph_var = as.matrix(forceSymmetric(alph_var))
        alph_mean = (1/sigmasq[i-1,1])*alph_var%*%t(kmat_dfs[[j]])%*%(y - mean(y)-omega*Reduce("+",(Map('%*%',kmat_dfs_gammad[-j],alpha_mats_k[-j]))))
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)

      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] = rinvgamma(1,length(y)/2+a/2,0.5*sum((y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)+b/2)

    } else if(distn == "ald"){

      # Gamma by K matrices
      gamma_lstr = as.list(matrix(prv.alpha(gamma),nrow = length(gamma[i-1,]), ncol=1))
      kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)


      for (j in 1:num_pwy){

        # sampling first path coefficients from spike slab prior
        ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(sqrt(sigmasq[i-1,])*ald_tau^2) + solve(alpha_prior_V[[j]])
        ald_bigB = solve(ald_bigB_inv[[j]])
        alph_mean = ald_bigB%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix((y- mean(y)-omega*Reduce("+",(Map('%*%',kmat_dfs_gammad[-j],alpha_mats_k[-j]))))-ald_theta*ald_z_vec)/(sqrt(sigmasq[i-1,])*ald_tau^2))
        alph_var = as.matrix(forceSymmetric(ald_bigB))
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

    if(distn == "mvn"){
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

    } else if(distn == "ald"){

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

        for(j in 1:num_pwy){

          if(rel_method == "gp"){
            Glstr = as.list(matrix(prv.alpha(gamma),nrow = length(gamma[i-1,]), ncol=1))
            KDFG = Map('*',kmat_dfs,Glstr)
            AMK = Map('*',alpha_mats_k,Glstr)
            y_tilde[[j]] = y - mean(y) - Reduce("+",(Map('%*%',KDFG,AMK))) + KDFG[[j]]%*%AMK[[j]]
          }

          p_gam = mean(gamma[(max((i-skipper),1)):i,j], na.rm=TRUE)

          elbo_vec = c()
          mu_list = list()

          for(r in 1:Restarts){
            outc = BIGM(y = y_tilde[[j]], pwy_dfs[[j]], corrmat = corr_mats[[j]], num_locs = ncol(pwy_dfs[[j]]), phi = (1/sigmasq[i,]), k_on = 1, N=300, rand = TRUE, prior_scaler = 1.0, lkli=distn)
            elbo_vec = c(elbo_vec, outc$elbo)
            mu_list = c(mu_list, list(outc$mu_r))
          }

          curr_xi_dfs[[j]][i,] = ((p_gam*mu_list[[which.max(elbo_vec)]]+(1-p_gam)*-1)>0)*1 #new v2

          # Rebuilding Pathways using only selected genes
          idx = replace(curr_xi_dfs[[j]][i,], is.na(curr_xi_dfs[[j]][i,]), 0)
          if(sum(idx)==0){idx=sample(1:ncol(pwy_dfs[[j]]),sample(1:ncol(pwy_dfs[[j]]),1))}
          kmat_dfs[[j]] = cov_gen(as.matrix(pwy_dfs[[j]][,idx]), theta=1 ,type="Gaussian")
        }

      } else {
        for(j in 1:num_pwy){
          curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
        }
      }
      # variational bayes implementation
    } else if(mthd == "VB"){
      if(i%%skipper==0){

        for (j in 1:num_pwy){

          print(paste("VB Iteration is: ", j))

          if(rel_method == "gp"){
            Glstr = as.list(matrix(prv.alpha(gamma),nrow = length(gamma[i-1,]), ncol=1))
            KDFG = Map('*',kmat_dfs,Glstr)
            AMK = Map('*',alpha_mats_k,Glstr)
            y_tilde[[j]] = y - mean(y) - Reduce("+",(Map('%*%',KDFG,AMK))) + KDFG[[j]]%*%AMK[[j]]
          }

          p_gam = mean(gamma[(max((i-skipper),1)):i,j], na.rm=TRUE)

          elbo_vec = c()
          mu_list = list()

          for(r in 1:Restarts){
            outc = VB(y = y_tilde[[j]], X = pwy_dfs[[j]], corrmat = corr_mats[[j]], num_locs = ncol(pwy_dfs[[j]]), Sigmat = sigmasq[i,], GP = 1, lambda = 1, iters = 30, rand = TRUE, prior_scaler = 1.0, lkli = distn)
            elbo_vec = c(elbo_vec, outc$elbo)
            mu_list = c(mu_list, list(outc$mu_r))
          }

          curr_xi_dfs[[j]][i,] = ((p_gam*mu_list[[which.max(elbo_vec)]]+(1-p_gam)*-1)>0)*1 #new v2

          print("VB done!")

          # Rebuilding Pathways using only selected genes
          idx = replace(curr_xi_dfs[[j]][i,], is.na(curr_xi_dfs[[j]][i,]), 0)
          if(sum(idx)==0){idx=sample(1:ncol(pwy_dfs[[j]]),sample(1:ncol(pwy_dfs[[j]]),1))}
          kmat_dfs[[j]] = cov_gen(as.matrix(pwy_dfs[[j]][,idx]), theta=1 ,type="Gaussian")
        }
      } else {
        for(j in 1:num_pwy){
          curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
        }
      }
    }
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

      for(j in 1:num_pwy){

        if(rel_method == "gp"){
          Glstr = gam_mod
          KDFG = Map('*',kmat_dfs,Glstr)
          AMK = Map('*',alpha_mats_k,Glstr)
          y_tilde[[j]] = y
        }

        p_gam = gam_mod[j]

        elbo_vec = c()
        mu_list = list()

        for(r in 1:Restarts){
          outc = BIGM(y = y_tilde[[j]], pwy_dfs[[j]], corrmat = corr_mats[[j]], num_locs = ncol(pwy_dfs[[j]]), phi = (1/sigmasq[i,]), k_on = 1, N=300, rand = TRUE, prior_scaler = 1.0, lkli=distn)
          elbo_vec = c(elbo_vec, outc$elbo)
          mu_list = c(mu_list, list(outc$mu_r))
        }

        curr_xi_dfs[[j]][i,] = ((p_gam*mu_list[[which.max(elbo_vec)]]+(1-p_gam)*-1)>0)*1 #new v2

        f_xi = c(f_xi,list(((p_gam*mu_list[[which.max(elbo_vec)]]+(1-p_gam)*-1)>0)*1)) #new v2
      }
    # variational bayes implementation
  } else if(mthd == "VB"){

      for (j in 1:num_pwy){

        if(rel_method == "gp"){
          Glstr = gam_mod
          KDFG = Map('*',kmat_dfs,Glstr)
          AMK = Map('*',alpha_mats_k,Glstr)
          y_tilde[[j]] = y
        }

        p_gam = gam_mod[j]

        elbo_vec = c()
        mu_list = list()

        for(r in 1:Restarts){
          outc = VB(y = y_tilde[[j]], X = pwy_dfs[[j]], corrmat = corr_mats[[j]], num_locs = ncol(pwy_dfs[[j]]), Sigmat = sigmasq[i,], GP = 1, lambda = 1, iters = 30, rand = TRUE, prior_scaler = 1.0, lkli = distn)
          elbo_vec = c(elbo_vec, outc$elbo)
          mu_list = c(mu_list, list(outc$mu_r))
        }

        f_xi = c(f_xi,list(((p_gam*mu_list[[which.max(elbo_vec)]]+(1-p_gam)*-1)>0)*1)) #new v2
      }
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
    kmat_dfs[[j]] = cov_gen(as.matrix(pwy_dfs[[j]][, xi_vec]), theta=1, type="Gaussian")
  }

  kmat_dfs_fin = kmat_dfs[gam_mod*seq(1,num_pwy)]
  alpha_mats = alpha_creator(sum(gam_mod), N_norm, kmat_dfs_fin)

  alpha_prior_V = list()
  for(i in 1:sum(gam_mod)){
    if(penalty == "function"){
      alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*kmat_dfs_fin[[i]]+sqrt(.Machine$double.eps)*diag(nrow = nx)))
    } else if(penalty == "weights"){
      alpha_prior_V = c(alpha_prior_V, list(sigmasq_alpha*diag(nrow=nx)))
    }
  }

  for(i in 2:N_norm){

    alpha_mats_k = lapply(alpha_mats, prv.alpha)

    if(distn == "mvn"){

      for (j in 1:sum(gam_mod)){

        if(sum(gam_mod)==1){
          Rem_term = 0
        } else {
          Rem_term = Reduce("+",(Map('%*%',kmat_dfs_fin[-j],alpha_mats_k[-j])))
        }

        # sampling first path coefficients from spike slab prior
        alph_var = sigmasq[i-1,1]*solve(sigmasq[i-1,1]*solve(alpha_prior_V[[j]])+t(kmat_dfs_fin[[j]])%*%(kmat_dfs_fin[[j]]))
        alph_var = as.matrix(forceSymmetric(alph_var))
        alph_mean = (1/sigmasq[i-1,1])*alph_var%*%t(kmat_dfs_fin[[j]])%*%(y - mean(y)-omega*Rem_term)
        alpha_mats[[j]][i,] = rmvnorm(1, alph_mean,alph_var)

      }

      # sampling new phi
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
      sigmasq[i,] = rinvgamma(1,length(y)/2+a/2,0.5*sum((y - mean(y) - Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))^2)+b/2)

    } else if(distn == "ald"){

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
        alph_var = as.matrix(forceSymmetric(ald_bigB))
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
  }

  alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[N_norm*0.7:N_norm, ]})
  alpha_mats_k = lapply(alpha_mats_k, function(x) {x <- colMeans(x, na.rm=TRUE)})

  MLN_mse = mean((y- mean(y)-Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k))))^2)

  if(mthd=="VB"){
    model_metric = outc$elbo
  } else if(mthd == "MCMC"){
    model_metric = sum(dnorm(y, mean=(mean(y)+Reduce("+",(Map('%*%',kmat_dfs_fin,alpha_mats_k)))), sd = 1, log = TRUE))
  }

  out_list = list()
  out_list[["gamma"]] = gam_mod
  out_list[["gamma_prob"]] = MLN_gamma_results
  for(j in 1:num_pwy){
    stringr_xi = paste0("xi.",j)
    out_list[[stringr_xi]] = f_xi[[j]]
    stringr_alpha = paste0("alpha.",j)
    out_list[[stringr_alpha]] = alpha_mats_k[[j]]
  }
  out_list[["mse"]] = MLN_mse
  out_list[["model.metric"]] = model_metric
  out_list[["kmats"]] = kmat_dfs_fin
  out_list[["alpha.mats"]] = alpha_mats_k

  return(out_list)

}

