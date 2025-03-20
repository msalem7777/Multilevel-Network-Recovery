#' VB
#'
#' Variational Bayes version of the MLNR algorithm.
#'
#' @param y A continuous response variable as a numeric vector.
#' @param X A dataframe or data matrix.
#' @param corrmat A matrix of prior interactions between elements. Can be set as the squared correlation.
#' @param num_locs A numeric value representing number of elements under examination.
#' @param prior_scaler A numeric value scaling the specified interaction matrix.
#' @param iters A numeric value for number of VB iterations.
#' @param k_on A logical value. FALSE (default). TRUE incorporates a matrix of marginal regressions as an external field effect -- see Fang & Kim 2016
#' @param rand A logical value. FALSE (default). TRUE randomizes the sampling order of the elements in the McMC exploration -- not recommended.
#' @param lkli A string. 'mvn' (default) uses a normal likelihood. 'ald' uses an asymmetric Laplace likelihood'
#' @param ald_skew A numeric value for the skew of the asymmetric LaPlace likelihood -- default 0.5, median.
#' @return A binary vector containing selection status for elements within sets.
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
VB = function(y, X, corrmat, num_locs, d_opt = 1, lambda = 1.0, iters=50, GP = 0, Sigmat = NULL, prior_scaler = 10, plot_on=FALSE, true_dist = 0, ker_wts = NULL, rand = TRUE, lambda_cv = 2.0, lkli = "mvn"){

  ELBO = numeric(iters)
  Hx_mean = numeric(iters)
  d_gp = c()

  if(GP == 0){
    library(ClusterR)
    clustr1 = ClusterR::GMM(as.matrix(y), gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10, em_iter = 50)
    mus = clustr1$centroids
    sds = clustr1$covariance_matrices

    emp_mean_num = max(mus)
    emp_sd_num = sqrt(sds[which.max(mus)])
    emp_sd_denom = sqrt(sds[which.min(mus)])

    if(true_dist==0){
      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(2*dnorm(as.matrix(y), mean = 0, sd = emp_sd_denom)))
    } else if(true_dist==1){

      logodds = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = emp_sd_num)/(dcorr(as.matrix(y), n = nrow(X), ncp = 0)))
    }

    p1 = sigmoid(logodds)
    mu = 2*p1-1

    a = mu + 0.5 * logodds

    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)

    if(true_dist==0){
      logp1 = log(dnorm(as.matrix(y), mean = emp_mean_num, sd = sqrt(sds[which.max(mus)])))
      logm1 = log(2*dnorm(as.matrix(y), mean = 0, sd = sqrt(sds[which.min(mus)])))
    } else if(true_dist==1){
      logp1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = emp_mean_num))
      logm1 = log(dcorr(as.matrix(y), n = nrow(X), ncp = 0))
    }

  } else if(GP == 1){

    logp1 = logm1 = c()
    logodds = fX_sig = c()
    fX = list()
    YGP = as.matrix(as.numeric(scale(y)), ncol=1)

    K_full_list = list()

    for(l in 1:num_locs){
      if((dim(X)[2]/num_locs)>1){
        mod_val = (dim(X)[2]/num_locs)
        XGP = as.matrix(X[,(mod_val*(l-1)+1):(mod_val*(l-1)+mod_val)])
      } else {
        XGP = as.matrix(X[,l])
      }

      Sigmat = 1
      K = plgp::covar.sep(XGP, d=d_opt, g=Sigmat)
      K_full_list[[l]] = K
      KX = plgp::covar.sep(XGP, d=d_opt, g=0)
      Ki = solve(K)

      fX = c(fX, list(KX %*% Ki %*% YGP))
      fX_sig = c(fX_sig, Sigmat)


      if(lkli == "mvn"){
        logm1 = c(logm1, sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE)))
        logp1 = c(logp1, sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE)))
        logodds = c(logodds, (sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE))-sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE))))
      } else if(lkli == "ald"){

        l1 = l2 = c()
        m1_ald = (mean(YGP)+fX[[l]])
        m2_ald = rep(mean(YGP), length(YGP))
        for(lp in 1:length(YGP)){
          l1 = c(l1, dALD(YGP[lp], m1_ald[lp]))
          l2 = c(l2, dALD(YGP[lp], m2_ald[lp]))
        }

        l1 = sum(log(l1))
        l2 = sum(log(l2))

        logm1 = c(logm1, l2)
        logp1 = c(logp1, l1)
        logodds = c(logodds, log_sum_exp((l1-l2)))
      }

    }

    clustr2 = ClusterR::GMM(as.matrix(logodds,ncol=1), gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10, em_iter = 50)
    clusters = ClusterR::predict_GMM(as.matrix(logodds), clustr2$centroids, clustr2$covariance_matrices, clustr2$weights)
    clusters = clusters$cluster_labels
    big = which.max(clustr2$centroids)
    small = which.min(clustr2$centroids)
    idx_big= min(logodds[clusters==big])
    idx_small= max(logodds[clusters==small])
    penalty = (idx_big-idx_small)/2
    penalty = idx_big-penalty

    logp1 = logp1 - penalty

    logodds = logodds - penalty


    p1 = sigmoid(logodds)
    mu = 2*p1-1

    a = mu + 0.5 * logodds

    qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)

  }

  if(rand == TRUE){
    rlistr = c(sample(seq(1,num_locs), replace = FALSE))
  } else if(rand == FALSE){
    rlistr = seq(1,num_locs)
  }

  diag(corrmat) = 0

  for (i in 1:iters){

    muNew = mu

    for (j in rlistr){

      Sbar = sum(prior_scaler*(corrmat[j,]^2)*muNew)

      muNew[j] = (1-lambda)*muNew[j] + lambda*tanh(Sbar + 0.5*logodds[j])

      ELBO[i] = ELBO[i] + 0.5*(Sbar * muNew[j])
    }

    mu = muNew

    a = mu + 0.5 * logodds
    qxp1 = sigmoid(+2*a) #q_i(x_i=+1)
    qxm1 = sigmoid(-2*a) #q_i(x_i=-1)
    Hx = -qxm1*log(qxm1+1e-10) - qxp1*log(qxp1+1e-10) #entropy

    ELBO[i] = ELBO[i] + sum(qxp1*logp1 + qxm1*logm1) + sum(Hx)
    Hx_mean[i] = mean(Hx)

  }

  if(plot_on == TRUE){
    plt_x = seq(1:(i-1))
    plot(plt_x,ELBO[1:(i-1)])
  }

  mu_r = ifelse(mu > 0, 1, 0)

  listr = list(mu_r = mu_r, elbo = ELBO[iters])

  return(listr)
}

