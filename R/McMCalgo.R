#' BIGM
#'
#' McMC version of the MLNR algorithm.
#'
#' @param y A continuous response variable as a numeric vector.
#' @param X A dataframe or data matrix.
#' @param corrmat A matrix of prior interactions between elements. Can be set as the squared correlation.
#' @param num_locs A numeric value representing number of elements under examination.
#' @param prior_scaler A numeric value scaling the specified interaction matrix.
#' @param N A numeric value for number of McMC iterations.
#' @param k_on A logical value. FALSE (default). TRUE incorporates a matrix of marginal regressions as an external field effect -- see Fang & Kim 2016
#' @param rand A logical value. FALSE (default). TRUE randomizes the sampling order of the elements in the McMC exploration -- not recommended.
#' @param lkli A string. 'mvn' (default) uses a normal likelihood. 'ald' uses an asymmetric Laplace likelihood'
#' @param ald_skew A numeric value for the skew of the asymmetric LaPlace likelihood -- default 0.5, median.
#' @return A binary vector containing selection status for elements within sets.
#' @import ClusterR
#' @import GGally
#' @import GIGrvg
#' @import LaplacesDemon
#' @import MASS
#' @import Matrix
#' @import ald
#' @import car
#' @import class
#' @import dplyr
#' @import expm
#' @import fMultivar
#' @import ggnetwork
#' @import ggplot2
#' @import grpreg
#' @import gtools
#' @import hetGP
#' @import infotheo
#' @import invgamma
#' @import laGP
#' @import mvtnorm
#' @import network
#' @import plgp
#' @import quantreg
#' @import robustbase
#' @import rpart
#' @import sna
#' @import splines
#' @import tidyr
#' @import tidyverse
#' @importFrom ClusterR GMM
#' @importFrom ClusterR predict_GMM
#' @importFrom ald rALD
#' @importFrom base sample
#' @importFrom plgp covar
#' @importFrom plgp covar.sep
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom utils combn
#' @export
BIGM = function(y, X, corrmat, num_locs, prior_scaler = 10, N=1000, k_on = 0, rand = FALSE, lkli = "mvn", ald_skew=0.5){
  y_s = scale(y)
  X_s = X
  MMR = matrix(0, nrow = dim(X_s)[1], ncol = dim(X_s)[2])
  h = matrix(0, nrow = 2, ncol = dim(X_s)[2])
  if(k_on==0){
    for(i in 1:(dim(X_s)[2])){
      MMR[,i] = (lm(y~X_s[,i]))$coefficients[2]*X_s[,i]
    }
  } else if(k_on==1) {

    logodds = fX_sig = c()
    fX = list()
    YGP = as.matrix(as.numeric(scale(y)), ncol=1)

    Sigmat = 1
    # d_opt = (4/(3*nrow(X_s)))^(0.2)*sqrt(Sigmat)
    d_opt = 1

    logodds = fX_sig = c()

    for(l in 1:(dim(X_s)[2])){

      XGP = as.matrix(X[,l])
      K = plgp::covar.sep(XGP, d=d_opt, g=Sigmat)
      KX = covar.sep(XGP, d=d_opt, g=0)
      Ki = solve(K)
      fX = c(fX, list(KX %*% Ki %*% YGP))
      fX_sig = c(fX_sig, Sigmat)

      if(lkli == "mvn"){
        h[1,l] = sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE))
        h[2,l] = sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE))
        logodds = c(logodds, (sum(dnorm(as.matrix(YGP), mean = fX[[l]], sd = sqrt(fX_sig[l]), log = TRUE))-sum(dnorm(as.matrix(YGP), mean = 0, sd = sqrt(fX_sig[l]), log = TRUE))))
      } else if(lkli == "ald"){

        l1 = l2 = c()
        m1_ald = (mean(YGP)+fX[[l]])
        m2_ald = rep(mean(YGP), length(YGP))
        for(lp in 1:length(YGP)){
          l1 = c(l1, dALD(YGP[lp], m1_ald[lp], p=ald_skew))
          l2 = c(l2, dALD(YGP[lp], m2_ald[lp], p=ald_skew))
        }

        l1 = sum(log(l1))
        l2 = sum(log(l2))

        h[1,l] = l2
        h[2,l] = l1
        logodds = c(logodds, log_sum_exp((l1-l2)))
      }

    }

    clustr2 = GMM(as.matrix(logodds,ncol=1), gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10, em_iter = 50)
    clusters = predict_GMM(as.matrix(logodds), clustr2$centroids, clustr2$covariance_matrices, clustr2$weights)
    clusters = clusters$cluster_labels
    big = which.max(clustr2$centroids)
    small = which.min(clustr2$centroids)
    idx_big= min(logodds[clusters==big])
    idx_small= max(logodds[clusters==small])
    penalty = (idx_big-idx_small)/2
    penalty = idx_big-penalty

  }
  if(rand == TRUE){
    rlistr = c(sample(seq(1,num_locs), replace = FALSE))
  } else if(rand == FALSE){
    rlistr = seq(1,num_locs)
  }

  diag(corrmat) = 0

  Sbar = prior_scaler*(corrmat^2)

  # Matrix to hold metropolis results
  cmat = matrix(NA, nrow = N+1, ncol = num_locs)
  cmat[1,] = sample(c(0,1), num_locs, replace = TRUE)
  curr_xi = cmat[1,]

  # Metropolis algorithm
  for (i in 1:N){
    for (j in rlistr){

      ref_curr_xi = replace(curr_xi, curr_xi==0, -1)
      xi_int_mat_curr = c(ref_curr_xi)%*%t(c(ref_curr_xi))

      rm(ref_curr_xi)

      U_current = sum(LaplacesDemon::upper.triangle(Sbar*xi_int_mat_curr, diag = FALSE))+sum(h[curr_xi[j]+1,j])-(penalty*curr_xi[j])

      prop_xi = curr_xi
      prop_xi[j] = abs(curr_xi[j]-1)

      ref_prop_xi = replace(prop_xi, prop_xi==0, -1)
      xi_int_mat_prop = c(ref_prop_xi)%*%t(c(ref_prop_xi))

      rm(ref_prop_xi)

      U_proposed = sum(LaplacesDemon::upper.triangle(Sbar*xi_int_mat_prop, diag = FALSE))+sum(h[prop_xi[j]+1,j])-(penalty*prop_xi[j])

      delta_U = U_proposed - U_current

      prob_flip = min(c(1, exp(delta_U)))

      flipper = rbinom(1,1,prob_flip)

      if (flipper == 1){
        curr_xi[j] = prop_xi[j]
      }

    }

    cmat[i+1,] = curr_xi
  }

  cmat_thin <- cmat[(round(0.75*N)+1):N, ]
  cmat_thin <- cmat_thin[(1:nrow(cmat_thin)) %% 5 == 0, ]

  curr_xi = colMeans(cmat_thin)

  mu_r = curr_xi
  ph_xi  = round(curr_xi)

  listr = list(mu_r = mu_r, elbo = 0)

  return(listr)

}

