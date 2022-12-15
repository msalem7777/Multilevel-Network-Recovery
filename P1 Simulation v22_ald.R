#######################################################
################     Libraries     ####################
#######################################################

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ald))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(stargazer))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(rpart))
suppressPackageStartupMessages(library(fMultivar))
suppressPackageStartupMessages(library(LaplacesDemon))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(invgamma))
suppressPackageStartupMessages(library(network))
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(grpreg))
suppressPackageStartupMessages(library(ggnetwork))
suppressPackageStartupMessages(library(GIGrvg))
suppressPackageStartupMessages(library(expm))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(infotheo)) # library for mutual information
suppressPackageStartupMessages(library(ClusterR)) # library for Gaussian Mixture Model

# competing methods
suppressPackageStartupMessages(library(glmnet)) # Lasso
suppressPackageStartupMessages(library(grplasso)) # Group Lasso
suppressPackageStartupMessages(library(SGL)) # Sparse Group Lasso
suppressPackageStartupMessages(library(BhGLM)) # Group spike-and-slab lasso
suppressPackageStartupMessages(library(BSGS)) # Bayesian Sparse Group Selection
suppressPackageStartupMessages(library(SSGL)) # Spike-and-Slab Group Lassos for Grouped Regression and Sparse Generalized Additive Models

#######################################################
################     BIGM          ####################
#######################################################

BIGM = function(path_features, y, curr_xi, N=1000, phi=1){
  y_s = scale(y)
  X_s = scale(path_features)
  MMR = matrix(0, nrow = dim(X_s)[1], ncol = dim(X_s)[2])
  for(i in 1:(dim(X_s)[2])){
    MMR[,i] = (lm(y~X_s[,i]))$coefficients[2]*X_s[,i]
  }
  
  J = -phi/2*t(MMR)%*%(MMR)
  h = phi*t(MMR)%*%(y-0.5*MMR%*%c(rep(1,dim(X_s)[2])))
  
  #curr_xi = sample(c(0,1), 10, replace = T)
  
  #Matrix to hold metropolis results
  cmat = matrix(NA, nrow = N+1, ncol = length(curr_xi))
  cmat[1,] = curr_xi
  
  # Metropolis algorithm
  for (i in 1:N){
    for (j in 1:length(curr_xi)){
      
      ref_curr_xi = replace(curr_xi, curr_xi==0, -1)
      xi_int_mat_curr = c(ref_curr_xi)%*%t(c(ref_curr_xi))
      xi_int_mat_curr = replace(xi_int_mat_curr, xi_int_mat_curr==-1, 0)
      rm(ref_curr_xi)
      
      U_current = -sum(J*upper.tri(xi_int_mat_curr, diag = FALSE))-sum(h*curr_xi)
      
      prop_xi = curr_xi
      prop_xi[j] = abs(curr_xi[j]-1)
      
      ref_prop_xi = replace(prop_xi, prop_xi==0, -1)
      xi_int_mat_prop = c(ref_prop_xi)%*%t(c(ref_prop_xi))
      xi_int_mat_prop = replace(xi_int_mat_prop, xi_int_mat_prop==-1, 0)
      rm(ref_prop_xi)
      
      U_proposed = -sum(J*upper.tri(xi_int_mat_prop, diag = FALSE))-sum(h*prop_xi)
      
      delta_U = U_proposed - U_current
      
      prob_flip = min(c(1, exp(-1*delta_U)))
      
      flipper = rbinom(1,1,prob_flip)  
      
      if (flipper == 1){
        curr_xi = prop_xi  
      }
      
    }
    
    cmat[i+1,] = curr_xi
    
  }
  
  #curr_xi = apply(cmat[(N+1-50+1):(N+1),], 2, modefunc)
  curr_xi = as.double(colSums(cmat[(N+1-50+1):(N+1),])/((N+1)-(N+1-50+1)+1)>0.5)
  
  return(curr_xi)
}
#######################################################
#######################################################

#######################################################
###########     Selector Function     #################
#######################################################

selector <- function(DF, iter_num){
  if(is.na(sum(DF[iter_num,])) & ((iter_num-skipper)==0)) return(DF[1, ])
  if(is.na(sum(DF[iter_num,]))) return(DF[iter_num-skipper, ])
  if(!is.na(sum(DF[iter_num,]))) return(DF[iter_num, ])
}

#######################################################
#########     Variational Inference     ###############
#######################################################
# Variational Inference

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

# Need to set lambda
# Need to set sigmasq: currently set at 0.25^2, to make sd 0.25 such that 95% of random noise correlation is less than 0.5
# Need to set h
# Need to check whether multiplying by 2 is necessary
# consider setting sigmasq to (max{kmeans{centers}}/3)^2


#X=second_path_features
VB = function(y, X, mus, full_dat, iter_num, corrmat, LB_mat, prior_scaler=1, gamma_status=1, LB=0.3, h=0.0, lambda = 1.0, iters=50, sigmasq=(max(mus$centers)/3)^2, sd_null = 0, plot_on=FALSE, dist_true = FALSE){
  
  ELBO = numeric(iters)
  Hx_mean = numeric(iters)
  
  diag(LB_mat) = 0
  
  #y = (cor(y, X))^2 # old  -- removed in v2
  #mus = kmeans(t(y), centers = 2) # old  -- removed in v2
  
  #y_full = (cor(y, full_dat[,2:dim(full_dat)[2]]))^2
  #ones = as.matrix(unlist(lapply(curr_xi_dfs,selector,iter_num))) # new
  #nonones = (ones==0)*1 # new
  #mus = c(t(ones)%*%t(as.matrix(y_full))/(sum(ones)),t(nonones)%*%t(as.matrix(y_full))/(sum(nonones))) # new
  #mus = as.data.frame(as.matrix(c(mus), nrow=1, ncol=2)) # new
  #colnames(mus) = c('centers') # new
  #y_ones = y_full[which(ones==1)] # new
  #y_nonones = y_full[which(nonones==1)] # new
  
  if(sd_null==0){
    sd_null = sqrt(sigmasq)
  }

  if(dist_true==FALSE){
    logodds = log(dnorm(y[,1], mean = max(mus), sd = sqrt(sigmasq))/(2*dnorm(y[,1], mean = 0, sd = sd_null)))
  } else if(dist_true==TRUE){
    logodds = log(dnorm(y[,1], mean = max(mus), sd = sqrt(sigmasq))/(df(y[,1], df1 = 1, df2 = length(y[,1])-2)))
  }
  
#  if(dist_true==FALSE){
#    logodds = log(dnorm(y, mean = max(mus$centers, na.rm = TRUE), sd = max(mus$centers, na.rm = TRUE)/3)*as.numeric(y<1)/(2*dnorm(y, mean = 0, sd = max(mus$centers, na.rm=TRUE)/3)*as.numeric(y>0)))
#  } else if(dist_true==TRUE){
#    logodds = log(dnorm(y, mean = max(mus$centers, na.rm = TRUE), sd = max(mus$centers, na.rm = TRUE)/3)*as.numeric(y<1)/(df(y, df1 = 1, df2 = length(y)-2)))
#  } # new
  

  
  p1 = sigmoid(logodds)
  mu = 2*p1-1
  #mu = sample(c(0,1),length(y), replace = T)
  
  a = mu + 0.5 * logodds
  
  qxp1 = sigmoid(+2*a)  #q_i(x_i=+1)
  qxm1 = sigmoid(-2*a)  #q_i(x_i=-1)
  
  logp1 = log(dnorm(y[,1], mean = max(mus), sd = sqrt(sigmasq)))
  logm1 = log(2*dnorm(y[,1], mean = 0, sd = sd_null))
  #j=1
  for (i in 1:iters){
    if(i>15){
      r = c(1,2,3,4,5,6,7,8,9,10)
      model = lm(ELBO[(i-10):(i-1)]~r)
      if(abs(coefficients(model)[2])<0.01){
        break
      }
    }
    
    muNew = mu
    
    for (j in 1:length(y)){
      
      Sbar = sum(prior_scaler*(corrmat[j,])*muNew)
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

  mu = replace(mu, mu == -1, 0)
  
  return(mu)
}  

#######################################################
#######################################################


#######################################################
#########     Gaussian Kernel Creation     ############
#######################################################

GaussianKernelize = function(path_features, curr_xi){
  
  holdr = as.matrix(path_features)%*%diag(curr_xi)
  if (sum(colSums(abs(holdr)))>0){
    holdr = (holdr[, colSums(abs(holdr)) != 0])
    kmat = exp(-1/2*as.matrix(dist(holdr, diag = TRUE, upper=TRUE)))
  } else {
    holdr = (as.matrix(path_features)[, sample(1:dim(path_features)[2],2)])
    kmat = exp(-1/2*as.matrix(dist(holdr, diag = TRUE, upper=TRUE)))
  }
  return(kmat)
}


#######################################################
#######################################################

#######################################################
pwy_s = 25
n_s = 100

for(pwy_s in c(25)){
  for(n_s in c(100)){
    #######################################################
    ##########     Simulation Settings     ################
    #######################################################
    
    nx=n_s # Number of observations
    num_pwy=pwy_s # Number of Pathways
    sim_iters = 100 # Number of simulation iterations 
    skipper = 100
    modar = 1
    

    # Placeholder matrices
    MLN_G_mat = as.data.frame(matrix(0, nrow = num_pwy, ncol=1))
    MLN_mat = as.data.frame(matrix(0, nrow = (1*num_pwy), ncol=1))
    Lasso_mat = as.data.frame(matrix(0, nrow = (1*num_pwy), ncol=1))
    GLasso_mat = as.data.frame(matrix(0, nrow = (1*num_pwy), ncol=1))
    SGL_mat = as.data.frame(matrix(0, nrow = (1*num_pwy), ncol=1))
    SSGL_mat = as.data.frame(matrix(0, nrow = (1*num_pwy), ncol=1))
    
    MLN_vec_mse = numeric(sim_iters)
    Lasso_vec_mse = numeric(sim_iters)
    Glasso_vec_mse = numeric(sim_iters)
    SGL_vec_mse = numeric(sim_iters)
    SSGL_vec_mse = numeric(sim_iters)
  
    
    s=1
    
    
    for(s in 1:sim_iters){
      print(paste("s is ",s))
      #######################################################
      #############     Building Data     ###################
      #######################################################
      
      # Need to change simulation setup such that means are gene specific
      # Need "high" intergene correlations between active genes
      #define x-values
      
      # simdat
      # Randomly generated 10 features from first path via standard normal
      sigma <- matrix(c(1,0.7,-0.3,0.6,0.7,1,-0.2,0.55,-0.3,-0.2,1,-0.1,0.6,0.55,-0.1,1), ncol=4)
      first_path_features = data.frame(matrix(rmvnorm(nx,rep(1.5,4),sigma), nrow = nx))
      
      #matrix(c(1,0.7,-0.3,0.6,0.7,1,-0.2,0.55,-0.3,-0.2,1,-0.1,0.6,0.55,-0.1,1), ncol=4)
      
      # Randomly generated 10 features from Second path via standard normal
      sigma2 <- matrix(c(1,-0.1,0.1,0.15,-0.1,1,0.65,0.55,0.1,0.65,1,0.6,0.15,0.55,0.6,1), ncol=4)
      second_path_features = data.frame(matrix(rmvnorm(nx,rep(-3.5,4),sigma2), nrow = nx))
      
      # Randomly generated 10 features from first path via standard normal
      #sigma3 <- matrix(c(1,0.7,-0.3,0.6,0.7,1,-0.2,0.55,-0.3,-0.2,1,-0.1,0.6,0.55,-0.1,1), ncol=4)
      #third_path_features = data.frame(matrix(rmvnorm(nx,rep(3.5,4),sigma2), nrow = nx))
      
      # Creating a y 
      y =  10*cos(first_path_features$X1)+3*(first_path_features$X2^2)+first_path_features$X1*first_path_features$X2*first_path_features$X4-sin(second_path_features$X3)*(second_path_features$X2)/4+2/5*(second_path_features$X3)*(second_path_features$X2)*(second_path_features$X4)+rALD(nx)#+rnorm(nx)
      
      first_path_features[dim(first_path_features)[1]+1,]=1
      second_path_features[dim(second_path_features)[1]+1,]=2
      #third_path_features[dim(third_path_features)[1]+1,]=3
      y[dim(first_path_features)[1]]=0
      
      dat = cbind(y,first_path_features,second_path_features)
      
      # Additional pathway creator 
      for(i in 3:num_pwy){
        lp = as.double(sample(4:50,1))
        #lp = rpois(1,20)
        #sigma = diag(rgamma(lp, shape=1, rate=1), lp)
        sigma = diag(lp)
        rho = matrix(runif(lp, 0.7, 0.9), lp)
        M1 = t(matrix(rho,nrow=lp*2-1,ncol=lp))[,1:lp]
        M1[lower.tri(M1)] = 0
        M1 = forceSymmetric(M1)
        diag(M1) = 1
        Sigma = sigma%*%M1%*%sigma
        nam = paste("path_features_", i, sep = "")
        path_features = data.frame(matrix(rmvnorm(nx,rep(rnorm(1),dim(sigma)[2]),sigma), nrow = nx))
        path_features[dim(path_features)[1]+1,] = i
        assign(nam, path_features)
        dat = cbind(dat, path_features)
      }
      
      colnames(dat) = c("y","X1", "X2","X3","X4","X11", "X12","X13","X14","X21", "X22","X23","X24")
      true_vec = numeric(dim(dat)[2]-1)
      true_vec[c(1,2,4,6,7,8)] = c(1,1,1,1,1,1)
      
      
      #######################################################
      #################     MLN     #########################
      #######################################################
      old <- Sys.time()
      
      N = 2000   # Number of iterations for SSVS
      #set.seed(1234)  # seed
      num_pw = num_pwy # any integer
      distn = "ald" # "mvn" or "ald"
      mthd = "VB" # "VB" or "MCMC"
      rel_method = "mi" # "mi" for mutual information or "cor" for correlation squared
      level_1_connected = 0 # 0 if no connection for supernodes, or 1 if connected
      
      # Creating a y 
      y = dat$y # simulation
      
      og_df = dat[-c((length(y))),]
      group_identifier_vector = c(dat[(length(y)),])
      group_list = paste(c(1:num_pw)) #simulation
      
      y = y[-c(length(y))]
      y =  as.numeric(y)
      
      # Standardizing y
      y = scale(y)
      y = as.numeric(y)
      
      # Extracting the grouped pathways from the dataset
      pathway_creator = function(dat, group_list, num_pw){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = as.data.frame(scale(dat[1:(dim(dat)[1]-1),c((dat[dim(dat)[1],])==group_list[i])]))
        }
        
        pwy_dfs = listr
        
        return(pwy_dfs)
      }
      # applying the above function
      pwy_dfs = pathway_creator(dat[2:dim(dat)[2]], group_list, num_pw)
      
      if(rel_method == "cor"){
        y_tilde = list()
        for(i in 1:num_pw){
          y_tilde[[i]] = as.data.frame(as.numeric(cor(y,pwy_dfs[[i]])^2))
        }
      } else if(rel_method == "mi"){
        y_tilde = list()
        for(i in 1:num_pw){
          y_tilde[[i]] = as.data.frame(matrix(NA, nrow = (dim(pwy_dfs[[i]])[2]), ncol = 1))
          for(j in 1:(dim(pwy_dfs[[i]])[2])){
            y_tilde[[i]][j,1] = mutinformation(discretize(y),discretize(pwy_dfs[[i]][j]), method="emp")
          }
        }
      }
      
      #mus = list()
      #for(i in 1:num_pw){
      #  mus[[i]] = kmeans(y_tilde[[i]], centers = 2)
      #}
      
      mus = list()
      sds = list()
      for(i in 1:num_pw){
        clustr1 = GMM(y_tilde[[i]],gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset",km_iter = 10,em_iter = 5)  
        mus[[i]] = clustr1$centroids
        sds[[i]] = clustr1$covariance_matrices
      }
      
      
      
      
      # Function to build correlation matrices for each of the pathways
      corr_mat_creator = function(df_list, num_pw){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = as.data.frame(cor(pwy_dfs[[i]]))
          diag(listr[[i]]) = 0
        }
        
        corr_mats = listr
        
        return(corr_mats)
      }
      # applying the above function
      corr_mats = corr_mat_creator(pwy_dfs, num_pw)
      
      # function to build lower bound matrices from correlation matrices using normal approximation
      corr_LBs_creator = function(corr_mats, n, num_pw, alpha=0.05){
        listr = list()
        for(i in 1:num_pw){
          J = corr_mats[[i]]
          z_crit = qnorm(1-alpha)
          listr[[i]] = as.data.frame(((1+J)/(1-J)*exp((-2*z_crit)/(sqrt(n-3)))-1)/((1+J)/(1-J)*exp((-2*z_crit)/(sqrt(n-3)))+1))
        }
        
        corr_LBs = listr
        
        return(corr_LBs)
      }
      # applying above function
      corr_LBs = corr_LBs_creator(corr_mats, n = length(y), num_pw)
      
      # Function to create kernelized matrices for each of the pathways
      kmat_creator = function(df_list, num_pw){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = exp(-1/2*as.matrix(dist(pwy_dfs[[i]], diag = TRUE, upper=TRUE)))
        }
        
        kmat_dfs = listr
        
        return(kmat_dfs)
      }
      
      # applying the above function
      kmat_dfs = kmat_creator(pwy_dfs, num_pw)
      
      # Function to create kernelized matrices for each of the pathways
      bigB_creator = function(kmat_dfs, num_pw){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = diag(dim(kmat_dfs[[i]])[1])
        }
        
        ald_bigB_inv = listr
        
        return(ald_bigB_inv)
      }
      
      # applying the above function
      ald_bigB_inv = bigB_creator(kmat_dfs, num_pw)
      
      # Setting up overall variance placeholder
      sigmasq = data.frame(rep(NA, N))
      sigmasq[1,] = 1
      
      # Setting up placeholder for coefficient variance
      sigmasq_alpha = data.frame(rep(NA, N))
      sigmasq_alpha[1,] = 1
      
      # Inverse Gamma parameters for overall variance prior
      a = 0.5
      b = 0.5
      
      # Inverse Gamma parameters for coefficient variance prior
      a_al = 0.5
      b_al = 0.5
      
      # total number of variables
      var_tot = 0
      for(i in 1:num_pw){
        var_tot = var_tot+dim(pwy_dfs[[i]])[2]
      }
      
      # Setting up placeholder for coefficients
      alpha_creator = function(num_pw, N, pwy_dfs, init = 0.5){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = matrix(NA, nrow = N, ncol = dim(pwy_dfs[[i]])[1])
          listr[[i]][1,] = init
        }
        
        alpha_mats = listr
        
        return(alpha_mats)
      }
      # applying the above function
      alpha_mats = alpha_creator(num_pw, N, pwy_dfs)
      
      # Setting up smoother for marginal coefficients
      smooth_creator = function(num_pw, n, reg=1){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = reg*diag(n)
        }
        
        smooth_mats = listr
        
        return(smooth_mats)
      }
      # applying the above function
      smooth_mats = smooth_creator(num_pw, nx, reg = 1)
      
      # Setting Binomial probability
      pi_0 = 0.5
      
      # Setting Asymmetric LaPlace quantile 
      ald_p = 0.5
      
      # Setting Asymmetric LaPlace Tau 
      ald_tau = sqrt(2/(ald_p*(1-ald_p)))
      
      # Setting Asymmetric LaPlace Theta 
      ald_theta = (1-2*ald_p)/(ald_p*(1-ald_p))
      
      # Setting z initial 
      #ald_z_vec = matrix(rexp(length(y), rate=1), nrow = length(y))
      ald_z_vec = rexp(length(y), rate=1)
      ald_z_mat = diag(ald_z_vec, nrow = length(y))
      
      # Setting up a placeholder for path selection variable (10 pathways)
      gamma = matrix(NA, nrow = N, ncol = num_pw)
      # Initializing path selection variable
      gamma[1,] = rep(1,num_pw)
      
      
      prob_notnull_creator = function(num_pw, N, pwy_dfs, init = 0){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = rep(init,N)
        }
        
        prob_not_null_vecs = listr
        
        return(prob_not_null_vecs)
      }
      # applying the above function
      prob_not_null_vecs = prob_notnull_creator(num_pw, N, pwy_dfs)
      
      # Setting up a placeholder for gene selection variable (10 pathways)
      curr_xi_creator = function(num_pw, N, pwy_dfs){
        listr = list()
        for(i in 1:num_pw){
          listr[[i]] = as.matrix(rbind(sample(c(1), dim(pwy_dfs[[i]])[2], replace = TRUE),matrix(NA, nrow=N-1, ncol=dim(pwy_dfs[[i]])[2])))
        }
        
        curr_xi_dfs = listr
        
        return(curr_xi_dfs)
      }
      # applying the above function
      curr_xi_dfs = curr_xi_creator(num_pw, N, pwy_dfs)
      
      if(level_1_connected == 1){
        kmi = Map('+',kmat_dfs, smooth_mats)
        kmi_s = lapply(kmi, solve) 
        marginal_alphas = Map('%*%',kmi_s, rep(list(y),num_pwy))
        kmat_alphad = as.data.frame(Map('%*%',kmat_dfs,marginal_alphas))
        colnames(kmat_alphad) = sprintf("pwy_%s",seq(1:num_pwy))
        W_mat = cor(kmat_alphad)
        W_mat = W_mat^2
        rm(kmi, kmi_s, marginal_alphas, kmat_alphad)
      }
      
      gam_mod = rep(1,num_pwy)
            
      i=2
      j=1
      for (iter_mod in 1:modar){
        for (i in 2:N){
          print(paste("i is ",i))
          alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i-1, ]})
          if(distn == "mvn"){
            for (j in 1:num_pw){
              
              # Gamma by K matrices
              gamma_lstr = as.list(matrix(gamma[i-1,],nrow = length(gamma[i-1,]), ncol=1)) 
              kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)
              
              # sampling first path coefficients from spike slab prior
              if((sum(colSums((kmat_dfs[[j]]))) == dim(kmat_dfs[[j]])[1]^2)||gam_mod[j]==0){
                alpha_mats[[j]][i,] = 0
              } else {
                alpha_mats[[j]][i,] = rmvnorm(1, as.matrix(forceSymmetric(as.matrix(ginv(sigmasq[i-1,]/sigmasq_alpha[i-1,]*kmat_dfs[[j]]+t(kmat_dfs[[j]])%*%kmat_dfs[[j]]))))%*%t(kmat_dfs[[j]])%*%(y-0.5*Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k)))+0.5*kmat_dfs[[j]]%*%alpha_mats[[j]][i-1,]),sigmasq[i-1,]*as.matrix(forceSymmetric(as.matrix(ginv(sigmasq[i-1,]/sigmasq_alpha[i-1,]*kmat_dfs[[j]]+t(kmat_dfs[[j]])%*%kmat_dfs[[j]])))))
              }
              
            }
            
            # sampling new phi
            #sigmasq[i,] = (length(y)/2*1/(0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)))
            #sigmasq[i,] = rgamma(1,length(y)/2,0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2))
            sigmasq[i,] = rinvgamma(1,length(y)/2+a,0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2))
            
          } else if(distn == "ald"){
            for (j in 1:num_pw){
              
              # Gamma by K matrices
              gamma_lstr = as.list(matrix(gamma[i-1,],nrow = length(gamma[i-1,]), ncol=1)) 
              kmat_dfs_gammad = Map('*',kmat_dfs,gamma_lstr)
              
              # sampling first path coefficients from spike slab prior
              if((sum(colSums((kmat_dfs[[j]]))) == dim(kmat_dfs[[j]])[1]^2)||gam_mod[j]==0){
                ald_bigB_inv[[j]] = diag(nx)
                alpha_mats[[j]][i,] = 0
              } else {
                ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(ald_tau^2) + diag(length(y)) #ald_bigB_inv[[j]]
                alpha_mats[[j]][i,] = solve(ald_bigB_inv[[j]])%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix(y-ald_theta*ald_z_vec)/(ald_tau^2))
              }
            }
            
            ald_delta_sq = ((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)/ald_tau^2
            ald_gamma_sq = 2 + (ald_theta^2)/(ald_tau^2)
            ald_delta = sqrt(ald_delta_sq)
            ald_gamma = sqrt(ald_gamma_sq)
            
            plc_hld_z = numeric(30)
            for (k in 1:length(ald_z_vec)){
              for(zico in 1:30){
                plc_hld_z[zico] = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
              }
              #ald_z_vec[k] = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
              ald_z_vec[k] = mean(plc_hld_z, na.rm=TRUE)
              ald_z_mat[k,k] = ald_z_vec[k]
            }
            
            # sampling new phi
            sigmasq[i,] = (length(y)/2*1/(0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)))#1/rgamma(1,length(y)/2,0.5*sum((y-kmat_1%*%alpha_1[i,]+kmat_2%*%alpha_2[i,])^2))
            
          }
          
          sigmasq_alpha[i,] = 1.0
          
          alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})
          
          if(distn == "mvn"){
            for (j in 1:num_pw){
              
              gamma_flipr = gamma_lstr
              gamma_flipr[[j]] = 1
              kmat_dfs_tstr = Map('*',kmat_dfs,gamma_flipr)
              
              if(level_1_connected == 1){  
                gamma_lstr_neg_coding = replace(unlist(gamma_lstr), unlist(gamma_lstr)==0, -1)
                delta_mat = as.matrix(gamma_lstr_neg_coding)%*%t(as.matrix(gamma_lstr_neg_coding))
                delta_mat = replace(delta_mat, delta_mat==-1, 0)
                upper_tri_delta_mat = upper.tri(delta_mat, diag = FALSE)
                prior_ising_e = exp(sum(W_mat*upper_tri_delta_mat))
                
                gamma_flipr_neg_coding = replace(unlist(gamma_flipr), unlist(gamma_flipr)==0, -1)
                delta_mat = as.matrix(gamma_flipr_neg_coding)%*%t(as.matrix(gamma_flipr_neg_coding))
                delta_mat = replace(delta_mat, delta_mat==-1, 0)
                upper_tri_delta_mat = upper.tri(delta_mat, diag = FALSE)
                prior_ising_d = exp(sum(W_mat*upper_tri_delta_mat))
                
              } else{
                prior_ising_d = 0.5
                prior_ising_e = 0.5
              }
              
              d = dmvnorm(y, (Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))),sigmasq[i,]*diag(length(y)))*prior_ising_d
              
              e = dmvnorm(y, (Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))-kmat_dfs[[j]]%*%alpha_mats[[j]][i,]),sigmasq[i,]*diag(length(y)))*prior_ising_e
              
              # Computing the probability that gamma j equals 1
              prob_not_null_vecs[[j]][i] = d/(d+e)
              
              if(is.na(prob_not_null_vecs[[j]][i]) == TRUE){
                prob_not_null_vecs[[j]][i] = 0
              } else {}
              
              # Updating our value for gamma 1
              if(prob_not_null_vecs[[j]][i]>=0.5){
                gamma[i,j] = 1
              } else {gamma[i,j] = 0}
            }
            
          } else if(distn == "ald"){
            
            for (j in 1:num_pw){
              
              gamma_flipr = gamma_lstr
              gamma_flipr[[j]] = 1
              kmat_dfs_tstr = Map('*',kmat_dfs,gamma_flipr)
              
              if(level_1_connected == 1){  
                gamma_lstr_neg_coding = replace(unlist(gamma_lstr), unlist(gamma_lstr)==0, -1)
                delta_mat = as.matrix(gamma_lstr_neg_coding)%*%t(as.matrix(gamma_lstr_neg_coding))
                delta_mat = replace(delta_mat, delta_mat==-1, 0)
                upper_tri_delta_mat = upper.tri(delta_mat, diag = FALSE)
                prior_ising_e = exp(sum(W_mat*upper_tri_delta_mat))
                
                gamma_flipr_neg_coding = replace(unlist(gamma_flipr), unlist(gamma_flipr)==0, -1)
                delta_mat = as.matrix(gamma_flipr_neg_coding)%*%t(as.matrix(gamma_flipr_neg_coding))
                delta_mat = replace(delta_mat, delta_mat==-1, 0)
                upper_tri_delta_mat = upper.tri(delta_mat, diag = FALSE)
                prior_ising_d = exp(sum(W_mat*upper_tri_delta_mat))
                
              } else{
                prior_ising_d = 0.5
                prior_ising_e = 0.5
              }
              
              Ld = sum(log(dnorm((y-ald_theta*ald_z_vec)/(ald_tau^2), Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))/(ald_tau^2),ald_z_vec)),na.rm=TRUE)+log(prior_ising_d)
              
              Le = sum(log(dnorm((y-ald_theta*ald_z_vec)/(ald_tau^2), (Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))-kmat_dfs[[j]]%*%alpha_mats[[j]][i,])/(ald_tau^2),ald_z_vec)),na.rm=TRUE)+log(prior_ising_e)
              
              
              d = prod(dnorm((y-ald_theta*ald_z_vec)/(ald_tau^2), Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))/(ald_tau^2),ald_z_vec))*prior_ising_d
              
              e = prod(dnorm((y-ald_theta*ald_z_vec)/(ald_tau^2), (Reduce("+",(Map('%*%',kmat_dfs_tstr,alpha_mats_k)))-kmat_dfs[[j]]%*%alpha_mats[[j]][i,])/(ald_tau^2),ald_z_vec))*prior_ising_e
              
              
              # Computing the probability that gamma j equals 1
              
              #print(1/(1+exp(Le-Ld)))
              #print(d/(d+e))
              
              prob_not_null_vecs[[j]][i] = (1/(1+exp(Le-Ld)))
              #prob_not_null_vecs[[j]][i] = (d/(d+e))
              
              if(prob_not_null_vecs[[j]][i] == "NaN"){
                prob_not_null_vecs[[j]][i] = 0
              } else {}
              
              #if(is.na(prob_not_null_vecs[[j]][i])){
              #  prob_not_null_vecs[[j]][i] = 0
              #} else {}
              
              # Updating our value for gamma 1
              if(prob_not_null_vecs[[j]][i]>=0.5){
                gamma[i,j] = 1
              } else {gamma[i,j] = 0}
            }
          }
          
          # Gene Selection
          if(mthd == 'MCMC'){
            if((i==2)|(i%%100==0)){
              for(j in 1:num_pw){
                curr_xi_dfs[[j]][i,] = BIGM(pwy_dfs[[j]], y, curr_xi_dfs[[j]][i-1,], phi = (1/sigmasq[i,]))  
                
                # Rebuilding Pathways using only selected genes
                kmat_dfs[[j]] = GaussianKernelize(pwy_dfs[[j]], curr_xi_dfs[[j]][i,])
              }
              
            } else {
              for(j in 1:num_pw){
                curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
              }
            }
            # variational bayes implementation
          } else if(mthd == "VB"){
            if((i%%skipper==0)){
              for (j in 1:num_pw){
                
                print(paste("VB Iteration is: ", j))
          
                p_gam = rbinom(dim(pwy_dfs[[j]])[2],1,mean(prob_not_null_vecs[[j]][(max((i-5),1)):i], na.rm=TRUE))
                
                curr_xi_dfs[[j]][i,] = ((p_gam*VB(y_tilde[[j]],pwy_dfs[[j]], mus=mus[[j]], sigmasq =max(sds[[j]]),  full_dat=og_df, iter_num = i, corrmat = corr_mats[[j]], LB_mat = corr_LBs[[j]], h=0, lambda = 1.0, iters = 15, dist_true = FALSE)+(1-p_gam)*-1)>0)*1 #new v2
                
                print("VB done!")
                
                # Rebuilding Pathways using only selected genes
                kmat_dfs[[j]] = GaussianKernelize(pwy_dfs[[j]], curr_xi_dfs[[j]][i,])
              }
            } 
          }
        }
        
        MLN_gamma_results = numeric(num_pw)
        for(j in 1:num_pw){
          MLN_gamma_results[j] = mean(gamma[N-500:N,j], na.rm=TRUE) 
        }
        
        gam_mod = MLN_gamma_results
        
        for(adjuster_switch in 1:num_pw){
          if(var(alpha_mats_k[[adjuster_switch]])==0){
            gam_mod[adjuster_switch]=0
          }
        }
        
        l_stack = length(gam_mod)-1
        while(sum((gam_mod>0)*1)>ceiling(((0.1^(1/modar))^(iter_mod))*length(gam_mod)) & l_stack>=1){
          gam_mod = (MLN_gamma_results>=nth(MLN_gamma_results,l_stack))*1
          l_stack = l_stack-1
        }
        
      }
      
      kmat_dfs = kmat_creator(pwy_dfs, num_pw) # recreating the matrices 
      
      curr_xi_dfs = curr_xi_creator(num_pw, N, pwy_dfs)
      
      rel_method = "mi"
      
      if(rel_method == "cor"){
        y_tilde = list()
        for(i in 1:num_pw){
          y_tilde[[i]] = as.data.frame(as.numeric(cor(y,pwy_dfs[[i]])^2))
        }
      } else if(rel_method == "mi"){
        y_tilde = list()
        for(i in 1:num_pw){
          y_tilde[[i]] = as.data.frame(matrix(NA, nrow = (dim(pwy_dfs[[i]])[2]), ncol = 1))
          for(j in 1:(dim(pwy_dfs[[i]])[2])){
            y_tilde[[i]][j,1] = mutinformation(discretize(y),discretize(pwy_dfs[[i]][j]), method="emp")
          }
        }
      }
      
      sd_null_vec = numeric()
      for(i in 1:num_pw){
        if(gam_mod[i]==0){
          sd_null_vec = append(sd_null_vec,y_tilde[[i]])  
        }
      }
      
      sd_null_vec = unlist(sd_null_vec)
      mad_null = mad(sd_null_vec)
      
      
      mus = list()
      sds = list()
      for(i in 1:num_pw){
        clustr1 = GMM(y_tilde[[i]],gaussian_comps = 2, dist_mode = "eucl_dist", seed_mode = "random_subset",km_iter = 10,em_iter = 5)  
        mus[[i]] = clustr1$centroids
        sds[[i]] = clustr1$covariance_matrices
      }
      
      
      for (i in 2:N){
        print(paste("i2 is ",i))
        alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i-1, ]})
        if(distn == "mvn"){
          for (j in 1:num_pw){
            
            # Gamma by K matrices
            kmat_dfs_gammad = Map('*',kmat_dfs,gam_mod)
            
            # sampling first path coefficients from spike slab prior
            if((sum(colSums((kmat_dfs[[j]]))) == dim(kmat_dfs[[j]])[1]^2)||gam_mod[j]==0){
              alpha_mats[[j]][i,] = 0
            } else {
              alpha_mats[[j]][i,] = rmvnorm(1, as.matrix(forceSymmetric(as.matrix(ginv(sigmasq[i-1,]/sigmasq_alpha[i-1,]*kmat_dfs[[j]]+t(kmat_dfs[[j]])%*%kmat_dfs[[j]]))))%*%t(kmat_dfs[[j]])%*%(y-0.5*Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k)))+0.5*kmat_dfs[[j]]%*%alpha_mats[[j]][i-1,]),sigmasq[i-1,]*as.matrix(forceSymmetric(as.matrix(ginv(sigmasq[i-1,]/sigmasq_alpha[i-1,]*kmat_dfs[[j]]+t(kmat_dfs[[j]])%*%kmat_dfs[[j]])))))
            }
            
          }
          
          # sampling new phi
         sigmasq[i,] = rinvgamma(1,length(y)/2+a,0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2))
          
        } else if(distn == "ald"){
          for (j in 1:num_pw){
            
            # Gamma by K matrices
            kmat_dfs_gammad = Map('*',kmat_dfs,gam_mod)
            
            # sampling first path coefficients from spike slab prior
            if((sum(colSums((kmat_dfs[[j]]))) == dim(kmat_dfs[[j]])[1]^2)||gam_mod[j]==0){
              ald_bigB_inv[[j]] = diag(nx)
              alpha_mats[[j]][i,] = 0
            } else {
              ald_bigB_inv[[j]] = t(kmat_dfs[[j]])%*%solve(ald_z_mat)%*%kmat_dfs[[j]]/(ald_tau^2) + diag(length(y)) #ald_bigB_inv[[j]]
              alpha_mats[[j]][i,] = solve(ald_bigB_inv[[j]])%*%(kmat_dfs[[j]]%*%(solve(ald_z_mat))%*%as.matrix(y-ald_theta*ald_z_vec)/(ald_tau^2))
            }
          }
          
          ald_delta_sq = ((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)/ald_tau^2
          ald_gamma_sq = 2 + (ald_theta^2)/(ald_tau^2)
          ald_delta = sqrt(ald_delta_sq)
          ald_gamma = sqrt(ald_gamma_sq)
          
          plc_hld_z = numeric(30)
          for (k in 1:length(ald_z_vec)){
            for(zico in 1:30){
              plc_hld_z[zico] = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
            }
            #ald_z_vec[k] = rgig(n=1, lambda = 0.5, chi = ald_delta_sq[k], psi = ald_gamma_sq)
            ald_z_vec[k] = mean(plc_hld_z, na.rm=TRUE)
            ald_z_mat[k,k] = ald_z_vec[k]
          }
          
          # sampling new phi
          sigmasq[i,] = (length(y)/2*1/(0.5*sum((y-Reduce("+",(Map('%*%',kmat_dfs_gammad,alpha_mats_k))))^2)))#1/rgamma(1,length(y)/2,0.5*sum((y-kmat_1%*%alpha_1[i,]+kmat_2%*%alpha_2[i,])^2))
          
        }
        
        sigmasq_alpha[i,] = 1.0
        
        alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[i, ]})

        length_vec = numeric()
        for (j in 1:num_pw){
          if (gam_mod[j]==0){
            next
          } else {
            length_vec = append(length_vec, length(curr_xi_dfs[[j]]))
          }
        }
        
        denom_prior_scaler = min(length_vec)
        
        # Gene Selection
        if(mthd == 'MCMC'){
          if((i==2)|(i%%100==0)){
            for(j in 1:num_pw){
              curr_xi_dfs[[j]][i,] = BIGM(pwy_dfs[[j]], y, curr_xi_dfs[[j]][i-1,], phi = (1/sigmasq[i,]))  
              
              # Rebuilding Pathways using only selected genes
              kmat_dfs[[j]] = GaussianKernelize(pwy_dfs[[j]], curr_xi_dfs[[j]][i,])
            }
            
          } else {
            for(j in 1:num_pw){
              curr_xi_dfs[[j]][i,] = curr_xi_dfs[[j]][i-1,]
            }
          }
          # variational bayes implementation
        } else if(mthd == "VB"){
          if((i%%skipper==0)){
            for (j in 1:num_pw){
              if (gam_mod[j]==0){
                curr_xi_dfs[[j]][i,]=0
              } else {
                
                fr <- function(x) {
                  stddev = x
                  abs(qnorm(0.975, mean=0, sd=stddev)-qnorm(0.05, mean = max(y_tilde[[j]]), sd=stddev))
                }
                
                var_est_opt_t1 = ((optimize(fr, c(0,1)))$minimum)^2
                
                print(paste("VB2 Iteration is: ", j))
                
                p_gam = rbinom(dim(pwy_dfs[[j]])[2],1,gam_mod[j])
                
                prior_scaler = denom_prior_scaler/length(curr_xi_dfs[[j]])
                
                curr_xi_dfs[[j]][i,] = ((p_gam*VB(y_tilde[[j]],pwy_dfs[[j]], mus=y_tilde[[j]], sigmasq = var_est_opt_t1, full_dat=og_df, sd_null = mad_null, prior_scaler = 1, iter_num = i, corrmat = corr_mats[[j]], LB_mat = corr_LBs[[j]], h=0, lambda = 1.0, iters = 15, dist_true = FALSE)+(1-p_gam)*-1)>0)*1 #new v2
              }

              print("VB2 done!")
              # Rebuilding Pathways using only selected genes
              kmat_dfs[[j]] = GaussianKernelize(pwy_dfs[[j]], curr_xi_dfs[[j]][i,])
            }
          } 
        }
      }
      
      mln_indic = 1
      MLN_results = numeric(dim(dat)[2]-1)
      for(j in 1:num_pw){
        for(k in 1:dim(pwy_dfs[[j]])[2]){
          MLN_results[mln_indic] = mean(curr_xi_dfs[[j]][(N-1000):N,k], na.rm=TRUE)
          mln_indic = mln_indic+1
        }
      }
      
      MLN_final = as.data.frame(cbind(MLN_results, MLN_results))
      
      alpha_mats_k = lapply(alpha_mats, function(x) {x <- x[(N-200):N, ]})
      alpha_mats_k = lapply(alpha_mats_k, function(x) {x <- colMeans(x, na.rm=TRUE)})
      #kmat_fin = Map('*',kmat_dfs,(colMeans(gamma[N:(N-200),])>=quantile(colMeans(gamma[N:(N-200),]),0.5))*1)
      kmat_fin = Map('*',kmat_dfs,gam_mod)
      MLN_mse = sum((Reduce("+",(Map('%*%',kmat_fin,alpha_mats_k)))-y)^2)
      print(paste("MSE MLN: ",MLN_mse))
      
      
      new <- Sys.time() - old # calculate difference
      print(paste("Elapsed time: ", new))
      #######################################################
      #######################################################
      
      #######################################################
      #################     Lasso     #######################
      #######################################################
      # Lasso Implementation
      
      ###########################   DATA   ############################
      # Creating a y 
      y = dat$y # simulation
      
      og_df = dat[-c((length(y))),]
      group_identifier_vector = c(dat[(length(y)),])
      #group_list = c("2","3","133") #data
      group_list = paste(c(1:num_pwy)) #simulation
      
      y = y[-c(length(y))]
      y =  as.numeric(y)
      
      # Standardizing y
      y = scale(y)
      y = as.numeric(y)
      
      y = y
      X = data.matrix(og_df[,2:(dim(og_df)[2])])
      #################################################################
      
      #perform k-fold cross-validation to find optimal lambda value
      cv_model <- cv.glmnet(X, y, alpha = 1)
      #find optimal lambda value that minimizes test MSE
      best_lambda <- cv_model$lambda.min
      best_model <- glmnet(X, y, alpha = 1, lambda = best_lambda)
      lasso_results = numeric(dim(X)[2])
      c = 1
      for(i in coef(best_model)@i){
        if(i==0){
          next
        }
        lasso_results[i] = coef(best_model)@x[c]
        c = c+1
      }
      lasso_results_indic = (lasso_results>0)*1
      lasso_final = as.data.frame(cbind(lasso_results, lasso_results_indic))
      #######################################################
      #######################################################
      print(paste("MSE LASSO: ",sum((y-predict(best_model, X))^2)))
      MSE_Lasso = sum((y-predict(best_model, X))^2)
      
      #######################################################
      #############     Group Lasso     #####################
      #######################################################
      # Group Lasso Implementation
      
      ###########################   DATA   ############################
      # Creating a y 
      y = dat$y # simulation
      
      og_df = dat[-c((length(y))),]
      group_identifier_vector = c(dat[(length(y)),])
      #group_list = c("2","3","133") #data
      group_list = paste(c(1:num_pwy)) #simulation
      
      y = y[-c(length(y))]
      y =  as.numeric(y)
      
      # Standardizing y
      y = scale(y)
      y = as.numeric(y)
      
      y = y
      X = data.matrix(og_df[,2:(dim(og_df)[2])])
      #################################################################
      
      group = unlist(group_identifier_vector[2:(dim(og_df)[2])])
      fit = grpreg(X, y, group, penalty="grLasso")
      cvfit <- cv.grpreg(X, y, group, penalty="grLasso")
      compare = as.numeric(abs(coef(fit, lambda=cvfit$lambda.min)))
      group_lasso_results = compare[2:length(compare)]
      group_lasso_results_indic = (group_lasso_results>0)*1
      group_lasso_final = as.data.frame(cbind(group_lasso_results, group_lasso_results_indic))
      #######################################################
      #######################################################
      print(paste("MSE GLASSO: ",sum((y-predict(fit, X))^2)))
      MSE_Glasso = sum((y-predict(fit, X))^2)
      
      #######################################################
      ##########     Sparse Group Lasso     #################
      #######################################################
      # Sparse Group Lasso Implementation
      
      ###########################   DATA   ############################
      # Creating a y 
      y = dat$y # simulation
      
      og_df = dat[-c((length(y))),]
      group_identifier_vector = c(dat[(length(y)),])
      #group_list = c("2","3","133") #data
      group_list = paste(c(1:num_pwy)) #simulation
      
      y = y[-c(length(y))]
      y =  as.numeric(y)
      
      # Standardizing y
      y = scale(y)
      y = as.numeric(y)
      
      y = y
      X = data.matrix(og_df[,2:(dim(og_df)[2])])
      #################################################################
      
      group = unlist(group_identifier_vector[2:(dim(og_df)[2])])
      data_SGL = list(x = X, y = y) 
      cvFit = cvSGL(data_SGL, group, type = "linear")
      #sglFit = SGL(data_SGL, group, type = "linear")
      SGL_coeffs = as.data.frame(cvFit$fit$beta[,which.min(cvFit$lldiff)])
      
      sparse_group_lasso_results = SGL_coeffs
      sparse_group_lasso_results_indic = (sparse_group_lasso_results>0)*1
      sparse_group_lasso_final = as.data.frame(cbind(sparse_group_lasso_results, sparse_group_lasso_results_indic))
      #######################################################
      #######################################################
      print(paste("MSE SGL: ",sum((y-X%*%as.matrix(SGL_coeffs))^2)))
      MSE_SGL = sum((y-X%*%as.matrix(SGL_coeffs))^2)
      
      
      
      #######################################################
      #################     SSGL     ########################
      ####################################################### 
      # Spike-and-Slab Group Lassos for Grouped Regression and Sparse Generalized Additive Models
      
      ###########################   DATA   ############################
      # Creating a y 
      y = dat$y # simulation
      
      og_df = dat[-c((length(y))),]
      group_identifier_vector = c(dat[(length(y)),])
      #group_list = c("2","3","133") #data
      group_list = paste(c(1:num_pwy)) #simulation
      
      y = y[-c(length(y))]
      y =  as.numeric(y)
      
      # Standardizing y
      y = scale(y)
      y = as.numeric(y)
      
      y = y
      X = data.matrix(og_df[,2:(dim(og_df)[2])])
      #################################################################
      
      group = unlist(group_identifier_vector[2:(dim(og_df)[2])])
      G = length(unique(group))
      lambda0seq = seq(1, 1000, by=2)
      #modSSGLcv = SSGLcv(Y=y, X=X, lambda1=1, lambda0seq = lambda0seq, groups = group, nFolds = 10)
      #modSSGL = SSGLpath(Y=y, X=X, lambda1=1, lambda0=modSSGLcv$lambda0, groups = group)
      
      modSSGLspr <- tryCatch(
        {
          SSGLspr(Y=y, x=X, lambda1=1, DF=2, nFolds = 5)
        },
        error = function(e){
          0
        }
      )
      
      #modSSGLspr = SSGLspr(Y=y, x=X, lambda1=1, DF=2)
      
      SSGL_coeffs = tryCatch(
        {
          as.data.frame(modSSGLspr$nonzero)
        },
        error = function(e){
          rep(0, dim(X)[2])
        }
      )
        
        
      
      SSGL_results = SSGL_coeffs
      SSGL_final = as.data.frame(cbind(SSGL_coeffs, SSGL_results))
      print(paste("MSE SSGL: ", sum((y-modSSGLspr$predY)^2)))
      MSE_SSGL = sum((y-modSSGLspr$predY)^2)
      #######################################################
      #######################################################
      
      MLN_G_mat = MLN_G_mat+gam_mod
      MLN_mat = MLN_mat+MLN_final[,2]
      Lasso_mat = Lasso_mat+lasso_final[,2]
      GLasso_mat = GLasso_mat+group_lasso_final[,2]
      SGL_mat = SGL_mat+sparse_group_lasso_final[,2]
      SSGL_mat = SSGL_mat+SSGL_final[,2]
      
      MLN_vec_mse[s] = MLN_mse
      Lasso_vec_mse[s] = MSE_Lasso
      Glasso_vec_mse[s] = MSE_Glasso
      SGL_vec_mse[s] = MSE_SGL
      SSGL_vec_mse[s] = MSE_SSGL
      
      
      string1s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_MLN_G_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      string2s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_MLN_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      string3s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_Lasso_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      string4s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_GLasso_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      string5s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_SGL_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      string6s = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_SSGL_mat", "_P", as.character(num_pwy),"_N", as.character(nx), "_S", as.character(s),".csv", sep = "")
      
      write.csv(MLN_G_mat, string1s, row.names = FALSE)
      write.csv(cbind(MLN_final[,2],group), string2s, row.names = FALSE)
      write.csv(cbind(lasso_final[,2],group), string3s, row.names = FALSE)
      write.csv(cbind(group_lasso_final[,2],group), string4s, row.names = FALSE)
      write.csv(cbind(sparse_group_lasso_final[,2],group), string5s, row.names = FALSE)
      write.csv(cbind(SSGL_final[,2],group), string6s, row.names = FALSE)
      
      write.csv(colnames(dat[,2:dim(dat)[2]]), "G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/actual_gene_names.csv", row.names = FALSE)
      write.csv(group, "G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/actual_gene_groups.csv", row.names = FALSE)
    }
    
    string1 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_MLN_G_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    string2 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_MLN_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    string3 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_Lasso_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    string4 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_GLasso_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    string5 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_SGL_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    string6 = paste("G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/RAP_SSGL_mat", "_P", as.character(num_pwy),"_N", as.character(nx),".csv", sep = "")
    
    write.csv(MLN_G_mat, string1, row.names = FALSE)
    write.csv(MLN_mat, string2, row.names = FALSE)
    write.csv(Lasso_mat, string3, row.names = FALSE)
    write.csv(GLasso_mat, string4, row.names = FALSE)
    write.csv(SGL_mat, string5, row.names = FALSE)
    write.csv(SSGL_mat, string6, row.names = FALSE)
    
    MSE_ph_mat = as.data.frame(cbind(MLN_vec_mse,Lasso_vec_mse,Glasso_vec_mse,SGL_vec_mse,SSGL_vec_mse))
    write.csv(MSE_ph_mat, "G:/My Drive/VTech/PhD Research/P1/Simulations_mod7_ald2/MSE_ph_mat.csv", row.names = FALSE)
  }
}
