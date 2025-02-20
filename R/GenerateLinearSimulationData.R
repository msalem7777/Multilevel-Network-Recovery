#' Generate Simulated Data
#'
#' This function generates simulated data based on specified distribution types
#' and a customizable number of features and paths.
#'
#' @param distn Character. Distribution type ("mvn" for Gaussian noise, "ald" for asymmetric Laplace noise).
#' @param num_samples Integer. Number of samples to generate.
#' @param num_sets Integer. Number of pathways to generate.
#' @param distn_mean Numeric. Mean of the distribution. Default is 0.
#' @param distn_sd Numeric. Standard deviation of the distribution. Default is 1.
#' @param skew Numeric. Skewness parameter for the asymmetric Laplace distribution. Default is 0.5.
#' @return A data frame containing the generated data.
#' @details
#' The function generates features for multiple pathways using multivariate normal distributions.
#' The first two pathways have predefined covariance structures, while additional pathways
#' are generated dynamically with random correlation matrices.
#'
#' @examples
#' # Generate data with Gaussian noise
#' data <- GenSimDat("mvn", num_samples = 100, num_sets = 3)
#'
#' # Generate data with asymmetric Laplace noise
#' data <- GenSimDat("ald", num_samples = 100, num_sets = 3, skew = 0.7)
#'
#' @import mvtnorm
#' @import Matrix
#' @importFrom ald rALD
#' @export

# Simulation Generation
# suppressPackageStartupMessages(library(mvtnorm))
# suppressPackageStartupMessages(library(ald))


GenSimDatLinear = function(distn, num_samples, num_sets, distn_mean=0, distn_sd=1, skew=0.5){

  library(mvtnorm)
  library(Matrix)
  library(ald)

  # Randomly generated 10 features from first path via standard normal
  sigma1 <- matrix(c(1,0.7,-0.1,0.6,0.7,1,-0.2,0.55,-0.1,-0.2,1,-0.1,0.6,0.55,-0.1,1), ncol=4)
  first_path_features = data.frame(matrix(rmvnorm(num_samples,rep(0,4),16*sigma1), nrow = num_samples))

  # Randomly generated 10 features from Second path via standard normal
  sigma2 <- matrix(c(1,-0.1,0.1,0.15,-0.1,1,0.65,0.55,0.1,0.65,1,0.6,0.15,0.55,0.6,1), ncol=4)
  second_path_features = scale(rowMeans(first_path_features[,c(1,2,4)]))+data.frame(matrix(rmvnorm(num_samples,rep(0,4),16*sigma2), nrow = num_samples))

  if(distn=="mvn"){

    y =  2*first_path_features$X1+3*first_path_features$X2+first_path_features$X1*first_path_features$X2*first_path_features$X4-second_path_features$X3+second_path_features$X2/4+2/5*second_path_features$X4+rnorm(num_samples, mean = distn_mean, sd = distn_sd)

  } else if(distn=="ald") {

    y =  2*first_path_features$X1+3*first_path_features$X2+first_path_features$X1*first_path_features$X2*first_path_features$X4-second_path_features$X3+second_path_features$X2/4+2/5*second_path_features$X4+rALD(num_samples, mu = distn_mean, sigma = distn_sd, p = skew)

  } else {

    stop("Invalid input for `distn`. Please choose a valid distribution: `mvn` for Gaussian noise or `ald` for asymmetric Laplace noise.")
  }

  first_path_features[dim(first_path_features)[1]+1,]=1
  second_path_features[dim(second_path_features)[1]+1,]=2
  y[dim(first_path_features)[1]]=0

  dat = cbind(y,first_path_features,second_path_features)

  # Additional pathway creator
  for(i in 3:num_sets){
    lp = 10
    sigma = diag(lp)
    rho = matrix(runif(lp, 0.7, 0.9), lp)
    M1 = t(matrix(rho,nrow=lp*2-1,ncol=lp))[,1:lp]
    M1[lower.tri(M1)] = 0
    M1 = Matrix::forceSymmetric(M1)
    diag(M1) = 1
    Sigma = sigma%*%M1%*%sigma
    nam = paste("path_features_", i, sep = "")
    path_features = data.frame(matrix(rmvnorm(num_samples,rep(rnorm(1),dim(sigma)[2]),sigma), nrow = num_samples))
    path_features[dim(path_features)[1]+1,] = i
    # assign(nam, path_features)
    dat = cbind(dat, path_features)
    rm(path_features)
  }

  # Start with the first set of column names
  namr = c("y",
           paste0("X1_", 1:4),  # X1 columns
           paste0("X2_", 1:4))  # X2 columns

  # Dynamically generate column names for the remaining sets (from 3 to num_sets)
  if (num_sets >= 3) {
    for (i in 3:num_sets) {
      namr = c(namr, paste0("X", i, "_", 1:10))
    }
  }

  colnames(dat) = namr

  return(dat)
}
