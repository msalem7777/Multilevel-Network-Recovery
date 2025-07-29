# **MLNR: Multilevel Network Recovery**  

## **Overview**  
MLNR is an R package designed for recovering multilevel networks from observed data. It provides methods for estimating network structures at multiple levels, integrating hierarchical dependencies, and performing statistical inference on complex networked systems.  

## **Features**  
- Multilevel network estimation from observational data  
- Support for hierarchical dependencies and latent structures  
- Statistical inference tools for network recovery  
- Visualization functions for network structures  

## **Installation**  
To install MLNR from GitHub, use:  

```r
# Install devtools if not already installed
install.packages("devtools")
library(devtools)

# Install MLNR from GitHub
devtools::install_github("msalem7777/Multilevel-Network-Recovery")
```

For a quick demo, run:

```r
library(MLNR)

demo_nonlinear = MLNR::demo_code_NonLinear()
demo_linear = MLNR::demo_code_Linear()
```

## **Usage**
To load the package and start recovering networks:

```r
library(MLNR) # load package

# df = # ... Insert appropriate command to load your data as a dataframe in R

# The algorithm assumes a specific format for your dataframe. The response variable is assumed to be in the first column. 
# The group identifiers are assumed to be in the last row. You can comment out the below lines for a sample dataframe.

# df = MLNR::GenSimDatLinear("mvn", n_samples=10, num_sets=3, distn_mean = 5, distn_sd = 1, skew = 0)
# View(df)

# Fit model
num_sets = #... Manually enter in the number of sets
model = MLNR(df, num_sets) # You can optimize the fit function's parameters for a better performance on your dataset 

# You can plot actual vs. model prediction
plot(model$yhat, model$y, main = "Actual vs. Predicted")

# You can also create a nested network visualization.
names_vector <- paste0("S_{", 1:num_sets, "}") # Change this line to create custom set names
CreateNetworkViz(dat, model[["gamma"]], model[["all_xi"]], num_sets, set_names = names_vector) # Creates a visualization of the network. Connected elements/sets are all relevant to the response

# For predicting on a test set, you need to pass a separate test dataset
# df_pred = # ... Insert appropriate command to load your data as a dataframe in R

# Below is a sample of what df_pred should look like
# df_pred = GenSimDatLinear("mvn", 3, num_sets, distn_mean = distn_mean, distn_sd = distn_sd, skew = skew)
# df_pred = df_pred[, 2:ncol(df_pred)]
# View(df_pred)

# The package has a separate prediction function which can be called as:
predictions = MLNR.predict(df_pred, model, cov_transform = "scale", scale_up = TRUE)

# The package also has a function to produce predictions with Credible Intervals
predictions = MLNR.predictCI(df_pred, model, cov_transform = "scale", interval = "credible", scale_up = TRUE)
```

## **Contributing**
Contributions are welcome! If you encounter issues or have feature requests, please open an issue or submit a pull request.

## **License**
This package is licensed under the MIT License.
