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

## **Usage**
Load the package and start recovering networks:

```r
library(MLNR)

demo_nonlinear = MLNR::demo_code_NonLinear()
demo_linear = MLNR::demo_code_Linear()
```

## **Contributing**
Contributions are welcome! If you encounter issues or have feature requests, please open an issue or submit a pull request.

## **License**
This package is licensed under the MIT License.
