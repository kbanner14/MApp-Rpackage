#  Model Averaged posteriors plot `{MAPP}`

The `{MAPP}` package provides plotting functions to help researchers digest results from the model averaging procedure when applied to partial regression coefficients. There are four major plotting functions:  

- `MApp_bms()` works with `bma` objects obtained from the `bms()` function @BMS (details about the `bms()` function can be found in `vignette("BMS")`.
- `MApp_MCMC` works with default output from the implementation of model averaging using the `OpenBUGS` @BUGS or self programmed RJMCMC samplers. 
- `MApp()` works with lists of posterior draws. 
- `MApp_IC()` works with approximate posterior model probabilities estimated with AIC or BIC and estimates and standard errors of partial regression coefficients from all individual models. 

# Load `{MAPP}` 

The package `{MAPP}` is contained in the directory `MApp\` within this repository. This package can be loaded by following these steps: 

1. Install package dependencies: `LearnBayes`, `beanplot`, `dplyr`, and `BMS`. Use `install.packages("packagename")` to install these packages.
2. Download all files in the `MApp-Rpackage` repository from the _Download Zip_ button on the right hand side of the screen 
3. Set your working directory to the location of the `MApp-Rpackage` repository on your computer. 
4. Run `devtools::load_all("MApp/.")` to load the package. Note that you must have the package `devtools` installed (to install `devtools`, run `install.packages("devtools")`).

# Use `{MAPP}`

Examples for all plotting functions are provided in the package vignette `MApp\vignettes\my-vignette.Rmd`. The vignette is viewable within this repository by navigating to `MApp\vignettes`. Documentation for all functions can be viewed by previewing the `.rd` files in the `MApp\man\` directory. The easiest way to do this is to open the files in RStudio and hit the preview button.




