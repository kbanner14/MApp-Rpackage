# `{MAPP}`

The `{MAPP}` package provides plotting functions to help researchers digest results from the model averaging procedure when applied to partial regression coefficients. There are four major plotting functions:  

- `MApp_bms()` works with `bma` objects obtained from the `bms()` function @BMS (details about the `bms()` function can be found in `vignette("BMS")`.
- `MApp_MCMC` works with default output from the implementation of model averaging using the `OpenBUGS` @BUGS or self programmed RJMCMC samplers. 
- `MApp()` works with lists of posterior draws. 
- `MApp_IC()` works with approximate posterior model probabilities estimated with AIC or BIC and estimates and standard errors of partial regression coefficients from all individual models. 

#Load `{MAPP}` 

`{MAPP}` can be loaded in three steps: 

1. Download all files in the `MApp/` repository from [here](https://github.com/kbanner14/MApp-Rpackage/tree/master/MApp)
2. Set your working directory to the location of the `MApp` repository on your computer. 
3. Run `devtools::load_all(".")` to load the package. Note that you must have the package `devtools` installed (to install run `install.packages("devtools")`).

# Use `{MAPP}`

Documentation is provided for all functions through `help()` and examples for all plotting functions are provided in the package vignette `MApp\vignettes\my-vignette.Rmd`. The vignette is viewable within this repository by navigating to `Mapp\vignettes`. 


