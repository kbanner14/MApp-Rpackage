# `{MAPP}`

The `{MAPP}` package provides plotting functions to help researchers digest results from the model averaging procedure when applied to partial regression coefficients. There are four major plotting functions:  

- `MApp_bms()` works with `bma` objects obtained from the `bms()` function in the `BMS` package (details about the `bms()` function can be found in the `BMS` tutorial [here: http://bms.zeugner.eu/tutorials/bms.pdf](http://bms.zeugner.eu/tutorials/bms.pdf)).
- `MApp_MCMC` works with default output from the implementation of model averaging using the `OpenBUGS` or self programmed RJMCMC samplers. 
- `MApp_AIC()` works with a data frame, conducts AIC or AICc based model averaging for all-subsets regression. Returns the MAP plot and Model averaged results. 
- `MApp_IC()` works with approximate posterior model probabilities estimated with AIC or BIC and estimates and standard errors of partial regression coefficients from all individual models. 

# Load `{MAPP}` 

`MAPP` is not currently availble on CRAN, so it must be installed from GitHub or loaded from your local disk. The package, `devtools` is necessary for both installs. To install `devtools`, run `install.packages("devtools")`.  

## GitHub Install
1. Install the latest version of `{MAPP}` from GitHub. _[maybe I should have a version number 0.0.3?]_

```{r install, echo = T, message = FALSE}
devtools::install_github("kbanner14/MApp-Rpackage", subdir = "MApp")
library(MApp)
```

2. Install package dependencies: `LearnBayes`, `beanplot`, `magritter`, `dplyr`, and `BMS`. Use `install.packages("packagename")` to install these packages.

## `{MAPP}` Local Disk

1. Install package dependencies: `LearnBayes`, `beanplot`, `dplyr`, and `BMS`. Use `install.packages("packagename")` to install these packages.
2. Download all files in the `MApp-Rpackage/` repository from the _Download Zip_ button on the right hand side of the screen [here](https://github.com/kbanner14/MApp-Rpackage).
3. Set your working directory to the location of the `MApp-Rpackage` repository on your computer. 
4. Run `devtools::load_all("MApp/.")` to load the package. Note that you must have the packages `devtools` and `roxygen2` installed (to install, run `install.packages(c("devtools", "roxygen2"))`).
5. Load the documentation by running `devtools::document(".")`

# Use `{MAPP}`

Documentation is provided for all functions through `help()` and examples for all plotting functions are provided in the package vignette `MApp\vignettes\my-vignette.Rmd`. The vignette is viewable within this repository by navigating to `Mapp\vignettes`. 


