## ----install, echo = TRUE, message = FALSE, eval = FALSE-----------------
#  # install the package
#  library(devtools)
#  install_github("kbanner14/MApp-Rpackage", subdir = "MApp")
#  # load the package
#  library(MApp)

## ----echo = TRUE, eval=TRUE, results='hide', message=FALSE---------------
# load the pckage
devtools::load_all(".")

# load the R documentation files
devtools::document(".")

## ----test, echo = TRUE, eval = FALSE-------------------------------------
#  ? MApp_bms

## ---- brain, results='asis', echo = TRUE, eval = FALSE-------------------
#  data(brainData)
#  head(brainData)

## ----brain2, echo = FALSE, results='asis'--------------------------------
knitr::kable(head(brainData, 5))

## ----bmsobj--------------------------------------------------------------
library(BMS)
brain <- bms(brainData, g = "UIP", 
                    mprior = "uniform", 
                    user.int = F)
class(brain)

## ----bmsfig, fig.show = 'hold', fig.width = 7, fig.height = 5, fig.cap="MAP plot for the partial regression coefficients associated with each explanatory variable in the brainData dataset", message=FALSE, results='hide'----
map_brain <- MApp_bms(x = brain, plot_wind = c(1,3),
                      mod_names = c("BGL", "BG", "BL", "B", "G", "GL", "L", "Null"))

## ----echo = TRUE, eval = FALSE-------------------------------------------
#  map_brain

## ----echo = FALSE, results = "asis"--------------------------------------
knitr::kable(map_brain, digits = 3)

## ---- results='asis', eval = FALSE---------------------------------------
#  data(bfat)
#  head(bfat)

## ---- echo = FALSE, results='asis'---------------------------------------
knitr::kable(head(bfat, 5))

## ----message=FALSE, results='hide'---------------------------------------
# create bma object using bfat data
bfat_bms <- bms(bfat, mprior = "uniform", g = "UIP", user.int = F)

## ------------------------------------------------------------------------
# see how much of the posterior model mass is accounted for by the top 500
# models
sum(pmp.bma(bfat_bms)[,1])

## ----bfatplot, fig.width=7, fig.height= 5, message= FALSE, results='hide'----
map_bfat <- MApp_bms(bfat_bms, plot_wind = c(1,3), 
                        include_coef = c(13, 2, 4), 
                        max_display = 10, 
                        num_sims = 1000)

## ----bfat-tab, echo = TRUE, eval = FALSE---------------------------------
#  map_bfat

## ----bfat-tab2, echo = FALSE, results = 'asis'---------------------------
knitr::kable(map_bfat, digits = 3)

## ---- echo = FALSE, results = 'hide',fig.width= 8, fig.height=5, eval = FALSE, message=FALSE----
#  MApp_IC(brainData, plot_wind = c(1,3),
#          type = "BIC", w_plus = FALSE)

## ---- echo = FALSE, results = 'asis',fig.width= 8, fig.height=5----------
brain_bic <- MApp_IC(brainData, plot_wind = c(1,3), 
                     type = "BIC", w_plus = FALSE)
knitr::kable(brain_bic[,c(1:5)], digits = 3)
knitr::kable(brain_bic[,c(1,c(6:9))], digits = 3)
knitr::kable(brain_bic[,c(1,c(10:13))], digits = 3)

## ----echo = TRUE, results= "asis",fig.width= 8, fig.height=5, eval = FALSE----
#  MApp_IC(brainData, plot_wind = c(1,3),
#          type = "BIC", w_plus = TRUE)

## ----echo = FALSE, results= "asis",fig.width= 8, fig.height=5------------
brain_bic2 <- MApp_IC(brainData, plot_wind = c(1,3), 
                     type = "BIC", w_plus = TRUE)

## ----ICbrain, echo = TRUE, results='asis'--------------------------------
post_means <- t(brain$topmod$betas())
post_means[post_means != 0] <- 1
inmat <- post_means

IC_approx <- approx_pmp(inmat, Xmat = brainData[,-1], Yvec = brainData[,1])
knitr::kable(IC_approx[,c(1:5)], digits = 3)
knitr::kable(IC_approx[,c(1,6:9)], digits = 3)
knitr::kable(IC_approx[,c(1,10:13)], digits = 3)

## ----echo = TRUE, fig.width= 8, fig.height=5-----------------------------
# BIC weights
x_coef <- IC_approx[c(2:10), c(7:9)]
x_se <- IC_approx[c(2:10), c(11:13)]
pmp <- as.numeric(as.character(IC_approx[c(3:10),2]))
MApp_IC_gen(x_coef = x_coef, x_se = x_se, pmp = pmp, inmat = inmat)

## ----mcmc-matrix, echo = FALSE, warning=FALSE, message=FALSE-------------
set.seed(12)
mcmc_mat <- data.frame(LearnBayes::rmnorm(700, mean = c(-1,3,6), varcov = diag(3)))
mcmc_mat$model <- rep(paste0("M",1:7), each = 100)
mcmc_mat <- mcmc_mat[,c(4,1,2,3)]
inmat <- matrix(c(1,1,1,
  			  1,1,0, 
				  1,0,1,
				  0,1,1,
				  1,0,0,
				  0,1,0,
				  0,0,1), 
				 byrow = T, 
				 nrow = 7, 
				 ncol = 3)

for(i in 1:7){
  idx <- which(mcmc_mat[,1] == paste0("M",i))
  mcmc_mat[idx, c(2:4)] <- t(apply(mcmc_mat[idx,c(2:4)], 1 , function(x){x*inmat[i,]}))
  }

# shuffle order
mcmc_mat <- mcmc_mat[sample(1:700, size = 700, rep = F),]

## ----MCMC, results='asis', eval = FALSE, echo = TRUE---------------------
#  head(mcmc_mat)

## ----MCMC2, results='asis', echo = FALSE---------------------------------
knitr::kable(head(mcmc_mat))

## ----MCMCmat, fig.width=7, fig.height=5, echo=TRUE-----------------------
MApp_MCMC(mcmc_mat, plot_wind = c(1,3))

## ----MCMCcommon, fig.width=7, fig.height=5, echo=TRUE--------------------
MApp_MCMC(mcmc_mat, plot_wind = c(1,3), max_display = "common3")

## ---- echo = FALSE, message = FALSE--------------------------------------
data_sim <- function(n = 60, p = 5, betas = c(0,0,0,1,1.2), 
                    sig_y = 2.5, sig_x = rep(1, 5), cor_x = 0,
                    cor_vars = c(3,5), tol = 0.05, truth4 = F, 
                    b42 = 0.5, b43 = 0.25, ...){
  # for now we are just going to assume that two of the variables are 
  # correlated X3 and X5, can change this with cor_vars
  try(if(length(betas) != p)  
    stop("dimensions of betas and number of parameters (p) must agree"))
  
    # draw iid normal x's, then use cor matrix (R) and Cholesky 
    # factorization to get desired correlation structure in x's
    
    var_x <- sig_x^2
    cov <- diag(var_x)
    R <- diag(p)
    R[cor_vars[1], cor_vars[2]] <- cor_x
    R[cor_vars[2], cor_vars[1]] <- cor_x
    U <- t(chol(R))
    
    
    # simulate the x matrix 
    Xmat <- LearnBayes::rmnorm(n, mean = betas, varcov = cov)
    Xmat <- t(U %*% t(Xmat))
    c <- abs(cor_x - cor(Xmat[,cor_vars[1]], Xmat[, cor_vars[2]]))
    
    while(c >= tol) {
      Xmat <- LearnBayes::rmnorm(n, mean = betas, varcov = cov)
      Xmat <- t(U %*% t(Xmat))
      c <- abs(cor_x - cor(Xmat[,cor_vars[1]], Xmat[, cor_vars[2]]))
    }
    
    dimnames(Xmat)[[2]] <- paste0("X", 1:p)
    
    # generate y 
    epsilon <- rnorm(n, 0, sig_y)
    # check for truth 4 scenario
    if(truth4 == F){
      Y <- Xmat%*%matrix(betas) + epsilon
    } else {
        X4sq <- Xmat[,4]^2
        X4cu <- Xmat[,4]^3
        Xmat <- cbind(Xmat, X4sq, X4cu)
        betas <- c(betas, b42, b43)
        Y <- Xmat%*%matrix(betas) + epsilon
        Xmat <- Xmat[,1:5]
      }
    # make the dat frame
    dat <- scale(data.frame(Y,Xmat), center = T, scale = T)
    dat <- data.frame(dat)
    return(dat)
    
}

fake_dat <- data_sim(cor_x = 0.8)
fd <- bms(fake_dat, g = "UIP", mprior = "uniform", user.int = F)
inmat <- t(fd$topmod$betas())
inmat[inmat !=0 ] <- 1
mcmc_list <- bms_post_sim(inmat, Xmat = fake_dat[,-1], Yvec = fake_dat[,1], num_sims = 1000)
names(mcmc_list) <- paste("M", 1:32, sep = "_")
samps <- round(pmp.bma(fd)[,2]*1000)
for(i in 1:32){
  mcmc_list[[i]] <- mcmc_list[[i]][1:samps[i],]
}
mcmc_list <- mcmc_list[1:23]

## ---- fig.width = 10, fig.height= 6, out.width= ".9\\linewidth"----------
MApp_list(mcmc_list, plot_wind = c(1,5), max_display = 10)

