#' @title Compute AICc for a model
#' @param x an object of the class \code{lm} or \code{glm}. 
#' @return AICc for the model fit. 

aicc <- function(x){
  n <- length(x$fitted) 
  k <- length(coef(x)) + 1
  out <- AIC(x) + (2*k*(k+1))/(n-k-1)
  return(out)
}

#' @title Compute approximate posterior model probabilities using the AIC
#'   approximation for a set of models.
#' @param AIC_vec A vector of AIC values for a set of models.
#' @return The approximate posterior model probability. (Formulae
#' according to Link & Barker, 2010; Burhnam and Anderson, 2002)

pmp_AIC <- function(AIC_vec){
  top <- min(AIC_vec)
  delta <- AIC_vec - top
  out <- exp(-.5*delta)/sum(exp(-.5*delta))
  return(out)
}

#' @title Compute approximate posterior model probabilities using the BIC
#'   approximation for a set of models.
#' @param BIC_vec A vector of BIC values for a set of models.
#' @param mod_prior A vector of prior probabilities for the set of models.
#' @return The approximate posterior model probability. (Formulae
#' according to Link & Barker, 2010; Burhnam and Anderson, 2002)

pmp_BIC <- function(BIC_vec, mod_prior){
  top <- min(BIC_vec)
  delta <- BIC_vec - top
  out <- exp(-.5*delta*mod_prior)/sum(exp(-.5*delta*mod_prior))
  out
}

#' @title Compute model averaged estimates for partial regression
#'   coefficients.
#' @description This function computes point estimates for model averaged
#'   partial regression coefficients using a weighted averaged of maximum
#'   likelihood estimates from individual models according to approximate
#'   posterior model probabilities.
#' @param weights A vector of approximate posterior model probabilities
#'   (commonly referred to as model weights in the IC literature).
#' @param coef_mat A \eqn{J} by \eqn{p} matrix of maximum likelihood
#'   estimates for the partial regression coefficients from all individual
#'   models considered.
#' @param inmat A matrix defining which variables are in each model.
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#'   Weights are normalized conditional on the models meeting this criterion. 
#' @return The point estimate of the model averaged partial regression
#'   coefficients.
est_MA <- function(weights, coef_mat, inmat, w_plus = FALSE){
  inmat <- cbind(rep(1,dim(inmat)[1]), inmat)
  if(w_plus == FALSE){
  out <- coef_mat*weights
  out <- apply(out, 2, sum)
  } else {
    out <- rep(NA, dim(coef_mat)[2])
    for(i in 1:dim(coef_mat)[2]){
      mods <- which(inmat[,i] == 1)
      w_norm <- weights[mods]/(sum(weights[mods]))
      out[i] <- coef_mat[mods, i] %*% w_norm
    }
  }
  return(out)
}

#' @title Compute the standard error for model averaged
#'   partial regression coefficients.
#' @description This function computes standard errors for model averaged
#'   partial regression coefficients using the formula specified in
#'   Link and Barker 2010 (see also Burhnam and Anderson (2002)).
#' @param weights A vector of approximate posterior model probabilities
#'   (commonly referred to as model weights in the IC literature).
#' @param coef_mat A \eqn{J} by \eqn{p} matrix of maximum likelihood
#'   estimates for the partial regression coefficients from all individual
#'   models considered.
#'@param se_mat A \eqn{J} by \eqn{p} matrix of maximum likelihood estimates
#'  of standard errors for the partial regression coefficients from all
#'  individual models considered.
#' @param inmat A matrix with dimensions \eqn{J \times p} specifying which
#'   variables are included in each model. The order of models must correspond
#'   to the order of \code{coef_mat} a \code{ses_mat} (it is recommended, but not
#'   required to sort these in decreasing posterior model probability).
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#'   Weights are normalized conditional on the models meeting this criterion.
#' @return The point estimate of the model averaged partial regression
#'   coefficients.
se_MA <- function(weights, coef_mat, se_mat, inmat, w_plus = FALSE){
  
  ests_ma <- est_MA(weights, coef_mat, inmat, w_plus = w_plus)
  inmat <- cbind(rep(1,dim(inmat)[1]), inmat)
  ests_mat <- matrix(rep(ests_ma, dim(coef_mat)[1]),
                       nrow = dim(coef_mat)[1], byrow = T)
  var_mat <- se_mat^2
  if(w_plus == TRUE){
    out <- numeric(dim(coef_mat)[2])
    dev_mat <- matrix(0, nrow = dim(coef_mat)[1], ncol = dim(coef_mat)[2])
    for(i in 1:dim(coef_mat)[2]){
      mods <- which(inmat[,i] == 1)
      w_norm <- weights[mods]/(sum(weights[mods]))
      devs <- (coef_mat[mods,i] - ests_ma[i])^2
      out[i] <- w_norm%*%sqrt(var_mat[mods,i] + devs)
    }
  } else {
    dev_mat <- (coef_mat - ests_mat)^2
    out <- weights%*%sqrt(var_mat + dev_mat)
    }
  
  return(out)
}

#' @title Compute individual and model averaged results using the AIC
#'   and BIC approximations.
#' @description Create a \code{data.frame} summarizing individual and
#'   model averaged results for partial regression coefficients from
#'   a set of models.
#' @param inmat A matrix with dimensions \eqn{K \times p} specifying which
#'   variables are included in each model.
#' @param Xmat A matrix of the potential covariates to consider
#'   (dimension \eqn{n \times p}).
#' @param Yvec A vector of responses (dimension \eqn{n \times 1}).
#' @param mod_names An optional vector specifying the names of the \eqn{J}
#'   models being considered.
#' @param mod_prior A vector of prior probabilities for the set of models
#'   for computing BIC model weights. If \code{NULL}, the prior defaults to
#'   \eqn{\frac{1}{J}} for all \eqn{J} models.
#' @param family An input specifying the \eqn{glm} family to use in
#'   estimation and computation of IC criterion.
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#'   Weights are normalized conditional on the models meeting this criterion.
#' @return A \code{data.frame} with information about the following variables
#'   for all models in the model set and the model averaged result: Model,
#'   approximate weights for AIC and BIC, AIC, BIC, and ML coefficient estimates
#'   and their associated standard errors. Pieces of this output can be used with
#'   the \code{MApp_IC} function.

approx_pmp <- function(inmat, Xmat, Yvec, mod_names = NULL,
                       mod_prior = NULL, family = "gaussian", 
                       w_plus = FALSE, aic_c = TRUE){
  #function to compute
  all_dat <- data.frame(Yvec,Xmat)
  num_x <- dim(Xmat)[2]
  num_mod <- dim(inmat)[1]
  if(is.null(mod_prior))
    mod_prior = rep(1/num_mod, num_mod)
  if(is.null(mod_names))
    mod_names = paste0("M", 1:num_mod)
  pmpA <- rep(NA, num_mod)
  pmpB <- rep(NA, num_mod)
  aic <- rep(NA, num_mod)
  bic <- rep(NA, num_mod)
  allModels <- as.list(1:num_mod)

  for(i in 1:num_mod){
    if (sum(inmat[i,]) != 0) {
      xs1 <- names(Xmat[which(inmat[i,] == 1)])
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1", xs, sep = "+")
      fit <- glm(form, family = family, data = all_dat)
    } else {
      fit <- glm(Yvec ~ 1, family = family, data = all_dat)
    }
    ests <- summary(fit)$coefficients[,1]
    ses <- summary(fit)$coefficients[,2]
    allModels[[i]] <- data.frame(ests, ses)
    if(aic_c == TRUE){
      aic[i] <- aicc(fit)
    } else {
      aic[i] <- AIC(fit)
    }
    bic[i] <- BIC(fit)
  }

  pmpA <- pmp_AIC(aic)
  pmpB <- pmp_BIC(bic, mod_prior)

  coef_mat <- matrix(0, nrow = num_mod, ncol = num_x + 1)
  se_mat<- matrix(0, nrow = num_mod, ncol = num_x + 1)
  
  for(ndx in 1:num_mod){
    coef_mat[ndx, c(1, which(inmat[ndx, ] == 1)+1)] <- allModels[[ndx]][,1]
    se_mat[ndx, c(1, which(inmat[ndx, ] == 1)+1)] <- allModels[[ndx]][,2]
  }
  
  ma_estsAIC <- est_MA(pmpA, coef_mat, inmat = inmat, w_plus = w_plus)
  ma_estsBIC <- est_MA(pmpB, coef_mat, inmat = inmat, w_plus = w_plus)
  
  ma_seAIC <- se_MA(pmpA, coef_mat, se_mat, inmat, w_plus)
  ma_seBIC <- se_MA(pmpB, coef_mat, se_mat, inmat, w_plus)
  
  
  coef_mat <- data.frame(round(coef_mat, digits = 4))
  idx_full <- which(apply(inmat, 1, sum) == num_x)
  names(coef_mat) <- paste("Est", c("Int",names(Xmat[which(inmat[idx_full,]==1)])),
                           sep = "_")

  se_mat <- data.frame(round(se_mat, digits = 4))
  names(se_mat) <- paste("SE", c("Int",names(Xmat[which(inmat[idx_full,]==1)])),
                         sep = "_")
  
  if(aic_c == TRUE){
    mod_names <- c("MA_AICc", "MA_BIC", mod_names)
  } else{
    mod_names <- c("MA_AIC", "MA_BIC", mod_names)
  }

  pmpA <- c("MA", "MA", round(pmpA, digits = 4))
  pmpB <- c( "MA", "MA", round(pmpB, digits = 4))
  aic <- c("--", "--", round(aic, digits = 4))
  bic <- c("--", "--", round(bic, digits = 4))

  coef_mat <- rbind(ma_estsAIC, ma_estsBIC, coef_mat)
  se_mat <- rbind(as.vector(ma_seAIC), as.vector(ma_seBIC), se_mat)
  
  # tell the user what type of MA they used
  disp_type <- ifelse(aic_c == TRUE & type == "AIC", "AICc", 
                      ifelse(type == "BIC", "BIC", "AIC"))
  disp_wplus <- ifelse(w_plus == TRUE, "Models where coefficient is not set to 0", 
                       "All models in the model set")
  message(paste("Type of information criteria used in MAP plot:", disp_type))
  message(paste("Subset of models used in MA:", disp_wplus))
  if(aic_c == TRUE){
    return(
      data.frame(Model = mod_names, AICc_pmp = pmpA,
                 BIC_pmp = pmpB,
                 AICc = aic,
                 BIC = bic,
                 coef_mat,
                 se_mat)
    )
    
  } else {
    return(
      data.frame(Model = mod_names, AIC_pmp = pmpA,
                 BIC_pmp = pmpB,
                 AIC = aic,
                 BIC = bic,
                 coef_mat,
                 se_mat)
    )
  }
  
}
