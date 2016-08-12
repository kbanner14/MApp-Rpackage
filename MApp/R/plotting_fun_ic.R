#' @title Model Averaged Posteriors Plot
#' @description Works with results from information theoretic approximations to
#'   posterior model probabilities. Creates a graphical display showing all
#'   relevant output from the model averaging procedure and facilitates comparisons
#'   between individual models and their model averaged counterpart by displaying
#'   confidence intervals from all models.
#' @param x_coef A matrix of coefficient estimates for all variables from each
#'   individual model considered.
#' @param x_se A matrix of standard errors for all coefficients from each
#'   individual model considered.
#' @param pmp A vector of approximate posterior model weights for all models
#'   considered
#' @param inmat A matrix with dimensions \eqn{J \times p} specifying which
#'   variables are included in each model. The order of models must correspond
#'   to the order of \code{x_coef} a \code{x_se} (it is recommended, but not
#'   required to sort these in decreasing posterior model probability).
#' @param mod_names A vector of length \eqn{J} specifying the names of the models.
#'   If left blank, models will be named \eqn{M_1, M_2,..., M_J}.
#' @param var_names A vector of length \eqn{p} specifying the names of the
#'   varibles, if left blank the names will be determined by the column names
#'   of \code{x_coef}.
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns for specifying the partition of the plotting window.
#' @return The \code{MApp} plot.
#'
MApp_IC <- function(x_coef, x_se, pmp, inmat, var_names = NULL,
                    mod_names = NULL, plot_wind = NULL){
  # assuming x is a matrix of estimates and standard errors
  # from each individual model (row) and the model averaged
  # result.

  k <- nrow(x_coef) -1
  p <- ncol(x_coef)
  pmp <- round(as.numeric(as.character(pmp)), 4)
  pep_text <- ifelse(pep == 0, "<0.0000", as.character(pep))
  pip <- t(inmat)%*%pmp
  pep <- 1-pip
  if ( is.null(mod_names) ) {
    mod_names <- c("MA",paste0("M", 1:k))
  }

  if ( is.null(plot_wind) ) {
    plot_wind <- c(1,p)
  }

  if ( is.null(var_names))
    var_names <- names(x_coef)

  par(mfrow = plot_wind, las = 1)
  for(i in 1:p){

    L68 <- x_coef[,i] - x_se[,i]
    U68 <- x_coef[,i] + x_se[,i]
    L95 <- x_coef[,i] - 2*x_se[,i]
    U95 <- x_coef[,i] + 2*x_se[,i]
    L99 <- x_coef[,i] - 3*x_se[,i]
    U99 <- x_coef[,i] + 3*x_se[,i]
    est <- x_coef[,i]

    col99 <- c("#e31a1c", rep("black", k))
    col95 <- c("#1f78b4", rep("#66c2a5", k))
    col68 <- c("#b2df8a", rep("#fc8d62", k))
    colest <- c("#e31a1c", rep("black", k))


    plot(c(1,1), xlim = c(min(L99), max(U99)), ylim = c(1, (k+1)),
         type = "n", xlab = "", yaxt = "n", ylab = "",
         main = var_names[i])

    axis(2, at = 1, labels = mod_names[1], tick = T, col.axis = col99[1])
    axis(2, at = 2:(k+1), labels = mod_names[2:(k+1)], tick = T,
         col.axis = col99[2])

    pmp_disp <- ifelse(pmp == 0, "<0.0000", as.character(pmp))
    text(min(L99)+ sd(L99)/5, 2:(k+1), pmp_disp, pos = 3)
    middle <- median(c(min(L99), max(U99)))
    text(middle, 1.1, paste0("Pr(beta_", i, ".MA = 0|y)= ",pep_text[i]),
         col = "#e31a1c", pos = 3)

    L68[L68 == 0] <- NA
    U68[U68 == 0] <- NA
    L95[L95 == 0] <- NA
    U95[U95 == 0] <- NA
    L99[L99 == 0] <- NA
    U99[U99 == 0] <- NA
    est[est == 0] <- NA
    segments(c(est[1], 0), c(1.5,0), c(est[1], 0), c(k+1.5, k+1.5),lty = c(4,4),
             col = c("#e31a1c", "gray"))
    abline(h = 1.5, col = "#e31a1c", lty = 2)
    abline(h = 2.5:(k+.5), col = "lightgray", lty = 2)
    segments(L99, 1:(k+1), U99, 1:(k+1), col = col99, lty = 1, lwd = 4)
    segments(L95, 1:(k+1), U95, 1:(k+1), col = col95, lty = 1, lwd = 6)
    segments(L68, 1:(k+1), U68, 1:(k+1), col = col68, lty = 1, lwd = 8)
    points(est, 1:(k+1), col = colest, pch = "|", cex = 1.5)


  }

}
#' @title Model Averaged Posteriors Plot for AIC weights
#' @description Generates results from information theoretic approximations to
#'   posterior model probabilities. Creates a graphical display showing all
#'   relevant output from the model averaging procedure and facilitates comparisons
#'   between individual models and their model averaged counterpart by displaying
#'   confidence intervals from all models.
#' @param x A data frame with the response variable in the first column 
#'   and x-variables in the remaining columns. 
#' @param mod_names A vector of length \eqn{J} specifying the names of the models.
#'   If left blank, models will be named \eqn{M_1, M_2,..., M_J}.
#' @param var_names A vector of length \eqn{p} specifying the names of the
#'   varibles, if left blank the names will be determined by the column names
#'   of \code{x}.
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns for specifying the partition of the plotting window.
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#' @param mod_prior A vector of length \eqn{K} specifying the prior model
#' probabilities for each regression model in \eqn{\mathcal{M}}. Default is 
#' uniform.
#' @param family A character defining which link function to use in the 
#'   linear regression models. Default is \code{"gaussian"}.
#' @return The \code{MApp} plot and a summary of the MA results for both AIC 
#' and BIC weights.

MApp_AIC <- function(x, var_names = NULL,
                    mod_names = NULL, plot_wind = NULL, w_plus = FALSE, 
                    family= "gaussian", mod_prior = NULL){
  # assuming x is a matrix of estimates and standard errors
  # from each individual model (row) and the model averaged
  # result.
  Yvec <- x[,1]
  Xmat <- x[,-1]
  inmat <- t(BMS::bms(x, mprior = "uniform", g = "UIP", user.int = F)$topmod$betas())
  inmat[inmat != 0] <- 1
  
  p <- ncol(inmat)
  k <- nrow(inmat)
  
  temp <- approx_pmp(inmat, Xmat, Yvec, mod_names = mod_names, 
                     mod_prior = mod_prior, family = family, w_plus = w_plus)
  pmp <- as.numeric(as.character(temp[-c(1,2),2]))
  pip <- t(inmat)%*%pmp
  pep <- 1-pip
  pep_text <- ifelse(pep == 0, "<0.0000", as.character(pep))
  
  x_coef <- temp[-2, c(7:(6+p))]
  x_se <- temp[-2,c((8+p):ncol(temp))]
  
  if ( is.null(mod_names) ) {
    mod_names <- c("MA",paste0("M", 1:k))
  }
  
  if ( is.null(plot_wind) ) {
    plot_wind <- c(1,p)
  }
  
  if ( is.null(var_names))
    var_names <- names(Xmat)
  
  par(mfrow = plot_wind, las = 1)
  for(i in 1:p){
    
    L68 <- x_coef[,i] - x_se[,i]
    U68 <- x_coef[,i] + x_se[,i]
    L95 <- x_coef[,i] - 2*x_se[,i]
    U95 <- x_coef[,i] + 2*x_se[,i]
    L99 <- x_coef[,i] - 3*x_se[,i]
    U99 <- x_coef[,i] + 3*x_se[,i]
    est <- x_coef[,i]
    
    col99 <- c("#e31a1c", rep("black", k))
    col95 <- c("#1f78b4", rep("#66c2a5", k))
    col68 <- c("#b2df8a", rep("#fc8d62", k))
    colest <- c("#e31a1c", rep("black", k))
    
    
    plot(c(1,1), xlim = c(min(L99), max(U99)), ylim = c(1, (k+1)),
         type = "n", xlab = "", yaxt = "n", ylab = "",
         main = var_names[i])
    
    axis(2, at = 1, labels = mod_names[1], tick = T, col.axis = col99[1])
    axis(2, at = 2:(k+1), labels = mod_names[2:(k+1)], tick = T,
         col.axis = col99[2])
    
    pmp_disp <- ifelse(pmp == 0, "<0.0000", as.character(pmp))
    text(min(L99)+ sd(L99)/5, 2:(k+1), pmp_disp, pos = 3)
    middle <- median(c(min(L99), max(U99)))
    text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~.(pep_text[i])),
         col = "#e31a1c", pos = 3)
    
    L68[L68 == 0] <- NA
    U68[U68 == 0] <- NA
    L95[L95 == 0] <- NA
    U95[U95 == 0] <- NA
    L99[L99 == 0] <- NA
    U99[U99 == 0] <- NA
    est[est == 0] <- NA
    segments(c(est[1], 0), c(1.5,0), c(est[1], 0), c(k+1.5, k+1.5),lty = c(4,4),
             col = c("#e31a1c", "gray"))
    abline(h = 1.5, col = "#e31a1c", lty = 2)
    abline(h = 2.5:(k+.5), col = "lightgray", lty = 2)
    segments(L99, 1:(k+1), U99, 1:(k+1), col = col99, lty = 1, lwd = 4)
    segments(L95, 1:(k+1), U95, 1:(k+1), col = col95, lty = 1, lwd = 6)
    segments(L68, 1:(k+1), U68, 1:(k+1), col = col68, lty = 1, lwd = 8)
    points(est, 1:(k+1), col = colest, pch = "|", cex = 1.5)
    
    
  }
  return(temp)
}


