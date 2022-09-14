#' @title Model Averaged Posteriors Plot for MA Using Information Criteria
#' @description Generates results from approximations to
#'   posterior model probabilities (PMPs) based on information criteria. 
#'   Creates a graphical display showing all relevant output from the model 
#'   averaging procedure and facilitates comparisons between individual models 
#'   and their model averaged counterpart by displaying
#'   confidence intervals from all models. Also displays the results in the form 
#'   of a data frame. \emph{This function works best with modest numbers of 
#'   explanatory variables because it uses} \code{bms} \emph{to define an all-subsets
#'   model set. 
#'   By default,} \code{bms} \emph{will save results from 500 models. When 
#'   there are 9 explanatory variables, there 512 models in the all-subsets model 
#'   set. When there are 9 or more explanatory variables, it is 
#'   important to check to make sure the sum of the PMPs for those 500 models 
#'   is large enough for them to be a reasonable subset to condition IC 
#'   approximate model averaging on. For more control with
#'   larger model sets, and use with other packages IC-based results use} 
#'   \code{MApp_IC_large()}.
#' @param x A data frame with the response variable in the first column 
#'   and x-variables in the remaining columns. 
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns for specifying the partition of the plotting window.
#' @param type An argument allowing the user to specify the type of information
#'   criterion to use in approximating the PMPs. Default is \code{BIC}. Will 
#'   also accommodate AIC weights if the user specifies \code{type = "AIC"}. 
#' @param aic_c A logical statement that goes with \code{type = "AIC"} to specify
#'   whether or not to use adjusted AIC (AICc) weights or raw AIC weights. The 
#'   default is \code{aic_c = TRUE}.
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#' @param var_names A vector of length \eqn{p} specifying the names of the
#'   variables, if left blank the names will be determined by the column names
#'   of \code{x}.
#' @param mod_names A vector of length \eqn{J} specifying the names of the models.
#'   If left blank, models will be named \eqn{M_1, M_2,..., M_J}.
#' @param mod_prior A vector of length \eqn{J} specifying the prior model
#'   probabilities for each regression model in \eqn{\mathcal{M}}. Default is 
#'   discrete uniform.
#' @param family A character defining which link function to use in the 
#'   linear regression models. Default is \code{"gaussian"}.
#' @return The \code{MApp} plot and a summary of the MA results for both types of 
#'  weights: AIC or AICc, and BIC.

MApp_IC <- function(x, plot_wind, type = "BIC", aic_c = TRUE, 
                    w_plus = FALSE, var_names = NULL, 
                    mod_names = NULL, family= "gaussian", 
                    mod_prior = NULL, max_display = NULL){
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
                     mod_prior = mod_prior, family = family, 
                     w_plus = w_plus, aic_c = aic_c)
  
  if(type == "AIC"){
    pmp <- as.numeric(as.character(temp[-c(1,2),2]))
    pip <- t(inmat)%*%pmp
    pep <- 1-pip
    pep <- round(pep, 3)
    pep_text <- ifelse(pep == 0, expression(" "%~~%0), as.character(pep))
    # sort by pmp 
    x_coef <- temp[-2, c(7:(6+p))]
    x_mods <- cbind(pmp, x_coef[2:nrow(x_coef),])
    x_mods <- x_mods[order(pmp, decreasing = T), 2:ncol(x_mods)]
    x_coef[2:nrow(x_coef),] <- x_mods
    x_se <- temp[-2,c((8+p):ncol(temp))]
    # sort by pmp
    x_mse <- cbind(pmp, x_se[2:nrow(x_se),])
    x_mse <- x_mse[order(pmp, decreasing = T), 2:ncol(x_mse)]
    x_se[2:nrow(x_se),] <- x_mse
    if(is.null(max_display) & nrow(x_coef) > 20){
      x_coef <- x_coef[1:21, ]
      x_se <- x_se[1:21, ]
      max_display <- 20
    } else {
      if(is.null(max_display)){
        max_display <- k
      } else {
        max_display <- max_display
        x_coef <- x_coef[1:(max_display+1), ]
        x_se <- x_se[1:(max_display+1), ]
      }
    }
  } else {
    pmp <- as.numeric(as.character(temp[-c(1,2),3]))
    pip <- t(inmat)%*%pmp
    pep <- 1-pip
    pep <- round(pep, 3)
    pep_text <- ifelse(pep == 0, expression(" "%~~%0), as.character(pep))
    
    x_coef <- temp[-1, c(7:(6+p))]
    # sort by pmp 
    x_mods <- cbind(pmp, x_coef[2:nrow(x_coef),])
    x_mods <- x_mods[order(pmp, decreasing = T), 2:ncol(x_mods)]
    x_coef[2:nrow(x_coef),] <- x_mods
    x_se <- temp[-1,c((8+p):ncol(temp))]
    # sort by pmp
    x_mse <- cbind(pmp, x_se[2:nrow(x_se),])
    x_mse <- x_mse[order(pmp, decreasing = T), 2:ncol(x_mse)]
    x_se[2:nrow(x_se),] <- x_mse
    if(is.null(max_display) & nrow(x_coef) > 20){
      x_coef <- x_coef[1:21, ]
      x_se <- x_se[1:21, ]
      max_display <- 20
    } else {
      if(is.null(max_display)){
        max_display <- k
      } else {
        max_display <- max_display
        x_coef <- x_coef[1:(max_display+1), ]
        x_se <- x_se[1:(max_display+1), ]
      }
    }
  }
  
  if ( is.null(mod_names) ) {
    mod_names <- c("MA",paste0("M", 1:k))
  } else {
    mod_names <- c("MA", mod_names)
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
    
    col99 <- c("#e31a1c", rep("black", max_display))
    col95 <- c("#1f78b4", rep("#66c2a5", max_display))
    col68 <- c("#b2df8a", rep("#fc8d62", max_display))
    colest <- c("#e31a1c", rep("black", max_display))
    
    
    plot(c(1,1), xlim = c(min(L99), max(U99)), ylim = c(1, (max_display+1)),
         type = "n", xlab = "", yaxt = "n", ylab = "",
         main = var_names[i])
    
    axis(2, at = 1, labels = mod_names[1], tick = T, col.axis = col99[1])
    axis(2, at = 2:(max_display+1), labels = mod_names[2:(max_display+1)], tick = T,
         col.axis = col99[2])
    
    pmp_disp <- round(pmp[order(pmp, decreasing = T)], 3)
    pmp_disp <- ifelse(pmp_disp == 0,expression(" "%~~%0), as.character(pmp_disp))
    text(min(L99)+ sd(L99)/5, 2:(max_display+1), pmp_disp[1:max_display], pos = 3)
    middle <- median(c(min(L99), max(U99)))
    if (as.character(pep_text[i]) == "\" \" %~~% 0"){
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~""%~~%0),
           col = "#e31a1c", pos = 3)
    } else {
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~.(as.character(pep_text[i]))),
           col = "#e31a1c", pos = 3)
    }
    
    L68[L68 == 0] <- NA
    U68[U68 == 0] <- NA
    L95[L95 == 0] <- NA
    U95[U95 == 0] <- NA
    L99[L99 == 0] <- NA
    U99[U99 == 0] <- NA
    est[est == 0] <- NA
    segments(c(est[1], 0), c(1.5,0), c(est[1], 0), c(max_display+1.5, max_display+1.5),lty = c(4,4),
             col = c("#e31a1c", "gray"))
    abline(h = 1.5, col = "#e31a1c", lty = 2)
    abline(h = 2.5:(max_display+.5), col = "lightgray", lty = 2)
    segments(L99, 1:(max_display+1), U99, 1:(max_display+1), col = col99, lty = 1, lwd = 4)
    segments(L95, 1:(max_display+1), U95, 1:(max_display+1), col = col95, lty = 1, lwd = 6)
    segments(L68, 1:(max_display+1), U68, 1:(max_display+1), col = col68, lty = 1, lwd = 8)
    points(est, 1:(max_display + 1), col = colest, pch = "|", cex = 1.5)
  }
  disp_type <- ifelse(aic_c == TRUE & type == "AIC", "AICc", 
                      ifelse(type == "BIC", "BIC", "AIC"))
  message(paste("Type of information criteria used in MAP plot:", disp_type))
  if(type == "AIC"){
    pmp_order <- order(as.numeric(as.character(temp[-c(1,2),2])), decreasing = T)
    temp_sub <- temp[-c(1,2),]
    temp_ord <- temp_sub[pmp_order, ]
    temp[3:nrow(temp), ] <- temp_ord
    return(temp[1:(max_display+2), ])
  } else {
    pmp_order <- order(as.numeric(as.character(temp[-c(1,2),3])), decreasing = T)
    temp_sub <- temp[-c(1,2),]
    temp_ord <- temp_sub[pmp_order, ]
    temp[3:nrow(temp), ] <- temp_ord
    return(temp[1:(max_display+2), ])
  }
}

#' @title Model Averaged Posteriors Plot for MA Using Information Criteria with dark theme
#' @description Generates results from approximations to
#'   posterior model probabilities (PMPs) based on information criteria. 
#'   Creates a graphical display showing all relevant output from the model 
#'   averaging procedure and facilitates comparisons between individual models 
#'   and their model averaged counterpart by displaying
#'   confidence intervals from all models. Also displays the results in the form 
#'   of a data frame. \emph{This function works best with modest numbers of 
#'   explanatory variables because it uses} \code{bms} \emph{to define an all-subsets
#'   model set. 
#'   By default,} \code{bms} \emph{will save results from 500 models. When 
#'   there are 9 explanatory variables, there 512 models in the all-subsets model 
#'   set. When there are 9 or more explanatory variables, it is 
#'   important to check to make sure the sum of the PMPs for those 500 models 
#'   is large enough for them to be a reasonable subset to condition IC 
#'   approximate model averaging on. For more control with
#'   larger model sets, and use with other packages IC-based results use} 
#'   \code{MApp_IC_large()}.
#' @param x A data frame with the response variable in the first column 
#'   and x-variables in the remaining columns. 
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns for specifying the partition of the plotting window.
#' @param type An argument allowing the user to specify the type of information
#'   criterion to use in approximating the PMPs. Default is \code{BIC}. Will 
#'   also accommodate AIC weights if the user specifies \code{type = "AIC"}. 
#' @param aic_c A logical statement that goes with \code{type = "AIC"} to specify
#'   whether or not to use adjusted AIC (AICc) weights or raw AIC weights. The 
#'   default is \code{aic_c = TRUE}.
#' @param w_plus A logical statement defining which type of MA to use. 
#'   Default is \code{FALSE}, which allows all models defined in \code{inmat} to 
#'   be used in the averaging for each variable. If set to \code{TRUE}, MA
#'   estimates will be made using only the models for which the coefficient
#'   associated with each variable is not set to exactly zero. 
#' @param var_names A vector of length \eqn{p} specifying the names of the
#'   variables, if left blank the names will be determined by the column names
#'   of \code{x}.
#' @param mod_names A vector of length \eqn{J} specifying the names of the models.
#'   If left blank, models will be named \eqn{M_1, M_2,..., M_J}.
#' @param mod_prior A vector of length \eqn{J} specifying the prior model
#'   probabilities for each regression model in \eqn{\mathcal{M}}. Default is 
#'   discrete uniform.
#' @param family A character defining which link function to use in the 
#'   linear regression models. Default is \code{"gaussian"}.
#' @return The \code{MApp} plot and a summary of the MA results for both types of 
#'  weights: AIC or AICc, and BIC.

MApp_IC_dark <- function(x, plot_wind, type = "BIC", aic_c = TRUE, 
                    w_plus = FALSE, var_names = NULL, 
                    mod_names = NULL, family= "gaussian", 
                    mod_prior = NULL, max_display = NULL){
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
                     mod_prior = mod_prior, family = family, 
                     w_plus = w_plus, aic_c = aic_c)
  
  if(type == "AIC"){
    pmp <- as.numeric(as.character(temp[-c(1,2),2]))
    pip <- t(inmat)%*%pmp
    pep <- 1-pip
    pep <- round(pep, 3)
    pep_text <- ifelse(pep == 0, expression(" "%~~%0), as.character(pep))
    # sort by pmp 
    x_coef <- temp[-2, c(7:(6+p))]
    x_mods <- cbind(pmp, x_coef[2:nrow(x_coef),])
    x_mods <- x_mods[order(pmp, decreasing = T), 2:ncol(x_mods)]
    x_coef[2:nrow(x_coef),] <- x_mods
    x_se <- temp[-2,c((8+p):ncol(temp))]
    # sort by pmp
    x_mse <- cbind(pmp, x_se[2:nrow(x_se),])
    x_mse <- x_mse[order(pmp, decreasing = T), 2:ncol(x_mse)]
    x_se[2:nrow(x_se),] <- x_mse
    if(is.null(max_display) & nrow(x_coef) > 20){
      x_coef <- x_coef[1:21, ]
      x_se <- x_se[1:21, ]
      max_display <- 20
    } else {
      if(is.null(max_display)){
        max_display <- k
      } else {
        max_display <- max_display
        x_coef <- x_coef[1:(max_display+1), ]
        x_se <- x_se[1:(max_display+1), ]
      }
    }
  } else {
    pmp <- as.numeric(as.character(temp[-c(1,2),3]))
    pip <- t(inmat)%*%pmp
    pep <- 1-pip
    pep <- round(pep, 3)
    pep_text <- ifelse(pep == 0, expression(" "%~~%0), as.character(pep))
    
    x_coef <- temp[-1, c(7:(6+p))]
    # sort by pmp 
    x_mods <- cbind(pmp, x_coef[2:nrow(x_coef),])
    x_mods <- x_mods[order(pmp, decreasing = T),2:ncol(x_mods)]
    x_coef[2:nrow(x_coef),] <- x_mods
    x_se <- temp[-1,c((8+p):ncol(temp))]
    # sort by pmp
    x_mse <- cbind(pmp, x_se[2:nrow(x_se),])
    x_mse <- x_mse[order(pmp, decreasing = T),2:ncol(x_mse)]
    x_se[2:nrow(x_se),] <- x_mse
    if(is.null(max_display) & nrow(x_coef) > 20){
      x_coef <- x_coef[1:21, ]
      x_se <- x_se[1:21, ]
      max_display <- 20
    } else {
      if(is.null(max_display)){
        max_display <- k
      } else {
        max_display <- max_display
        x_coef <- x_coef[1:(max_display+1), ]
        x_se <- x_se[1:(max_display+1), ]
      }
    }
  }
  
  if ( is.null(mod_names) ) {
    mod_names <- c("MA",paste0("M", 1:k))
  } else {
    mod_names <- c("MA", mod_names)  
  }
  
  
  if ( is.null(plot_wind) ) {
    plot_wind <- c(1,p)
  }
  
  if ( is.null(var_names))
    var_names <- names(Xmat)
  
  par(mfrow = plot_wind, las = 1, col.axis = "#f0f0f0", 
      col.main = "#f0f0f0", col.lab = "#f0f0f0", fg = "#f0f0f0")
  for(i in 1:p){
    
    L68 <- x_coef[,i] - x_se[,i]
    U68 <- x_coef[,i] + x_se[,i]
    L95 <- x_coef[,i] - 2*x_se[,i]
    U95 <- x_coef[,i] + 2*x_se[,i]
    L99 <- x_coef[,i] - 3*x_se[,i]
    U99 <- x_coef[,i] + 3*x_se[,i]
    est <- x_coef[,i]
    
    col99 <- c("#e6f598", rep("#dfc27d", max_display))
    col95 <- c("#1f78b4", rep("#66c2a5", max_display))
    col68 <- c("#b2df8a", rep("#fc8d62", max_display))
    colest <- c("#e6f598", rep("#dfc27d", max_display))
    
    
    plot(c(1,1), xlim = c(min(L99), max(U99)), ylim = c(1, (max_display+1)),
         type = "n", xlab = "", yaxt = "n", ylab = "",
         main = var_names[i])
    
    axis(2, at = 1, labels = mod_names[1], tick = T, col.axis = col99[1])
    axis(2, at = 2:(max_display+1), labels = mod_names[2:(max_display+1)], tick = T,
         col.axis = col99[2])
    
    pmp_disp <- round(pmp[order(pmp, decreasing = T)], 3)
    pmp_disp <- ifelse(pmp_disp == 0,expression(" "%~~%0), as.character(pmp_disp))
    text(min(L99)+ sd(L99)/5, 2:(max_display+1), pmp_disp[1:max_display], pos = 3)
    middle <- median(c(min(L99), max(U99)))
    if (as.character(pep_text[i]) == "\" \" %~~% 0"){
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~""%~~%0),
           col = "#e6f598", pos = 3)
    } else {
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~.(as.character(pep_text[i]))),
           col = "#e6f598", pos = 3)
    }
    
    L68[L68 == 0] <- NA
    U68[U68 == 0] <- NA
    L95[L95 == 0] <- NA
    U95[U95 == 0] <- NA
    L99[L99 == 0] <- NA
    U99[U99 == 0] <- NA
    est[est == 0] <- NA
    segments(c(est[1], 0), c(1.5,0), c(est[1], 0), c(max_display+1.5, max_display+1.5),lty = c(4,4),
             col = c("#e6f598", "gray"))
    abline(h = 1.5, col = "#e6f598", lty = 2)
    abline(h = 2.5:(max_display+.5), col = "lightgray", lty = 2)
    segments(L99, 1:(max_display+1), U99, 1:(max_display+1), col = col99, lty = 1, lwd = 4)
    segments(L95, 1:(max_display+1), U95, 1:(max_display+1), col = col95, lty = 1, lwd = 6)
    segments(L68, 1:(max_display+1), U68, 1:(max_display+1), col = col68, lty = 1, lwd = 8)
    points(est, 1:(max_display + 1), col = colest, pch = "|", cex = 1.5)
  }
  disp_type <- ifelse(aic_c == TRUE & type == "AIC", "AICc", 
                      ifelse(type == "BIC", "BIC", "AIC"))
  message(paste("Type of information criteria used in MAP plot:", disp_type))
  if(type == "AIC"){
    pmp_order <- order(as.numeric(as.character(temp[-c(1,2),2])), decreasing = T)
    temp_sub <- temp[-c(1,2),]
    temp_ord <- temp_sub[pmp_order, ]
    temp[3:nrow(temp), ] <- temp_ord
    return(temp[1:(max_display+2), ])
  } else {
    pmp_order <- order(as.numeric(as.character(temp[-c(1,2),3])), decreasing = T)
    temp_sub <- temp[-c(1,2),]
    temp_ord <- temp_sub[pmp_order, ]
    temp[3:nrow(temp), ] <- temp_ord
    return(temp[1:(max_display+2), ])
  }
}

#' @description Works with results from information theoretic approximations to
#'   posterior model probabilities. Creates a graphical display showing all
#'   relevant output from the model averaging procedure and facilitates comparisons
#'   between individual models and their model averaged counterpart by displaying
#'   confidence intervals from all models
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

MApp_IC_gen <- function(x_coef, x_se, pmp, inmat, var_names = NULL,
                    mod_names = NULL, plot_wind = NULL){
  # assuming x is a matrix of estimates and standard errors
  # from each individual model (row) and the model averaged
  # result.
  
  k <- nrow(x_coef) -1
  p <- ncol(x_coef)
  pip <- t(inmat)%*%pmp
  pep <- round(1-pip, 4)
  pmp <- round(as.numeric(as.character(pmp)), 4)
  pep_text <- ifelse(pep == 0, expression(" "%~~%0), as.character(pep))
  
  if ( is.null(mod_names) ) {
    mod_names <- c("MA",paste0("M", 1:k))
  } else {
    mod_names <- c("MA", mod_names)
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
    
    pmp_disp <- ifelse(pmp == 0, expression(" "%~~%0), as.character(pmp))
    text(min(L99)+ sd(L99)/5, 2:(k+1), pmp_disp, pos = 3)
    middle <- median(c(min(L99), max(U99)))
    if (as.character(pep_text[i]) == "\" \" %~~% 0"){
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~""%~~%0),
           col = "#e31a1c", pos = 3)
    } else {
      text(middle, 1.1, bquote("Pr("~beta[.(i)~", MA"]~"= 0 | y ) ="~.(as.character(pep_text[i]))),
           col = "#e31a1c", pos = 3)
    }
    
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

