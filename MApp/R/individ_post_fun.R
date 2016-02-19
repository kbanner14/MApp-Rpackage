#' @title Compute analytical posterior variance covariance matrix of partial
#'   regression coefficients
#' @description Compute posterior variance covariance matrix for an individual
#'   regression model after model averaging. Posterior variance-covariance
#'   matrices are generated assuming model averaging with Zellner's
#'   \eqn{g}-prior [add citation here]. If no value of \eqn{g} is specified, the
#'   unit information prior (g = N) is assumed. placed on the regression
#'   parameters.
#'
#' @param lmObject An object of class \code{lm}
#' @return The posterior variance for the individual \code{lm} object.
cov_fun <- function(lmObject, g = NULL, ...) {

  # function to specify posterior vcov structure under g-prior
  # Input a lm object and g. Output a vcov matrix.
  # if no g is specified default is UIP (g = N)

  R2 <- summary(lmObject)$r.squared
  XtXinv <- summary(lmObject)$cov.unscaled
  N <- dim(lmObject$model)[1]
  if(is.null(g))
    g <- N
  shrink <- g/(1+g)
  y <- as.matrix(lmObject$model[, 1])
  ybar <- mean(y)
  ymybar <- y - ybar
  constant <- as.numeric((t(ymybar) %*% (ymybar)/(N - 3)) * shrink * (1 - (shrink * R2)))
  cov <- constant * XtXinv
  return(cov[-1, -1])
}



#' @title Convert a list of posterior samples into a long data frame.
#' @description This function will create a long data frame to be used as an
#'   input in \code{\link{MApost.fun}}.
vec_fun <- function(post.list, coef) {
  # Function to convert from wide to long format by model for each coef

  num.vars <- dim(post.list[[1]])[2]
  num.mods <- length(post.list)
  num.draws <- dim(post.list[[1]])[1]
  len <- num.draws * num.mods
  vector <- rep(NA, len)
  for (i in c(0, seq(1:(num.mods - 1)))) {
    vector[(i * (num.draws) + 1):((i + 1) * num.draws)] <- post.list[[i + 1]][, coef]
  }
  return(vector)
}

#' @title Compute analytical posterior means
#'
#' @description Compute posterior means for an individual regression model after
#'   model averaging. Posterior means are generated assuming model averaging
#'   with Zellner's \eqn{g}-prior [add citation here]. If no value of \eqn{g} is
#'   specified, the unit information prior (g = N) is assumed. placed on the
#'   regression parameters.
#'
#' @param lmObject A lm object
#' @return The posterior mean for the individual\code{lm} object.
mean_fun <- function(lmObject, g = NULL,...){
  # function to compute posterior mean under g-prior
  # Input a lm object and g. Output a mean vector
  # if no g is specified default is UIP (g = N)
  N <- dim(lmObject$model)[1]
  if(is.null(g))
    g <- N
  shrink <- g/(g + 1)
  betas <- shrink * coef(lmObject)
  return(betas[-1])
}

#' @title Simulate draws from posterior distributions of partial regression
#'   coefficients.
#' @family MApp functions
#' @description Simulate draws from individual posterior distributions for all
#'   models in the model set. Posteriors are simulated assuming a
#'    g-prior was used on the coefficints.
#' @param input_mat A matrix defining the model set.
#' @param Xmat The numeric p by n X matrix for the regression using all p
#'   covariates.
#' @param Yvec The numeric vector of responses.
#' @param num_sims The number of draws to simulate.
#' @return A list of simulated posterior distributions for all partial
#'   regression coefficients for each model considered.
bms_post_sim <- function(input_mat, Xmat, Yvec, num_sims, g = NULL) {

  # Function to simulate draws from MVt for individual model
  # posteriors, cov_fun and mean_fun must be loaded.

  all_dat <- data.frame(Yvec, Xmat)
  num_x <- dim(Xmat)[2]
  num_models <- dim(input_mat)[1]
  allMods <- as.list(1:(2 * num_models))
  posts <- as.list(1:num_models)
  names_posts <- rep(NA, num_models)
  N <- dim(all_dat)[1]
  for (j in 1:num_models) {
    names_posts[j] <- paste("M", j, sep = "_")
  }

  names(posts) <- names_posts

  # see if the full null model is included
  # if not obtain posteriors according to g-prior for all models
  if(sum(apply(input_mat,1, sum) == 0) == 0){

    idx <- which(apply(input_mat, 1, sum) != 0)
    # draws for coefficients from all other models are mvt's
    # obtain posterior moments
    for (i in idx){
      xs1 <- colnames(Xmat)[which(input_mat[i, ] == 1)]
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

    # simulate posterior draws according to moments
    for (i in idx) {
      posts[[i]] <- data.frame(matrix(rep(input_mat[i, ], num_sims),
                                      nrow = num_sims, 
                                      ncol = num_x, 
                                      byrow = T))
      names(posts[[i]]) <- colnames(Xmat)
      posts[[i]][posts[[i]] == 0] <- NA
      num_vars <- length(which(input_mat[i, ] == 1))
      var_idx <- which(input_mat[i, ] == 1)
      posts[[i]][, var_idx] <- LearnBayes::rmt(num_sims, df = N - num_vars - 1,
                                                S = allMods[[i + num_models]],
                                                mean = allMods[[i]])  
         }

  } else {

    idx <- which(apply(input_mat,1, sum) !=0)
    idx_null <- which(apply(input_mat, 1, sum) == 0)

    # obtain posterior moments for each of the models

    allMods[[idx_null + (num_models)]] <- 0
    allMods[[idx_null]] <- 0

    # Draws from coefficients from null model are 0

    posts[[idx_null]] <- data.frame(matrix(rep(input_mat[idx_null, ], num_sims),
                                           nrow = num_sims, 
                                           ncol = num_x, 
                                           byrow = T))
    posts[[idx_null]][posts[[idx_null]] == 0] <- NA

    # draws for coefficients from all other models are mvt's
    # obtain posterior moments
    for (i in idx){
      xs1 <- colnames(Xmat)[which(input_mat[i, ] == 1)]
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

    # simulate posterior draws according to moments
    for (i in idx) {
      posts[[i]] <- data.frame(matrix(rep(input_mat[i, ], num_sims),
                                      nrow = num_sims, 
                                      ncol = num_x, 
                                      byrow = T))
      names(posts[[i]]) <- colnames(Xmat)
      posts[[i]][posts[[i]] == 0] <- NA
      num_vars <- length(which(input_mat[i, ] == 1))
      var_idx <- which(input_mat[i, ] == 1)
      posts[[i]][, var_idx] <- LearnBayes::rmt(num_sims, df = N - num_vars - 1,
                                                S = allMods[[i + num_models]],
                                                mean = allMods[[i]])
    }

  }
  
  return(posts)
}

#' @title Generate individual and MA posterior distirbuitons 
#' @family MApp functions
#' @description Simulate draws from individual posterior distributions for all
#'   models in the model set. Posteriors are simulated assuming a
#'    g-prior was used on the coefficients. Use individual posterior draws and
#' 	 posterior model probabilities to sample from the MA posterior distribuiton.
#' @param x a \code{bma} object generated by \code{bms()}
#' @param num_draws The size of the model averaged posterior distribution.
#' @param mod_names a vector of the model names
#' @return \code{list}. The first element is the model averaged posterior 
#' 		distribuitons for the partial regression coefficients 
#'		associated with each explanatory variable and the number 
#' 		of samples taken from each model \eqn{\mathcal{M}}. 
MApost_bms_sim <- function(x, num_draws = 5000, 
                           mod_names = NULL, ...) {
  # function to genereate individual and MA posterior distributions 
  # for coefficients 
  
  # pull off informaiton from bms
  results <- data.frame(coef(x))
  Yvec <- x$X.data[,1]
  Xmat <- x$X.data[,-1]
  all_dat <- data.frame(Yvec, Xmat)
  N <- dim(all_dat)[1]
  
  # Create input matrix to work in simulation function. Rows ordered w.r.t
  # posterior model probability
  post_means <- x$topmod$betas()
  inmat <- post_means
  inmat[inmat != 0] <- 1
  inmat <- t(inmat)
  weights <- sort(pmp.bma(x)[,1], decreasing = T)
  g <- x$gprior.info$g
  K <- dim(inmat)[1]
  if (is.null(mod_names))
    mod_names <- paste("M", seq(1:K), sep = "")
     
  var_names <- colnames(Xmat)
  
  # create MA posterior for all coefficients simultaneously
  # generate samps per model using multinomial distribution
  
  samps <- rmultinom(num_draws, 1, weights)
  samps_mod <- apply(samps, 1, sum)
  names(samps_mod) <- mod_names
  
  index_samp <-which(samps_mod != 0 ) 
  num_models <- length(index_samp)
  allMods <- as.list(1:(2 * num_models))
  
  # check if null model is one of the sampled models
  inmat <- inmat[index_samp,]
  if(sum(apply(inmat,1, sum) == 0) == 0){
  	idx <- which(apply(inmat, 1, sum) != 0)
    # draws for coefficients from all other models are mvt's
    # obtain posterior moments
    for (i in idx){
      xs1 <- colnames(Xmat)[which(inmat[i, ] == 1)]
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

    MA_posts <- data.frame()    
    for(m in index_samp){ 
    	draws <- data.frame(matrix(rep(inmat[m, ], samps_mod[m]),
                                      nrow = samps_mod[m], 
                                      ncol = ncol(inmat), 
                                      byrow = T))
        num_vars <- length(which(inmat[m, ] == 1))
        size_vcov <- length(allMods[[(m + num_models)]])
        var_idx <- which(inmat[m, ] == 1)
        draws[,var_idx] <- LearnBayes::rmt(samps_mod[m], df = N - num_vars - 1,
        										S = allMods[[(m + num_models)]],
        										mean = allMods[[m]])
        MA_posts <- rbind(MA_posts, draws)
          }
    
    names(MA_posts) <- var_names
    

  } else {
	
    idx <- which(apply(inmat,1, sum) !=0)
    idx_null <- which(apply(inmat, 1, sum) == 0)
    
    # obtain posterior moments for each of the models
    
    allMods[[idx_null + (num_models)]] <- 0
    allMods[[idx_null]] <- 0

    # draws for coefficients from all other models are mvt's
    # obtain posterior moments
    for (i in idx){
      xs1 <- colnames(Xmat)[which(inmat[i, ] == 1)]
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

   MA_posts <- data.frame(matrix(rep(0, samps_mod[idx_null]*length(var_names)),
   									  nrow = samps_mod[idx_null], 
   									  ncol = length(var_names)))    
    for(m in index_samp[-idx_null]){ 
    	draws <- data.frame(matrix(rep(inmat[m, ], samps_mod[m]),
                                      nrow = samps_mod[m], 
                                      ncol = ncol(inmat), 
                                      byrow = T))
        num_vars <- length(which(inmat[m, ] == 1))
        size_vcov <- length(allMods[[(m + num_models)]])
        var_idx <- which(inmat[m, ] == 1)
        draws[,var_idx] <- LearnBayes::rmt(samps_mod[m], df = N - num_vars - 1,
        										S = allMods[[(m + num_models)]],
        										mean = allMods[[m]])
        MA_posts <- rbind(MA_posts, draws)
          }
    
    names(MA_posts) <- var_names
 
  }

  return(list(MA_posts, samps_mod))
}


