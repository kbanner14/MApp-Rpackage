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
      xs1 <- names(Xmat[which(input_mat[i, ] == 1)])
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

    # simulate posterior draws according to moments
    for (i in idx) {
      posts[[i]] <- data.frame(matrix(rep(input_mat[i, ], num_sims),
                                      nrow = num_sims, ncol = num_x, byrow = T))
      names(posts[[i]]) <- names(Xmat)
      posts[[i]][posts[[i]] == 0] <- NA
      if (length(allMods[[(i + num_models)]]) == 1) {
        posts[[i]][, which(input_mat[i, ] == 1)] <- LearnBayes::rmt(num_sims,
                                                                    df = N - length(which(input_mat[i, ] == 1)) - 1,
                                                                    S = allMods[[i + num_models]],
                                                                    mean = allMods[[i]])
      } else {
        posts[[i]][, which(input_mat[i, ] == 1)] <- LearnBayes::rmt(num_sims,
                                                                    df = N - length(which(input_mat[i, ] == 1)) - 1,
                                                                    S = allMods[[i + num_models]],
                                                                    mean = allMods[[i]])
      }
    }

  } else {

    idx <- which(apply(input_mat,1, sum) !=0)
    idx_null <- which(apply(input_mat, 1, sum) == 0)

    # obtain posterior moments for each of the models

    allMods[[idx_null + (num_models)]] <- 0
    allMods[[idx_null]] <- 0

    # Draws from coefficients from null model are 0

    posts[[idx_null]] <- data.frame(matrix(rep(input_mat[idx_null, ], num_sims),
                                           nrow = num_sims, ncol = num_x, byrow = T))
    posts[[idx_null]][posts[[idx_null]] == 0] <- NA

    # draws for coefficients from all other models are mvt's
    # obtain posterior moments
    for (i in idx){
      xs1 <- names(Xmat[which(input_mat[i, ] == 1)])
      xs <- paste(xs1, collapse = "+")
      form <- paste("Yvec ~ 1 +", xs, sep = "")
      fit <- lm(form, data = all_dat)
      allMods[[i + (num_models)]] <- cov_fun(fit, g)
      allMods[[i]] <- mean_fun(fit, g)
    }

    # simulate posterior draws according to moments
    for (i in idx) {
      posts[[i]] <- data.frame(matrix(rep(input_mat[i, ], num_sims),
                                      nrow = num_sims, ncol = num_x, byrow = T))
      names(posts[[i]]) <- names(Xmat)
      posts[[i]][posts[[i]] == 0] <- NA
      if (length(allMods[[(i + num_models)]]) == 1) {
        posts[[i]][, which(input_mat[i, ] == 1)] <- LearnBayes::rmt(num_sims,
                                                                    df = N - length(which(input_mat[i, ] == 1)) - 1,
                                                                    S = allMods[[i + num_models]],
                                                                    mean = allMods[[i]])
      } else {
        posts[[i]][, which(input_mat[i, ] == 1)] <- LearnBayes::rmt(num_sims,
                                                                    df = N - length(which(input_mat[i, ] == 1)) - 1,
                                                                    S = allMods[[i + num_models]],
                                                                    mean = allMods[[i]])
      }
    }

  }


  return(posts)
}


cov.fun <- function(lmObject, g = NULL, ...) {
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
mean.fun <- function(lmObject, g = NULL,...){
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
#'   models in the model set. Posteriors are simulated assuming the unit
#'   information formulation (g = N) of the g-prior was used in the model
#'   averaging.
#' @param input.mat A matrix defining the model set.
#' @param Xmat The numeric p by n X matrix for the regression using all p
#'   covariates.
#' @param Yvec The numeric vector of responses.
#' @param num.sims The number of draws to simulate.
#' @return A list of simulated posterior distributions for all partial
#'   regression coefficients for each model considered.
sim.post.fun <- function(input.mat, Xmat, Yvec, num.sims) {
    # Function to simulate draws from MVt for individual model
    # posteriors, cov.fun and mean.fun must be loaded.
    all.dat <- data.frame(Yvec, Xmat)
    num.x <- dim(Xmat)[2]
    num.models <- dim(input.mat)[1]
    allMods <- as.list(1:(2 * num.models))
    posts <- as.list(1:num.models)
    names.posts <- rep(NA, num.models)
    N <- dim(all.dat)[1]
    for (j in 1:num.models) {
        names.posts[j] <- paste("M", j, sep = "_")
    }
    names(posts) <- names.posts
    for (i in 1:num.models) {
        xs1 <- names(Xmat[which(input.mat[i, ] == 1)])
        xs <- paste(xs1, collapse = "+")
        form <- paste("Yvec ~", xs, sep = "")
        fit <- lm(form, data = all.dat)
        allMods[[i + (num.models)]] <- cov.fun(fit)
        allMods[[i]] <- mean.fun(fit)
    }

    for (i in 1:num.models) {
        posts[[i]] <- data.frame(matrix(rep(input.mat[i, ], num.sims),
                                        nrow = num.sims, ncol = num.x, byrow = T))
        names(posts[[i]]) <- names(Xmat)
        posts[[i]][posts[[i]] == 0] <- NA
        if (length(allMods[[(i + num.models)]]) == 1) {
            posts[[i]][, which(input.mat[i, ] == 1)] <- LearnBayes::rmt(num.sims,
            	df = N - length(which(input.mat[i, ] == 1)) - 1,
              S = allMods[[i + num.models]],
              mean = allMods[[i]])
        } else {
            posts[[i]][, which(input.mat[i, ] == 1)] <- LearnBayes::rmt(num.sims,
            	df = N - length(which(input.mat[i, ] == 1)) - 1,
              S = allMods[[i + num.models]],
              mean = allMods[[i]])
        }
    }
    return(posts)
}

#' @title Convert a list of posterior samples into a long data frame.
#' @description This function will create a long data frame to be used as an
#'   input in \code{\link{MApost.fun}}.
vec.fun2 <- function(post.list, coef) {
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

#' @title Create MA distribution using draws from individual models
#'
#' @description Sample from individual posteriors to generate an approximate
#'   model averaged posterior distribution according to the unit information
#'   prior.
#' @param post.dataframe A data frame in long format with posterior draws of a
#'   parameter of interest from all models in the model set.
#'   \code{\link{vec.fun2}} will convert a list of posterior samples into a long
#'   data frame.
#' @param models An integer specifying the number of models in the model set
#' @param weight A vector of posterior model weights of length \code{models}
#' @param num.draws An integer specifying the number of draws from each
#'   posterior vector.
MApost.fun <- function(post.dataframe, models, weight, num.draws) {

    post.dataframe$Include <- as.numeric(is.na(post.dataframe$Post.Vec) == FALSE)
    keep <- with(subset(post.dataframe, Include == 1), unique(Model))
    idx.in <- rep(NA, length(keep))
    for (i in 1:length(keep)) {
        idx.in[i] <- which(unique(post.dataframe$Model) == keep[i])
    }
    tot.w <- sum(weight[idx.in])
    weight[-idx.in] <- 0
    if (tot.w == 0) {
        weight <- rep(0, models)
    } else {
        weight <- weight / tot.w
    }
    samps.per.model <- round(weight * num.draws, 0)

    # Keep only those that have nonzero samps associated with them.

    idx.in <- idx.in[which(samps.per.model[idx.in] != 0)]
    post.samps <- sum(samps.per.model)
    if (post.samps == 0) {
        MApost <- rep(NA, num.draws)
    } else {
        MApost <- rep(NA, post.samps)
    }
    if (post.samps == 0) {
        post.dataframe[, 1] <- 0
        MApost <- data.frame(MApost, rep("MA", num.draws))
        names(MApost) <- c("Post.Vec", "Model")
        Allpost <- rbind(MApost, post.dataframe[, 1:2])
        return(Allpost)
    } else {
        end <- 0
        for (i in idx.in) {
            MApost[(end + 1):(sum(samps.per.model[1:i]))] <- sample(
            		post.dataframe[(((i - 1) * num.draws) + 1):(i * num.draws), 1],
            		samps.per.model[i],
            		rep = F)
            end <- sum(samps.per.model[1:i])
        }
        MApost <- data.frame(MApost, rep("MA", post.samps))
        names(MApost) <- c("Post.Vec", "Model")
        Allpost <- rbind(MApost, post.dataframe[, 1:2])
        return(Allpost)
    }
}
