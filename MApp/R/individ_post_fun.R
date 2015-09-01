#' @title Compute analytical posterior variance covariance matrix of
#'   partial regression coefficients
#' @description Compute posterior variance covariance matrix
#' for an individual
#' regression model after model averaging.
#' Posterior means are generated assuming
#' model averaging with the \code{bms()} function,
#' with a unit information prior (g = N) placed
#' on the regression parameters.
#'
#' @param lmObject An object of class \code{lm}
#' @return The posterior variance for the individual \code{lm} object.
cov.fun <- function(lmObject) {
  R2 <- summary(lmObject)$r.squared
  X <- model.matrix(lmObject)
  N <- dim(lmObject$model)[1]
  y <- as.matrix(lmObject$model[, 1])
  shrink <- N/(N + 1)
  ybar <- mean(y)
  ymybar <- y - ybar
  constant <- as.numeric((t(ymybar) %*% (ymybar)/(N - 3)) * shrink * (1 - (shrink * R2)))
  XtXinv <- solve(t(X) %*% X)
  cov <- constant * XtXinv
  return(cov[-1, -1])
}

#' @title Compute analytical posterior means
#'
#' @description Compute posterior means for an individual
#' regression model after model averaging.
#' Posterior means are generated assuming
#' model averaging with the \code{bms()} function,
#' with a unit information prior (g = N) placed
#' on the regression parameters.
#'
#' @param lmObject A lm object
#' @return The posterior mean for the individual\code{lm} object.
mean.fun <- function(lmObject) {
    N <- dim(lmObject$model)[1]
    shrink <- N/(N + 1)
    betas <- shrink * coef(lmObject)
    return(betas[-1])
}

#' @title Analytical posterior simulation function for partial regression coefficients
#' @description Simulate draws from individual posterior distributions for
#' all models in the model set. Posteriors are simulated assuming
#' the unit information formulation  (g = N) of the g-prior was used in
#' the model aeraging.
#'
#' @param input.mat A matrix defining the model set.
#' @param Xmat The numeric p x n X matrix for the regression
#'  using all p covariates.
#' @param Yvec The numeric vector of responses.
#' @param num.sims The number of draws to simulate.
#' @return A list of simulated posterior distribuitons for all
#'   partial regression coefficients for each model considered.
sim.post.fun <- function(input.mat, Xmat, Yvec, num.sims) {
    # Function to simulate draws from MVt for individual model
    # posteriors, cov.fun and mean.fun must be loaded.
    data <- data.frame(Yvec, Xmat)
    num.x <- dim(Xmat)[2]
    num.models <- dim(input.mat)[1]
    allMods <- as.list(1:(2 * num.models))
    posts <- as.list(1:num.models)
    names.posts <- rep(NA, num.models)
    N <- dim(data)[1]
    for (j in 1:num.models) {
        names.posts[j] <- paste("M", j, sep = "_")
    }
    names(posts) <- names.posts
    for (i in 1:num.models) {
        xs1 <- names(Xmat[which(input.mat[i, ] == 1)])
        xs <- paste(xs1, collapse = "+")
        form <- paste("Yvec ~", xs, sep = "")
        fit <- lm(form, data = data)
        allMods[[i + (num.models)]] <- cov.fun(fit)
        allMods[[i]] <- mean.fun(fit)
    }

    for (i in 1:num.models) {
        posts[[i]] <- data.frame(matrix(rep(input.mat[i, ], num.sims), nrow = num.sims,
            ncol = num.x, byrow = T))
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
#' @description This function will create a long data frame to be used as an input
#'   in \code{MApost.fun}.
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
#' model averaged posterior distribution according to the unit information prior.
#' @param post.dataframe A data frame in long format with posterior draws of a
#'   parameter of interest from all models in the model set. \code{vec.fun2} will
#'   convert a list of posterior samples into a long data frame.
#' @param models An integer specifying the number of models in the model set
#' @param weight A vector of posterior model weights of length \code{models}
#' @param num.draws An integer specifying the number of draws from each posterior
#'   vector.
MApost.fun <- function(post.dataframe, models, weight, num.draws) {

    post.dataframe$Include <- (as.numeric(is.na(post.dataframe$Post.Vec)) - 1) * (-1)
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
