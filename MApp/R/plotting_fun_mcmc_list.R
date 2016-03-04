#' @title Model Averaged Posteriors Plot
#' @description Works with lists of posterior draws of partial regression
#'   coefficients from multiple models. Creates a graphical display showing all
#'   relevant output from the model averaging procedure. Used to graphically
#'   compare posterior distributions of partial regression coefficients from
#'   individual models in the model set to the posterior distributions resulting
#'   from the model averaged partial regression coefficients. Also used to
#'   compare posterior standard deviations of the continuous piece of the model
#'   averaged distribution to corresponding standard deviations from individual
#'   posterior distributions in tabular form.
#' @param mcmc.list A list of posterior samples of partial regression
#'   coefficients of interest for each model in the model set.
#' @param g A numeric value specifying the value of the hyperparameter \eqn{g}.
#' @param weights A vector of posterior model weights for all models in the
#'   model set.
#' @param PIP A vector of posterior inclusion probabilities for all potential
#'   covariates.
#' @param plot.wind A vector of length 2 specifying the number of rows and
#'   columns to partition the plotting window into.
#' @param max.display An integer specifying the number of models to display in
#'   the plotting window.
#' @param mod.names A vector specifying the names of the models. If left blank,
#'   models will be named \eqn{M_1, M_2,..., M_K}, where \eqn{M_1} will have the
#'   largest posterior model probability and \eqn{M_K} has the smallest
#'   posterior model probability.
#' @param include.coef A vector specifying which partial regression coefficients
#'   to display.
#' @return The \code{MApp} plot and a table of posterior standard deviations for
#'   the top \code{max.display} individual models specified, along with the
#'   posterior standard deviation of the model averaged parameter.
MApp <- function(mcmc.list, g, weights, PIP, plot.wind, max.display = NULL,
                   mod.names = NULL, include.coef = NULL, ...) {

    K <- length(mcmc.list)
    p <- dim(mcmc.list[[1]])[2]  # Number of parameters
    if (is.null(mod.names))
        mod.names <- paste("M", seq(1:K), sep = "")
    names(mcmc.list) <- paste("M", seq(1:K), sep = "")  # Create names for models
    num.draws <- dim(mcmc.list[[1]])[1]  # draws from each posterior beta vector
    var.names <- names(mcmc.list[[1]])

    # Create a data frame for each coefficient estimate and track which model it
    # came from

    coef.frames <- list(1:p)

    for (i in 1:p) {
        coef.frames[[i]] <- vec.fun2(mcmc.list, i)
        coef.frames[[i]] <- data.frame(coef.frames[[i]], rep(mod.names, each = num.draws))
        names(coef.frames[[i]]) <- c("Post.Vec", "Model")
        coef.frames[[i]]$Model <- factor(coef.frames[[i]]$Model, levels = mod.names)
    }

    for (i in 1:p) {
        coef.frames[[i]] <- MApost.fun(coef.frames[[i]], K, weights, num.draws)
    }

    # Only plot up to the max.disply model (after creating MA distributions)

    if (is.null(max.display)) {
        for (i in 1:p) {
            coef.frames[[i]] <- coef.frames[[i]]
        }
        max.display <- K
    } else {
        for (i in 1:p) {
            coef.frames[[i]]$mod <- as.numeric(coef.frames[[i]]$Model)
            coef.frames[[i]] <- subset(coef.frames[[i]], mod <= (max.display + 1))
            coef.frames[[i]]$Model <- factor(coef.frames[[i]]$Model)
            coef.frames[[i]] <- coef.frames[[i]][, 1:2]
        }
    }

    # Find MaxMin of each beanplot

    MaxMin <- matrix(NA, nrow = p, ncol = 2)

    for (i in 1:p) {
        MaxMin[i, 1] <- min(coef.frames[[i]][, 1], na.rm = T)
        MaxMin[i, 2] <- max(coef.frames[[i]][, 1], na.rm = T)
    }

    # Only make plots for coefficients that appeared

    include.beans <- c(1:p)
    num.null <- sum((MaxMin[, 1] & MaxMin[, 2] == 0))
    if (num.null == 0) {
        include.beans <- include.beans
    } else {
        include.beans <- include.beans[-which((MaxMin[, 1] & MaxMin[, 2] == 0))]
    }


    if (length(include.beans) < p) {
        SD <- matrix(NA, nrow = p, ncol = (max.display + 1))
        for (i in include.beans) {
            for (j in 1:(max.display + 1)) {
                SD[i, j] <- round(sd(
                					subset(coef.frames[[i]],
                					Model == levels(coef.frames[[i]]$Model)[j])[, 1]),
                					3)
            }
        }
        SD[which(is.na(SD))] <- 0
        test.in <- apply(SD, 1, sum)
        SD <- SD[-which(test.in == 0), ]
        SD[which(SD == 0)] <- "---"
        SD <- as.table(SD)
        dimnames(SD)[[1]] <- var.names[include.beans]
        dimnames(SD)[[2]] <- levels(coef.frames[[i]]$Model)
    } else {

        SD <- matrix(NA, nrow = p, ncol = (max.display + 1))
        for (i in include.beans) {
            for (j in 1:(max.display + 1)) {
                SD[i, j] <- round(sd(
                					subset(coef.frames[[i]],
                					Model == levels(coef.frames[[i]]$Model)[j])[, 1]),
                					3)
            }
        }
        SD[which(is.na(SD))] <- "---"
        SD <- as.table(SD)
        dimnames(SD)[[1]] <- var.names[include.beans]
        dimnames(SD)[[2]] <- levels(coef.frames[[i]]$Model)
    }


    par(mfrow = c(plot.wind), las = 1)

    # Plot the beans
    if (is.null(include.coef)) {
        include.coef <- include.beans
    }

    for (i in include.coef) {

        # Initilize plot
        allplot <- beanplot::beanplot(Post.Vec ~ Model, names = levels(coef.frames[[i]]$Model),
        				      data = coef.frames[[i]],
        					  what = c(0, 0, 0, 0),
        					  border = 1,
        					  main = paste("Coefficient for ",
        								     var.names[i],
        									 sep = ""),
        									 horizontal = T)

        # Add top model beans
        beanplot::beanplot(Post.Vec ~ Model, side = "second", show.names = F,
        		  data = coef.frames[[i]],
        		  add = T,
        		  subset = Model != "MA",
        		  what = c(0, 1, 1, 0),
        		  border = 1,
        		  ll = 0.12,
                  horizontal = T,
            	  col = c("lightgray"))

        # Add MA posterior mean line
        MA.mean <- mean(subset(coef.frames[[i]], Model == "MA")$Post.Vec, na.rm = T)
        segments(MA.mean, 1, MA.mean, (max.display + 1.35), col = "black", lty = 4)

        # Add MA bean
        beanplot::beanplot(Post.Vec ~ Model, side = "second", show.names = F,
        		   data = coef.frames[[i]],
        		   add = T,
        		   subset = Model == "MA",
        		   what = c(0, 1, 1, 1),
        	       border = 1,
        		   ll = 0.12,
            	   horizontal = T,
                   col = c("black", "gray", "red", "darkgrey"))

        # Add horizontal lines
        modLines <- seq(1:(max.display + 1))
        abline(h = c(1, modLines), lty = 2, col = "gray")

        # Add model weights. display 0.000 for ~0's to convey rounding
        disp.weights <- round(weights[1:max.display], 3)
        disp.weights <- ifelse(disp.weights == 0, "< 0.000", as.character(disp.weights))

        # Add text to the plot
        text(rep(MaxMin[i, 1], max.display + 1),
        				c(1:(max.display + 1) + 0.25),
        				c("Avg", disp.weights),
        				col = "red",
        				cex = 1,
        				pos = 1)

        # Add and make sure PIP's are not rounded to exactly 0.
        disp.PEP <- round(1 - PIP, 3)
        disp.PEP <- ifelse(disp.PEP == 0, "<0.000", as.character(disp.PEP))
        vadjust <- ifelse(max.display <= 10, 0.35, 0)
        text(mean(MaxMin[i, ]),
        		    vadjust,
        			  paste("Pr(beta_",
        			    paste(include.beans[i],
        				  paste(".MA = 0",
            				paste("|y)",
            				  disp.PEP[i],
            				  sep = "="))),
            			      sep = ""),
            				  col = "black",
            				  cex = 1,
            				  pos = 3)
        text(MaxMin[i, 2], (max.display + 1.5),
        	 paste("Sum PMP", round(sum(weights), 3),
        	   sep = "="),
        	   col = "black",
        	   cex = 1,
        	   pos = 2)
    }

    return(SD)
}
