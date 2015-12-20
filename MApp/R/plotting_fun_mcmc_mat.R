#' @title Model Averaged Posteriors Plot
#' @description Works with posterior draws of partial regression
#'   coefficients from multiple models. Creates a graphical display showing all
#'   relevant output from the model averaging procedure. Used to graphically
#'   compare posterior distributions of partial regression coefficients from
#'   individual models in the model set to the posterior distributions resulting
#'   from the model averaged partial regression coefficients. Also used to
#'   compare posterior standard deviations of the continuous piece of the model
#'   averaged distribution to corresponding standard deviations from individual
#'   posterior distributions in tabular form.
#' @param mcmc_draws A matrix of posterior samples of partial regression
#'   coefficients of interest with a column indicating which model the draws
#'   came from
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns to partition the plotting window into.
#' @param max_display An integer specifying the number of models to display in
#'   the plotting window.
#' @param mod_names A vector specifying the names of the models. If left blank,
#'   models will be named \eqn{M_1, M_2,..., M_K}, where \eqn{M_1} will have the
#'   largest posterior model probability and \eqn{M_K} has the smallest
#'   posterior model probability.
#' @param include_coef A vector specifying which partial regression coefficients
#'   to display.
#' @return The \code{MApp} plot and a table of posterior standard deviations for
#'   the top \code{max.display} individual models specified, along with the
#'   posterior standard deviation of the model averaged parameter.
#'
MApp_MCMC <- function(mcmc_draws, plot_wind = c(1,1),
                      max_display = NULL,
                      mod_names = NULL,
                      include_coef = NULL, ...) {
  # number of models
  p <- dim(mcmc_draws)[2] - 1
  K <- length(unique(mcmc_draws[,(p+1)]))
  var_names <- names(mcmc_draws)[1:p]

  coef_frames <- as.list(1:p)
  for(i in 1:p){
    coef_frames[[i]] <- mcmc_draws[,c(i,ncol(mcmc_test))]
    ma <- data.frame(coef_frames[[i]][,1], rep("MA", nrow(coef_frames[[i]])))
    names(coef_frames[[i]]) <- c("post_vec", "model")
    names(ma) <- names(coef_frames[[i]])
    coef_frames[[i]] <- rbind(ma, coef_frames[[i]])
  }

  names(coef_frames) <- var_names

  if (is.null(mod_names))
    mod_names <- paste("M", seq(1:K), sep = "")

  if (is.null(max_display)) {
    for (i in 1:p) {
      coef_frames[[i]] <- coef_frames[[i]]
    }
    max_display <- K

  } else {
    for (i in 1:p) {
      coef_frames[[i]]$mod <- as.numeric(coef_frames[[i]]$model)
      coef_frames[[i]] <- subset(coef_frames[[i]], mod <= (max_display + 1))
      coef_frames[[i]]$model <- factor(coef_frames[[i]]$model)
      coef_frames[[i]] <- coef_frames[[i]][, 1:2]
    }

  }

  # Find MaxMin of each beanplot

  MaxMin <- matrix(NA, nrow = p, ncol = 2)

  for (i in 1:p) {
    MaxMin[i, 1] <- min(coef_frames[[i]][, 1], na.rm = T)
    MaxMin[i, 2] <- max(coef_frames[[i]][, 1], na.rm = T)
  }

  # Only make plots for coefficients that appeared

  include_beans <- c(1:p)
  num_null <- sum((MaxMin[, 1] & MaxMin[, 2]) == 0)
  if (num_null == 0) {
    include_beans <- include_beans
  } else {
    include_beans <- include_beans[-which((MaxMin[, 1] & MaxMin[, 2]) == 0)]
  }

  # only display included coefficients in SD table
  if (is.null(include_coef)) {
    include_coef <- include_beans
  }
  # find cts piece of MA posterior SD for each coefficient and compare to individual models
  SD <- matrix(NA, nrow = p, ncol = (max_display + 1))

  if (length(include_coef) == p) {
    for (i in include_coef){
      sds <- coef_frames[[i]] %>% group_by(model) %>%
        summarise(SD = round(sd(post_vec, na.rm = T), 4))
      SD[i, ] <- t(sds[,"SD"])
    }

  } else {
    for (i in include_coef) {
      sds <- coef_frames[[i]] %>% group_by(model) %>%
        summarise(SD = round(sd(post_vec, na.rm = T), 4))
      SD[i, ] <- t(sds[,"SD"])
    }
    test.in <- apply(SD, 1, sum, na.rm = T)
    SD <- SD[-which(test.in ==0),]
  }

  SD[which(is.na(SD))] <- "---"
  SD <- as.table(SD)
  dimnames(SD)[[1]] <- var_names[include_coef]
  dimnames(SD)[[2]] <- levels(coef_frames[[i]]$model)


  # make the beanplot

  par(mfrow = c(plot_wind), las = 1, ask = T)
  for (i in include_coef) {

    # Initilize plot
    allplot <- beanplot::beanplot(post_vec ~ model,
                                  names = levels(coef_frames[[i]]$model),
                                  data = coef_frames[[i]],
                                  what = c(0, 0, 0, 0),
                                  border = 1,
                                  main = paste("Coefficient for ",
                                               var_names[i],
                                               sep = ""),
                                  horizontal = T, bw = "nrd0")

    # Add top model beans
    beanplot::beanplot(post_vec ~ model, side = "second",
                       show.names = F,
                       data = coef_frames[[i]],
                       add = T,
                       subset = model != "MA",
                       what = c(0, 1, 1, 0),
                       border = 1,
                       ll = 0.12,
                       horizontal = T,
                       col = c("lightgray"), bw = "nrd0")

    # Add MA posterior mean line
    MA.mean <- mean(subset(coef_frames[[i]], model == "MA")$post_vec, na.rm = T)
    segments(MA.mean, 1, MA.mean, (max_display + 1.35), col = "black", lty = 4)

    # Add MA bean
    beanplot::beanplot(post_vec ~ model, side = "second",
                       show.names = F,
                       data = coef_frames[[i]],
                       add = T,
                       subset = model == "MA",
                       what = c(0, 1, 1, 1),
                       border = 1,
                       ll = 0.12,
                       horizontal = T,
                       col = c("black", "gray", "red", "darkgrey"), bw = "nrd0")

    # Add horizontal lines
    modLines <- seq(1:(max_display + 1))
    abline(h = c(1, modLines), lty = 2, col = "gray")

    # Add model weights. display 0.000 for ~0's to convey rounding
    disp.weights <- round(weights[1:max_display], 3)
    disp.weights <- ifelse(disp.weights == 0, "< 0.000", as.character(disp.weights))

    # Add text to the plot
    text(rep(MaxMin[i, 1], max_display + 1),
         c(1:(max_display + 1) + 0.25),
         c("Avg", disp.weights),
         col = "red",
         cex = 1,
         pos = 1)

    # Add and make sure PIP's are not rounded to exactly 0.
    disp.PEP <- round(1 - PIP, 3)
    disp.PEP <- ifelse(disp.PEP == 0, "<0.000", as.character(disp.PEP))
    vadjust <- ifelse(max_display <= 10, 0.35, 0)
    text(mean(MaxMin[i, ]),
         vadjust,
         paste("Pr(beta_",
               paste(include_beans[i],
                     paste(".MA = 0",
                           paste("|y)",
                                 disp.PEP[i],
                                 sep = "="))),
               sep = ""),
         col = "black",
         cex = 1,
         pos = 3)
    text(MaxMin[i, 2], (max_display + 1.5),
         paste("Sum PMP", round(sum(weights), 3),
               sep = "="),
         col = "black",
         cex = 1,
         pos = 2)
  }

  readline("Table of posterior standard deviations: compare individual models to the continuous component of the MA posterior.\n\n (press enter to display table)\n")
  return(SD)
}

