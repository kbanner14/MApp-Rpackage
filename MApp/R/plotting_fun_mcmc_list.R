#' @title Model Averaged Posteriors Plot
#' @description Works with lists of posterior draws of partial regression
#'   coefficients from multiple models. Creates a graphical display showing all
#'   relevant output from the model averaging procedure. Used to graphically
#'   compare posterior distributions of partial regression coefficients from
#'   individual models in the model set to the posterior distributions resulting
#'   from the model averaged partial regression coefficients. Also used to
#'   compare posterior standard deviations of the model
#'   averaged distribution to corresponding standard deviations from individual
#'   posterior distributions in tabular form.
#' @param mcmc_list A list of posterior samples of partial regression
#'   coefficients of interest for each model in the model set.
#' @param plot_wind A vector of length 2 specifying the number of rows and
#'   columns to partition the plotting window into.
#' @param max_display An integer specifying the number of models to display in
#'   the plotting window.
#' @param mod_names A vector specifying the names of the models. If left blank,
#'   models will be named \eqn{M_1, M_2,..., M_J}. 
#' @param include_coef A vector specifying which partial regression coefficients
#'   to display.
#' @return MAP plot and a table of posterior standard deviations for
#'   the top \code{max_display} individual models specified, along with the
#'   posterior standard deviation from the model averaged result.
MApp_list <- function(mcmc_list, plot_wind, max_display = NULL,
                       mod_names = NULL, include_coef = NULL, ...) {
  
    `%>%` <- magrittr::`%>%`
    K <- length(mcmc_list)
    p <- dim(mcmc_list[[1]])[2]  # Number of parameters
    if (is.null(mod_names))
        mod_names <- paste("M", seq(1:K), sep = "")
    names(mcmc_list) <- paste("M", seq(1:K), sep = "")  
    # Create names for models
    var_names <- names(mcmc_list[[1]])

    # Create a data frame for each coefficient estimate 
    # and track which model it came from
    
    draws_mod <- sapply(mcmc_list, function(x){ dim(x)[1]})
    pmps <- draws_mod/sum(draws_mod)
    
    all_draws <- do.call(rbind, mcmc_list)
    mod_vec <- character()
    for(i in 1:length(mod_names)){
      vec <- rep(mod_names[i], draws_mod[i])
      mod_vec <- c(mod_vec, vec)
    }
    all_draws$Model <- factor(mod_vec)
    coef_frames <- list(1:p)
    
    for (i in 1:p) {
      
      coef_frames[[i]] <- data.frame(all_draws[,i], all_draws[,"Model"])
        names(coef_frames[[i]]) <- c("Post.Vec", "Model")
        coef_frames[[i]]$Model <- factor(coef_frames[[i]]$Model, levels = mod_names)
    }

    for (i in 1:p) {
        ma_post <- coef_frames[[i]]$Post.Vec
        ma_post <- data.frame(Post.Vec = ma_post, Model = "MA")
        coef_frames[[i]] <- rbind(ma_post, coef_frames[[i]])
    }

    # index of three commonly chosen models
    full <- which(apply(inmat, 1, sum) == p)
    choiceIdx <- ifelse(full == 1 , 1, full)
    choiceIdx <- c(1,choiceIdx)
    choiceIdx <- unique(choiceIdx)
    
    if (is.null(max_display)) {
      for (i in 1:p) {
        coef_frames[[i]] <- coef_frames[[i]]
      }
      max_display <- K
    } else {
      if (max_display == "common3") {
        for (i in 1:p) {
          coef_frames[[i]]$mod <- as.numeric(coef_frames[[i]]$Model)
          coef_frames[[i]] <- subset(coef_frames[[i]], mod == 1 |
                                       mod == choiceIdx[1]+1 |
                                       mod == choiceIdx[2]+1)
          coef_frames[[i]]$Model <- factor(coef_frames[[i]]$Model)
          coef_frames[[i]] <- coef_frames[[i]][, 1:2]
        }
        max_display <- length(choiceIdx)
      } else {
        for (i in 1:p) {
          coef_frames[[i]]$mod <- as.numeric(coef_frames[[i]]$Model)
          coef_frames[[i]] <- subset(coef_frames[[i]], mod <= (max_display + 1))
          coef_frames[[i]]$Model <- factor(coef_frames[[i]]$Model)
          coef_frames[[i]] <- coef_frames[[i]][, 1:2]
        }
        
      }}
    
    # Find MaxMin of each beanplot
    
    MaxMin <- matrix(NA, nrow = p, ncol = 2)
    
    for (i in 1:p) {
      if(sum(coef_frames[[i]][,1], na.rm = T) == 0){
        MaxMin[i,1] <- 0
        MaxMin[i,2] <- 0
      } else {
        MaxMin[i, 1] <- min(coef_frames[[i]][, 1], na.rm = T)
        MaxMin[i, 2] <- max(coef_frames[[i]][, 1], na.rm = T)
      }
    }
    
    # Only make plots for coefficients that appeared
    
    include_beans <- c(1:p)
    num_null <- sum(MaxMin[, 1] == 0 & MaxMin[, 2] == 0)
    if (num_null == 0) {
      include_beans <- include_beans
    } else {
      include_beans <- include_beans[-which(MaxMin[,2] == 0 & MaxMin[,1] == 0)]
    }
    
    # only display included coefficients in SD table
    if (is.null(include_coef)) {
      include_coef <- include_beans
    }
    # find cts piece of MA posterior SD for each coefficient and compare 
    # to individual models. Also find MA post from draws - to compare to 
    # other SDs and also analytical result from BMS. 
    SD <- matrix(NA, nrow = length(include_coef), ncol = (max_display + 2))
    
    if (length(include_coef) == p) {
      for (i in include_coef){
        sds <- coef_frames[[i]] %>% dplyr::group_by(Model) %>%
          dplyr::summarise(SD = round(sd(Post.Vec, na.rm = T), 4))
        coef_full <- coef_frames[[i]]
        coef_full[is.na(coef_full[,1]), 1] <- 0
        sd_full <- coef_full %>% dplyr::group_by(Model)%>% 
          dplyr::summarise(SD_full = round(sd(Post.Vec), 4))
        sds <- unlist(c(sd_full[1,2], sds[,2]))
        SD[i, ] <- sds
        
      }
    } else {
      j = 1
      for (i in include_coef) {
        sds <- coef_frames[[i]] %>% dplyr::group_by(Model) %>%
          dplyr::summarise(SD = round(sd(Post.Vec, na.rm = T), 4))
        coef_full <- coef_frames[[i]]
        coef_full[is.na(coef_full[,1]), 1] <- 0
        sd_full <- coef_full %>% dplyr::group_by(Model)%>% 
          dplyr::summarise(SD_full = round(sd(Post.Vec), 4))
        sds <- unlist(c(sd_full[1,2], sds[,2]))
        SD[j, ] <- sds
        j = j+1
      }
    }
    
    SD[which(is.na(SD))] <- "---"
    SD <- as.table(SD)
    
    SD <- rbind(SD,c("---","---",round(pmps[1:max_display],4)))
    dimnames(SD)[[2]] <- c("MA", "MA_cts", levels(coef_frames[[i]]$Model)[-1])
    
    dimnames(SD)[[1]] <- c(var_names[include_coef], "PMP")
    # Add and make sure PIP's are not rounded to exactly 0.
    PEP <- sapply(coef_frames, function(x){
      sum(is.na(subset(x, Model == "MA")$Post.Vec))/dim(subset(x, Model == "MA"))[1]}
)
    disp.PEP <- round(PEP, 3)
    disp.PEP <- ifelse(disp.PEP == 0, "<0.000", as.character(disp.PEP))
    # make vertical adjust for larger plots
    vadjust <- ifelse(max_display <= 8, 0.35, 0.05)
    
    # make the beanplot
    par(mfrow = c(plot_wind), las = 1, mar = c(2.1,3.1,3.1,.5))
    for (i in include_coef) {
      
      # Initilize plot
      allplot <- beanplot::beanplot(Post.Vec ~ Model,
                                    names = levels(coef_frames[[i]]$Model),
                                    data = coef_frames[[i]],
                                    what = c(0, 0, 0, 0),
                                    border = 1,
                                    main = paste("Coefficient for ",
                                                 var_names[i],
                                                 sep = ""),
                                    horizontal = T, bw = "nrd0", log = "")
      
      # Add top model beans
      beanplot::beanplot(Post.Vec ~ Model, side = "second",
                         show.names = F,
                         data = coef_frames[[i]],
                         add = T,
                         subset = Model != "MA",
                         what = c(0, 1, 1, 0),
                         border = 1,
                         ll = 0.12,
                         horizontal = T,
                         col = c("lightgray"), bw = "nrd0",
                         log = "")
      
      # Add MA posterior mean line
      ma.mean <- subset(coef_frames[[i]], Model =="MA")$Post.Vec
      ma.mean[is.na(ma.mean)] <- 0
      MA.mean <- mean(ma.mean)
      segments(MA.mean, 1, MA.mean, (max_display + 1.35), col = "black", lty = 4)
      
      # Add MA bean
      beanplot::beanplot(Post.Vec ~ Model, side = "second",
                         show.names = F,
                         data = coef_frames[[i]],
                         add = T,
                         subset = Model == "MA",
                         what = c(0, 1, 1, 1),
                         border = 1,
                         ll = 0.12,
                         horizontal = T,
                         col = c("black", "gray", "red", "darkgrey"),
                         bw = "nrd0", log = "")
      
      # Add horizontal lines
      modLines <- seq(1:(max_display + 1))
      abline(h = c(1, modLines), lty = 2, col = "gray")
      
      # Add model weights. display 0.000 for ~0's to convey rounding
      disp.weights <- round(pmps[1:max_display], 3)
      disp.weights <- ifelse(disp.weights == 0, "< 0.000", as.character(disp.weights))
      
      # Add text to the plot
      text(rep(MaxMin[i, 1], max_display + 1),
           c(1:(max_display + 1) + 0.25),
           c("Avg", disp.weights),
           col = "red",
           cex = 1,
           pos = 1)
      
      # Add PEP's
      if (disp.PEP[i] == "< 0.000") {
        text(mean(MaxMin[i, ]),
             vadjust,
             bquote("Pr("~beta[.(include_beans[i])~", MA"]~"= 0 | y )"~.(disp.PEP[i])), 
             col = "black",
             cex = 1,
             pos = 3)
      } else {
        text(mean(MaxMin[i, ]),
             vadjust,
             bquote("Pr("~beta[.(include_beans[i])~", MA"]~"= 0 | y ) ="~.(disp.PEP[i])), 
             col = "black",
             cex = 1,
             pos = 3)
      }
      text(MaxMin[i, 2], (max_display + 1.5),
           bquote(sum(w[j], j = "j=1", .(nrow(inmat))) == .(round(sum(pmps), 3))),
           col = "black",
           cex = 1,
           pos = 2)
    }
    
    
    # provide message to user before table is printed
    readline("Table of posterior standard deviations: 
             compare individual models to the MA posterior.\n\n 
             (press enter to display table)\n")
    # print table
    return(SD)
}
