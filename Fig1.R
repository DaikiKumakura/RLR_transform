## -------------------------------------------------------------------------------
library(rEDM)
library(deSolve)
library(rgr)
library(tidyverse)
library(tseries)


## -------------------------------------------------------------------------------
BestEmbed <- function(data = data, E = 1:50) {
    simp <- simplex(time_series = data, E = E, silent = T)
    opt_embed <- simp$E[which.min(simp$rmse)]
    return(opt_embed)
}

DoCCM <- function(data = data, E = 3,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(10, 71, by = 1),
                  num_samples = 100)
{
    ccm(data, E = E,
        lib_column = lib_column,
        target_column = target_column,
        lib_sizes = lib_sizes,
        num_samples = num_samples,
        random_libs = TRUE, replace = FALSE, silent = TRUE)
}

raw2comp <- function(data) {
        dataout <- data
        dataout_comp <- dataout[, c(2:3)]
        dataout_comp$sum <- c(apply(dataout_comp, 1, sum))
        dataout_comp$x_re <- c(dataout$x / dataout_comp$sum)
        dataout_comp$y_re <- c(dataout$y / dataout_comp$sum)
        dataout$x_re <- c(dataout_comp$x_re)
        dataout$y_re <- c(dataout_comp$y_re)
        dataout_composition <- dataout[, c(1, 4:5)]
        names(dataout_composition) <- c("time", "x", "y")
        dataout <- dataout_composition
        d <- dataout
        return(d)
}

comp2rlr <- function(data) {
    d <- data
    d$lx <- c(log(d$y / d$x))
    d$ly <- c(log(d$x / d$y))
    d$x <- c(- 1 / (1 + d$lx))
    d$y <- c(- 1 / (1 + d$ly))
    df <- d[, c(1:3)]
    d <- df
    return(d)
}


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) 
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
raw_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(raw_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
raw_res <- mutate(raw_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
raw_res <- mutate(raw_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
raw_res <- mutate(raw_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
raw_res <- mutate(raw_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
raw_res <- mutate(raw_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
raw_res <- mutate(raw_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
raw_plot <- ggplot(raw_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
raw_plot <- raw_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
raw_plot <- raw_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "Raw")
raw_plot


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) 
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- raw2comp(d)


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
comp_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(comp_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
comp_res <- mutate(comp_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
comp_res <- mutate(comp_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
comp_res <- mutate(comp_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
comp_res <- mutate(comp_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
comp_res <- mutate(comp_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
comp_res <- mutate(comp_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
comp_plot <- ggplot(comp_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
comp_plot <- comp_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
comp_plot <- comp_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "Composition")
comp_plot


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) 
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- raw2comp(d)
d <- comp2rlr(d)


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
rlr_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(rlr_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
rlr_res <- mutate(rlr_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
rlr_res <- mutate(rlr_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
rlr_res <- mutate(rlr_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
rlr_res <- mutate(rlr_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
rlr_res <- mutate(rlr_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
rlr_res <- mutate(rlr_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
rlr_plot <- ggplot(rlr_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
rlr_plot <- rlr_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
rlr_plot <- rlr_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "RLR-transformed")
rlr_plot


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) 
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
dataout <- d
dataout_logAbs <- dataout[, c(2:3)]
dataout_logAbs$x_log <- c(log(dataout_logAbs$x + 1))
dataout_logAbs$y_log <- c(log(dataout_logAbs$y + 1))
dataout$x_re <- c(dataout_logAbs$x_log)
dataout$y_re <- c(dataout_logAbs$y_log)
dataout_logAbsolute <- dataout[, c(1, 4:5)]
names(dataout_logAbsolute) <- c("time", "x", "y")
dataout <- dataout_logAbsolute
d <- dataout


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
lograw_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(lograw_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
lograw_res <- mutate(lograw_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
lograw_res <- mutate(lograw_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
lograw_res <- mutate(lograw_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
lograw_res <- mutate(lograw_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
lograw_res <- mutate(lograw_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
lograw_res <- mutate(lograw_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
lograw_plot <- ggplot(lograw_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
lograw_plot <- lograw_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
lograw_plot <- lograw_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "Log(Raw+1)")
lograw_plot


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1)
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- raw2comp(d)
dataout_logComp <- d[, c(2:3)]
dataout_logComp$x_log <- c(log(dataout_logComp$x + 1))
dataout_logComp$y_log <- c(log(dataout_logComp$y + 1))
d$x_re <- c(dataout_logComp$x_log)
d$y_re <- c(dataout_logComp$y_log)
dataout_logComposition <- d[, c(1, 4:5)]
names(dataout_logComposition) <- c("time", "x", "y")
d <- dataout_logComposition


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
logcomp_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(logcomp_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
logcomp_res <- mutate(logcomp_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
logcomp_res <- mutate(logcomp_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
logcomp_res <- mutate(logcomp_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
logcomp_res <- mutate(logcomp_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
logcomp_res <- mutate(logcomp_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
logcomp_res <- mutate(logcomp_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
logcomp_plot <- ggplot(logcomp_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
logcomp_plot <- logcomp_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
logcomp_plot <- logcomp_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "Log(Composition+1)")
logcomp_plot


## -------------------------------------------------------------------------------
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) 
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
dataout <- d
dataout_clr <- clr(dataout[2:3])
dataout_clr <- as.data.frame(dataout_clr)
dataout$x_re <- c(dataout_clr$x)
dataout$y_re <- c(dataout_clr$y)
dataout_clrTransform <- dataout[, c(1, 4:5)]
names(dataout_clrTransform) <- c("time", "x", "y")
dataout <- dataout_clrTransform
d <- dataout


## -------------------------------------------------------------------------------
adf.test(d$y)


## -------------------------------------------------------------------------------
y_optembed <- BestEmbed(data = d$y, E = 1:50)
y_xmap_x <- DoCCM(data = d, E = y_optembed,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(50, NROW(d$y), by = 50),
                  num_samples = 100)
y_xmap_x_means <- ccm_means(y_xmap_x)


## -------------------------------------------------------------------------------
y_surr <- make_surrogate_ebisuzaki(d$y)
y_surr_out <- lapply(seq_len(NCOL(y_surr)),
                     function(i) {
                       surr_ts <- y_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$x),
                                         E = y_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$y), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
y_surr_means <- do.call(rbind, y_surr_out)


## -------------------------------------------------------------------------------
y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
y_surr_lower_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
y_xmap_x_res <- data.frame(y_xmap_x_means$lib_size, 
                           y_xmap_x_means$rho, 
                           y_surr_upper_CI$rho, 
                           y_surr_lower_CI$rho)
names(y_xmap_x_res) <- c("lib_size", "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower")


## -------------------------------------------------------------------------------
y_xmap_x_plot <- ggplot(y_xmap_x_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual"))
y_xmap_x_plot <- y_xmap_x_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI"))
y_xmap_x_plot <- y_xmap_x_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
y_xmap_x_plot <- y_xmap_x_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
y_xmap_x_plot


## -------------------------------------------------------------------------------
adf.test(d$x)


## -------------------------------------------------------------------------------
x_optembed <- BestEmbed(data = d$x, E = 1:50)
x_xmap_y <- DoCCM(data = d, E = x_optembed,
                  lib_column = "x",
                  target_column = "y",
                  lib_sizes = seq(50, NROW(d$x), by = 50),
                  num_samples = 100)
x_xmap_y_means <- ccm_means(x_xmap_y)


## -------------------------------------------------------------------------------
x_surr <- make_surrogate_ebisuzaki(d$x)
x_surr_out <- lapply(seq_len(NCOL(x_surr)),
                     function(i) {
                       surr_ts <- x_surr[, i]
                       xmap_out <- DoCCM(data = cbind(surr_ts, d$y),
                                         E = x_optembed,
                                         lib_column = 1,
                                         target_column = 2,
                                         lib_sizes = seq(50, NROW(d$x), by = 50),
                                         num_samples = 100)
                       xmap_means <- ccm_means(xmap_out)
                     })
x_surr_means <- do.call(rbind, x_surr_out)


## -------------------------------------------------------------------------------
x_surr_upper_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.975)
x_surr_lower_CI <- ccm_means(x_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
x_xmap_y_res <- data.frame(x_xmap_y_means$lib_size, 
                           x_xmap_y_means$rho, 
                           x_surr_upper_CI$rho, 
                           x_surr_lower_CI$rho)
names(x_xmap_y_res) <- c("lib_size", "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")


## -------------------------------------------------------------------------------
x_xmap_y_plot <- ggplot(x_xmap_y_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
x_xmap_y_plot <- x_xmap_y_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_x_95%CI"))
x_xmap_y_plot <- x_xmap_y_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
x_xmap_y_plot <- x_xmap_y_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
x_xmap_y_plot


## -------------------------------------------------------------------------------
clr_res <- data.frame(y_xmap_x_res$lib_size, 
                      y_xmap_x_res$y_Xmap_x_actual,
                      y_xmap_x_res$y_Xmap_x_upper,
                      y_xmap_x_res$y_Xmap_x_lower,
                      x_xmap_y_res$x_Xmap_y_actual,
                      x_xmap_y_res$x_Xmap_y_upper,
                      x_xmap_y_res$x_Xmap_y_lower)
names(clr_res) <- c("lib_size", 
                    "y_Xmap_x_actual", "y_Xmap_x_upper", "y_Xmap_x_lower",
                    "x_Xmap_y_actual", "x_Xmap_y_upper", "x_Xmap_y_lower")
clr_res <- mutate(clr_res, y_Xmap_x_lower = if_else(condition = y_Xmap_x_lower < 0, true = 0, false = y_Xmap_x_lower))
clr_res <- mutate(clr_res, x_Xmap_y_lower = if_else(condition = x_Xmap_y_lower < 0, true = 0, false = x_Xmap_y_lower))
clr_res <- mutate(clr_res, y_Xmap_x_upper = if_else(condition = y_Xmap_x_upper < 0, true = 0, false = y_Xmap_x_upper))
clr_res <- mutate(clr_res, x_Xmap_y_upper = if_else(condition = x_Xmap_y_upper < 0, true = 0, false = x_Xmap_y_upper))
clr_res <- mutate(clr_res, y_Xmap_x_actual = if_else(condition = y_Xmap_x_actual < 0, true = 0, false = y_Xmap_x_actual))
clr_res <- mutate(clr_res, x_Xmap_y_actual = if_else(condition = x_Xmap_y_actual < 0, true = 0, false = x_Xmap_y_actual))


## -------------------------------------------------------------------------------
clr_plot <- ggplot(clr_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = y_Xmap_x_actual, colour = "y_Xmap_x_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = x_Xmap_y_actual, colour = "x_Xmap_y_actual"))
clr_plot <- clr_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = y_Xmap_x_lower, ymax = y_Xmap_x_upper, fill = "y_Xmap_x_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = x_Xmap_y_lower, ymax = x_Xmap_y_upper, fill = "x_Xmap_y_95%CI"))
clr_plot <- clr_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))), title = "clr-Transformed")
clr_plot


## -------------------------------------------------------------------------------
sessionInfo()

