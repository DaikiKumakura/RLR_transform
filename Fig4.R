## -------------------------------------------------------------------------------
library(rEDM)
library(tidyverse)
library(plotly)
library(tseries)


## -------------------------------------------------------------------------------
BestEmbed <- function(data = data, E = 1:10) {
    simp <- simplex(time_series = data, E = E, silent = T)
    opt_embed <- simp$E[which.min(simp$rmse)]
    return(opt_embed)
}

DoCCM <- function(data = data, E = 3,
				  lib_column = "paramecium",
				  target_column = "didinium",
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


## -------------------------------------------------------------------------------
data(paramecium_didinium)
str(paramecium_didinium)


## -------------------------------------------------------------------------------
para_adf <- adf.test(paramecium_didinium$paramecium)
para_adf


## -------------------------------------------------------------------------------
para_optembed <- BestEmbed(data = paramecium_didinium$paramecium, E = 1:10)
para_xmap_didi <- DoCCM(data = paramecium_didinium, E = para_optembed,
						lib_column = "paramecium",
						target_column = "didinium",
						lib_sizes = seq(10, 70, by = 1),
						num_samples = 100)
para_xmap_didi_means <- ccm_means(para_xmap_didi)


## -------------------------------------------------------------------------------
para_surr <- make_surrogate_ebisuzaki(paramecium_didinium$paramecium)
para_surr_out <- lapply(seq_len(NCOL(para_surr)), 
                           function(i) {
                             surr_ts <- para_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium$didinium),
                             					  E = para_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
para_surr_means <- do.call(rbind, para_surr_out)


## -------------------------------------------------------------------------------
para_surr_upper_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.975)
para_surr_lower_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
para_xmap_didi_res <- data.frame(para_xmap_didi_means$lib_size, 
				  				 para_xmap_didi_means$rho, 
				                 para_surr_upper_CI$rho, 
				                 para_surr_lower_CI$rho)
names(para_xmap_didi_res) <- c("lib_size", "ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower")


## -------------------------------------------------------------------------------
para_xmap_didi_plot <- ggplot(para_xmap_didi_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual"))
para_xmap_didi_plot <- para_xmap_didi_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI"))
para_xmap_didi_plot <- para_xmap_didi_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
para_xmap_didi_plot <- para_xmap_didi_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
para_xmap_didi_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
didi_adf <- adf.test(paramecium_didinium$didinium)
didi_adf


## -------------------------------------------------------------------------------
didi_optembed <- BestEmbed(data = paramecium_didinium$didinium, E = 1:10)
didi_xmap_para <- DoCCM(data = paramecium_didinium, E = didi_optembed,
						lib_column = "didinium",
						target_column = "paramecium",
						lib_sizes = seq(10, 71, by = 1),
						num_samples = 100)
didi_xmap_para_means <- ccm_means(didi_xmap_para)


## -------------------------------------------------------------------------------
didi_surr <- make_surrogate_ebisuzaki(paramecium_didinium$didinium)
didi_surr_out <- lapply(seq_len(NCOL(didi_surr)), 
                           function(i) {
                             surr_ts <- didi_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium$paramecium),
                             					  E = didi_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
didi_surr_means <- do.call(rbind, didi_surr_out)


## -------------------------------------------------------------------------------
didi_surr_upper_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.975)
didi_surr_lower_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
didi_xmap_para_res <- data.frame(didi_xmap_para_means$lib_size, 
				  				 didi_xmap_para_means$rho, 
				                 didi_surr_upper_CI$rho, 
				                 didi_surr_lower_CI$rho)
names(didi_xmap_para_res) <- c("lib_size", "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")


## -------------------------------------------------------------------------------
didi_xmap_para_plot <- ggplot(didi_xmap_para_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
didi_xmap_para_plot <- didi_xmap_para_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
didi_xmap_para_plot <- didi_xmap_para_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
didi_xmap_para_plot <- didi_xmap_para_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
didi_xmap_para_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
para_xmap_didi_res <- para_xmap_didi_res[-61, ]
abs_res <- data.frame(para_xmap_didi_res$lib_size, 
				      para_xmap_didi_res$ParaXmapDidi_actual,
				      para_xmap_didi_res$ParaXmapDidi_upper,
				      para_xmap_didi_res$ParaXmapDidi_lower,
				      didi_xmap_para_res$DidiXmapPara_actual,
				      didi_xmap_para_res$DidiXmapPara_upper,
				      didi_xmap_para_res$DidiXmapPara_lower)
names(abs_res) <- c("lib_size", 
					"ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower",
			        "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")
abs_res <- mutate(abs_res, ParaXmapDidi_lower = if_else(condition = ParaXmapDidi_lower < 0, true = 0, false = ParaXmapDidi_lower))
abs_res <- mutate(abs_res, DidiXmapPara_lower = if_else(condition = DidiXmapPara_lower < 0, true = 0, false = DidiXmapPara_lower))


## -------------------------------------------------------------------------------
abs_plot <- ggplot(abs_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
abs_plot <- abs_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
abs_plot <- abs_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
abs_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
paramecium_didinium_comp <- data.frame(paramecium_didinium$paramecium, paramecium_didinium$didinium)
paramecium_didinium_comp$sum <- c(apply(paramecium_didinium, 1, sum))
names(paramecium_didinium_comp) <- c("paramecium", "didinium", "sum")
paramecium_didinium_comp$paramecium_comp <- c(paramecium_didinium_comp$paramecium / paramecium_didinium_comp$sum)
paramecium_didinium_comp$didinium_comp <- c(paramecium_didinium_comp$didinium / paramecium_didinium_comp$sum)
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp$time <- paramecium_didinium$time
paramecium_didinium_comp <- paramecium_didinium_comp[, c(3, 1, 2)]
names(paramecium_didinium_comp) <- c("time", "paramecium", "didinium")


## -------------------------------------------------------------------------------
para_adf <- adf.test(paramecium_didinium_comp$paramecium)
para_adf


## -------------------------------------------------------------------------------
para_optembed <- BestEmbed(data = paramecium_didinium_comp$paramecium, E = 1:10)
para_xmap_didi <- DoCCM(data = paramecium_didinium_comp, E = para_optembed,
						lib_column = "paramecium",
						target_column = "didinium",
						lib_sizes = seq(10, 71, by = 1),
						num_samples = 100)
para_xmap_didi_means <- ccm_means(para_xmap_didi)


## -------------------------------------------------------------------------------
para_surr <- make_surrogate_ebisuzaki(paramecium_didinium_comp$paramecium)
para_surr_out <- lapply(seq_len(NCOL(para_surr)), 
                           function(i) {
                             surr_ts <- para_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium_comp$didinium),
                             					  E = para_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
para_surr_means <- do.call(rbind, para_surr_out)


## -------------------------------------------------------------------------------
para_surr_upper_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.975)
para_surr_lower_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
para_xmap_didi_res <- data.frame(para_xmap_didi_means$lib_size, 
				  				 para_xmap_didi_means$rho, 
				                 para_surr_upper_CI$rho, 
				                 para_surr_lower_CI$rho)
names(para_xmap_didi_res) <- c("lib_size", "ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower")


## -------------------------------------------------------------------------------
para_xmap_didi_plot <- ggplot(para_xmap_didi_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual"))
para_xmap_didi_plot <- para_xmap_didi_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI"))
para_xmap_didi_plot <- para_xmap_didi_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
para_xmap_didi_plot <- para_xmap_didi_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
para_xmap_didi_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
didi_adf <- adf.test(paramecium_didinium_comp$didinium)
didi_adf


## -------------------------------------------------------------------------------
didi_optembed <- BestEmbed(data = paramecium_didinium_comp$didinium, E = 1:10)
didi_xmap_para <- DoCCM(data = paramecium_didinium_comp, E = didi_optembed,
						lib_column = "didinium",
						target_column = "paramecium",
						lib_sizes = seq(10, 71, by = 1),
						num_samples = 100)
didi_xmap_para_means <- ccm_means(didi_xmap_para)


## -------------------------------------------------------------------------------
didi_surr <- make_surrogate_ebisuzaki(paramecium_didinium_comp$didinium)
didi_surr_out <- lapply(seq_len(NCOL(didi_surr)), 
                           function(i) {
                             surr_ts <- didi_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium_comp$paramecium),
                             					  E = didi_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
didi_surr_means <- do.call(rbind, didi_surr_out)


## -------------------------------------------------------------------------------
didi_surr_upper_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.975)
didi_surr_lower_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
didi_xmap_para_res <- data.frame(didi_xmap_para_means$lib_size, 
				  				 didi_xmap_para_means$rho, 
				                 didi_surr_upper_CI$rho, 
				                 didi_surr_lower_CI$rho)
names(didi_xmap_para_res) <- c("lib_size", "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")


## -------------------------------------------------------------------------------
didi_xmap_para_plot <- ggplot(didi_xmap_para_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
didi_xmap_para_plot <- didi_xmap_para_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
didi_xmap_para_plot <- didi_xmap_para_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
didi_xmap_para_plot <- didi_xmap_para_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
didi_xmap_para_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
comp_res <- data.frame(para_xmap_didi_res$lib_size, 
				      para_xmap_didi_res$ParaXmapDidi_actual,
				      para_xmap_didi_res$ParaXmapDidi_upper,
				      para_xmap_didi_res$ParaXmapDidi_lower,
				      didi_xmap_para_res$DidiXmapPara_actual,
				      didi_xmap_para_res$DidiXmapPara_upper,
				      didi_xmap_para_res$DidiXmapPara_lower)
names(comp_res) <- c("lib_size", 
					 "ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower",
			         "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")
comp_res <- mutate(comp_res, ParaXmapDidi_lower = if_else(condition = ParaXmapDidi_lower < 0, true = 0, false = ParaXmapDidi_lower))
comp_res <- mutate(comp_res, DidiXmapPara_lower = if_else(condition = DidiXmapPara_lower < 0, true = 0, false = DidiXmapPara_lower))


## -------------------------------------------------------------------------------
comp_plot <- ggplot(comp_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
comp_plot <- comp_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
comp_plot <- comp_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
comp_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
paramecium_didinium_comp <- data.frame(paramecium_didinium$paramecium, paramecium_didinium$didinium)
paramecium_didinium_comp$sum <- c(apply(paramecium_didinium, 1, sum))
names(paramecium_didinium_comp) <- c("paramecium", "didinium", "sum")
paramecium_didinium_comp$paramecium_comp <- c(paramecium_didinium_comp$paramecium / paramecium_didinium_comp$sum)
paramecium_didinium_comp$didinium_comp <- c(paramecium_didinium_comp$didinium / paramecium_didinium_comp$sum)
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp <- paramecium_didinium_comp[, -1]
paramecium_didinium_comp$time <- paramecium_didinium$time
paramecium_didinium_comp <- paramecium_didinium_comp[, c(3, 1, 2)]
names(paramecium_didinium_comp) <- c("time", "paramecium", "didinium")

# Convert composition to RLR-transformed
paramecium_didinium_rlr <- paramecium_didinium_comp
paramecium_didinium_rlr$ln_paramecium <- c(log(paramecium_didinium_rlr$didinium / paramecium_didinium_rlr$paramecium))
paramecium_didinium_rlr$ln_didinium <- c(log(paramecium_didinium_rlr$paramecium / paramecium_didinium_rlr$didinium))
paramecium_didinium_rlr$paramecium_rlr <- c(-1 / (1 + paramecium_didinium_rlr$ln_paramecium))
paramecium_didinium_rlr$didinium_rlr <- c(-1 / (1 + paramecium_didinium_rlr$ln_didinium))
paramecium_didinium_rlr <- paramecium_didinium_rlr[, -2]
paramecium_didinium_rlr <- paramecium_didinium_rlr[, -2]
paramecium_didinium_rlr <- paramecium_didinium_rlr[, -2]
paramecium_didinium_rlr <- paramecium_didinium_rlr[, -2]
names(paramecium_didinium_rlr) <- c("time", "paramecium", "didinium")


## -------------------------------------------------------------------------------
para_adf <- adf.test(paramecium_didinium_rlr$paramecium)
para_adf


## -------------------------------------------------------------------------------
para_optembed <- BestEmbed(data = paramecium_didinium_rlr$paramecium, E = 1:10)
para_xmap_didi <- DoCCM(data = paramecium_didinium_rlr, E = para_optembed,
						lib_column = "paramecium",
						target_column = "didinium",
						lib_sizes = seq(10, 71, by = 1),
						num_samples = 100)
para_xmap_didi_means <- ccm_means(para_xmap_didi)


## -------------------------------------------------------------------------------
para_surr <- make_surrogate_ebisuzaki(paramecium_didinium_rlr$paramecium)
para_surr_out <- lapply(seq_len(NCOL(para_surr)), 
                           function(i) {
                             surr_ts <- para_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium_rlr$didinium),
                             					  E = para_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
para_surr_means <- do.call(rbind, para_surr_out)


## -------------------------------------------------------------------------------
para_surr_upper_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.975)
para_surr_lower_CI <- ccm_means(para_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
para_xmap_didi_res <- data.frame(para_xmap_didi_means$lib_size, 
				  				 para_xmap_didi_means$rho, 
				                 para_surr_upper_CI$rho, 
				                 para_surr_lower_CI$rho)
names(para_xmap_didi_res) <- c("lib_size", "ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower")


## -------------------------------------------------------------------------------
para_xmap_didi_plot <- ggplot(para_xmap_didi_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual"))
para_xmap_didi_plot <- para_xmap_didi_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI"))
para_xmap_didi_plot <- para_xmap_didi_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
para_xmap_didi_plot <- para_xmap_didi_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
para_xmap_didi_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
didi_adf <- adf.test(paramecium_didinium_rlr$didinium)
didi_adf


## -------------------------------------------------------------------------------
didi_optembed <- BestEmbed(data = paramecium_didinium_rlr$didinium, E = 1:10)
didi_xmap_para <- DoCCM(data = paramecium_didinium_rlr, E = didi_optembed,
						lib_column = "didinium",
						target_column = "paramecium",
						lib_sizes = seq(10, 71, by = 1),
						num_samples = 100)
didi_xmap_para_means <- ccm_means(didi_xmap_para)


## -------------------------------------------------------------------------------
didi_surr <- make_surrogate_ebisuzaki(paramecium_didinium_rlr$didinium)
didi_surr_out <- lapply(seq_len(NCOL(didi_surr)), 
                           function(i) {
                             surr_ts <- didi_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, paramecium_didinium_rlr$paramecium),
                             					  E = didi_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, 71, by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
didi_surr_means <- do.call(rbind, didi_surr_out)


## -------------------------------------------------------------------------------
didi_surr_upper_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.975)
didi_surr_lower_CI <- ccm_means(didi_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
didi_xmap_para_res <- data.frame(didi_xmap_para_means$lib_size, 
				  				 didi_xmap_para_means$rho, 
				                 didi_surr_upper_CI$rho, 
				                 didi_surr_lower_CI$rho)
names(didi_xmap_para_res) <- c("lib_size", "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")


## -------------------------------------------------------------------------------
didi_xmap_para_plot <- ggplot(didi_xmap_para_res) + geom_line(size = 1, 
															  alpha = 1.5, 
															  aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
didi_xmap_para_plot <- didi_xmap_para_plot + geom_ribbon(alpha = 0.5, 
														 aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
didi_xmap_para_plot <- didi_xmap_para_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
didi_xmap_para_plot <- didi_xmap_para_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
didi_xmap_para_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
didi_xmap_para_res <- didi_xmap_para_res[-57, ]
didi_xmap_para_res <- didi_xmap_para_res[-56, ]
didi_xmap_para_res <- didi_xmap_para_res[-55, ]
didi_xmap_para_res <- didi_xmap_para_res[-54, ]
rlr_res <- data.frame(para_xmap_didi_res$lib_size, 
				      para_xmap_didi_res$ParaXmapDidi_actual,
				      para_xmap_didi_res$ParaXmapDidi_upper,
				      para_xmap_didi_res$ParaXmapDidi_lower,
				      didi_xmap_para_res$DidiXmapPara_actual,
				      didi_xmap_para_res$DidiXmapPara_upper,
				      didi_xmap_para_res$DidiXmapPara_lower)
names(rlr_res) <- c("lib_size", 
					 "ParaXmapDidi_actual", "ParaXmapDidi_upper", "ParaXmapDidi_lower",
			         "DidiXmapPara_actual", "DidiXmapPara_upper", "DidiXmapPara_lower")
rlr_res <- mutate(rlr_res, DidiXmapPara_actual = if_else(condition = DidiXmapPara_actual < 0, true = 0, false = DidiXmapPara_actual))
rlr_res <- mutate(rlr_res, ParaXmapDidi_actual = if_else(condition = ParaXmapDidi_actual < 0, true = 0, false = ParaXmapDidi_actual))
rlr_res <- mutate(rlr_res, ParaXmapDidi_lower = if_else(condition = ParaXmapDidi_lower < 0, true = 0, false = ParaXmapDidi_lower))
rlr_res <- mutate(rlr_res, DidiXmapPara_lower = if_else(condition = DidiXmapPara_lower < 0, true = 0, false = DidiXmapPara_lower))


## -------------------------------------------------------------------------------
rlr_plot <- ggplot(rlr_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = ParaXmapDidi_actual, colour = "ParaXmapDidi_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = DidiXmapPara_actual, colour = "DidiXmapPara_actual"))
rlr_plot <- rlr_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = ParaXmapDidi_lower, ymax = ParaXmapDidi_upper, fill = "ParaXmapDidi_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = DidiXmapPara_lower, ymax = DidiXmapPara_upper, fill = "DidiXmapPara_95%CI"))
rlr_plot <- rlr_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
rlr_plot 
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
sessionInfo()

