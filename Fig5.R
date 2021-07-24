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
d <- read.csv("Fig5.csv")
d$clo <- d$Clostridium
d$gna <- d$gnavus


## -------------------------------------------------------------------------------
adf.test(d$gna)


## -------------------------------------------------------------------------------
gna_optembed <- BestEmbed(data = d$gna, E = 1:10)
gna_xmap_clo <- DoCCM(data = d, E = gna_optembed,
					  lib_column = "gna",
					  target_column = "clo",
					  lib_sizes = seq(10, NROW(d$gna), by = 1),
					  num_samples = 100)
gna_xmap_clo_means <- ccm_means(gna_xmap_clo)


## -------------------------------------------------------------------------------
gna_surr <- make_surrogate_ebisuzaki(d$gna)
gna_surr_out <- lapply(seq_len(NCOL(gna_surr)), 
                           function(i) {
                             surr_ts <- gna_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$clo),
                             					  E = gna_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$gna), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
gna_surr_means <- do.call(rbind, gna_surr_out)


## -------------------------------------------------------------------------------
gna_surr_upper_CI <- ccm_means(gna_surr_means, FUN = quantile, probs = 0.975)
gna_surr_lower_CI <- ccm_means(gna_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
gna_xmap_clo_res <- data.frame(gna_xmap_clo_means$lib_size, 
				  			   gna_xmap_clo_means$rho, 
				               gna_surr_upper_CI$rho, 
				               gna_surr_lower_CI$rho)
names(gna_xmap_clo_res) <- c("lib_size", "GnaXmapClo_actual", "GnaXmapClo_upper", "GnaXmapClo_lower")


## -------------------------------------------------------------------------------
gna_xmap_clo_plot <- ggplot(gna_xmap_clo_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = GnaXmapClo_actual, colour = "GnaXmapClo_actual"))
gna_xmap_clo_plot <- gna_xmap_clo_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = GnaXmapClo_lower, ymax = GnaXmapClo_upper, fill = "GnaXmapClo_95%CI"))
gna_xmap_clo_plot <- gna_xmap_clo_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
gna_xmap_clo_plot <- gna_xmap_clo_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
gna_xmap_clo_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
adf.test(d$clo)


## -------------------------------------------------------------------------------
clo_optembed <- BestEmbed(data = d$clo, E = 1:10)
clo_xmap_gna <- DoCCM(data = d, E = clo_optembed,
					  lib_column = "clo",
					  target_column = "gna",
					  lib_sizes = seq(10, NROW(d$clo), by = 1),
					  num_samples = 100)
clo_xmap_gna_means <- ccm_means(clo_xmap_gna)


## -------------------------------------------------------------------------------
clo_surr <- make_surrogate_ebisuzaki(d$clo)
clo_surr_out <- lapply(seq_len(NCOL(clo_surr)), 
                           function(i) {
                             surr_ts <- clo_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$gna),
                             					  E = clo_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$clo), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
clo_surr_means <- do.call(rbind, clo_surr_out)


## -------------------------------------------------------------------------------
clo_surr_upper_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.975)
clo_surr_lower_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
clo_xmap_gna_res <- data.frame(clo_xmap_gna_means$lib_size, 
				  			   clo_xmap_gna_means$rho, 
				               clo_surr_upper_CI$rho, 
				               clo_surr_lower_CI$rho)
names(clo_xmap_gna_res) <- c("lib_size", "CloXmapGna_actual", "CloXmapGna_upper", "CloXmapGna_lower")


## -------------------------------------------------------------------------------
clo_xmap_gna_plot <- ggplot(clo_xmap_gna_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = CloXmapGna_actual, colour = "CloXmapGna_actual"))
clo_xmap_gna_plot <- clo_xmap_gna_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = CloXmapGna_lower, ymax = CloXmapGna_upper, fill = "CloXmapGna_95%CI"))
clo_xmap_gna_plot <- clo_xmap_gna_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
clo_xmap_gna_plot <- clo_xmap_gna_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
clo_xmap_gna_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
gna_xmap_clo_res <- gna_xmap_clo_res[-40, ]
gna_xmap_clo_res <- gna_xmap_clo_res[-39, ]
gna_xmap_clo_res <- gna_xmap_clo_res[-38, ]
gna_xmap_clo_res <- gna_xmap_clo_res[-37, ]
gna_xmap_clo_res <- gna_xmap_clo_res[-36, ]
gna_xmap_clo_res <- gna_xmap_clo_res[-35, ]
raw_res <- data.frame(gna_xmap_clo_res$lib_size, 
				      gna_xmap_clo_res$GnaXmapClo_actual,
				      gna_xmap_clo_res$GnaXmapClo_upper,
				      gna_xmap_clo_res$GnaXmapClo_lower,
				      clo_xmap_gna_res$CloXmapGna_actual,
				      clo_xmap_gna_res$CloXmapGna_upper,
				      clo_xmap_gna_res$CloXmapGna_lower)
names(raw_res) <- c("lib_size", 
					"GnaXmapClo_actual", "GnaXmapClo_upper", "GnaXmapClo_lower",
			        "CloXmapGna_actual", "CloXmapGna_upper", "CloXmapGna_lower")
raw_res <- mutate(raw_res, GnaXmapClo_lower = if_else(condition = GnaXmapClo_lower < 0, true = 0, false = GnaXmapClo_lower))
raw_res <- mutate(raw_res, CloXmapGna_lower = if_else(condition = CloXmapGna_lower < 0, true = 0, false = CloXmapGna_lower))


## -------------------------------------------------------------------------------
raw_plot <- ggplot(raw_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = GnaXmapClo_actual, colour = "GnaXmapClo_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = CloXmapGna_actual, colour = "CloXmapGna_actual"))
raw_plot <- raw_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = GnaXmapClo_lower, ymax = GnaXmapClo_upper, fill = "GnaXmapClo_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = CloXmapGna_lower, ymax = CloXmapGna_upper, fill = "CloXmapGna_95%CI"))
raw_plot <- raw_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
raw_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
d <- read.csv("Fig5.csv")
d$c1Clo <- c(log(d$gnavus / d$Clostridium))
d$c1Gna <- c(log(d$Clostridium / d$gnavus))
d$clo <- c(-1 / (1 + d$c1Clo))
d$gna <- c(-1 / (1 + d$c1Gna))


## -------------------------------------------------------------------------------
adf.test(d$gna)


## -------------------------------------------------------------------------------
gna_optembed <- BestEmbed(data = d$gna, E = 1:10)
gna_xmap_clo <- DoCCM(data = d, E = gna_optembed,
					  lib_column = "gna",
					  target_column = "clo",
					  lib_sizes = seq(10, NROW(d$gna), by = 1),
					  num_samples = 100)
gna_xmap_clo_means <- ccm_means(gna_xmap_clo)


## -------------------------------------------------------------------------------
gna_surr <- make_surrogate_ebisuzaki(d$gna)
gna_surr_out <- lapply(seq_len(NCOL(gna_surr)), 
                           function(i) {
                             surr_ts <- gna_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$clo),
                             					  E = gna_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$gna), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
gna_surr_means <- do.call(rbind, gna_surr_out)


## -------------------------------------------------------------------------------
gna_surr_upper_CI <- ccm_means(gna_surr_means, FUN = quantile, probs = 0.975)
gna_surr_lower_CI <- ccm_means(gna_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
gna_xmap_clo_res <- data.frame(gna_xmap_clo_means$lib_size, 
				  			   gna_xmap_clo_means$rho, 
				               gna_surr_upper_CI$rho, 
				               gna_surr_lower_CI$rho)
names(gna_xmap_clo_res) <- c("lib_size", "GnaXmapClo_actual", "GnaXmapClo_upper", "GnaXmapClo_lower")


## -------------------------------------------------------------------------------
gna_xmap_clo_plot <- ggplot(gna_xmap_clo_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = GnaXmapClo_actual, colour = "GnaXmapClo_actual"))
gna_xmap_clo_plot <- gna_xmap_clo_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = GnaXmapClo_lower, ymax = GnaXmapClo_upper, fill = "GnaXmapClo_95%CI"))
gna_xmap_clo_plot <- gna_xmap_clo_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
gna_xmap_clo_plot <- gna_xmap_clo_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
gna_xmap_clo_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
adf.test(d$clo)


## -------------------------------------------------------------------------------
clo_optembed <- BestEmbed(data = d$clo, E = 1:10)
clo_xmap_gna <- DoCCM(data = d, E = clo_optembed,
					  lib_column = "clo",
					  target_column = "gna",
					  lib_sizes = seq(10, NROW(d$clo), by = 1),
					  num_samples = 100)
clo_xmap_gna_means <- ccm_means(clo_xmap_gna)


## -------------------------------------------------------------------------------
clo_surr <- make_surrogate_ebisuzaki(d$clo)
clo_surr_out <- lapply(seq_len(NCOL(clo_surr)), 
                           function(i) {
                             surr_ts <- clo_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$gna),
                             					  E = clo_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$clo), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
clo_surr_means <- do.call(rbind, clo_surr_out)


## -------------------------------------------------------------------------------
clo_surr_upper_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.975)
clo_surr_lower_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
clo_xmap_gna_res <- data.frame(clo_xmap_gna_means$lib_size, 
				  			   clo_xmap_gna_means$rho, 
				               clo_surr_upper_CI$rho, 
				               clo_surr_lower_CI$rho)
names(clo_xmap_gna_res) <- c("lib_size", "CloXmapGna_actual", "CloXmapGna_upper", "CloXmapGna_lower")


## -------------------------------------------------------------------------------
clo_xmap_gna_plot <- ggplot(clo_xmap_gna_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = CloXmapGna_actual, colour = "CloXmapGna_actual"))
clo_xmap_gna_plot <- clo_xmap_gna_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = CloXmapGna_lower, ymax = CloXmapGna_upper, fill = "CloXmapGna_95%CI"))
clo_xmap_gna_plot <- clo_xmap_gna_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
clo_xmap_gna_plot <- clo_xmap_gna_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
clo_xmap_gna_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
clo_xmap_gna_res <- clo_xmap_gna_res[-41, ]
clo_xmap_gna_res <- clo_xmap_gna_res[-40, ]
clo_xmap_gna_res <- clo_xmap_gna_res[-39, ]
clo_xmap_gna_res <- clo_xmap_gna_res[-38, ]
clo_xmap_gna_res <- clo_xmap_gna_res[-37, ]
rlr_res <- data.frame(gna_xmap_clo_res$lib_size, 
				      gna_xmap_clo_res$GnaXmapClo_actual,
				      gna_xmap_clo_res$GnaXmapClo_upper,
				      gna_xmap_clo_res$GnaXmapClo_lower,
				      clo_xmap_gna_res$CloXmapGna_actual,
				      clo_xmap_gna_res$CloXmapGna_upper,
				      clo_xmap_gna_res$CloXmapGna_lower)
names(rlr_res) <- c("lib_size", 
					"GnaXmapClo_actual", "GnaXmapClo_upper", "GnaXmapClo_lower",
			        "CloXmapGna_actual", "CloXmapGna_upper", "CloXmapGna_lower")
rlr_res <- mutate(rlr_res, GnaXmapClo_lower = if_else(condition = GnaXmapClo_lower < 0, true = 0, false = GnaXmapClo_lower))
rlr_res <- mutate(rlr_res, CloXmapGna_lower = if_else(condition = CloXmapGna_lower < 0, true = 0, false = CloXmapGna_lower))


## -------------------------------------------------------------------------------
rlr_plot <- ggplot(rlr_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = GnaXmapClo_actual, colour = "GnaXmapClo_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = CloXmapGna_actual, colour = "CloXmapGna_actual"))
rlr_plot <- rlr_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = GnaXmapClo_lower, ymax = GnaXmapClo_upper, fill = "GnaXmapClo_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = CloXmapGna_lower, ymax = CloXmapGna_upper, fill = "CloXmapGna_95%CI"))
rlr_plot <- rlr_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
rlr_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
d <- read.csv("Fig5.csv")
d$clo <- d$Clostridium
d$bac <- d$bacteroidetes


## -------------------------------------------------------------------------------
adf.test(d$bac)


## -------------------------------------------------------------------------------
bac_optembed <- BestEmbed(data = d$bac, E = 1:10)
bac_xmap_clo <- DoCCM(data = d, E = bac_optembed,
					  lib_column = "bac",
					  target_column = "clo",
					  lib_sizes = seq(10, NROW(d$bac), by = 1),
					  num_samples = 100)
bac_xmap_clo_means <- ccm_means(bac_xmap_clo)


## -------------------------------------------------------------------------------
bac_surr <- make_surrogate_ebisuzaki(d$bac)
bac_surr_out <- lapply(seq_len(NCOL(bac_surr)), 
                           function(i) {
                             surr_ts <- bac_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$clo),
                             					  E = bac_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$bac), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
bac_surr_means <- do.call(rbind, bac_surr_out)


## -------------------------------------------------------------------------------
bac_surr_upper_CI <- ccm_means(bac_surr_means, FUN = quantile, probs = 0.975)
bac_surr_lower_CI <- ccm_means(bac_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
bac_xmap_clo_res <- data.frame(bac_xmap_clo_means$lib_size, 
				  			   bac_xmap_clo_means$rho, 
				               bac_surr_upper_CI$rho, 
				               bac_surr_lower_CI$rho)
names(bac_xmap_clo_res) <- c("lib_size", "BacXmapClo_actual", "BacXmapClo_upper", "BacXmapClo_lower")


## -------------------------------------------------------------------------------
bac_xmap_clo_plot <- ggplot(bac_xmap_clo_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = BacXmapClo_actual, colour = "BacXmapClo_actual"))
bac_xmap_clo_plot <- bac_xmap_clo_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = BacXmapClo_lower, ymax = BacXmapClo_upper, fill = "BacXmapClo_95%CI"))
bac_xmap_clo_plot <- bac_xmap_clo_plot + scale_colour_manual("", values = "blue")  + scale_fill_manual("", values = "grey12")
bac_xmap_clo_plot <- bac_xmap_clo_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
bac_xmap_clo_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
adf.test(d$clo)


## -------------------------------------------------------------------------------
clo_optembed <- BestEmbed(data = d$clo, E = 1:10)
clo_xmap_bac <- DoCCM(data = d, E = clo_optembed,
					  lib_column = "clo",
					  target_column = "bac",
					  lib_sizes = seq(10, NROW(d$clo), by = 1),
					  num_samples = 100)
clo_xmap_bac_means <- ccm_means(clo_xmap_bac)


## -------------------------------------------------------------------------------
clo_surr <- make_surrogate_ebisuzaki(d$clo)
clo_surr_out <- lapply(seq_len(NCOL(clo_surr)), 
                           function(i) {
                             surr_ts <- clo_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$bac),
                             					  E = clo_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$clo), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
clo_surr_means <- do.call(rbind, clo_surr_out)


## -------------------------------------------------------------------------------
clo_surr_upper_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.975)
clo_surr_lower_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
clo_xmap_bac_res <- data.frame(clo_xmap_bac_means$lib_size, 
				  			   clo_xmap_bac_means$rho, 
				               clo_surr_upper_CI$rho, 
				               clo_surr_lower_CI$rho)
names(clo_xmap_bac_res) <- c("lib_size", "CloXmapBac_actual", "CloXmapBac_upper", "CloXmapBac_lower")


## -------------------------------------------------------------------------------
clo_xmap_bac_plot <- ggplot(clo_xmap_bac_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = CloXmapBac_actual, colour = "CloXmapBac_actual"))
clo_xmap_bac_plot <- clo_xmap_bac_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = CloXmapBac_lower, ymax = CloXmapBac_upper, fill = "CloXmapBac_95%CI")) + scale_fill_manual("", values = "grey12")
clo_xmap_bac_plot <- clo_xmap_bac_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
clo_xmap_bac_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
bac_xmap_clo_res <- bac_xmap_clo_res[-37, ]
bac_xmap_clo_res <- bac_xmap_clo_res[-36, ]
bac_xmap_clo_res <- bac_xmap_clo_res[-35, ]
raw_res <- data.frame(bac_xmap_clo_res$lib_size, 
				      bac_xmap_clo_res$BacXmapClo_actual,
				      bac_xmap_clo_res$BacXmapClo_upper,
				      bac_xmap_clo_res$BacXmapClo_lower,
				      clo_xmap_bac_res$CloXmapBac_actual,
				      clo_xmap_bac_res$CloXmapBac_upper,
				      clo_xmap_bac_res$CloXmapBac_lower)
names(raw_res) <- c("lib_size", 
					"BacXmapClo_actual", "BacXmapClo_upper", "BacXmapClo_lower",
			        "CloXmapBac_actual", "CloXmapBac_upper", "CloXmapBac_lower")
raw_res <- mutate(raw_res, BacXmapClo_lower = if_else(condition = BacXmapClo_lower < 0, true = 0, false = BacXmapClo_lower))
raw_res <- mutate(raw_res, CloXmapBac_lower = if_else(condition = CloXmapBac_lower < 0, true = 0, false = CloXmapBac_lower))


## -------------------------------------------------------------------------------
raw_plot <- ggplot(raw_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = BacXmapClo_actual, colour = "BacXmapClo_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = CloXmapBac_actual, colour = "CloXmapBac_actual"))
raw_plot <- raw_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = BacXmapClo_lower, ymax = BacXmapClo_upper, fill = "BacXmapClo_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = CloXmapBac_lower, ymax = CloXmapBac_upper, fill = "CloXmapBac_95%CI"))
raw_plot <- raw_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
raw_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
d <- read.csv("Fig5.csv")
d$c1Clo <- c(log(d$bacteroidetes / d$Clostridium))
d$c1Bac <- c(log(d$Clostridium / d$bacteroidetes))
d$clo <- c(-1 / (1 + d$c1Clo))
d$bac <- c(-1 / (1 + d$c1Bac))


## -------------------------------------------------------------------------------
adf.test(d$bac)


## -------------------------------------------------------------------------------
bac_optembed <- BestEmbed(data = d$bac, E = 1:10)
bac_xmap_clo <- DoCCM(data = d, E = bac_optembed,
					  lib_column = "bac",
					  target_column = "clo",
					  lib_sizes = seq(10, NROW(d$bac), by = 1),
					  num_samples = 100)
bac_xmap_clo_means <- ccm_means(bac_xmap_clo)


## -------------------------------------------------------------------------------
bac_surr <- make_surrogate_ebisuzaki(d$bac)
bac_surr_out <- lapply(seq_len(NCOL(bac_surr)), 
                           function(i) {
                             surr_ts <- bac_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$clo),
                             					  E = bac_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$bac), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
bac_surr_means <- do.call(rbind, bac_surr_out)


## -------------------------------------------------------------------------------
bac_surr_upper_CI <- ccm_means(bac_surr_means, FUN = quantile, probs = 0.975)
bac_surr_lower_CI <- ccm_means(bac_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
bac_xmap_clo_res <- data.frame(bac_xmap_clo_means$lib_size, 
				  			   bac_xmap_clo_means$rho, 
				               bac_surr_upper_CI$rho, 
				               bac_surr_lower_CI$rho)
names(bac_xmap_clo_res) <- c("lib_size", "BacXmapClo_actual", "BacXmapClo_upper", "BacXmapClo_lower")


## -------------------------------------------------------------------------------
bac_xmap_clo_plot <- ggplot(bac_xmap_clo_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = BacXmapClo_actual, colour = "BacXmapClo_actual"))
bac_xmap_clo_plot <- bac_xmap_clo_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = BacXmapClo_lower, ymax = BacXmapClo_upper, fill = "BacXmapClo_95%CI"))
bac_xmap_clo_plot <- bac_xmap_clo_plot + scale_colour_manual("", values = "blue")  + scale_fill_manual("", values = "grey12")
bac_xmap_clo_plot <- bac_xmap_clo_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
bac_xmap_clo_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
adf.test(d$clo)


## -------------------------------------------------------------------------------
clo_optembed <- BestEmbed(data = d$clo, E = 1:10)
clo_xmap_bac <- DoCCM(data = d, E = clo_optembed,
					  lib_column = "clo",
					  target_column = "bac",
					  lib_sizes = seq(10, NROW(d$clo), by = 1),
					  num_samples = 100)
clo_xmap_bac_means <- ccm_means(clo_xmap_bac)


## -------------------------------------------------------------------------------
clo_surr <- make_surrogate_ebisuzaki(d$clo)
clo_surr_out <- lapply(seq_len(NCOL(clo_surr)), 
                           function(i) {
                             surr_ts <- clo_surr[, i]
                             xmap_out <- DoCCM(data = cbind(surr_ts, d$bac),
                             					  E = clo_optembed, 
                                                  lib_column = 1,
                                                  target_column = 2,
                                                  lib_sizes = seq(10, NROW(d$clo), by = 1),
                                                  num_samples = 100)
                             xmap_means <- ccm_means(xmap_out)
                           }
                    )
clo_surr_means <- do.call(rbind, clo_surr_out)


## -------------------------------------------------------------------------------
clo_surr_upper_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.975)
clo_surr_lower_CI <- ccm_means(clo_surr_means, FUN = quantile, probs = 0.025)


## -------------------------------------------------------------------------------
clo_xmap_bac_res <- data.frame(clo_xmap_bac_means$lib_size, 
				  			   clo_xmap_bac_means$rho, 
				               clo_surr_upper_CI$rho, 
				               clo_surr_lower_CI$rho)
names(clo_xmap_bac_res) <- c("lib_size", "CloXmapBac_actual", "CloXmapBac_upper", "CloXmapBac_lower")


## -------------------------------------------------------------------------------
clo_xmap_bac_plot <- ggplot(clo_xmap_bac_res) + geom_line(size = 1, 
														  alpha = 1.5, 
														  aes(x = lib_size, y = CloXmapBac_actual, colour = "CloXmapBac_actual"))
clo_xmap_bac_plot <- clo_xmap_bac_plot + geom_ribbon(alpha = 0.5, 
													 aes(x = lib_size, ymin = CloXmapBac_lower, ymax = CloXmapBac_upper, fill = "CloXmapBac_95%CI"))
clo_xmap_bac_plot <- clo_xmap_bac_plot + scale_colour_manual("", values = "blue") + scale_fill_manual("", values = "grey12")
clo_xmap_bac_plot <- clo_xmap_bac_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
clo_xmap_bac_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
clo_xmap_bac_res <- clo_xmap_bac_res[-41, ]
rlr_res <- data.frame(bac_xmap_clo_res$lib_size, 
				      bac_xmap_clo_res$BacXmapClo_actual,
				      bac_xmap_clo_res$BacXmapClo_upper,
				      bac_xmap_clo_res$BacXmapClo_lower,
				      clo_xmap_bac_res$CloXmapBac_actual,
				      clo_xmap_bac_res$CloXmapBac_upper,
				      clo_xmap_bac_res$CloXmapBac_lower)
names(rlr_res) <- c("lib_size", 
					"BacXmapClo_actual", "BacXmapClo_upper", "BacXmapClo_lower",
			        "CloXmapBac_actual", "CloXmapBac_upper", "CloXmapBac_lower")
rlr_res <- mutate(rlr_res, BacXmapClo_lower = if_else(condition = BacXmapClo_lower < 0, true = 0, false = BacXmapClo_lower))
rlr_res <- mutate(rlr_res, CloXmapBac_lower = if_else(condition = CloXmapBac_lower < 0, true = 0, false = CloXmapBac_lower))


## -------------------------------------------------------------------------------
rlr_plot <- ggplot(rlr_res) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = BacXmapClo_actual, colour = "BacXmapClo_actual")) + geom_line(size = 1, alpha = 1.5, aes(x = lib_size, y = CloXmapBac_actual, colour = "CloXmapBac_actual")) 
rlr_plot <- rlr_plot + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = BacXmapClo_lower, ymax = BacXmapClo_upper, fill = "BacXmapClo_95%CI")) + geom_ribbon(alpha = 0.5, aes(x = lib_size, ymin = CloXmapBac_lower, ymax = CloXmapBac_upper, fill = "CloXmapBac_95%CI"))
rlr_plot <- rlr_plot + labs(x = "Library size", y = expression(paste("CMS"(rho))))
rlr_plot
#fig <- ggplotly()
#fig


## -------------------------------------------------------------------------------
sessionInfo()

