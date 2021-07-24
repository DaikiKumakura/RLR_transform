## -------------------------------------------------------------------------------
library(rEDM)
library(deSolve)
library(tidyverse)
library(tseries)
library(DT)


## -------------------------------------------------------------------------------
BestEmbed <- function(data = data, E = 1:10) {
    simp <- simplex(time_series = data, E = E, silent = T)
    opt_embed <- simp$E[which.min(simp$rmse)]
    return(opt_embed)
}

DoCCM <- function(data = data, E = 3,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(10, 1000, by = 1),
                  num_samples = 100)
{
    ccm(data, E = E,
        lib_column = lib_column,
        target_column = target_column,
        lib_sizes = lib_sizes,
        num_samples = num_samples,
        random_libs = TRUE, replace = FALSE, silent = TRUE)
}

raw <- function(a_y_to_x, a_x_to_y) {
        d0 <- as.data.frame(matrix(rep(NA,4*1000),ncol=4))
        colnames(d0) <- c("time","x","y","z")
        d0[,1] <- 1:1000
        d0[1,2:4] <- c(0.1, 0.1, 0.1) # initial value
        for(i in 1:999){

                d0[i+1,2] <- d0[i,2]*(3.80 - 3.80*d0[i,2] - a_y_to_x*d0[i,3] - 0.02*d0[i,4])
                d0[i+1,3] <- d0[i,3]*(3.50 - a_x_to_y*d0[i,2] - 3.50*d0[i,3] - 0.08*d0[i,4])
                d0[i+1,4] <- d0[i,4]*(3.70 - 0.025*d0[i,2] - 0.02*d0[i,3] - 3.70*d0[i,4])
        }
        d <- d0
        return(d)
}

raw2comp <- function(data) {
        dataout <- data
        dataout_comp <- dataout[, c(2:4)]
        dataout_comp$sum <- c(apply(dataout_comp, 1, sum))
        dataout_comp$x_re <- c(dataout$x / dataout_comp$sum)
        dataout_comp$y_re <- c(dataout$y / dataout_comp$sum)
        dataout_comp$z_re <- c(dataout$z / dataout_comp$sum)
        dataout$x_re <- c(dataout_comp$x_re)
        dataout$y_re <- c(dataout_comp$y_re)
        dataout$z_re <- c(dataout_comp$z_re)
        dataout_composition <- dataout[, c(1, 5:7)]
        names(dataout_composition) <- c("time", "x", "y", "z")
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
strength_x_to_y_raw <- c()
correct_x_to_y_raw <- c()
result_correct_raw <- data.frame(strength_x_to_y_raw, correct_x_to_y_raw)


## -------------------------------------------------------------------------------
XtoY <- c(2, 4, 6)
Stg <- c(0.01, 0.1)
for (i in XtoY) {
	for (j in Stg) {
		for (k in seq(1, 3, by = 1)) {
			# dataset
			alpha_y_to_x <- k * j
			alpha_x_to_y <- k * j * i
			d <- raw(alpha_y_to_x, alpha_x_to_y)

			# CCM on actual time series (x -> y)
			y_optembed <- BestEmbed(data = d$y, E = 1:50)
			y_xmap_x <- DoCCM(data = d, E = y_optembed,
							lib_column = "y",
							target_column = "x",
							lib_sizes = seq(50, NROW(d$y), by = 50),
							num_samples = 100)
			y_xmap_x_means <- ccm_means(y_xmap_x)
			act_x_to_y <- y_xmap_x_means[20, 9]

			# CCM on surrogate time series
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

			# 95% confidence interval
			y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
			surr_upper <- y_surr_upper_CI[20, 9]

			# CCM on actual time series (y -> x)
			x_optembed <- BestEmbed(data = d$x, E = 1:50)
			x_xmap_y <- DoCCM(data = d, E = x_optembed,
							lib_column = "x",
							target_column = "y",
							lib_sizes = seq(50, NROW(d$x), by = 50),
							num_samples = 100)
			x_xmap_y_means <- ccm_means(x_xmap_y)
			act_y_to_x <- x_xmap_y_means[20, 9]

			# evaluate act and surr_upper
			if (act_x_to_y > surr_upper) {
				if (act_x_to_y > act_y_to_x) {
					corr <- 1
				} else {
					corr <- 0
				}
			} else {
				corr <- 0
			}

			# record
			result_correct_raw <- rbind(result_correct_raw, c(i, corr))
		}
	}
}


## -------------------------------------------------------------------------------
result_correct_raw <- result_correct_raw[-1, ]
write.csv(result_correct_raw, "result_correct_raw.csv")
DT::datatable(result_correct_raw)


## -------------------------------------------------------------------------------
strength_x_to_y_com <- c()
correct_x_to_y_com <- c()
result_correct_com <- data.frame(strength_x_to_y_com, correct_x_to_y_com)


## -------------------------------------------------------------------------------
XtoY <- c(2, 4, 6)
Stg <- c(0.01, 0.1)
for (i in XtoY) {
	for (j in Stg) {
		for (k in seq(1, 3, by = 1)) {
			# dataset
			alpha_y_to_x <- k * j
			alpha_x_to_y <- k * j * i
			d <- raw(alpha_y_to_x, alpha_x_to_y)
			d <- raw2comp(d)

			# CCM on actual time series (x -> y)
			y_optembed <- BestEmbed(data = d$y, E = 1:50)
			y_xmap_x <- DoCCM(data = d, E = y_optembed,
							lib_column = "y",
							target_column = "x",
							lib_sizes = seq(50, NROW(d$y), by = 50),
							num_samples = 100)
			y_xmap_x_means <- ccm_means(y_xmap_x)
			act_x_to_y <- y_xmap_x_means[20, 9]

			# CCM on surrogate time series
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

			# 95% confidence interval
			y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
			surr_upper <- y_surr_upper_CI[20, 9]

			# CCM on actual time series (y -> x)
			x_optembed <- BestEmbed(data = d$x, E = 1:50)
			x_xmap_y <- DoCCM(data = d, E = x_optembed,
							lib_column = "x",
							target_column = "y",
							lib_sizes = seq(50, NROW(d$x), by = 50),
							num_samples = 100)
			x_xmap_y_means <- ccm_means(x_xmap_y)
			act_y_to_x <- x_xmap_y_means[20, 9]

			# evaluate act and surr_upper
			if (act_x_to_y > surr_upper) {
				if (act_x_to_y > act_y_to_x) {
					corr <- 1
				} else {
					corr <- 0
				}
			} else {
				corr <- 0
			}

			# record
			result_correct_com <- rbind(result_correct_com, c(i, corr))
		}
	}
}


## -------------------------------------------------------------------------------
result_correct_com <- result_correct_com[-1, ]
write.csv(result_correct_com, "result_correct_com.csv")
DT::datatable(result_correct_raw)


## -------------------------------------------------------------------------------
strength_x_to_y_rlr <- c()
correct_x_to_y_rlr <- c()
result_correct_rlr <- data.frame(strength_x_to_y_rlr, correct_x_to_y_rlr)


## -------------------------------------------------------------------------------
XtoY <- c(2, 4, 6)
Stg <- c(0.01, 0.1)
for (i in XtoY) {
	for (j in Stg) {
		for (k in seq(1, 3, by = 1)) {
			# dataset
			alpha_y_to_x <- k * j
			alpha_x_to_y <- k * j * i
			d <- raw(alpha_y_to_x, alpha_x_to_y)
			d <- raw2comp(d)
			d <- comp2rlr(d)

			# CCM on actual time series (x -> y)
			y_optembed <- BestEmbed(data = d$y, E = 1:50)
			y_xmap_x <- DoCCM(data = d, E = y_optembed,
							lib_column = "y",
							target_column = "x",
							lib_sizes = seq(50, NROW(d$y), by = 50),
							num_samples = 100)
			y_xmap_x_means <- ccm_means(y_xmap_x)
			act_x_to_y <- y_xmap_x_means[20, 9]

			# CCM on surrogate time series
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

			# 95% confidence interval
			y_surr_upper_CI <- ccm_means(y_surr_means, FUN = quantile, probs = 0.975)
			surr_upper <- y_surr_upper_CI[20, 9]

			# CCM on actual time series (y -> x)
			x_optembed <- BestEmbed(data = d$x, E = 1:50)
			x_xmap_y <- DoCCM(data = d, E = x_optembed,
							lib_column = "x",
							target_column = "y",
							lib_sizes = seq(50, NROW(d$x), by = 50),
							num_samples = 100)
			x_xmap_y_means <- ccm_means(x_xmap_y)
			act_y_to_x <- x_xmap_y_means[20, 9]

			# evaluate act and surr_upper
			if (act_x_to_y > surr_upper) {
				if (act_x_to_y > act_y_to_x) {
					corr <- 1
				} else {
					corr <- 0
				}
			} else {
				corr <- 0
			}

			# record
			result_correct_rlr <- rbind(result_correct_rlr, c(i, corr))
		}
	}
}


## -------------------------------------------------------------------------------
result_correct_rlr <- result_correct_rlr[-1, ]
write.csv(result_correct_rlr, "result_correct_rlr.csv")
DT::datatable(result_correct_raw)


## -------------------------------------------------------------------------------
correct_raw <- result_correct_raw$X1
correct_com <- result_correct_com$X0
correct_rlr <- result_correct_rlr$X1
strength <- result_correct_raw$X2
result_integrate <- data.frame(strength, correct_raw, correct_com, correct_rlr)
DT::datatable(result_integrate)
write.csv(result_integrate, "result_correct_integrate.csv")


## -------------------------------------------------------------------------------
p_abs <- sum(result_integrate$correct_raw) / 18
p_com <- sum(result_integrate$correct_com) / 18
p_rlr <- sum(result_integrate$correct_rlr) / 18
p_data <- data.frame(
  Data = c("Absolute", "Composition", "RLR"),
  Percent = c(p_abs, p_com, p_rlr)
)
DT::datatable(p_data)


## -------------------------------------------------------------------------------
library(ggsci)
plt_correct <- p_data %>% 
  ggplot(aes(x = Data, y = Percent, fill = Data)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_npg() +
  theme(axis.text = element_text(angle = 0), text=element_text(size = 18, family = "Arial")) +
  ylim(0, 1) +
  labs(x = "Dataset", y = "Correct")
plt_correct
ggsave(file = "Fig3.eps", plot = plt_correct)


## -------------------------------------------------------------------------------
sessionInfo()

