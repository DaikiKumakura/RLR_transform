
##### Parameter dependence ######
##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)
library(plotly)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:50) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Converts the absolute to composition
# 2species
comp_trans_2 <- function(data) {
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


### Converts the absolute to relative
# 2species
rel_trans_2 <- function(data) {
        d <- data
        d$lx <- c(log((d$y / d$x)))
        d$ly <- c(log((d$x / d$y)))
        d$x <- c(-1 / (1 + d$lx))
        d$y <- c(-1 / (1 + d$ly))
        df <- d[, c(1:3)]
        d <- df
        return(d)
}



# Absolute
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list <- data.frame(Byx, Ry, RhoYy)
names(list) <- c("Byx", "Ry", "RhoYy")
for (k in seq(0, 0.8, by = 0.005)) {
	for (j in seq(3.5, 3.9, by = 0.005)) {
		d0 <- as.data.frame(matrix(rep(NA,3*1400),ncol=3))
		colnames(d0) <- c("time","x","y")
		d0[,1] <- 1:1400
		d0[1,2:3] <- c(0.1, 0.1) 
		for(i in 1:1399){
        		d0[i+1,2] <- d0[i,2]*((7.3984 - j) - (7.3984 - j)*d0[i,2])
        		d0[i+1,3] <- d0[i,3]*(j - k*d0[i,2] - j*d0[i,3])
		}
		d <- d0[1000:1400, ]
		#d <- comp_trans_2(d)
		#d <- rel_trans_2(d)
		e.x <- BestEmbed(d$x)
		# CCM
		ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
		raw.x.model <- ccm.raw.x$model_output[[1]]
		cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
		list <- rbind(list, c(k, j, cor_x))
	}
}

list <- list[-1, ]

write.csv(list, "rho_abs.csv")




# Composition
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list <- data.frame(Byx, Ry, RhoYy)
names(list) <- c("Byx", "Ry", "RhoYy")
for (k in seq(0, 0.8, by = 0.005)) {
        for (j in seq(3.5, 3.9, by = 0.005)) {
                d0 <- as.data.frame(matrix(rep(NA,3*1400),ncol=3))
                colnames(d0) <- c("time","x","y")
                d0[,1] <- 1:1400
                d0[1,2:3] <- c(0.1, 0.1) 
                for(i in 1:1399){
                        d0[i+1,2] <- d0[i,2]*((7.3984 - j) - (7.3984 - j)*d0[i,2])
                        d0[i+1,3] <- d0[i,3]*(j - k*d0[i,2] - j*d0[i,3])
                }
                d <- d0[1000:1400, ]
                d <- comp_trans_2(d)
                #d <- rel_trans_2(d)
                e.x <- BestEmbed(d$x)
                # CCM
                ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
                raw.x.model <- ccm.raw.x$model_output[[1]]
                cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
                list <- rbind(list, c(k, j, cor_x))
        }
}

list <- list[-1, ]

write.csv(list, "rho_comp.csv")





# RLR
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list <- data.frame(Byx, Ry, RhoYy)
names(list) <- c("Byx", "Ry", "RhoYy")
for (k in seq(0, 0.8, by = 0.005)) {
        for (j in seq(3.5, 3.9, by = 0.005)) {
                d0 <- as.data.frame(matrix(rep(NA,3*1400),ncol=3))
                colnames(d0) <- c("time","x","y")
                d0[,1] <- 1:1400
                d0[1,2:3] <- c(0.1, 0.1) 
                for(i in 1:1399){
                        d0[i+1,2] <- d0[i,2]*((7.3984 - j) - (7.3984 - j)*d0[i,2])
                        d0[i+1,3] <- d0[i,3]*(j - k*d0[i,2] - j*d0[i,3])
                }
                d <- d0[1000:1400, ]
                d <- comp_trans_2(d)
                d <- rel_trans_2(d)
                e.x <- BestEmbed(d$x)
                # CCM
                ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
                raw.x.model <- ccm.raw.x$model_output[[1]]
                cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
                list <- rbind(list, c(k, j, cor_x))
        }
}

list <- list[-1, ]

write.csv(list, "rho_de.csv")


