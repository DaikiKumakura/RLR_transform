##### Code #####


### Packages ###
library(rEDM)
library(ggplot2)
library(plyr)
library(reshape2)


### Functions ###
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(3, NROW(x)), RNGseed = 123, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(3, NROW(y)), RNGseed = 123, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 9)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rho)
	# rhoYX is X->Y, and rhoXY is Y->X
	names(dataout.rho) <- c("lib_size", "rhoYX", "rhoXY")
	return(dataout.rho)	
}


plotccm <- function(dataout.rho) {
	dataset.rho <- melt(dataout.rho, id = "lib_size", variable.name = "rho")
	dataset.rho$rho <- as.factor(dataset.rho$rho)
	pCCM <- ggplot(data = dataset.rho, aes(x = lib_size, y = value, color = rho))
	pCCM <-	pCCM + geom_line(alpha = 1, size = 1) 
	pCCM <-	pCCM + labs(x = "Library Size", y = expression(rho)) 
	pCCM <-	pCCM + scale_color_hue(name = "Direction", labels = c(rhoYX = "X -> Y", rhoXY = "Y -> X"))
	return(pCCM)
}


doccm_rmse <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = 3, lib_sizes = seq(3, NROW(x)), RNGseed = 123, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = 3, lib_sizes = seq(3, NROW(x)), RNGseed = 123, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rmse <- ccm.tmp.y[, c(7, 11)]
	dataout.rmse$rmse2 <- c(ccm.tmp.x$rmse)
	# RMSE_YX is X->Y, and RMSE_XY is Y->X
	names(dataout.rmse) <- c("lib_size", "RMSE_YX", "RMSE_XY")
	return(dataout.rmse)	
}


plotccm_rmse <- function(dataout.rho) {
	dataset.rho <- melt(dataout.rho, id = "lib_size", variable.name = "rho")
	dataset.rho$rho <- as.factor(dataset.rho$rho)
	pCCM <- ggplot(data = dataset.rho, aes(x = lib_size, y = value, color = rho))
	pCCM <-	pCCM + geom_line(alpha = 1, size = 1) 
	pCCM <-	pCCM + labs(x = "Library Size", y = "RMSE") 
	pCCM <-	pCCM + scale_color_hue(name = "Direction", labels = c(RMSE_YX = "X -> Y", RMSE_XY = "Y -> X"))
	return(pCCM)
}


BestEmbed <- function(time.series, E = 1:10) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


comp_trans_2 <- function(data) {
        dataout <- data
        dataout_comp <- dataout
        dataout_comp$sum <- c(apply(dataout_comp, 1, sum))
        dataout_comp$x_re <- c(dataout$x / dataout_comp$sum)
        dataout_comp$y_re <- c(dataout$y / dataout_comp$sum)
        dataout$x_re <- c(dataout_comp$x_re)
        dataout$y_re <- c(dataout_comp$y_re)
        dataout_composition <- dataout[, c(3:4)]
        names(dataout_composition) <- c("x", "y")
        dataout <- dataout_composition
        d <- dataout
        return(d)
}



rel_trans_2 <- function(data) {
        d <- data
        d$lx <- c(log((d$y / d$x)))
        d$ly <- c(log((d$x / d$y)))
        d$x <- c(-1 / (1 + d$lx))
        d$y <- c(-1 / (1 + d$ly))
        df <- d[, c(1:2)]
        d <- df
        return(d)
}


### Datasets ###
data(paramecium_didinium)
str(paramecium_didinium)
d <- paramecium_didinium
d <- d[, c(2:3)]
names(d) <- c("x", "y")
# paramecium is "x", and didinium is "y"


### Analysis(rho) ###
# Raw
res_raw <- doccm(d$x, d$y)
plotccm(res_raw)
res_raw_rmse <- doccm_rmse(d$x, d$y)
plotccm_rmse(res_raw_rmse)


# Absolute
d$x <- c((d$x - min(d$x)) / (max(d$x) - min(d$x)))
d$y <- c((d$y - min(d$y)) / (max(d$y) - min(d$y)))
res_abs <- doccm(d$x, d$y)
plotccm(res_abs)
res_abs_rmse <- doccm_rmse(d$x, d$y)
plotccm_rmse(res_abs_rmse)


# Composition
d <- comp_trans_2(d)
res_comp <- doccm(d$x, d$y)
plotccm(res_comp)
res_comp_rmse <- doccm_rmse(d$x, d$y)
plotccm_rmse(res_comp_rmse)


# RLR
d <- rel_trans_2(d)
res_de <- doccm(d$x, d$y)
plotccm(res_de)
res_de_rmse <- doccm_rmse(d$x, d$y)
plotccm_rmse(res_de_rmse)


