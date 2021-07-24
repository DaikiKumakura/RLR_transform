## -------------------------------------------------------------------------------
library(rEDM)
library(tidyverse)


## -------------------------------------------------------------------------------
BestEmbed <- function(data = data, E = 1:50) {
    simp <- simplex(time_series = data, E = E, silent = T)
    opt_embed <- simp$E[which.min(simp$rmse)]
    return(opt_embed)
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
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list_raw <- data.frame(Byx, Ry, RhoYy)


## -------------------------------------------------------------------------------
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
		e.x <- BestEmbed(d$x)
		# CCM
		ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
		raw.x.model <- ccm.raw.x$model_output[[1]]
		cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
		list_raw <- rbind(list_raw, c(k, j, cor_x))
	}
}
list_raw <- list_raw[-1, ]
list_raw <- mutate(list_raw, RhoYy = if_else(condition = RhoYy < 0, true = 0, false = RhoYy))
write.csv(list_raw, "rho_raw.csv")


## -------------------------------------------------------------------------------
list_raw %>% 
  ggplot(aes(x = Byx, y = Ry, fill = RhoYy)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("darkcyan", "grey90","darkorenge"),limits=c(0,5)) +
  theme_bw() +
  labs(x = expression(paste(beta[21])), y = expression(paste({r[2]}))) +
  theme(text=element_text(size = 20, family = "Arial"),
        legend.position = "none"
        ) 


## -------------------------------------------------------------------------------
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list_com <- data.frame(Byx, Ry, RhoYy)


## -------------------------------------------------------------------------------
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
		d <- raw2comp(d)
		e.x <- BestEmbed(d$x)
		# CCM
		ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
		raw.x.model <- ccm.raw.x$model_output[[1]]
		cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
		list_com <- rbind(list_com, c(k, j, cor_x))
	}
}
list_com <- list[-1, ]
list_com <- mutate(list_com, RhoYy = if_else(condition = RhoYy < 0, true = 0, false = RhoYy))
write.csv(list_com, "rho_com.csv")


## -------------------------------------------------------------------------------
list_com %>% 
  ggplot(aes(x = Byx, y = Ry, fill = RhoYy)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("darkcyan", "grey90","darkorenge"),limits=c(0,5)) +
  theme_bw() +
  labs(x = expression(paste(beta[21])), y = expression(paste({r[2]}))) +
  theme(text=element_text(size = 20, family = "Arial"),
        legend.position = "none"
        ) 


## -------------------------------------------------------------------------------
Byx <- c(0)
Ry <- c(0)
RhoYy <- c(0)
list_rlr <- data.frame(Byx, Ry, RhoYy)


## -------------------------------------------------------------------------------
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
		d <- raw2comp(d)
		d <- comp2rlr(d)
		e.x <- BestEmbed(d$x)
		# CCM
		ccm.raw.x <- ccm(cbind(d$x, d$y), E = e.x, lib_sizes = seq(10, NROW(d$x)), RNGseed = 123, silent = TRUE, stats_only = FALSE)
		raw.x.model <- ccm.raw.x$model_output[[1]]
		cor_x <- cor(raw.x.model$obs, raw.x.model$pred, use = "complete.obs")
		list_rlr <- rbind(list_rlr, c(k, j, cor_x))
	}
}
list_rlr <- list_rlr[-1, ]
list_rlr <- mutate(list_rlr, RhoYy = if_else(condition = RhoYy < 0, true = 0, false = RhoYy))
write.csv(list_rlr, "rho_rlr.csv")


## -------------------------------------------------------------------------------
list_rlr %>% 
  ggplot(aes(x = Byx, y = Ry, fill = RhoYy)) +
  geom_tile() +
  scale_fill_gradientn(colours=c("darkcyan", "grey90","darkorenge"),limits=c(0,5)) +
  theme_bw() +
  labs(x = expression(paste(beta[21])), y = expression(paste({r[2]}))) +
  theme(text=element_text(size = 20, family = "Arial"),
        legend.position = "none"
        ) 

