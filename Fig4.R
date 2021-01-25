##### Rho_Raw #####
##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:10) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Do CCM (For 2 species)
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(2, NROW(x)), num_samples = 50000, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(2, NROW(y)), num_samples = 50000, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 9)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rho)
	# rhoxx is X->Y, and rhoyy is Y->X
	names(dataout.rho) <- c("lib_size", "rhoxx", "rhoyy")
	return(dataout.rho)	
}


### Case.1, Clostridium Vs gnavus
d <- read.csv("Fig4.csv")

d$clo <- d$Clostridium
d$gna <- d$gnavus

c1 <- doccm(d$clo, d$gna)
names(c1) <- c("lib_size", "rhoCC", "rhoGG")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = expression(rho)) 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Ruminococcus_gnavus", rhoGG = "Ruminococcus_gnavus -> Clostridium"))

write.csv(c1, "CloVsgna.csv")

pdf("CloVsgna.pdf", height=10, width=10)
p1
dev.off()


### Case.2, Clostridium Vs bacteroidetes
d <- read.csv("Fig4.csv")

d$clo <- d$Clostridium
d$bif <- d$bacteroidetes

c1 <- doccm(d$clo, d$bac)
names(c1) <- c("lib_size", "rhoCC", "rhoBB")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = expression(rho)) 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Bacteroidetes", rhoBB = "Bacteroidetes -> Clostridium"))

write.csv(c1, "CloVsbac.csv")

pdf("CloVsbac.pdf", height=10, width=10)
p1
dev.off()




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

##### Rho_RLR-Transform #####
##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:10) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Do CCM (For 2 species)
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(2, NROW(x)), num_samples = 50000, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(2, NROW(y)), num_samples = 50000, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 9)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rho)
	# rhoxx is X->Y, and rhoyy is Y->X
	names(dataout.rho) <- c("lib_size", "rhoxx", "rhoyy")
	return(dataout.rho)	
}


### Case.1, Clostridium Vs gnavus
d <- read.csv("Fig4.csv")

d$c1Clo <- c(log(d$gnavus / d$Clostridium))
d$c1Gna <- c(log(d$Clostridium / d$gnavus))
d$clo <- c(-1 / (1 + d$c1Clo))
d$gna <- c(-1 / (1 + d$c1Gna))

c1 <- doccm(d$clo, d$gna)
names(c1) <- c("lib_size", "rhoCC", "rhoGG")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = expression(rho)) 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Ruminococcus_gnavus", rhoGG = "Ruminococcus_gnavus -> Clostridium"))

write.csv(c1, "CloVsgna.csv")

pdf("CloVsgna.pdf", height=10, width=10)
p1
dev.off()


### Case.2, Clostridium Vs bacteroidetes
d <- read.csv("Fig4.csv")

d$c1Clo <- c(log(d$bacteroidetes / d$Clostridium))
d$c1Bac <- c(log(d$Clostridium / d$bacteroidetes))
d$clo <- c(-1 / (1 + d$c1Clo))
d$bac <- c(-1 / (1 + d$c1Bac))

c1 <- doccm(d$clo, d$bac)
names(c1) <- c("lib_size", "rhoCC", "rhoBB")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = expression(rho)) 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Bacteroidetes", rhoBB = "Bacteroidetes -> Clostridium"))

write.csv(c1, "CloVsbac.csv")

pdf("CloVsbac.pdf", height=10, width=10)
p1
dev.off()




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

##### RMSE_Raw #####
##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:10) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Do CCM (For 2 species)
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(2, NROW(x)), num_samples = 50000, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(2, NROW(y)), num_samples = 50000, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 11)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rmse)
	# rhoxx is X->Y, and rhoyy is Y->X
	names(dataout.rho) <- c("lib_size", "rhoxx", "rhoyy")
	return(dataout.rho)	
}


### Case.1, Clostridium Vs gnavus
d <- read.csv("Fig4.csv")

d$clo <- d$Clostridium
d$gna <- d$gnavus

c1 <- doccm(d$clo, d$gna)
names(c1) <- c("lib_size", "rhoCC", "rhoGG")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = "RMSE") 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Ruminococcus_gnavus", rhoGG = "Ruminococcus_gnavus -> Clostridium"))

write.csv(c1, "CloVsgna.csv")

pdf("CloVsgna.pdf", height=10, width=10)
p1
dev.off()


### Case.2, Clostridium Vs bacteroidetes
d <- read.csv("Fig4.csv")

d$clo <- d$Clostridium
d$bif <- d$bacteroidetes

c1 <- doccm(d$clo, d$bac)
names(c1) <- c("lib_size", "rhoCC", "rhoBB")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = "RMSE") 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Bacteroidetes", rhoBB = "Bacteroidetes -> Clostridium"))

write.csv(c1, "CloVsbac.csv")

pdf("CloVsbac.pdf", height=10, width=10)
p1
dev.off()




# --------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

##### RMSE_RLR-Transform #####
##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:10) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Do CCM (For 2 species)
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(2, NROW(x)), num_samples = 50000, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(2, NROW(y)), num_samples = 50000, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 11)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rmse)
	# rhoxx is X->Y, and rhoyy is Y->X
	names(dataout.rho) <- c("lib_size", "rhoxx", "rhoyy")
	return(dataout.rho)	
}


### Case.1, Clostridium Vs gnavus
d <- read.csv("Fig4.csv")

d$c1Clo <- c(log(d$gnavus / d$Clostridium))
d$c1Gna <- c(log(d$Clostridium / d$gnavus))
d$clo <- c(-1 / (1 + d$c1Clo))
d$gna <- c(-1 / (1 + d$c1Gna))

c1 <- doccm(d$clo, d$gna)
names(c1) <- c("lib_size", "rhoCC", "rhoGG")

c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1) 
p1 <- p1 + labs(x = "Library Size", y = "RMSE") 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Ruminococcus_gnavus", rhoGG = "Ruminococcus_gnavus -> Clostridium"))

write.csv(c1, "CloVsgna.csv")

pdf("CloVsgna.pdf", height=10, width=10)
p1
dev.off()


### Case.2, Clostridium Vs bacteroidetes
d <- read.csv("Fig4.csv")

d$c1Clo <- c(log(d$bacteroidetes / d$Clostridium))
d$c1Bac <- c(log(d$Clostridium / d$bacteroidetes))
d$clo <- c(-1 / (1 + d$c1Clo))
d$bac <- c(-1 / (1 + d$c1Bac))

c1 <- doccm(d$clo, d$bac)
names(c1) <- c("lib_size", "rhoCC", "rhoBB")

corr <- cor(d$clo, d$bac)
c1.melt <- melt(c1, id = "lib_size", variable.name = "rho")
c1.melt$rho <- as.factor(c1.melt$rho)
p1 <- ggplot(data = c1.melt, aes(x = lib_size, y = value, color = rho))
p1 <- p1 + geom_line(alpha = 1, size = 1)  
p1 <- p1 + labs(x = "Library Size", y = "RMSE") 
p1 <- p1 + scale_color_hue(name = "Direction", labels = c(rhoCC = "Clostridium -> Bacteroidetes", rhoBB = "Bacteroidetes -> Clostridium"))

write.csv(c1, "CloVsbac.csv")

pdf("CloVsbac.pdf", height=10, width=10)
p1
dev.off()



