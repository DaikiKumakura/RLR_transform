####### 2 species; competition ########


##### Package #####
library(rEDM)
library(deSolve)
library(ggplot2)
library(plyr)
library(reshape2)
library(rgr)


##### Useful Functions #####


### Check best embedding dimension
BestEmbed <- function(time.series, E = 1:50) {
	simplex.res <- simplex(time.series, E = E, silent = T)
	bestE <- simplex.res$E[which.min(simplex.res$rmse)]
	return(bestE)
}


### Visualizes Population Dynamics (For 2 species)
plotdyn <- function(dataout) {
	dataset <- melt(dataout, id = "time", variable.name = "Species")
	pDyn <-  ggplot(data = dataset, aes(x = time, y = value, color = Species)) + geom_line(alpha = 1, 
    	size = 0.5) + labs(x = "Time", y = "Population size") + scale_color_hue(name = "Species", labels = c(x = "X", y = "Y")) 
	pDyn <- pDyn + xlim(50, 150)
	return(pDyn)
}


### Do CCM (For 2 species)
doccm <- function(x, y) {
	e.x <- BestEmbed(x)
	e.y <- BestEmbed(y)
	ccm.raw.x <- ccm(cbind(x, y), E = e.x, lib_sizes = seq(10, NROW(x)), RNGseed = 123, silent = TRUE)
	ccm.raw.y <- ccm(cbind(y, x), E = e.y, lib_sizes = seq(10, NROW(y)), RNGseed = 123, silent = TRUE)
	ccm.tmp.x <- ccm_means(ccm.raw.x)
	ccm.tmp.y <- ccm_means(ccm.raw.y)
	dataout.rho <- ccm.tmp.y[, c(7, 9)]
	dataout.rho$rho2 <- c(ccm.tmp.x$rho)
	# rhoYX is X->Y, and rhoXY is Y->X
	names(dataout.rho) <- c("lib_size", "rhoYX", "rhoXY")
	return(dataout.rho)	
}


### Visualizes Cross-Map (For 2 species)
plotccm <- function(dataout.rho) {
	dataset.rho <- melt(dataout.rho, id = "lib_size", variable.name = "rho")
	dataset.rho$rho <- as.factor(dataset.rho$rho)
	pCCM <- ggplot(data = dataset.rho, aes(x = lib_size, y = value, color = rho))
	pCCM <-	pCCM + geom_line(alpha = 1, size = 1) 
	pCCM <-	pCCM + labs(x = "Library Size", y = expression(rho)) 
	pCCM <-	pCCM + scale_color_hue(name = "Direction", labels = c(rhoYX = "X -> Y", rhoXY = "Y -> X"))
	pCCM <- pCCM + ylim(-1, 1)
	return(pCCM)
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


### Visualizes the Omega Dynamics
plotWdyn <- function(data) {
	set <- melt(data, id = "lib_size", variable.name = "Data")
	set$Data <- as.factor(set$Data)
	pWDy <- ggplot(data = set, aes(x = lib_size, y = value, color = Data))
	pWDy <- pWDy + geom_line(alpha = 1, size = 1) 
	pWDy <- pWDy + labs(x = "Library Size", y = expression(omega))
	return(pWDy) 
}


### Converts the absolute to relative (RLR)
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


##### Creates the Omega list for some data #####
num <- c(0)
data <- c(0)
w <- c(0)
list_w <- data.frame(num, data, w)
names(list_w) <- c("num", "data", "w")


##### Absolute data #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_Abs.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_Abs.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_Abs.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(2, "Absolute", sum_w))


##### Log(Absolute + 1) #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
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


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_LogAbs.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_LogAbs.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_LogAbs.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(3, "LogAbsolute", sum_w))


##### Composition #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- comp_trans_2(d)


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_Comp.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_Comp.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_Comp.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(4, "Composition", sum_w))


##### Log(Composition + 1) #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- comp_trans_2(d)
dataout_logComp <- d[, c(2:3)]
dataout_logComp$x_log <- c(log(dataout_logComp$x + 1))
dataout_logComp$y_log <- c(log(dataout_logComp$y + 1))
d$x_re <- c(dataout_logComp$x_log)
d$y_re <- c(dataout_logComp$y_log)
dataout_logComposition <- d[, c(1, 4:5)]
names(dataout_logComposition) <- c("time", "x", "y")
d <- dataout_logComposition


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_LogComp.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_LogComp.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_LogComp.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(5, "LogComposition", sum_w))


##### clrTransform #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
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


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_clr.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_clr.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_clr.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(6, "clrTransform", sum_w))


##### RelativeConversion (RLR) #####


### Creates the datasets
d0 <- as.data.frame(matrix(rep(NA,3*1000),ncol=3))
colnames(d0) <- c("time","x","y")
d0[,1] <- 1:1000
d0[1,2:3] <- c(0.1, 0.1) # 初期値
for(i in 1:999){
        d0[i+1,2] <- d0[i,2]*(3.8 - 3.80*d0[i,2] - 0.02*d0[i,3])
        d0[i+1,3] <- d0[i,3]*(3.5 - 0.10*d0[i,2] - 3.50*d0[i,3])
}
d <- d0
d <- comp_trans_2(d)
d <- rel_trans_2(d)


### Visualizes the Population Dynamics
pop_dyn <- "PopDyn_Rel.png"
png(pop_dyn, height = 493, width = 835)
plotdyn(d)
dev.off()


### Calculates the CCM, Visulizes the Cross-Map and Write Omega for csv file
out <- doccm(d$x, d$y)
plt_ccm <- "CCM_Rel.png"
png(plt_ccm, height = 493, width = 835)
plotccm(out)
dev.off()
out$w <- c(out$rhoYX - out$rhoXY)
write.csv(out, "Rho_Rel.csv")
mer_w <- out[, 4]
mer_w <- as.data.frame(mer_w)
sum_w <- colSums(mer_w)
list_w <- rbind(list_w, c(7, "RelativeConversion", sum_w))


##### Write the Omega list for some data #####
write.csv(list_w, "omega_list.csv")


##### Visualize the Omega dynamics #####
# the row name for "df" are 2 (lib_size, and w)
df <- read.csv("omega_dynamics.csv")
plotWdyn(df)




