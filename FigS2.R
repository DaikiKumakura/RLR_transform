## -------------------------------------------------------------------------------
library(tseries)
library(rEDM)
library(DT)
library(tidyverse)


## -------------------------------------------------------------------------------
data(paramecium_didinium)
df <- paramecium_didinium
DT::datatable(df)


## -------------------------------------------------------------------------------
BestEmbed <- function(data = data, E = 1:10) {
    simp <- simplex(time_series = data, E = E, silent = T)
    opt_embed <- simp$E[which.min(simp$rmse)]
    return(opt_embed)
}

DoCCM <- function(data = data, E = 2,
                  lib_column = "y",
                  target_column = "x",
                  lib_sizes = seq(10, 72, by = 1),
                  num_samples = 100,
                  tp = tp)
{
    ccm(data, E = E,
        lib_column = lib_column,
        target_column = target_column,
        lib_sizes = lib_sizes,
        num_samples = num_samples,
        tp = tp,
        random_libs = TRUE, replace = FALSE, silent = TRUE)
}


## -------------------------------------------------------------------------------
para_adf <- adf.test(df$paramecium)
para_adf
didi_adf <- adf.test(df$didinium)
didi_adf


## -------------------------------------------------------------------------------
vars <- names(df)[c(2, 3)]
pars <- expand.grid(lib_column = vars, target_column = vars, tp = -10:10)
pars <- pars[pars$lib_column !=pars$target_column, ]
#pars$E_didi <- BestEmbed(df$didinium)
#pars$E_para <- BestEmbed(df$paramecium)
pars$E_didi <- 3
pars$E_para <- 3

out_didi <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_didi[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_didi$direction <- paste(out_didi$lib_column, "xmap to\n", out_didi$target_column)

out_para <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_para[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_para$direction <- paste(out_para$lib_column, "xmap to\n", out_para$target_column)
out_abs <- out_para


## -------------------------------------------------------------------------------
out_abs %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
  ylim(0, 1.0) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "Raw; E = 3")


## -------------------------------------------------------------------------------
out_abs %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
  ylim(0.4, 0.9) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "Raw; E = 3")


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
df <- paramecium_didinium_comp


## -------------------------------------------------------------------------------
para_adf <- adf.test(df$paramecium)
para_adf
didi_adf <- adf.test(df$didinium)
didi_adf


## -------------------------------------------------------------------------------
vars <- names(df)[c(2, 3)]
pars <- expand.grid(lib_column = vars, target_column = vars, tp = -10:10)
pars <- pars[pars$lib_column !=pars$target_column, ]
#pars$E_didi <- BestEmbed(df$didinium)
#pars$E_para <- BestEmbed(df$paramecium)
pars$E_didi <- 3
pars$E_para <- 3

out_didi <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_didi[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_didi$direction <- paste(out_didi$lib_column, "xmap to\n", out_didi$target_column)

out_para <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_para[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_para$direction <- paste(out_para$lib_column, "xmap to\n", out_para$target_column)
out_com <- out_para


## -------------------------------------------------------------------------------
out_com %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "Com; E = 3")


## -------------------------------------------------------------------------------
out_com %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
  ylim(0.5, 1) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "Com; E = 3")


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

df <- paramecium_didinium_rlr


## -------------------------------------------------------------------------------
para_adf <- adf.test(df$paramecium)
para_adf
didi_adf <- adf.test(df$didinium)
didi_adf


## -------------------------------------------------------------------------------
vars <- names(df)[c(2, 3)]
pars <- expand.grid(lib_column = vars, target_column = vars, tp = -10:10)
pars <- pars[pars$lib_column !=pars$target_column, ]
#pars$E_didi <- BestEmbed(df$didinium)
#pars$E_para <- BestEmbed(df$paramecium)
pars$E_didi <- 3
pars$E_para <- 3

out_didi <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_didi[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_didi$direction <- paste(out_didi$lib_column, "xmap to\n", out_didi$target_column)

out_para <- do.call(rbind,
                        lapply(seq_len(NROW(pars)),
                        function(i) {
                          ccm(df,
                              E = pars$E_para[i],
                              lib_sizes = NROW(df),
                              random_libs = FALSE,
                              lib_column = pars$lib_column[i],
                              target_column = pars$target_column[i],
                              tp = pars$tp[i],
                              silent = TRUE)
}))
# Add an additional column
out_para$direction <- paste(out_para$lib_column, "xmap to\n", out_para$target_column)


## -------------------------------------------------------------------------------
out_rlr <- mutate(out_para, rho = if_else(condition = rho < 0, true = 0, false = rho))


## -------------------------------------------------------------------------------
out_rlr %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
  ylim(0, 1) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "RLR; E = 3")


## -------------------------------------------------------------------------------
out_rlr %>%
    ggplot(aes(x = tp, y = rho, color = lib_column)) +
    geom_line(size = 1, alpha = 1) +
  geom_point(size = 2, alpha = 1, color = "black") +
    theme_bw() +
  xlim(-5, 5) +
  ylim(0, 0.1) +
    scale_color_hue(name = "Direction", labels = c("didinium" = "paramecium -> didinium", "paramecium" = "didinium -> paramecium")) +
  labs(x = "Time delay", y = expression(paste("CMS")), title = "RLR; E = 3")


## -------------------------------------------------------------------------------
sessionInfo()

