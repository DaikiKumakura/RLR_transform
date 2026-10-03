# Run from the repository root. Pin the legacy API used by the figure scripts.
stopifnot(as.character(packageVersion('rEDM')) == '0.7.3')
library(rEDM)
for (f in list.files(pattern='^Fig.*[.]R$')) parse(f)
expr <- parse('Fig1.R')
for (e in expr) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) &&
      is.symbol(e[[2]]) && as.character(e[[2]]) %in%
      c('BestEmbed','DoCCM','raw2comp','comp2rlr')) eval(e)
}
# Retain the original DoCCM parameters, with an explicit verification seed.
body(DoCCM)[[2]]$RNGseed <- 123L
d <- data.frame(time=1:500, x=NA_real_, y=NA_real_)
d[1,2:3] <- c(0.1,0.1)
for (i in 1:499) {
  d$x[i+1] <- d$x[i]*(3.8-3.80*d$x[i]-0.02*d$y[i])
  d$y[i+1] <- d$y[i]*(3.5-0.10*d$x[i]-3.50*d$y[i])
}
d <- tail(d,120)
row.names(d) <- NULL
comp <- raw2comp(d)
stopifnot(all(abs(comp$x+comp$y-1)<1e-12))
rlr <- comp2rlr(comp)
stopifnot(all(is.finite(as.matrix(rlr))))
results <- do.call(rbind, lapply(list(raw=d,composition=comp,rlr=rlr), function(x) {
  E <- BestEmbed(x$x, E=1:3)
  r <- ccm_means(DoCCM(x,E=E,lib_sizes=c(20,40,60),num_samples=10))
  stopifnot(nrow(r)==3,all(is.finite(r$rho)))
  r$E <- E
  r
}))
results$representation <- rep(c('raw','composition','rlr'),each=3)
print(results[,c('representation','lib_size','rho','E')])
if (nzchar(Sys.getenv('RLR_SMOKE_OUTPUT'))) {
  write.csv(results,Sys.getenv('RLR_SMOKE_OUTPUT'),row.names=FALSE)
}
cat('Small original-API CCM smoke check passed; full figures are not validated.\n')
