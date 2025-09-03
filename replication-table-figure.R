# title: "Replication codes for Table 1 and Figure 1"
# author: "Andreï V. Kostyrka (andrei.kostyrka@gmail.com)"
# date: '4th of September 2025'
# Tested with R v. 4.5.1 and pnd v. 0.1.1

rm(list = ls())
tryCatch(setwd("~/Dropbox/HSE/14/pnd/wydss/"), error = function(e) return(NULL))

# For EXACT replication of the results, use the SAME package version
# Warning: this option is aggressive and may overwrite your existing installed pnd
library(devtools)
devtools::install_version("pnd", "0.1.1")

library(pnd)
library(parallel)

changeFont <- function() {
  a <- tryCatch(par(family = "Fira Sans"), error = function(e) return(NULL))
  if (is.null(a)) a <- tryCatch(par(family = "Noto Sans"), error = function(e) return(NULL))
  if (is.null(a)) a <- tryCatch(par(family = "Verdana"), error = function(e) return(NULL))
  if (is.null(a)) a <- tryCatch(par(family = "Gill Sans"), error = function(e) return(NULL))
  if (is.null(a)) a <- tryCatch(par(family = "Geneva"), error = function(e) return(NULL))
}

##########################################
# Figure 1
##########################################
a  <- step.K(sin, 1)
aa <- step.K(function(x) x, 1)
cairo_pdf("figure1.pdf", 6, 3)
changeFont()
par(mar = c(4, 4, 0, 0.5) + .1, mfrow = c(1, 2))
plot(a$iterations$h, a$iterations$est.error[, 1], log = "xy", bty = 'n', pch = 16, cex = 0.7,
     xlab = "Step size", ylab = "Combined error")
legend("topleft", "f(x) = sin(x)", bty = "n")
plot(aa$iterations$h, aa$iterations$est.error[, 1], log = "xy", bty = 'n', pch = 16, cex = 0.7,
     xlab = "Step size", ylab = "Combined error")
legend("top", "f(x) = x", bty = "n")
dev.off()


###########################################
# Table 1
###########################################
n <- 10000
imgh <- 720
imgv <- 320
cores <- detectCores() / 2

linear <- TRUE
k <- 1

# 1:8
for (k in c(1:8)) {
  print(k)
  
  if (k == 1) {
    f  <- sin
    fp <- cos
    fppp <- function(x) -cos(x)
    # xsl <- c(-pi, pi)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-1, 1e+5)
    flab <- "sin"
  } else if (k == 2) {
    f  <- exp
    fp <- exp
    fppp <- exp
    # xsl <- c(-10, 1)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e1)
    flab <- "exp"
  } else if (k == 3) {
    f  <- log
    fp <- function(x) 1/x
    fppp <- function(x) 2/x^3
    # xsl <- c(1e-2, 10)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+5)
    flab <- "log"
  } else if (k == 4) {
    f  <- sqrt
    fp <- function(x) 0.5/sqrt(x)
    fppp <- function(x) 3/8*x^(-5/2)
    # xsl <- c(1e-2, 10)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+5)
    flab <- "sqrt"
  } else if (k == 5) {
    f  <- atan
    fp <- function(x) 1/(x^2 + 1)
    fppp <- function(x) (6*x^2 - 2)/(1+x^2)^3  # Simplify[D[ArcTan[x], {x, 3}]]
    # xsl <- c(-10, 10)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+5)
    flab <- "atan"
  } else if (k == 6) {
    f  <- function(x) pi*x + 2
    fp <- function(x) rep(pi, length(x))
    fppp <- function(x) x*0
    # xsl <- c(-1e2, 1e2)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+5)
    flab <- "poly1"
  } else if (k == 7) {
    f  <- function(x) x^2
    fp <- function(x) 2*x
    fppp <- function(x) x*0
    # xsl <- c(-1e2, 1e2)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+4)
    flab <- "poly2"
  } else if (k == 8) {
    f  <- function(x) sin(x^2 + 10^6 * x)
    fp <- function(x) (1e6 + 2*x) * cos(x^2 + 10^6 * x)
    # D[Sin[x^2 + 10^6*x], {x, 3}]
    fppp <- function(x) -(1e6 + 2*x)^3 * cos(1e6*x + x^2) -  6*(1e6 + 2*x)*sin(1e6*x + x^2)
    # xsl <- c(-10, 10)
    xsl <- c(0.1, 12.5)
    xse <- c(1e-2, 1e+1)
    flab <- "illsin"
  }
  
  set.seed(1)
  xgrid <- sort(runif(n, min = xsl[1], max = xsl[2]))
  xax.vals <- pretty(xsl)
  xl <- xsl
  
  # Analytical minimiser of etrunc + eround under certain assumptions
  fHopt <- function(x) {
    f1 <- f(x)
    f3 <- fppp(x)
    if (abs(f1) < 16 * .Machine$double.eps) f1 <- 16 * .Machine$double.eps  # Safeguarding against zero steps
    if (abs(f3) < .Machine$double.eps) list(par = NA, value = NA)  # Division by zero
    h <- (1.5 * abs(f1 / f3) * .Machine$double.eps)^(1/3)
    v <- (-0.5*f(x-h) + 0.5*f(x+h)) / h
    return(list(par = h, value = v))
  }
  
  fHRoT  <- function(x) {
    ltt <- x < .Machine$double.eps^(1/3)
    h <- .Machine$double.eps^(1/3) * (max(1, abs(x)))
    v <- (-0.5*f(x-h) + 0.5*f(x+h)) / h
    return(list(par = h, value = v))
  }
  
  me <- .Machine$double.eps / 2
  # Set cores <- 1 to get overhead-free timings
  t.opt    <- system.time(getHopt    <- mclapply(xgrid, function(x) fHopt(x = x), mc.cores = 1))
  t.RoT    <- system.time(getHRoT    <- mclapply(xgrid, function(x) fHRoT(x = x), mc.cores = 1))
  t.plugin <- system.time(getHplugin <- mclapply(xgrid, function(x) step.plugin(x = x, FUN = f, max.rel.error = me), mc.cores = cores))
  t.CR     <- system.time(getHCR     <- mclapply(xgrid, function(x) step.CR(x = x, FUN = f, max.rel.error = me, acc.order = 1), mc.cores = cores))
  t.CRm    <- system.time(getHCRm    <- mclapply(xgrid, function(x) step.CR(x = x, FUN = f, max.rel.error = me, acc.order = 2), mc.cores = cores))
  t.DV     <- system.time(getHDV     <- mclapply(xgrid, function(x) step.DV(x = x, FUN = f, max.rel.error = me), mc.cores = cores))
  t.SW     <- system.time(getHSW     <- mclapply(xgrid, function(x) step.SW(x = x, FUN = f, max.rel.error = me), mc.cores = cores))
  t.M      <- system.time(getHM      <- mclapply(xgrid, function(x) step.M(x = x, FUN = f, max.rel.error = me), mc.cores = cores))
  t.K      <- system.time(getHK      <- mclapply(xgrid, function(x) step.K(x = x, FUN = f, max.rel.error = me), mc.cores = cores))
  
  tms <- c(t.opt[3], t.RoT[3], t.CR[3], t.CRm[3], t.DV[3], t.SW[3], t.M[3], t.K[3])
  ms <- tms / length(xgrid) * 1000
  
  res <- list(getHopt, getHRoT, getHplugin, getHCR, getHCRm, getHDV, getHSW, getHM, getHK)
  names(res) <- c("opt", "RoT", "plugin", "CR", "CRm", "DV", "SW", "M", "K")
  
  h <- do.call(cbind, lapply(res, function(x) unlist(lapply(x, "[[", "par"))))
  
  # d <- apply(h, 2, function(x) tryCatch(density(log10(x), bw = "bcv"), error = function(e) density(log10(x), bw = "nrd")))
  # plot(NULL, NULL, xlim = range(unlist(lapply(d, "[[", "x"))), ylim = range(unlist(lapply(d, "[[", "y"))), bty = "n")
  # for (i in 1:ncol(h)) lines(d[[i]], col = i %% 6 + 1, lty = i %% 5 + 1, lwd = 1.5)
  
  for (i in 1:ncol(h)) {
    if (all(!is.finite(h[, i]))) next
    plot(xgrid, log10(h[, i]), log = if (linear) "" else "x", col = i %% 6 + 1, main = colnames(h)[i], pch = ".", ylim = quantile(log10(h[is.finite(h)]), c(0.01, 0.99)))
    l <- predict(loess(log10(h[, i]) ~ xgrid, span = 0.1, degree = 1))
    lines(xgrid, l, col = "#FFFFFFEE", lwd = 5)
    lines(xgrid, l, col = i %% 6 + 1, lwd = 2)
    # Sys.sleep(1)
  }
  
  g.true <- fp(xgrid)
  g <- do.call(cbind, lapply(res, function(x) unlist(lapply(x, "[[", "value"))))
  
  
  err <- g.true - g
  abserr  <- as.data.frame(abs(err))
  absrerr <- as.data.frame(abs(err / g.true))
  
  jj <- c(6, 7, 8, 9)
  
  cat(paste0(sprintf("%1.2e", apply(absrerr, 2, median)[jj]), collapse = " & "), "\n")
  
  ranks <- t(apply(absrerr[, jj], 1, rank))
  cat("New = best in", 100*mean(ranks[, "K"] == apply(ranks, 1, min)), "%\n")
  
  lxgrid <- if (!linear) log10(xgrid) else xgrid
  labserr <- log10(abserr)
  labserr[labserr == -Inf] <- min(unlist(labserr)[is.finite(unlist(labserr))])
  xnice <- seq(quantile(lxgrid, 0.001), quantile(lxgrid, 0.999), length.out = 201)
  l <- do.call(cbind, parallel::mclapply(labserr, function(y){
    if (all(is.nan(y))) return(xnice * NA)
    predict(loess(y ~ lxgrid, degree = 1, family = "symmetric",
                  control = loess.control(surface = "direct"), span = 0.2), newdata = xnice)
  }, mc.cores = cores))
  
  mycols <- c("#ebce2b", "#702c8c", "#db6917", "#96cde6", "#ba1c30", "#c0bd7f", "#7f7e80", "#5fa641", "#d485b2", "#4277b6")
  
  # One plot with everything
  cairo_pdf(paste0("compare-methods-", flab, "-", if (linear) "lin" else "exp", ".pdf"), imgh/90, imgv/60)
  par(mar = c(4, 4, 0, 0) + .1)
  changeFont()
  matplot(if (linear) xnice else 10^xnice, l, type = "l", log = if (linear) "" else "x", lty = 1:5, col = mycols, lwd = 2, xlab = "Evaluation point", ylab = "log10(|error|)", bty = "n")
  legend("top", legend = colnames(abserr), ncol = 4, lwd = 2, col = mycols, lty = 1:5, bg = "#FFFFFFCC", box.col = "#FFFFFF00")
  dev.off()
  
}
