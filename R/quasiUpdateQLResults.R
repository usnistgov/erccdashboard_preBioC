QL.results <- function (fit, Dispersion = "Deviance", spline.df = NULL, Plot = TRUE) 
{
  if (!Dispersion %in% c("Deviance", "Pearson")) 
    stop("Unidentified Dispersion: Dispersion must be either 'Deviance' or 'Pearson'.")
  LRT <- fit$LRT
  phi.hat <- fit$phi.hat.dev
  mn.cnt <- fit$mn.cnt
  den.df <- fit$den.df
  num.df <- fit$num.df
  Model = fit$Model
  if (Dispersion == "Pearson") 
    phi.hat <- fit$phi.hat.pearson
  if (length(num.df) == 1) 
    num.df <- rep(num.df, ncol(LRT))
  shrink.phi <- function(phi.hat, den.df) {
    phi.hat[phi.hat <= 0] <- min(phi.hat[phi.hat > 0])
    z <- log(phi.hat)
    z[z == Inf] <- max(z[z != Inf])
    z[z == -Inf] <- min(z[z != -Inf])
    mnz <- mean(z)
    d0arg <- var(z) - trigamma((den.df)/2)
    if (d0arg > 0) {
      dif <- function(x, y) abs(trigamma(x) - y)
      inverse.trigamma <- function(y) optimize(dif, interval = c(0, 
                                                                 10000), y = y)$minimum
      d0 <- 2 * inverse.trigamma(d0arg)
      phi0 <- exp(mnz - digamma((den.df)/2) + digamma(d0/2) - 
                    log(d0/(den.df)))
      phi.shrink <- ((den.df) * phi.hat + d0 * phi0)/(den.df + 
                                                        d0)
    }
    else {
      phi.shrink <- rep(exp(mnz), length(z))
      d0 <- Inf
      phi0 <- exp(mnz)
    }
    return(list(phi.shrink = phi.shrink, d0 = d0, phi0 = phi0))
  }
  phi.hat[phi.hat < 0] <- min(phi.hat[phi.hat > 0])
  phi.hat2 <- phi.hat
  if (Model == "Poisson") 
    phi.hat2[phi.hat < 1] <- 1
  shrink <- shrink.phi(phi.hat, den.df)
  phi.shrink <- shrink[[1]]
  est.d0 <- shrink[[2]]
  if (Model == "Poisson") 
    phi.shrink[phi.shrink < 1] <- 1
  y <- log(phi.hat)
  y[y == -Inf] <- min(y[y != -Inf])
  y[y == Inf] <- max(y[y != Inf])
  spline.fit <- if (is.null(spline.df)) 
    smooth.spline(x = log(mn.cnt), y = y)
  else smooth.spline(x = log(mn.cnt), y = y, df = spline.df)
  spline.pred <- predict(spline.fit, x = log(mn.cnt))$y
  fit.method <- "spline"
  y2 <- phi.hat/exp(spline.pred)
  shrink <- shrink.phi(y2, den.df)
  D0 <- shrink[[2]]
  phi0 <- shrink[[3]]
  print(paste("Spline scaling factor:", phi0))
  if (Plot) {
    dev.new(height = 9)
    nf <- layout(matrix(1:2, 2, 1), heights = c(7, 2))
    par(mai = c(1, 1.2, 1, 0.2))
    suppressWarnings(plot(log(mn.cnt), y, xlab = expression(log(bar(y)[phantom() %.% 
                                                                         phantom()] * phantom()[phantom() %.% phantom()])), 
                          ylab = expression(log(hat(Phi))), main = "Estimated Dispersion\n versus Average Count", 
                          pch = 1, cex.lab = 2, cex.axis = 2, cex.main = 2))
    lines(sort(log(mn.cnt)), spline.pred[order(mn.cnt)], 
          col = 2, lwd = 3)
    RR <- NULL
    sort.mn <- log(sort(mn.cnt))
    qq <- c(0.05, 0.95)
    ord.F <- y[order(mn.cnt)]
    bins <- c(1 + round(length(y) * (0:19)/20), length(y))
    for (ii in 1:(length(bins) - 1)) {
      ind <- bins[ii]:bins[ii + 1]
      RR <- rbind(RR, c(quantile(ord.F[ind], qq), median(sort.mn[ind])))
    }
    lines(sort(log(mn.cnt)), spline.pred[order(mn.cnt)] + 
            log(phi0 * qf(0.95, den.df, D0)), col = 4, lwd = 3)
    lines(sort(log(mn.cnt)), spline.pred[order(mn.cnt)] + 
            log(phi0 * qf(0.05, den.df, D0)), col = 4, lwd = 3)
    lines(RR[, 3], RR[, 2], col = 3, lwd = 3)
    lines(RR[, 3], RR[, 1], col = 3, lwd = 3)
    par(mai = rep(0, 4))
    plot(19, 19, col = "white", axes = FALSE)
    suppressWarnings(legend("top", legend = c(paste("Fitted", 
                                                    fit.method, "with", signif(spline.fit$df, 2), "df"), 
                                              "0.05 & 0.95 Quantiles from Empirical Distribution", 
                                              "0.05 & 0.95 Quantiles from Scaled F-Distribution"), 
                            lwd = 3, lty = 1, col = 2:4, cex = 1.5))
  }
  phi.spline <- (D0 * exp(spline.pred) * phi0 + (den.df) * 
                   phi.hat)/(D0 + den.df)
  if (D0 == Inf) {
    warning(cat("D0 estimate is infinity for QLSpline (there's little scatter in original dispersion estimates around fitted spline).\n QLSpline dispersions set to fitted cubic spline (use 'Plot=TRUE' to view) with no uncertainty."))
    #phi.spline <- exp(spline.fit$fitted.values)
    phi.spline <- exp(spline.pred)
  }
  if (Model == "Poisson") 
    phi.spline[phi.spline < 1] <- 1
  log.pval <- F.stat <- list(QL = NULL, QLShrink = NULL, QLSpline = NULL)
  for (i in 1:ncol(LRT)) {
    log.pval[[1]] <- cbind(log.pval[[1]], pf(LRT[, i]/phi.hat2, 
                                             num.df[i], den.df, lower.tail = FALSE, log.p = TRUE))
    log.pval[[2]] <- cbind(log.pval[[2]], pf(LRT[, i]/phi.shrink, 
                                             num.df[i], est.d0 + den.df, lower.tail = FALSE, log.p = TRUE))
    log.pval[[3]] <- cbind(log.pval[[3]], pf(LRT[, i]/phi.spline, 
                                             num.df[i], D0 + den.df, lower.tail = FALSE, log.p = TRUE))
    F.stat[[1]] <- cbind(F.stat[[1]], LRT[, i]/phi.hat2)
    F.stat[[2]] <- cbind(F.stat[[2]], LRT[, i]/phi.shrink)
    F.stat[[3]] <- cbind(F.stat[[3]], LRT[, i]/phi.spline)
  }
  pval <- log.pval
  for (ii in 1:3) {
    pval[[ii]] <- exp(log.pval[[ii]])
    colnames(F.stat[[ii]]) <- colnames(log.pval[[ii]]) <- colnames(pval[[ii]]) <- colnames(LRT)
  }
  d0 <- c(QLShrink = est.d0, QLSpline = D0)
  estimate.m0 <- function(p, B = 20) {
    m <- length(p)
    m0 <- m
    bin <- c(-0.1, (1:B)/B)
    bin.counts = rep(0, B)
    for (i in 1:B) {
      bin.counts[i] = sum((p > bin[i]) & (p <= bin[i + 
                                                     1]))
    }
    tail.means <- rev(cumsum(rev(bin.counts))/(1:B))
    temp <- bin.counts - tail.means
    index <- min((1:B)[temp <= 0])
    m0 <- B * tail.means[index]
    return(m0)
  }
  jabes.q <- function(p, B = 20) {
    m = length(p)
    m0 = estimate.m0(p, B)
    k = 1:m
    ord = order(p)
    p[ord] = (p[ord] * m0)/(1:m)
    qval = p
    qval[ord] = rev(cummin(rev(qval[ord])))
    return(qval)
  }
  qval <- pval
  m0 <- NULL
  for (ii in 1:length(qval)) {
    M0 <- NULL
    Qval <- qval[[ii]]
    for (jj in 1:ncol(Qval)) {
      M0 <- c(M0, estimate.m0(Qval[!is.na(Qval[, jj]), 
                                   jj]))
      qval[[ii]][!is.na(Qval[, jj]), jj] = jabes.q(Qval[!is.na(Qval[, 
                                                                    jj]), jj])
    }
    m0 <- rbind(m0, M0)
  }
  colnames(m0) <- colnames(qval[[1]])
  rownames(m0) <- names(log.pval)
  return(list(P.values = pval, log.P.values = log.pval, Q.values = qval, 
              F.stat = F.stat, m0 = m0, d0 = d0))
}