# ---------------------------------------------------------------------------
# Gaussian-mixture RI verification. Copied (definition only) from
# Skripte_und_Pakete/mclust_functions.R, with one fix: `y.max` is now
# initialised before the plotting block so plot.it = FALSE no longer errors
# (the original returned y.max unconditionally). See SPEC.md sections 7.3, 10.
# requireNamespace("mclust")
# ---------------------------------------------------------------------------

lab_mclust <- function(x, lognormal, remove.extremes,
                       targets = NULL,
                       plot.it = TRUE, add.boxplot = TRUE,
                       plot.legend = TRUE, pos.legend = "topright",
                       plot.bic = FALSE,
                       xlim = NULL, ylim = NULL,
                       main = "", xlab = "",
                       hist.bins = 50, model = NULL, n.cluster = NULL,
                       apply.rounding = TRUE, digits = NULL) {

  x <- na.omit(as.numeric(x[x > 0]))
  lx <- length(x)
  if (lx < 200)
    warning(paste("The data has", lx, "numeric elements, where 200 are recommended."))

  if (remove.extremes) x <- x[x < (median(x) + 6 * IQR(x))]
  xx <- if (lognormal) log(x) else x

  if (length(xx) < 200 & lx >= 200)
    warning(paste("After removal of extremes, the data has", length(xx),
                  "numeric elements, where 200 are recommended."))

  mc <- if (is.null(model)) mclust::Mclust(xx, G = n.cluster)
        else mclust::Mclust(xx, modelNames = model, G = n.cluster)

  for (i in 1:mc$G)
    if (is.na(mc$parameters$variance$sigmasq[i]))
      mc$parameters$variance$sigmasq[i] <- mc$parameters$variance$sigmasq[1]

  res.tab <- data.frame(matrix(NA, mc$G, 5))
  for (i in 1:mc$G) {
    res.tab[i, 3] <- round(mc$parameters$pro[i] * 100, 1)
    res.tab[i, 4] <- round(mc$parameters$mean[i], 3)
    res.tab[i, 5] <- round(sqrt(mc$parameters$variance$sigmasq[i]), 3)
  }
  if (lognormal) {
    colnames(res.tab) <- c("ll", "ul", "percent", "meanlog", "sdlog")
    for (i in 1:mc$G) {
      res.tab[i, 1] <- qlnorm(0.025, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
      res.tab[i, 2] <- qlnorm(0.975, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
    }
  } else {
    colnames(res.tab) <- c("ll", "ul", "percent", "mean", "sd")
    for (i in 1:mc$G) {
      res.tab[i, 1] <- qnorm(0.025, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
      res.tab[i, 2] <- qnorm(0.975, mc$parameters$mean[i], sqrt(mc$parameters$variance$sigmasq[i]))
    }
  }
  if (is.null(digits)) { digits <- 2 - floor(log10(median(xx))); if (digits < 0) digits <- 0 }
  if (apply.rounding) {
    res.tab[, 1] <- round(res.tab[, 1], digits)
    res.tab[, 2] <- round(res.tab[, 2], digits)
  }

  y.max <- NA   # <-- fix: defined unconditionally (was only set inside plot block)

  if (plot.it) {
    col <- rainbow(9)
    d <- density(x)
    y.max <- max(d$y) * 1.1
    if (add.boxplot) y.max <- y.max * 1.4
    breaks <- if (length(x) > 200) seq(0.9 * min(x), 1.1 * max(x), length.out = hist.bins) else "Sturges"
    if (is.null(xlim)) xlim <- c(min(x), max(x))
    if (is.null(ylim)) ylim <- c(0, y.max * 1.1)
    hist(x, freq = FALSE, breaks = breaks, col = "white", border = "grey",
         xlim = xlim, ylim = ylim, yaxt = "n", main = main, xlab = xlab, ylab = "")
    box(); lines(d, lty = 2)
    if (add.boxplot)
      boxplot(x, horizontal = TRUE, at = y.max * 0.9, boxwex = y.max * 0.1, add = TRUE)
    if (lognormal) {
      for (i in 1:nrow(res.tab))
        curve(dlnorm(x, meanlog = mc$parameters$mean[i], sdlog = sqrt(mc$parameters$variance$sigmasq[i])) * mc$parameters$pro[i],
              from = min(x), to = max(x), lwd = 2, col = col[i], add = TRUE)
    } else {
      for (i in 1:nrow(res.tab))
        curve(dnorm(x, mean = mc$parameters$mean[i], sd = sqrt(mc$parameters$variance$sigmasq[i])) * mc$parameters$pro[i],
              from = min(x), to = max(x), lwd = 2, col = col[i], add = TRUE)
    }
    if (!is.null(targets)) {
      lines(rep(targets[1], 2), c(0, y.max * 0.8), lty = 2)
      lines(rep(targets[2], 2), c(0, y.max * 0.8), lty = 2)
      text(targets, rep(y.max * 0.85, 2), targets)
    }
    if (plot.legend)
      legend(pos.legend,
             paste0(round(res.tab$ll, digits), "-", round(res.tab$ul, digits),
                    " (", res.tab$percent, "%)"),
             lwd = 2, col = col[1:nrow(res.tab)], cex = 0.8, bty = "n")
    if (plot.bic) plot(mc$BIC)
  }

  n.c <- if (is.null(n.cluster)) mc$G else paste(mc$G, "from", deparse(n.cluster))
  list(n.cluster = noquote(n.c), stats = res.tab, y.max = y.max, BIC = mc$BIC)
}
