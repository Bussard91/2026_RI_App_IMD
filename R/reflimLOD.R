# ---------------------------------------------------------------------------
# reflimLOD: reference-limit estimation in the presence of values below the
# limit of detection (LOD). Functions copied verbatim (definitions only) from
# Skripte_und_Pakete/iboxplot95_lod_gh.r, the example/demo block removed.
# The CI bootstrap's internal call to the stale name `lodiboxplot95` has been
# corrected to `lod.iboxplot95` (SPEC.md section 10, item 1).
# ---------------------------------------------------------------------------

# requireNamespace("truncnorm") — provides ptruncnorm / dtruncnorm

box.cox.trans <- function(x, lambda = 1) {
  if (lambda == 0) return(log(x)) else return((x^lambda - 1) / lambda)
}

box.cox.inv.trans <- function(x, lambda = 1) {
  if (lambda == 0) return(exp(x)) else return((x * lambda + 1)^(1 / lambda))
}

lod.boxplot.truncation <- function(measured.values, lod, n.lod, right.quantile = 0.75) {
  n.min.values <- 2 * n.lod + 1
  quantile.factor <- qnorm(0.975) / qnorm(right.quantile)

  x <- c(rep(lod - 1, n.lod), measured.values)
  if (length(x) < n.min.values)
    return(list(selected.values = NA, upper.truncation = NA, median = NA, lod = lod, n.lod = n.lod))

  x.qmr <- quantile(x, probs = c(0.5, right.quantile))
  upper.truncation <- x.qmr[1] + (x.qmr[2] - x.qmr[1]) * quantile.factor

  x.length.old <- length(x)
  x <- x[x <= upper.truncation]
  if (length(x) < n.min.values)
    return(list(selected.values = NA, upper.truncation = NA, median = NA, lod = lod, n.lod = n.lod))
  if (length(x) == x.length.old)
    return(list(selected.values = x[-(1:n.lod)], upper.truncation = upper.truncation,
                median = x.qmr[1], lod = lod, n.lod = n.lod))

  median.mod <- 0.5 / 0.975
  right.quantile.mod <- right.quantile / 0.975

  while (length(x) < x.length.old) {
    x.length.old <- length(x)
    x.qmr <- quantile(x, probs = c(median.mod, right.quantile.mod))
    upper.truncation <- x.qmr[1] + (x.qmr[2] - x.qmr[1]) * quantile.factor
    x <- x[x <= upper.truncation]
    if (length(x) < n.min.values)
      return(list(selected.values = NA, upper.truncation = NA, median = NA, lod = lod, n.lod = n.lod))
  }

  list(selected.values = x[-(1:n.lod)], upper.truncation = upper.truncation,
       median = x.qmr[1], lod = lod, n.lod = n.lod)
}

lod.iboxplot95 <- function(measured.values, lod, n.lod, lambda = 0, right.quantile = 0.75) {
  transformed.measured.values <- box.cox.trans(measured.values, lambda = lambda)
  transformed.lod <- box.cox.trans(lod, lambda = lambda)

  normal.result <- lod.boxplot.truncation(transformed.measured.values, transformed.lod,
                                          n.lod, right.quantile = right.quantile)

  if (is.na(normal.result$upper.truncation))
    return(list(lower.limit = NA, upper.limit = NA, mu.log = NA, sigma.log = NA,
                upper.truncation = NA, selected.values = NA, lod = lod, n.lod = n.lod))

  obj.fun <- function(pars) {
    pnorm.lod <- truncnorm::ptruncnorm(transformed.lod, b = normal.result$upper.truncation,
                                       mean = pars[1], sd = pars[2])
    -n.lod * log(pnorm.lod) -
      length(normal.result$selected.values) * log(1 - pnorm.lod) -
      sum(log(truncnorm::dtruncnorm(normal.result$selected.values,
                                    a = transformed.lod, b = normal.result$upper.truncation,
                                    mean = pars[1], sd = pars[2])))
  }

  pars.initial <- c(normal.result$median,
                    (normal.result$upper.truncation - normal.result$median) / qnorm(0.975))
  optim.result <- optim(pars.initial, obj.fun)

  lims <- box.cox.inv.trans(qnorm(c(0.025, 0.975), mean = optim.result$par[1],
                                  sd = optim.result$par[2]), lambda = lambda)
  lims1perc <- box.cox.inv.trans(qnorm(c(0.01, 0.99), mean = optim.result$par[1],
                                       sd = optim.result$par[2]), lambda = lambda)

  list(lower.limit = lims[1], upper.limit = lims[2],
       percentile1 = lims1perc[1], percentile99 = lims1perc[2],
       mu.log = optim.result$par[1], sigma.log = optim.result$par[2],
       upper.truncation = box.cox.inv.trans(normal.result$upper.truncation, lambda = lambda),
       selected.values = box.cox.inv.trans(normal.result$selected.values, lambda = lambda),
       lod = lod, n.lod = n.lod)
}

ci.lodiboxplot95 <- function(measured.values, lod, n.lod, lambda = 0, right.quantile = 0.75,
                             conf.level = 0.95, n.bootstrap = 1000) {
  lower.limits <- upper.limits <- percentile1s <- percentile99s <-
    mu.logs <- sigma.logs <- upper.truncations <- rep(NA, n.bootstrap)

  x <- c(rep(lod / 2, n.lod), measured.values)
  for (i in 1:n.bootstrap) {
    xi <- sample(x, length(x), replace = TRUE)
    mvi <- subset(xi, xi >= lod)
    n.lodi <- sum(xi < lod)
    resi <- lod.iboxplot95(mvi, lod, n.lodi, lambda = lambda, right.quantile = right.quantile)
    if (!is.na(resi$lower.limit)) {
      lower.limits[i] <- resi$lower.limit; upper.limits[i] <- resi$upper.limit
      percentile1s[i] <- resi$percentile1; percentile99s[i] <- resi$percentile99
      mu.logs[i] <- resi$mu.log; sigma.logs[i] <- resi$sigma.log
      upper.truncations[i] <- resi$upper.truncation
    }
  }
  conf.lims <- c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2)
  list(lower.limit.ci = quantile(lower.limits, conf.lims, na.rm = TRUE),
       upper.limit.ci = quantile(upper.limits, conf.lims, na.rm = TRUE),
       percentile1.ci = quantile(percentile1s, conf.lims, na.rm = TRUE),
       percentile99.ci = quantile(percentile99s, conf.lims, na.rm = TRUE),
       mu.log.ci = quantile(mu.logs, conf.lims, na.rm = TRUE),
       sigma.log.ci = quantile(sigma.logs, conf.lims, na.rm = TRUE),
       upper.truncation.ci = quantile(upper.truncations, conf.lims, na.rm = TRUE))
}
