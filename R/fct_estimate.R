# ---------------------------------------------------------------------------
# Per-group reference-interval estimation. See SPEC.md sections 5, 6, 7.
# reflimR is the primary estimator; reflimLOD takes over when the lower limit
# sits at/below the LOD; groups with > 20% sub-LOD values are blocked.
# Pure logic (depends on reflimR + reflimLOD.R), no Shiny.
# ---------------------------------------------------------------------------

LOD_BLOCK_SHARE <- 0.20   # > 20% sub-LOD  ->  estimation not possible (SPEC 6.2)

#' Estimate one group's reference interval, applying the LOD decision rule.
#' @param values numeric value vector for the group (NA for sub-LOD rows)
#' @param n_lod  number of sub-LOD (< LOD) rows in the group
#' @param lod_value analyte LOD (NA if none)
#' @param n_min reflimR minimum n
#' @return list with status ("ok"/"blocked"/"failed"), method, n, n_lod, share,
#'         lower, upper, lower.ci, upper.ci, x (measurable values), note
estimate_group <- function(values, n_lod, lod_value, n_min = 200) {
  measurable <- values[!is.na(values) & values > 0]
  n_total <- length(measurable) + n_lod
  share   <- if (n_total > 0) n_lod / n_total else 0

  base <- list(n = n_total, n_lod = n_lod, share = share,
               lower = NA_real_, upper = NA_real_,
               lower.ci = c(NA_real_, NA_real_), upper.ci = c(NA_real_, NA_real_),
               x = measurable, lod_value = lod_value, lognormal = FALSE, note = NULL)

  # --- block if too many sub-LOD values -----------------------------------
  if (!is.na(lod_value) && share > LOD_BLOCK_SHARE) {
    return(modifyList(base, list(
      status = "blocked", method = NA_character_,
      note = sprintf("Anteil < LOD = %.0f %% (> %.0f %%): Schätzung mit diesen Algorithmen nicht möglich.",
                     100 * share, 100 * LOD_BLOCK_SHARE))))
  }

  # --- primary: reflimR ---------------------------------------------------
  r <- tryCatch(reflimR::reflim(measurable, plot.it = FALSE, n.min = n_min),
                error = function(e) NULL, warning = function(w) {
                  suppressWarnings(reflimR::reflim(measurable, plot.it = FALSE, n.min = n_min))
                })
  if (is.null(r)) {
    return(modifyList(base, list(status = "failed", method = NA_character_,
                                 note = "reflimR-Schätzung fehlgeschlagen.")))
  }
  lim <- r$limits
  lower <- unname(lim["lower.lim"]); upper <- unname(lim["upper.lim"])
  lower.ci <- unname(c(lim["lower.lim.low"], lim["lower.lim.upp"]))
  upper.ci <- unname(c(lim["upper.lim.low"], lim["upper.lim.upp"]))
  method <- "reflimR"; note <- NULL; lognormal <- isTRUE(r$lognormal)

  # --- LOD switch: lower limit at/below LOD -------------------------------
  if (!is.na(lod_value) && n_lod > 0 &&
      !is.na(lower) && (lower <= lod_value || lower.ci[1] <= lod_value)) {
    lodres <- tryCatch(lod.iboxplot95(measurable, lod_value, n_lod),
                       error = function(e) NULL)
    if (!is.null(lodres) && !is.na(lodres$lower.limit)) {
      lower <- lodres$lower.limit; upper <- lodres$upper.limit
      lower.ci <- c(NA_real_, NA_real_); upper.ci <- c(NA_real_, NA_real_)
      method <- "reflimLOD"; lognormal <- TRUE
      note <- sprintf("Algorithmus auf reflimLOD umgestellt (%d Werte < LOD = %s, %.0f %%).",
                      n_lod, lod_value, 100 * share)
    }
  }

  modifyList(base, list(status = "ok", method = method, lognormal = lognormal,
                        lower = round_lim(lower), upper = round_lim(upper),
                        lower.ci = round_lim(lower.ci), upper.ci = round_lim(upper.ci),
                        note = note))
}

#' Sensible digit rounding for a limit (mirrors reflimR's scale-based rounding).
round_lim <- function(v) {
  m <- suppressWarnings(stats::median(abs(v[is.finite(v)]), na.rm = TRUE))
  if (!is.finite(m) || m == 0) return(round(v, 2))
  d <- 2 - floor(log10(m)); if (d < 0) d <- 0
  round(v, d)
}

#' Map reflimR interpretation text -> traffic-light colour.
interpret_to_color <- function(txt) {
  if (is.na(txt)) return("grey")
  if (grepl("within", txt))   return("green")
  if (grepl("slightly", txt)) return("yellow")
  if (grepl("markedly", txt)) return("red")
  "grey"
}

#' Traffic-light verdict comparing an estimate against the lab's current limits.
#' reflimR groups use reflimR's own interpretation (authoritative);
#' reflimLOD groups use permissible_uncertainty band overlap (green/red).
#' @return list(lower = colour, upper = colour) or NULL when inputs incomplete
group_verdict <- function(est, current_lower, current_upper, n_min = 200) {
  if (est$status != "ok") return(NULL)
  if (is.null(current_lower) || is.null(current_upper) ||
      is.na(current_lower) || is.na(current_upper)) return(NULL)

  if (est$method == "reflimR") {
    r <- tryCatch(suppressWarnings(
      reflimR::reflim(est$x, targets = c(current_lower, current_upper),
                      plot.it = FALSE, n.min = n_min)),
      error = function(e) NULL)
    if (is.null(r)) return(NULL)
    return(list(lower = interpret_to_color(r$interpretation[["lower.limit"]]),
                upper = interpret_to_color(r$interpretation[["upper.limit"]])))
  }

  # reflimLOD: band-overlap rule via permissible_uncertainty
  pe <- reflimR::permissible_uncertainty(est$lower, est$upper)
  pt <- reflimR::permissible_uncertainty(current_lower, current_upper)
  ov <- function(a, b) max(a[1], b[1]) <= min(a[2], b[2])
  low <- if (ov(c(pe["lower.lim.low"], pe["lower.lim.upp"]),
               c(pt["lower.lim.low"], pt["lower.lim.upp"]))) "green" else "red"
  up  <- if (ov(c(pe["upper.lim.low"], pe["upper.lim.upp"]),
               c(pt["upper.lim.low"], pt["upper.lim.upp"]))) "green" else "red"
  list(lower = low, upper = up)
}

COLOR_HEX  <- c(green = "#1a9850", yellow = "#d6a400", red = "#d73027", grey = "#888888")
COLOR_TEXT <- c(green = "konsistent", yellow = "Prüfung empfohlen",
                red = "Prüfung dringend empfohlen", grey = "—")

#' Compare two reference-interval estimates (e.g. reflimR vs refineR) per limit
#' using permissible_uncertainty band overlap. Green = bands overlap (agree),
#' red = bands separated (differ beyond permissible uncertainty). See SPEC 7.1.
compare_methods <- function(lo1, up1, lo2, up2) {
  if (any(is.na(c(lo1, up1, lo2, up2)))) return(NULL)
  p1 <- reflimR::permissible_uncertainty(lo1, up1)
  p2 <- reflimR::permissible_uncertainty(lo2, up2)
  ov <- function(a, b) max(a[1], b[1]) <= min(a[2], b[2])
  list(
    lower = if (ov(c(p1["lower.lim.low"], p1["lower.lim.upp"]),
                   c(p2["lower.lim.low"], p2["lower.lim.upp"]))) "green" else "red",
    upper = if (ov(c(p1["upper.lim.low"], p1["upper.lim.upp"]),
                   c(p2["upper.lim.low"], p2["upper.lim.upp"]))) "green" else "red")
}

#' Plain-language interpretation of an mclust cluster result (German). SPEC 7.3.
#' @param stats the $stats data.frame from lab_mclust (columns ll, ul, percent)
mclust_sentence <- function(stats) {
  g <- nrow(stats)
  if (g <= 1)
    return("Die Daten enthalten scheinbar nur eine Population – keine Hinweise auf versteckte Subpopulationen.")
  o <- order(stats$percent, decreasing = TRUE)
  main <- o[1]
  main_pct <- stats$percent[main]
  others <- o[-1]
  descr <- vapply(others, function(i) {
    dir <- if (stats$ll[i] > stats$ul[main]) "höheren"
           else if (stats$ul[i] < stats$ll[main]) "niedrigeren" else "überlappenden"
    sprintf("eine Gruppe mit %s Werten (%.0f %%)", dir, stats$percent[i])
  }, character(1))
  sprintf(paste0("Die Daten enthalten scheinbar %d Subpopulationen – die Hauptgruppe (%.0f %%) und %s. ",
                 "Das kann darauf hindeuten, dass das Referenzintervall durch nicht-gesunde Proben beeinflusst wird."),
          g, main_pct, paste(descr, collapse = " sowie "))
}
