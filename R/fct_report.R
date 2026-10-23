# ---------------------------------------------------------------------------
# Self-contained HTML report + CSV assembly. See SPEC.md section 8.
# Plots are rendered to PNG and embedded as base64 data URIs, so the report
# needs neither pandoc nor LaTeX and works identically locally and on
# shinyapps.io. Depends on htmltools + knitr (image_uri). No Shiny.
# ---------------------------------------------------------------------------

IMD_TEAL_R     <- "#005357"
IMD_CHARCOAL_R <- "#1D1D1B"

#' X-axis label from analyte + unit, e.g. "Hämoglobin [g/dL]".
axis_label <- function(analyte, unit) {
  a <- if (!is.null(analyte) && nzchar(analyte)) analyte else "Messwert"
  if (!is.null(unit) && nzchar(unit)) paste0(a, " [", unit, "]") else a
}

#' Plot title from analyte + group label, e.g. "Hämoglobin – M 0-56".
plot_title <- function(analyte, grp) {
  if (!is.null(analyte) && nzchar(analyte)) paste0(analyte, " – ", grp) else grp
}

#' Age-trend scatter for one sex: grey points (age vs value) with, per age
#' segment, the segment MEAN drawn as a horizontal line (= the decision-tree
#' leaf mean when groups match the tree), labelled, and red vertical lines
#' separating the groups. y-axis is scaled to the full data range so nothing
#' is clipped.
plot_age_trend <- function(d_sex, groups_sex, main = "",
                           xlab = "Alter (Jahre)", ylab = "Messwert") {
  x <- d_sex$age; y <- d_sex$value
  ok <- !is.na(x) & !is.na(y) & y > 0
  x <- x[ok]; y <- y[ok]
  if (length(y) == 0) { plot.new(); title(main); return(invisible()) }

  g <- groups_sex[order(groups_sex$from), , drop = FALSE]
  n <- nrow(g)
  means <- vapply(seq_len(n), function(i) {
    seg <- y[x >= g$from[i] & x <= g$to[i]]
    if (length(seg)) mean(seg) else NA_real_
  }, numeric(1))

  # y-range covers all points AND all mean lines, with a small margin
  rng <- range(c(y, means), na.rm = TRUE)
  pad <- diff(rng) * 0.06; if (pad == 0) pad <- max(abs(rng[1]) * 0.05, 1)
  ylim <- c(rng[1] - pad, rng[2] + pad)

  plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim,
       pch = 20, col = grDevices::rgb(0, 0, 0, 0.22), cex = 0.6)
  grid()
  # red vertical group separators
  if (n > 1) abline(v = g$to[-n] + 0.5, col = "red", lwd = 1.5, lty = 2)
  # segment mean lines + value labels
  for (i in seq_len(n)) {
    if (is.na(means[i])) next
    lines(c(g$from[i], g$to[i]), rep(means[i], 2), lwd = 3, col = "#005357")
    text((g$from[i] + g$to[i]) / 2, means[i], labels = round_lim(means[i]),
         pos = 3, col = "#005357", font = 2, cex = 0.9)
  }
}

#' Decision-tree image for one sex. Uses rpart.plot when available (nicer),
#' else falls back to base rpart plotting. Placeholder if the tree has no split.
plot_age_tree <- function(model, main = "") {
  if (is.null(model) || is.null(model$splits) || nrow(model$splits) == 0) {
    plot.new(); title(main); text(0.5, 0.5, "Kein Alterssplit\n(eine Altersgruppe)")
    return(invisible())
  }
  if (requireNamespace("rpart.plot", quietly = TRUE)) {
    rpart.plot::rpart.plot(model, type = 5, main = main,
                           box.palette = "Greens", shadow.col = "grey85")
  } else {
    plot(model, uniform = TRUE, margin = 0.1, main = main)
    text(model, use.n = TRUE, cex = 0.8)
  }
}

#' Diagnostic plot for a reflimLOD group (LOD line + estimated limits).
plot_lod_group <- function(e, tg, main = e$label, xlab = "") {
  hist(e$x, breaks = 40, col = "white", border = "grey",
       main = main, xlab = xlab, freq = FALSE)
  abline(v = e$lod_value, col = "red", lwd = 2, lty = 3)
  abline(v = c(e$lower, e$upper), lwd = 2)
  if (!is.null(tg)) abline(v = tg, col = "blue", lwd = 2, lty = 2)
  legend("topright",
         c("LOD", "reflimLOD-Grenzen", if (!is.null(tg)) "aktuell"),
         col = c("red", "black", if (!is.null(tg)) "blue"),
         lwd = 2, lty = c(3, 1, if (!is.null(tg)) 2), bty = "n", cex = 0.8)
}

#' Render a plotting expression to a base64 PNG data URI.
plot_to_uri <- function(expr, width = 580, height = 320) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = width, height = height, res = 92)
  tryCatch(force(expr), error = function(e) { plot.new(); text(0.5, 0.5, "Plot-Fehler") })
  grDevices::dev.off()
  uri <- knitr::image_uri(f)
  unlink(f)
  uri
}

#' Flatten the per-group results into a tidy data.frame (used for CSV export).
#' @param results list of list(est, current_lower, current_upper, verdict,
#'        refiner, mclust) as assembled in the server.
build_results_df <- function(results, analyte = "", unit = "") {
  do.call(rbind, lapply(results, function(g) {
    e <- g$est
    v <- g$verdict
    rr <- g$refiner
    data.frame(
      Analyt = analyte, Einheit = unit,
      Gruppe = e$label, Geschlecht = e$sex, von = e$from, bis = e$to,
      n = e$n, Methode = ifelse(is.na(e$method), "—", e$method),
      RI_unten = e$lower, RI_oben = e$upper,
      KI_unten_low = e$lower.ci[1], KI_unten_upp = e$lower.ci[2],
      KI_oben_low = e$upper.ci[1], KI_oben_upp = e$upper.ci[2],
      aktuell_unten = g$current_lower %||% NA, aktuell_oben = g$current_upper %||% NA,
      Ampel_unten = if (is.null(v)) "—" else unname(COLOR_TEXT[[v$lower]]),
      Ampel_oben  = if (is.null(v)) "—" else unname(COLOR_TEXT[[v$upper]]),
      refineR_unten = if (is.null(rr) || !isTRUE(rr$ok)) NA else rr$lower,
      refineR_oben  = if (is.null(rr) || !isTRUE(rr$ok)) NA else rr$upper,
      stringsAsFactors = FALSE)
  }))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Verification summary per group: reflimR vs refineR vs mclust + recommendation.
#' Shared by the Step-4 results table and the report. See SPEC.md section 7.
build_verification_df <- function(results) {
  do.call(rbind, lapply(results, function(g) {
    e <- g$est
    if (e$status != "ok") {
      return(data.frame(Gruppe = e$label, reflimR = "—", refineR = "—",
        Vergleich = "—", mclust_Cluster = "—", mclust_Hauptgruppe = "—",
        Empfehlung = ifelse(is.null(e$note), "—", e$note), stringsAsFactors = FALSE))
    }
    rr <- g$refiner; mm <- g$mclust; v <- g$verdict
    refineR_txt <- if (!is.null(rr) && isTRUE(rr$ok)) sprintf("%s – %s", rr$lower, rr$upper) else "—"
    cmp <- if (!is.null(rr) && isTRUE(rr$ok))
      compare_methods(e$lower, e$upper, rr$lower, rr$upper) else NULL
    vergleich <- if (is.null(cmp)) "—"
      else if (cmp$lower == "red" || cmp$upper == "red") "abweichend" else "übereinstimmend"
    mcl_cluster <- "—"; mcl_haupt <- "—"
    if (!is.null(mm) && isTRUE(mm$ok)) {
      st <- mm$stats; mcl_cluster <- as.character(nrow(st))
      k <- which.max(st$percent)
      mcl_haupt <- sprintf("%s – %s (%.0f %%)", st$ll[k], st$ul[k], st$percent[k])
    }
    mcl_multi <- !is.null(mm) && isTRUE(mm$ok) && nrow(mm$stats) > 1
    emp <- if (mcl_multi && identical(vergleich, "abweichend"))
             "Subpopulationen – mclust beachten"
           else if (identical(vergleich, "übereinstimmend")) "reflimR beibehalten"
           else if (identical(vergleich, "abweichend"))      "Abweichung – mclust prüfen"
           else if (mcl_multi)                               "mehrere Subpopulationen erkannt"
           else if (!is.null(mm) && isTRUE(mm$ok))           "eine Population (mclust)"
           else if (!is.null(v) && (v$lower == "red" || v$upper == "red")) "refineR-Verifikation empfohlen"
           else if (!is.null(v)) "reflimR konsistent"
           else "—"
    data.frame(Gruppe = e$label, reflimR = sprintf("%s – %s", e$lower, e$upper),
      refineR = refineR_txt, Vergleich = vergleich,
      mclust_Cluster = mcl_cluster, mclust_Hauptgruppe = mcl_haupt,
      Empfehlung = emp, stringsAsFactors = FALSE)
  }))
}

# colour dot for verdicts in the HTML table
.dot <- function(col) htmltools::span(style = sprintf(
  "display:inline-block;width:10px;height:10px;border-radius:50%%;background:%s;margin-right:6px;",
  COLOR_HEX[[col]]))

#' Build the full self-contained report as an htmltools tag structure.
build_report_html <- function(analyte, unit, results, min_bucket,
                              generated = Sys.time(), logo_path = NULL) {
  h <- htmltools::tags
  logo_uri <- if (!is.null(logo_path) && file.exists(logo_path))
    knitr::image_uri(logo_path) else NULL

  # ---- summary table ----
  sum_rows <- lapply(results, function(g) {
    e <- g$est; v <- g$verdict; rr <- g$refiner
    ri_txt  <- if (e$status == "ok") sprintf("%s – %s", e$lower, e$upper) else "—"
    cur_txt <- if (!is.null(g$current_lower) && !is.na(g$current_lower))
      sprintf("%s – %s", g$current_lower, g$current_upper) else "—"
    ref_txt <- if (!is.null(rr) && isTRUE(rr$ok)) sprintf("%s – %s", rr$lower, rr$upper) else "—"
    h$tr(
      h$td(e$label), h$td(e$n), h$td(ifelse(is.na(e$method), "—", e$method)),
      h$td(ri_txt), h$td(cur_txt), h$td(ref_txt),
      h$td(if (is.null(v)) "—" else htmltools::tagList(.dot(v$lower), COLOR_TEXT[[v$lower]])),
      h$td(if (is.null(v)) "—" else htmltools::tagList(.dot(v$upper), COLOR_TEXT[[v$upper]]))
    )
  })

  # ---- per-group detail sections ----
  detail <- lapply(results, function(g) {
    e <- g$est
    if (e$status != "ok") {
      return(h$div(h$h3(e$label), h$p(style = "color:#b02a37", e$note)))
    }
    tg <- if (!is.null(g$current_lower) && !is.na(g$current_lower) &&
              !is.null(g$current_upper) && !is.na(g$current_upper))
      c(g$current_lower, g$current_upper) else NULL

    xlb <- axis_label(analyte, unit); ttl <- plot_title(analyte, e$label)
    diag_uri <- plot_to_uri(if (e$method == "reflimLOD")
        plot_lod_group(e, tg, main = ttl, xlab = xlb)
      else
        suppressWarnings(reflimR::reflim(e$x, targets = tg, plot.it = TRUE,
          n.min = min_bucket, main = ttl, xlab = xlb)))

    # refineR plot beside reflimR when the model was computed
    ref_uri <- if (!is.null(g$refiner) && isTRUE(g$refiner$ok) && !is.null(g$refiner$model))
      plot_to_uri(plot(g$refiner$model, title = ttl, xlab = xlb)) else NULL
    plots <- if (is.null(ref_uri))
      h$img(src = diag_uri, style = "max-width:100%;")
    else
      h$div(style = "display:flex;gap:12px;flex-wrap:wrap;",
        h$div(style = "flex:1;min-width:280px;",
              h$strong("reflimR"), h$br(), h$img(src = diag_uri, style = "max-width:100%;")),
        h$div(style = "flex:1;min-width:280px;",
              h$strong("refineR"), h$br(), h$img(src = ref_uri, style = "max-width:100%;")))

    ref_block <- NULL
    if (!is.null(g$refiner) && isTRUE(g$refiner$ok)) {
      cmp <- compare_methods(e$lower, e$upper, g$refiner$lower, g$refiner$upper)
      cmp_txt <- function(c) c(green = "übereinstimmend", red = "abweichend")[[c]]
      ref_block <- h$p(h$strong("refineR: "),
        sprintf("%s – %s", g$refiner$lower, g$refiner$upper),
        if (!is.null(cmp)) sprintf("  (Vergleich: unten %s, oben %s)",
                                   cmp_txt(cmp$lower), cmp_txt(cmp$upper)))
    }

    mcl_block <- NULL
    if (!is.null(g$mclust) && isTRUE(g$mclust$ok)) {
      st <- g$mclust$stats
      mcl_uri <- plot_to_uri(suppressWarnings(lab_mclust(e$x, lognormal = e$lognormal,
                   remove.extremes = TRUE, targets = g$mclust$tg, plot.it = TRUE,
                   n.cluster = g$mclust$ncluster,
                   main = paste0(ttl, " (mclust)"), xlab = xlb)))
      crows <- lapply(seq_len(nrow(st)), function(k)
        h$tr(h$td(k), h$td(sprintf("%s – %s", st$ll[k], st$ul[k])),
             h$td(sprintf("%.0f %%", st$percent[k]))))
      mcl_block <- h$div(
        h$h4("mclust-Verifikation"),
        h$img(src = mcl_uri, style = "max-width:100%;"),
        h$table(class = "tbl", h$thead(h$tr(h$th("Cluster"), h$th("RI"), h$th("Anteil"))),
                h$tbody(crows)),
        h$p(h$em(mclust_sentence(st))))
    }

    h$div(class = "group",
      h$h3(e$label, h$span(class = "badge", e$method)),
      h$p(h$strong("reflimR-RI: "), sprintf("%s – %s", e$lower, e$upper),
          sprintf("  (90%%-KI unten %s–%s | oben %s–%s)",
                  e$lower.ci[1], e$lower.ci[2], e$upper.ci[1], e$upper.ci[2])),
      if (!is.null(e$note)) h$p(style = "color:#666", e$note),
      plots,
      ref_block, mcl_block)
  })

  css <- sprintf("
    body{font-family:Segoe UI,Arial,sans-serif;color:%s;max-width:900px;margin:24px auto;padding:0 16px;}
    h1{color:%s;border-bottom:3px solid %s;padding-bottom:8px;}
    h2,h3,h4{color:%s;} .group{border-top:1px solid #ddd;padding-top:12px;margin-top:24px;}
    table.tbl{border-collapse:collapse;margin:10px 0;font-size:.9rem;}
    table.tbl th,table.tbl td{border:1px solid #ccc;padding:4px 10px;text-align:left;}
    table.tbl th{background:%s;color:#fff;}
    .badge{background:%s;color:#fff;border-radius:12px;padding:2px 10px;margin-left:10px;font-size:.7rem;}
    .muted{color:#666;font-size:.85rem;}",
    IMD_CHARCOAL_R, IMD_TEAL_R, IMD_TEAL_R, IMD_TEAL_R, IMD_TEAL_R, IMD_TEAL_R)

  htmltools::tagList(
    h$head(h$meta(charset = "utf-8"), h$style(htmltools::HTML(css))),
    h$body(
      if (!is.null(logo_uri)) h$img(src = logo_uri, style = "height:42px;"),
      h$h1("Referenzintervall-Report"),
      h$p(h$strong("Analyt: "), ifelse(nzchar(analyte), analyte, "—"),
          h$strong("  Einheit: "), ifelse(nzchar(unit), unit, "—")),
      h$p(class = "muted", "Erstellt: ", format(generated, "%Y-%m-%d %H:%M"),
          " · reflimR ist die primäre Schätzung; refineR und mclust dienen der Verifikation."),
      h$h2("Übersicht"),
      h$table(class = "tbl",
        h$thead(h$tr(h$th("Gruppe"), h$th("n"), h$th("Methode"), h$th("reflimR-RI"),
                     h$th("aktuell"), h$th("refineR"), h$th("Ampel unten"), h$th("Ampel oben"))),
        h$tbody(sum_rows)),
      h$h2("Verifikations-Übersicht"),
      local({
        vdf <- build_verification_df(results)
        hdr <- c("Gruppe", "reflimR", "refineR", "reflimR↔refineR",
                 "mclust Cluster", "mclust Hauptgruppe", "Empfehlung")
        h$table(class = "tbl",
          h$thead(h$tr(lapply(hdr, h$th))),
          h$tbody(lapply(seq_len(nrow(vdf)), function(r)
            h$tr(lapply(vdf[r, ], function(cell) h$td(as.character(cell)))))))
      }),
      h$h2("Details je Gruppe"), detail,
      h$hr(),
      h$p(class = "muted",
        "Methoden: reflimR (Hoffmann et al.), refineR, mclust (Gaussian mixture). ",
        "Diese Auswertung ist eine indirekte Schätzung aus Routinedaten.")
    )
  )
}
