# ===========================================================================
# IMD Reference-Interval App  -  app.R
# Stepwise wizard front-end for indirect RI estimation. See ../SPEC.md
# v1 scaffold: Step 1 (upload/map/validate) + Step 2 (partitioning) functional;
# Steps 3-4 are placeholders to be implemented next.
# ===========================================================================

library(shiny)
library(bslib)
library(DT)
library(rpart)
library(reflimR)
library(truncnorm)
library(refineR)
library(mclust)

source("R/fct_validation.R", local = TRUE)
source("R/fct_partition.R",  local = TRUE)
source("R/reflimLOD.R",      local = TRUE)
source("R/lab_mclust.R",     local = TRUE)
source("R/fct_estimate.R",   local = TRUE)
source("R/fct_report.R",     local = TRUE)

IMD_TEAL     <- "#005357"
IMD_CHARCOAL <- "#1D1D1B"
MIN_BUCKET   <- 200L   # hard-coded minimum n per age group (SPEC 4.1)
MAXG         <- 16L    # max groups for which we pre-register dynamic outputs

imd_theme <- bs_theme(version = 5, primary = IMD_TEAL, secondary = IMD_CHARCOAL)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_fluid(
  theme = imd_theme,
  tags$head(includeCSS("www/styles.css")),

  div(class = "app-header",
      img(src = "IMD.jpg", alt = "IMD Labor Berlin"),
      div(
        div(class = "title", "Referenzintervall-Tool"),
        div(class = "subtitle",
            "Indirekte Schätzung & Verifikation von Referenzintervallen")
      )
  ),

  tabsetPanel(
    id = "wizard", type = "pills",

    # ---- Step 1 -----------------------------------------------------------
    nav_panel(
      "1 · Upload", value = "step1",
      div(class = "step-card",
        div(class = "disclaimer-box",
          strong("Hinweis zum Datenschutz / Note on data privacy"), br(),
          "Diese öffentliche Version dient ausschließlich zu ",
          strong("Schulungszwecken"),
          ". Laden Sie nur ", strong("anonymisierte, aggregierte Daten"),
          " hoch (Alter, Geschlecht, Messwert – keine Namen, IDs oder ",
          "Probendaten). Daten werden nur im Arbeitsspeicher verarbeitet und ",
          "nicht gespeichert.",
          br(), br(),
          checkboxInput("accept_terms",
                        "Ich bestätige, nur anonymisierte Daten zu verwenden.",
                        value = FALSE)
        ),
        fileInput("file", "CSV-Datei (Spalten: Alter; Geschlecht; Messwert)",
                  accept = c(".csv", ".txt"), buttonLabel = "Durchsuchen",
                  placeholder = "Keine Datei ausgewählt"),
        conditionalPanel(
          "output.file_loaded == true",
          layout_columns(
            col_widths = c(4, 4, 4),
            selectInput("map_age", "Spalte: Alter", choices = NULL),
            selectInput("map_sex", "Spalte: Geschlecht", choices = NULL),
            selectInput("map_val", "Spalte: Messwert", choices = NULL)
          ),
          layout_columns(
            col_widths = c(6, 6),
            textInput("analyte", "Analyt / Parameter", placeholder = "z. B. Hämoglobin"),
            textInput("unit", "Einheit", placeholder = "z. B. g/dL")
          ),
          actionButton("validate", "Prüfen", class = "btn-primary"),
          hr(),
          uiOutput("validation_summary")
        )
      )
    ),

    # ---- Step 2 -----------------------------------------------------------
    nav_panel(
      "2 · Partitionierung", value = "step2",
      div(class = "step-card",
        uiOutput("step2_body")
      )
    ),

    # ---- Step 3 (placeholder) --------------------------------------------
    nav_panel(
      "3 · reflimR", value = "step3",
      div(class = "step-card",
        uiOutput("step3_body"))
    ),

    # ---- Step 4 ----------------------------------------------------------
    nav_panel(
      "4 · Verifikation", value = "step4",
      div(class = "step-card",
        uiOutput("step4_body"))
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    raw     = NULL,   # raw uploaded data (character)
    valid   = NULL,   # result of validate_data()
    groups  = NULL,   # suggested partition table (read-only)
    trees   = NULL    # per-sex rpart age trees (for scatter + tree image)
  )
  rv4 <- reactiveValues(refiner = list(), mclust = list())  # Step 4 results by group id

  # ---- Upload & column mapping ----
  observeEvent(input$file, {
    req(input$file)
    # read as character to preserve "<x" entries; sniff separator
    raw <- tryCatch(
      utils::read.csv(input$file$datapath, sep = ";", header = TRUE,
                      colClasses = "character", check.names = FALSE,
                      fileEncoding = "UTF-8"),
      error = function(e) NULL)
    if (is.null(raw) || ncol(raw) < 2) {
      raw <- utils::read.csv(input$file$datapath, sep = ",", header = TRUE,
                             colClasses = "character", check.names = FALSE)
    }
    rv$raw <- raw
    rv$valid <- NULL; rv$groups <- NULL
    cols <- names(raw)
    guess <- function(pat) { hit <- grep(pat, cols, ignore.case = TRUE)[1]
                             if (is.na(hit)) cols[1] else cols[hit] }
    updateSelectInput(session, "map_age", choices = cols, selected = guess("alter|age"))
    updateSelectInput(session, "map_sex", choices = cols, selected = guess("geschl|sex"))
    updateSelectInput(session, "map_val", choices = cols, selected = guess("mess|wert|value|result"))
  })

  output$file_loaded <- reactive({ !is.null(rv$raw) })
  outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)

  # ---- Validate & continue ----
  observeEvent(input$validate, {
    if (!isTRUE(input$accept_terms)) {
      showNotification("Bitte bestätigen Sie zuerst den Datenschutzhinweis.",
                       type = "error"); return()
    }
    req(rv$raw)
    cols <- c(input$map_age, input$map_sex, input$map_val)
    if (anyDuplicated(cols)) {
      showNotification("Bitte drei unterschiedliche Spalten zuordnen.",
                       type = "error"); return()
    }
    res <- tryCatch(
      validate_data(rv$raw, input$map_age, input$map_sex, input$map_val),
      error = function(e) { showNotification(paste("Fehler:", conditionMessage(e)),
                                             type = "error"); NULL })
    req(res)
    rv$valid <- res
    rv$trees <- fit_age_trees(res$data, min_bucket = MIN_BUCKET)
    rv$groups <- build_partition(res$data, min_bucket = MIN_BUCKET)
  })

  observeEvent(input$to_step2, {
    req(rv$valid)
    updateTabsetPanel(session, "wizard", selected = "step2")
  })

  output$validation_summary <- renderUI({
    res <- rv$valid; req(res)
    lod_txt <- if (is.na(res$lod_value)) "keine <LOD-Werte erkannt"
      else sprintf("LOD = %s erkannt, %d Werte < LOD (%.1f%%)",
                   res$lod_value, res$n_lod, 100 * res$lod_share)
    tagList(
      h5("Validierung"),
      span(class = "summary-pill", sprintf("%d Zeilen gesamt", res$n_total)),
      span(class = "summary-pill", sprintf("%d gültig", res$n_kept)),
      span(class = "summary-pill",
           sprintf("M: %s | W: %s",
                   res$per_sex$n[res$per_sex$sex == "M"],
                   res$per_sex$n[res$per_sex$sex == "W"])),
      span(class = if (!is.na(res$lod_value)) "summary-pill warn-pill" else "summary-pill",
           lod_txt),
      if (sum(res$dropped$n) > 0) tagList(
        p(style = "margin-top:10px", sprintf("%d Zeilen entfernt:", sum(res$dropped$n))),
        tags$ul(lapply(which(res$dropped$n > 0), function(i)
          tags$li(sprintf("%d – %s", res$dropped$n[i], res$dropped$reason[i]))))),
      hr(),
      actionButton("to_step2", "Weiter →", class = "btn-primary")
    )
  })

  # ---- Step 2: partitioning ----
  output$step2_body <- renderUI({
    if (is.null(rv$groups)) {
      return(div(class = "locked-note",
                 "Bitte zuerst Schritt 1 abschließen (Daten hochladen & prüfen)."))
    }
    sexes <- intersect(c("M", "W"), unique(rv$groups$sex))
    plot_sections <- lapply(sexes, function(s) {
      tagList(
        h6(paste0("Alterstrend & Entscheidungsbaum – ",
                  ifelse(s == "M", "Männer", "Frauen"))),
        layout_columns(col_widths = c(7, 5),
          plotOutput(paste0("trend_", s), height = "300px"),
          plotOutput(paste0("tree_", s), height = "300px"))
      )
    })
    tagList(
      h5("Vorgeschlagene Gruppen"),
      p("Geschlecht wird immer getrennt; je Geschlecht wird automatisch ein ",
        "Alters-Entscheidungsbaum berechnet. Die vorgeschlagenen Altersgruppen:"),
      DTOutput("group_table"),
      hr(),
      p(class = "text-muted",
        "Streudiagramm: Alter vs. Messwert (Punkte). Grün = Mittelwert je ",
        "Altersgruppe (entspricht den Werten im Entscheidungsbaum), getrennt ",
        "durch rote Gruppengrenzen. Rechts der zugehörige Entscheidungsbaum."),
      do.call(tagList, plot_sections),
      hr(),
      actionButton("to_step3", "Weiter zu reflimR →", class = "btn-primary")
    )
  })

  # Static, read-only table of the suggested partitions.
  output$group_table <- renderDT({
    req(rv$groups)
    datatable(
      rv$groups, rownames = FALSE, selection = "none",
      options = list(dom = "t", paging = FALSE, ordering = FALSE),
      colnames = c("Geschlecht", "von (Jahre)", "bis (Jahre)", "n")
    )
  }, server = FALSE)

  # per-sex age-trend scatter + decision-tree image (Step 2)
  for (s in c("M", "W")) local({
    ss <- s
    output[[paste0("trend_", ss)]] <- renderPlot({
      req(rv$valid, rv$groups)
      gs <- rv$groups[rv$groups$sex == ss, , drop = FALSE]
      req(nrow(gs) > 0)
      d <- rv$valid$data
      plot_age_trend(d[d$sex == ss, , drop = FALSE], gs,
        main = plot_title(input$analyte, paste0("Alterstrend ", ss)),
        xlab = "Alter (Jahre)", ylab = axis_label(input$analyte, input$unit))
    })
    output[[paste0("tree_", ss)]] <- renderPlot({
      req(rv$trees)
      plot_age_tree(rv$trees[[ss]], main = paste0("Entscheidungsbaum ", ss))
    })
  })

  observeEvent(input$to_step3, {
    req(rv$groups)
    updateTabsetPanel(session, "wizard", selected = "step3")
  })

  # ---- Step 3: reflimR estimation & comparison ----
  estimates <- reactive({
    req(rv$groups, rv$valid)
    d <- rv$valid$data; lodv <- rv$valid$lod_value
    lapply(seq_len(nrow(rv$groups)), function(i) {
      g <- rv$groups[i, ]
      amax <- suppressWarnings(max(d$age[d$sex == g$sex], na.rm = TRUE))
      sub <- d[d$sex == g$sex & d$age >= g$from & d$age <= g$to, ]
      est <- estimate_group(sub$value, sum(sub$is_lod, na.rm = TRUE),
                            lodv, n_min = MIN_BUCKET)
      est$id  <- paste0("g", i)
      est$sex <- g$sex; est$from <- g$from; est$to <- g$to
      est$label <- paste0(g$sex, " ", group_label(g$from, g$to, amax))
      est
    })
  })

  output$step3_body <- renderUI({
    if (is.null(rv$groups)) {
      return(div(class = "locked-note", "Bitte zuerst Schritte 1–2 abschließen."))
    }
    ests <- estimates()
    cards <- lapply(ests, function(e) {
      if (e$status != "ok") {
        return(card(card_header(strong(e$label)),
                    div(style = "color:#b02a37", e$note)))
      }
      pid <- paste0("plot_", e$id);   vid <- paste0("verdict_", e$id)
      clo <- paste0("curr_lo_", e$id); cup <- paste0("curr_up_", e$id)
      method_badge <- span(
        class = if (e$method == "reflimLOD") "summary-pill warn-pill" else "summary-pill",
        e$method)
      warn <- if (e$n < MIN_BUCKET)
        span(class = "summary-pill warn-pill", sprintf("n=%d < %d", e$n, MIN_BUCKET))
      card(
        card_header(strong(e$label), " ", method_badge, " ", warn),
        layout_columns(
          col_widths = c(7, 5),
          plotOutput(pid, height = "240px"),
          div(
            p(strong("Geschätztes RI: "), sprintf("%s – %s", e$lower, e$upper)),
            p(class = "text-muted",
              sprintf("90%%-KI unten: %s–%s | oben: %s–%s",
                      e$lower.ci[1], e$lower.ci[2], e$upper.ci[1], e$upper.ci[2])),
            if (!is.null(e$note)) p(class = "text-muted", e$note),
            hr(),
            p(strong("Aktuelle Laborgrenzen:")),
            layout_columns(col_widths = c(6, 6),
                           numericInput(clo, "untere", value = NA),
                           numericInput(cup, "obere", value = NA)),
            uiOutput(vid)
          )
        )
      )
    })
    tagList(
      h5("reflimR-Schätzung & Vergleich"),
      p("Geben Sie je Gruppe Ihre aktuellen Referenzgrenzen ein, um die ",
        "Ampel-Bewertung zu erhalten. reflimR ist der primäre Schätzer."),
      do.call(tagList, cards),
      hr(),
      h6("Übersicht"),
      DTOutput("comparison_table"),
      br(),
      actionButton("to_step4", "Weiter zur Verifikation →", class = "btn-primary")
    )
  })

  # Fixed pool of per-group outputs/observers (avoids dynamic re-registration).
  # Each slot i reads the i-th estimate at runtime; unused slots stay inert.
  for (i in seq_len(MAXG)) local({
    ii  <- i
    gid <- paste0("g", ii)
    pid <- paste0("plot_g", ii);    vid <- paste0("verdict_g", ii)
    clo <- paste0("curr_lo_g", ii); cup <- paste0("curr_up_g", ii)
    rid <- paste0("run_ref_g", ii); mid <- paste0("run_mcl_g", ii)
    refout <- paste0("ref_out_g", ii); mout <- paste0("mcl_out_g", ii)
    mplot  <- paste0("mcl_plot_g", ii)
    s4r    <- paste0("s4_reflim_g", ii);  s4ref <- paste0("s4_refiner_g", ii)

    cur_est <- function() { es <- estimates(); if (ii <= length(es)) es[[ii]] else NULL }
    targets_now <- function() {
      cl <- input[[clo]]; cu <- input[[cup]]
      if (!is.null(cl) && !is.null(cu) && !is.na(cl) && !is.na(cu)) c(cl, cu) else NULL
    }

    # --- Step 3: diagnostic plot + verdict ---
    output[[pid]] <- renderPlot({
      e <- cur_est(); req(e, e$status == "ok")
      ttl <- plot_title(input$analyte, e$label); xlb <- axis_label(input$analyte, input$unit)
      if (e$method == "reflimLOD") plot_lod_group(e, targets_now(), main = ttl, xlab = xlb)
      else suppressWarnings(reflimR::reflim(e$x, targets = targets_now(),
             plot.it = TRUE, n.min = MIN_BUCKET, main = ttl, xlab = xlb))
    })
    output[[vid]] <- renderUI({
      e <- cur_est(); req(e, e$status == "ok")
      v <- group_verdict(e, input[[clo]], input[[cup]], n_min = MIN_BUCKET)
      if (is.null(v)) return(span(class = "text-muted", "Grenzen eingeben für Bewertung."))
      verdict_badges(v)
    })

    # --- Step 4: refineR ---
    observeEvent(input[[rid]], {
      e <- cur_est(); req(e, e$status == "ok")
      withProgress(message = "refineR läuft – das kann einige Minuten dauern …",
                   value = NULL, {
        mod <- tryCatch(refineR::findRI(Data = e$x, model = "BoxCox"),
                        error = function(err) NULL)
      })
      if (is.null(mod)) { rv4$refiner[[gid]] <- list(ok = FALSE); return() }
      ri <- refineR::getRI(mod, RIperc = c(0.025, 0.975))
      rv4$refiner[[gid]] <- list(ok = TRUE,
                                 lower = round_lim(ri$PointEst[1]),
                                 upper = round_lim(ri$PointEst[2]),
                                 model = mod)
    })
    # side-by-side diagnostic plots on the verification page
    output[[s4r]] <- renderPlot({
      e <- cur_est(); req(e, e$status == "ok")
      ttl <- plot_title(input$analyte, e$label); xlb <- axis_label(input$analyte, input$unit)
      if (e$method == "reflimLOD") plot_lod_group(e, targets_now(), main = ttl, xlab = xlb)
      else suppressWarnings(reflimR::reflim(e$x, targets = targets_now(),
             plot.it = TRUE, n.min = MIN_BUCKET, main = ttl, xlab = xlb))
    })
    output[[s4ref]] <- renderPlot({
      rr <- rv4$refiner[[gid]]; req(rr, isTRUE(rr$ok), !is.null(rr$model))
      e <- cur_est()
      plot(rr$model, title = plot_title(input$analyte, e$label),
           xlab = axis_label(input$analyte, input$unit))
    })
    output[[refout]] <- renderUI({
      rr <- rv4$refiner[[gid]]; if (is.null(rr)) return(NULL)
      if (!isTRUE(rr$ok)) return(div(style = "color:#b02a37",
                                     "refineR-Schätzung fehlgeschlagen."))
      e <- cur_est()
      cmp <- compare_methods(e$lower, e$upper, rr$lower, rr$upper)
      differ <- !is.null(cmp) && (cmp$lower == "red" || cmp$upper == "red")
      tagList(
        hr(),
        layout_columns(col_widths = c(6, 6),
          div(strong("reflimR"), plotOutput(s4r, height = "240px")),
          div(strong("refineR"), plotOutput(s4ref, height = "240px"))),
        p(strong("refineR-RI: "), sprintf("%s – %s", rr$lower, rr$upper)),
        if (!is.null(cmp)) tagList(p("Vergleich reflimR ↔ refineR:"), cmp_badges(cmp)),
        br(),
        if (differ) div(class = "summary-pill warn-pill",
              "Methoden weichen ab – mclust-Verifikation empfohlen (Button unten)")
        else div(class = "summary-pill",
                 "reflimR und refineR stimmen überein – reflimR wird beibehalten.")
      )
    })

    # --- Step 4: mclust ---
    observeEvent(input[[mid]], {
      e <- cur_est(); req(e, e$status == "ok")
      gsel <- input[[paste0("mcl_g_g", ii)]]
      ncl <- if (is.null(gsel) || is.na(gsel) || gsel < 1) NULL else as.integer(gsel)
      res <- tryCatch(suppressWarnings(lab_mclust(e$x, lognormal = e$lognormal,
               remove.extremes = TRUE, targets = targets_now(), plot.it = FALSE,
               n.cluster = ncl)),
               error = function(err) NULL)
      rv4$mclust[[gid]] <- if (is.null(res)) list(ok = FALSE)
        else list(ok = TRUE, stats = res$stats, tg = targets_now(), ncluster = ncl)
    })
    output[[mplot]] <- renderPlot({
      mm <- rv4$mclust[[gid]]; req(mm, isTRUE(mm$ok))
      e <- cur_est()
      suppressWarnings(lab_mclust(e$x, lognormal = e$lognormal, remove.extremes = TRUE,
                                  targets = mm$tg, plot.it = TRUE,
                                  main = paste0(plot_title(input$analyte, e$label), " (mclust)"),
                                  xlab = axis_label(input$analyte, input$unit),
                                  n.cluster = mm$ncluster))
    })
    output[[mout]] <- renderUI({
      mm <- rv4$mclust[[gid]]; if (is.null(mm)) return(NULL)
      if (!isTRUE(mm$ok)) return(div(style = "color:#b02a37",
                                     "mclust-Schätzung fehlgeschlagen."))
      st <- mm$stats
      rows <- lapply(seq_len(nrow(st)), function(k)
        tags$tr(tags$td(k), tags$td(sprintf("%s – %s", st$ll[k], st$ul[k])),
                tags$td(sprintf("%.0f %%", st$percent[k]))))
      tagList(
        hr(),
        plotOutput(mplot, height = "240px"),
        tags$table(class = "table table-sm", style = "max-width:360px",
          tags$thead(tags$tr(tags$th("Cluster"), tags$th("RI"), tags$th("Anteil"))),
          tags$tbody(rows)),
        p(tags$em(mclust_sentence(st)))
      )
    })
  })

  output$comparison_table <- renderDT({
    ests <- estimates()
    rows <- lapply(ests, function(e) {
      cl <- input[[paste0("curr_lo_", e$id)]]
      cu <- input[[paste0("curr_up_", e$id)]]
      v  <- if (e$status == "ok") group_verdict(e, cl, cu, MIN_BUCKET) else NULL
      data.frame(
        Gruppe = e$label, n = e$n,
        Methode = ifelse(is.na(e$method), "—", e$method),
        `RI unten` = e$lower, `RI oben` = e$upper,
        `aktuell unten` = ifelse(is.null(cl) || is.na(cl), NA, cl),
        `aktuell oben`  = ifelse(is.null(cu) || is.na(cu), NA, cu),
        `Ampel unten` = if (is.null(v)) "—" else unname(COLOR_TEXT[[v$lower]]),
        `Ampel oben`  = if (is.null(v)) "—" else unname(COLOR_TEXT[[v$upper]]),
        check.names = FALSE, stringsAsFactors = FALSE)
    })
    datatable(do.call(rbind, rows), rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE))
  })

  observeEvent(input$to_step4, {
    updateTabsetPanel(session, "wizard", selected = "step4")
  })

  # ---- Step 4: verification UI ----
  output$step4_body <- renderUI({
    if (is.null(rv$groups)) {
      return(div(class = "locked-note", "Bitte zuerst Schritte 1–3 abschließen."))
    }
    ests <- estimates()
    cards <- lapply(seq_along(ests), function(i) {
      e <- ests[[i]]
      if (e$status != "ok") {
        return(card(card_header(strong(e$label)), div(style = "color:#b02a37", e$note)))
      }
      clo <- paste0("curr_lo_g", i); cup <- paste0("curr_up_g", i)
      rid <- paste0("run_ref_g", i)
      v   <- group_verdict(e, input[[clo]], input[[cup]], MIN_BUCKET)
      red <- !is.null(v) && (v$lower == "red" || v$upper == "red")
      rec <- if (is.null(v))
               span(class = "text-muted", "Aktuelle Grenzen in Schritt 3 eingeben.")
             else if (red)
               div(class = "summary-pill warn-pill",
                   "reflimR weicht von den aktuellen Grenzen ab – refineR-Verifikation empfohlen")
             else
               div(class = "summary-pill",
                   "reflimR ist konsistent – Verifikation optional")
      card(
        card_header(strong(e$label), " ", span(class = "summary-pill", e$method)),
        p(strong("reflimR-RI: "), sprintf("%s – %s", e$lower, e$upper)),
        rec, br(), br(),
        actionButton(rid, "refineR ausführen",
                     class = if (red) "btn-primary" else "btn-outline-secondary"),
        uiOutput(paste0("ref_out_g", i)),
        hr(),
        p(class = "text-muted",
          "mclust-Verifikation jederzeit möglich (prüft auf versteckte Subpopulationen):"),
        layout_columns(col_widths = c(5, 7),
          numericInput(paste0("mcl_g_g", i), "Cluster (leer = automatisch)",
                       value = NA, min = 1, max = 9, step = 1),
          div(style = "margin-top:28px;",
              actionButton(paste0("run_mcl_g", i), "mclust ausführen",
                           class = "btn-outline-secondary"))),
        uiOutput(paste0("mcl_out_g", i))
      )
    })
    tagList(
      h5("Verifikation mit refineR & mclust"),
      p("Bei Abweichungen zwischen reflimR und Ihren aktuellen Grenzen empfiehlt ",
        "sich eine Verifikation mit refineR; weichen reflimR und refineR ",
        "voneinander ab, dient mclust als Kontrolle. Diese Berechnungen starten ",
        "Sie bewusst per Knopfdruck (refineR kann einige Minuten dauern)."),
      do.call(tagList, cards),
      hr(),
      h6("Verifikations-Ergebnisse"),
      p(class = "text-muted",
        "Zusammenfassung je Gruppe: reflimR, refineR (falls ausgeführt), Vergleich, ",
        "mclust-Cluster und Empfehlung."),
      DTOutput("verification_table"),
      hr(),
      h6("Export"),
      p(class = "text-muted",
        "Der Report enthält alle Gruppen, Tabellen und Diagnostik-Plots ",
        "(eigenständige HTML-Datei; PDF per Browser-Druck)."),
      downloadButton("dl_report", "Report (HTML) herunterladen", class = "btn-primary"),
      downloadButton("dl_csv", "Grenzen als CSV")
    )
  })

  # ---- assembled results (shared by report + CSV) ----
  build_results <- reactive({
    ests <- estimates()
    lapply(seq_along(ests), function(i) {
      e  <- ests[[i]]
      cl <- input[[paste0("curr_lo_g", i)]]
      cu <- input[[paste0("curr_up_g", i)]]
      list(est = e, current_lower = cl, current_upper = cu,
           verdict = if (e$status == "ok") group_verdict(e, cl, cu, MIN_BUCKET) else NULL,
           refiner = rv4$refiner[[paste0("g", i)]],
           mclust  = rv4$mclust[[paste0("g", i)]])
    })
  })

  output$verification_table <- renderDT({
    df <- build_verification_df(build_results())
    colnames(df) <- c("Gruppe", "reflimR", "refineR", "reflimR↔refineR",
                      "mclust Cluster", "mclust Hauptgruppe", "Empfehlung")
    datatable(df, rownames = FALSE, selection = "none",
              options = list(dom = "t", paging = FALSE, ordering = FALSE))
  })

  safe_name <- function() {
    a <- input$analyte; if (is.null(a) || !nzchar(a)) a <- "Analyt"
    gsub("[^A-Za-z0-9_-]", "_", a)
  }

  output$dl_report <- downloadHandler(
    filename = function() sprintf("RI_Report_%s_%s.html", safe_name(),
                                  format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      html <- build_report_html(
        analyte = input$analyte %||% "", unit = input$unit %||% "",
        results = build_results(), min_bucket = MIN_BUCKET,
        generated = Sys.time(), logo_path = "www/IMD.jpg")
      htmltools::save_html(html, file, background = "white")
    })

  output$dl_csv <- downloadHandler(
    filename = function() sprintf("RI_Grenzen_%s_%s.csv", safe_name(),
                                  format(Sys.Date(), "%Y%m%d")),
    content = function(file) {
      df <- build_results_df(build_results(), analyte = input$analyte %||% "",
                             unit = input$unit %||% "")
      utils::write.csv2(df, file, row.names = FALSE, fileEncoding = "UTF-8")
    })
}

`%||%` <- function(a, b) if (is.null(a)) b else a

verdict_badges <- function(v) {
  badge <- function(lab, col) span(style = sprintf(
    "display:inline-block;background:%s;color:#fff;border-radius:14px;padding:3px 10px;margin-right:6px;font-size:.8rem;",
    COLOR_HEX[[col]]), sprintf("%s: %s", lab, COLOR_TEXT[[col]]))
  tagList(badge("untere", v$lower), badge("obere", v$upper))
}

# badges for reflimR-vs-refineR method comparison (übereinstimmend / abweichend)
cmp_badges <- function(cmp) {
  txt <- c(green = "übereinstimmend", red = "abweichend")
  badge <- function(lab, col) span(style = sprintf(
    "display:inline-block;background:%s;color:#fff;border-radius:14px;padding:3px 10px;margin-right:6px;font-size:.8rem;",
    COLOR_HEX[[col]]), sprintf("%s: %s", lab, txt[[col]]))
  tagList(badge("untere", cmp$lower), badge("obere", cmp$upper))
}

shinyApp(ui, server)
