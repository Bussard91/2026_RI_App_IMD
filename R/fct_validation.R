# ---------------------------------------------------------------------------
# Pure data-validation & parsing helpers (no Shiny dependencies, unit-testable)
# Part of the IMD reference-interval app. See SPEC.md sections 3.4 and 6.
# ---------------------------------------------------------------------------

#' Parse a vector of raw measurement strings.
#'
#' Handles three forms found in IMD exports:
#'   - plain numeric        "13.3"  or German  "13,3"
#'   - below-LOD censored   "< 4.5" (LOD value parsed from the number after "<")
#'   - junk / empty         -> NA
#'
#' @return list with vectors: value (numeric, NA for sub-LOD), is_lod (logical),
#'         lod_value (numeric LOD per row, NA when not a "<x" entry)
parse_measurement <- function(raw) {
  raw <- trimws(as.character(raw))
  is_lod <- grepl("^<", raw)
  # number following "<"  ->  the analyte LOD
  lod_chr <- sub("^<\\s*", "", raw)
  lod_value <- suppressWarnings(as.numeric(gsub(",", ".", lod_chr)))
  lod_value[!is_lod] <- NA_real_
  # plain numeric value (NA for "<x" strings, which are not numeric)
  value <- suppressWarnings(as.numeric(gsub(",", ".", raw)))
  list(value = value, is_lod = is_lod, lod_value = lod_value)
}

#' Normalize a sex column to "M" / "W" (everything else -> NA).
#' IMD convention is M/W only; we accept simple case/lowercase variants.
normalize_sex <- function(x) {
  s <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(s))
  out[s %in% c("M", "MAN", "MANN", "MALE", "1")] <- "M"
  out[s %in% c("W", "F", "FRAU", "WOMAN", "WEIBLICH", "FEMALE", "2")] <- "W"
  out
}

#' Validate & clean an uploaded data frame given the user's column mapping.
#'
#' @param df        raw data frame (all columns read as character)
#' @param col_age   name of the age column
#' @param col_sex   name of the sex column
#' @param col_val   name of the value column
#' @return list(data, dropped, n_total, n_kept, per_sex, lod_value, n_lod, messages)
validate_data <- function(df, col_age, col_sex, col_val) {
  stopifnot(all(c(col_age, col_sex, col_val) %in% names(df)))

  age <- suppressWarnings(as.numeric(gsub(",", ".", as.character(df[[col_age]]))))
  age_int <- suppressWarnings(as.integer(round(age)))
  sex <- normalize_sex(df[[col_sex]])
  parsed <- parse_measurement(df[[col_val]])

  d <- data.frame(
    age       = age_int,
    sex       = sex,
    value     = parsed$value,
    is_lod    = parsed$is_lod,
    lod_value = parsed$lod_value,
    stringsAsFactors = FALSE
  )
  n_total <- nrow(d)

  drop_age <- is.na(d$age) | d$age < 0
  drop_sex <- is.na(d$sex)
  # a value is invalid only if it is NOT a sub-LOD entry and is non-numeric / <= 0
  drop_val <- (!d$is_lod) & (is.na(d$value) | d$value <= 0)
  keep <- !(drop_age | drop_sex | drop_val)

  dropped <- data.frame(
    reason = c("invalid/negative age", "unrecognized sex (not M/W)",
               "non-numeric or <= 0 value"),
    n = c(sum(drop_age), sum(drop_sex & !drop_age),
          sum(drop_val & !drop_age & !drop_sex)),
    stringsAsFactors = FALSE
  )

  clean <- d[keep, , drop = FALSE]
  rownames(clean) <- NULL

  # LOD detection: one analyte per run -> expect a single LOD value
  lod_vals <- unique(clean$lod_value[clean$is_lod])
  lod_vals <- lod_vals[!is.na(lod_vals)]
  lod_value <- if (length(lod_vals) >= 1) max(lod_vals) else NA_real_
  n_lod <- sum(clean$is_lod, na.rm = TRUE)

  per_sex <- as.data.frame(table(factor(clean$sex, levels = c("M", "W"))),
                           stringsAsFactors = FALSE)
  names(per_sex) <- c("sex", "n")

  list(
    data      = clean,
    dropped   = dropped,
    n_total   = n_total,
    n_kept    = nrow(clean),
    per_sex   = per_sex,
    lod_value = lod_value,
    n_lod     = n_lod,
    lod_share = if (nrow(clean) > 0) n_lod / nrow(clean) else NA_real_
  )
}
