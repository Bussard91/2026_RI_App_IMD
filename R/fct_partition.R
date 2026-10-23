# ---------------------------------------------------------------------------
# Sex + age partitioning helpers (no Shiny dependencies). See SPEC.md section 4.
# Always split by sex first, then build one rpart age tree per sex.
# rpart parameters are hard-coded; min n per group defaults to 200.
# ---------------------------------------------------------------------------

#' Fit the anova rpart age tree for one sex. @return rpart model or NULL.
fit_age_tree <- function(ages, values, min_bucket = 200) {
  ok <- !is.na(ages) & !is.na(values)
  ages <- ages[ok]; values <- values[ok]
  if (length(ages) < 2 * min_bucket) return(NULL)
  ctrl <- rpart::rpart.control(minbucket = min_bucket,
                               minsplit  = 2 * min_bucket, cp = 0.01)
  rpart::rpart(values ~ ages, method = "anova", control = ctrl)
}

#' Fit one age tree per sex. @return named list(M = model/NULL, W = model/NULL).
fit_age_trees <- function(d, min_bucket = 200) {
  out <- list()
  for (s in c("M", "W")) {
    sub <- d[d$sex == s, , drop = FALSE]
    out[[s]] <- if (nrow(sub) > 0) fit_age_tree(sub$age, sub$value, min_bucket) else NULL
  }
  out
}

#' Extract age split points for one sex via an anova rpart tree.
#' @return numeric vector of split points (continuous age cut values); empty if
#'         too few observations or the tree found no split.
partition_one_sex <- function(ages, values, min_bucket = 200) {
  rp <- fit_age_tree(ages, values, min_bucket)
  if (is.null(rp) || is.null(rp$splits) || nrow(rp$splits) == 0) return(numeric(0))
  sort(unique(as.numeric(rp$splits[, "index"])))
}

#' Convert continuous split points into integer-year age groups.
#' A split at 14.5 yields ... to=14 | from=15 ... (ages are integer years).
#' @return data.frame(from, to) covering [min(ages), max(ages)]
groups_from_splits <- function(ages, splits) {
  amin <- min(ages, na.rm = TRUE)
  amax <- max(ages, na.rm = TRUE)
  if (length(splits) == 0) {
    return(data.frame(from = amin, to = amax))
  }
  cuts  <- sort(unique(as.integer(ceiling(splits))))  # lower bound of upper group
  cuts  <- cuts[cuts > amin & cuts <= amax]
  froms <- c(amin, cuts)
  tos   <- c(cuts - 1L, amax)
  data.frame(from = as.integer(froms), to = as.integer(tos))
}

#' Build the full suggested partition table over both sexes.
#' Tree is fit on measurable (non-LOD) values; n counts ALL rows in the group.
#' @param d cleaned data frame with columns age, sex, value, is_lod
#' @return data.frame(sex, from, to, n)
build_partition <- function(d, min_bucket = 200) {
  out <- list()
  for (s in c("M", "W")) {
    sub <- d[d$sex == s, , drop = FALSE]
    if (nrow(sub) == 0) next
    sp <- partition_one_sex(sub$age, sub$value, min_bucket)
    g  <- groups_from_splits(sub$age, sp)
    g$sex <- s
    out[[s]] <- g
  }
  if (length(out) == 0) {
    return(data.frame(sex = character(0), from = integer(0),
                      to = integer(0), n = integer(0)))
  }
  res <- do.call(rbind, out)
  res <- res[, c("sex", "from", "to")]
  res$n <- mapply(function(s, f, t) sum(d$sex == s & d$age >= f & d$age <= t,
                                        na.rm = TRUE),
                  res$sex, res$from, res$to)
  rownames(res) <- NULL
  res
}

#' Validate an (edited) group table against the data range.
#' Rules per sex: ascending, contiguous (no gaps / no overlaps), within range.
#' @return list(ok = logical, messages = character vector of problems)
validate_groups <- function(groups, d) {
  msg <- character(0)
  if (nrow(groups) == 0) return(list(ok = FALSE, messages = "No groups defined."))
  if (any(is.na(groups$from)) || any(is.na(groups$to)))
    msg <- c(msg, "All boundaries must be whole numbers.")
  if (any(groups$to < groups$from, na.rm = TRUE))
    msg <- c(msg, "Each group's 'to' must be >= its 'from'.")

  for (s in unique(groups$sex)) {
    g <- groups[groups$sex == s, , drop = FALSE]
    g <- g[order(g$from), , drop = FALSE]
    rng <- range(d$age[d$sex == s], na.rm = TRUE)
    if (nrow(g) >= 1) {
      if (min(g$from) > rng[1])
        msg <- c(msg, sprintf("%s: youngest age %d is not covered.", s, rng[1]))
      if (max(g$to) < rng[2])
        msg <- c(msg, sprintf("%s: oldest age %d is not covered.", s, rng[2]))
    }
    if (nrow(g) >= 2) {
      gaps <- g$from[-1] - g$to[-nrow(g)]
      if (any(gaps != 1))
        msg <- c(msg, sprintf("%s: groups must be contiguous (gap/overlap at %s).",
                              s, paste(g$to[-nrow(g)][gaps != 1], collapse = ", ")))
    }
  }
  list(ok = length(msg) == 0, messages = msg)
}

#' Pretty age-range label, e.g. "0-14", "61+".
group_label <- function(from, to, amax) {
  ifelse(to >= amax, paste0(from, "+"), paste0(from, "-", to))
}
