# =============================================================================
# import_legacy.R -- Family A importer: "OPIHI data YYYY.xlsx" (2004-2020)
#
# LAYOUT OF THESE FILES (one tab per file)
#   Each field trip is a vertical BLOCK of columns, blocks sit side by side:
#     rows 1-11  : metadata labels + values (Date, Start Time, ..., School)
#     row 12     : running quadrat number (1, 2, 3 ...)
#     row 13     : "Organisms" | "transect-1 quadrat-1 at point 0.00" | ...
#     rows 14+   : one row per taxon, one column per quadrat, cell = points
#                  hit out of 25 (blank = not detected)
#     below that : a "TOTALS:" section the spreadsheet computed itself
#
# KNOWN QUIRKS HANDLED (all found by inspecting the raw files)
#   * Blocks are found by the "Organisms" header, NOT the "Date" label, because
#     one 2017 trip's Date label does not match exactly and would merge into
#     its neighbour.
#   * In some 2019 blocks the Location label/value sit 2 columns right of the
#     others, so each label is searched for by NAME within its block.
#   * Two header styles: "transect-1 quadrat-1 at point 0.00" / "T-1 Q-1 at ...".
#   * One 2017 quadrat column has data but a blank header. Its transect/quadrat
#     are inferred from row 12 + neighbouring headers and flagged "inferred".
#   * One 2020 block has no metadata at all (flagged; excluded downstream).
# =============================================================================

LEGACY_META_LABELS <- c(
  date_raw            = "^date",
  start_time_raw      = "^start",
  stop_time_raw       = "^stop",
  n_transects_raw     = "^number of tran",
  transect_length_raw = "^length",
  point_spacing_raw   = "^point",
  n_quad_per_t_raw    = "^number of quad",
  site_raw            = "^location",
  latitude_raw        = "^lat",
  longitude_raw       = "^lon",
  school_raw          = "^school"
)

LEGACY_HEADER_RE <- stringr::regex(
  "^t(?:ransect)?\\s*-?\\s*(\\d+)\\s*q(?:uadrat)?\\s*-?\\s*(\\d+)\\s*(?:at\\s*point\\s*(-?[\\d.]+))?",
  ignore_case = TRUE
)

# -----------------------------------------------------------------------------
# parse_legacy_file()
# Input : path to one legacy .xlsx
# Output: list of tibbles
#   surveys  - one row per field-trip block (raw metadata kept as text)
#   quadrats - one row per quadrat column (sampling effort, incl. empty ones)
#   counts   - one row per non-empty cell: block x quadrat x taxon x count
#   totals   - the sheet's own TOTALS section (for verification only)
#   issues   - anything odd, with location, never silently dropped
# -----------------------------------------------------------------------------
parse_legacy_file <- function(path) {
  fn <- basename(path)
  m  <- read_sheet_matrix(path, 1)

  org_row <- which(apply(m, 1, function(r) any(tolower(r) == "organisms", na.rm = TRUE)))[1]
  if (is.na(org_row)) stop("No 'Organisms' header row found in ", fn)
  meta_rows <- (org_row - 13) + seq_along(LEGACY_META_LABELS)   # the 11 rows above the quadrat-number row (Excel rows 1-11)
  num_row   <- m[org_row - 1, ]
  hdr       <- m[org_row, ]

  starts <- which(tolower(hdr) == "organisms")
  ends   <- c(starts[-1] - 1, ncol(m))

  surveys <- list(); quadrats <- list(); counts <- list()
  totals  <- list(); issues   <- list()
  add_issue <- function(id, type, detail) {
    issues[[length(issues) + 1]] <<- tibble::tibble(block_id = id, issue = type, detail = detail)
  }

  for (k in seq_along(starts)) {
    j0 <- starts[k]; j1 <- ends[k]
    block_id <- sprintf("%s|block%02d", fn, k)

    # ---- metadata: find each label by name inside the block -----------------
    search_cols <- max(1, j0 - 3):j1
    meta <- purrr::imap_chr(LEGACY_META_LABELS, function(re, nm) {
      vals <- m[meta_rows[which(names(LEGACY_META_LABELS) == nm)], search_cols]
      hit  <- which(stringr::str_detect(tolower(vals), re))[1]
      if (is.na(hit) || hit >= length(vals)) return(NA_character_)
      after <- vals[(hit + 1):length(vals)]
      after <- after[!is.na(after)]
      if (length(after) == 0) NA_character_ else after[1]
    })

    # ---- quadrat columns ----------------------------------------------------
    if (j1 <= j0) { add_issue(block_id, "no_quadrat_columns", ""); next }
    q <- tibble::tibble(col = (j0 + 1):j1, header = hdr[(j0 + 1):j1], qnum = num_row[(j0 + 1):j1])
    mt <- stringr::str_match(q$header, LEGACY_HEADER_RE)
    q$transect <- as.integer(mt[, 2])
    q$quadrat  <- as.integer(mt[, 3])
    q$point_m  <- suppressWarnings(as.numeric(mt[, 4]))
    q$header_status <- dplyr::case_when(
      !is.na(q$transect) ~ "ok",
      is.na(q$header)    ~ "blank",
      TRUE               ~ "unparsed"
    )
    for (r in which(q$header_status == "unparsed")) {
      add_issue(block_id, "unparsed_header", paste0("col ", q$col[r], ": ", q$header[r]))
    }

    # ---- taxon rows: from below "Organisms" down to the first blank ---------
    below   <- m[(org_row + 1):nrow(m), j0]
    first_na <- which(is.na(below))[1]
    n_sp    <- if (is.na(first_na)) length(below) else first_na - 1
    sp_rows <- if (n_sp > 0) (org_row + 1):(org_row + n_sp) else integer(0)
    if (n_sp == 0) add_issue(block_id, "no_taxa_listed", "")

    # ---- recover blank headers that still contain data ----------------------
    if (n_sp > 0) {
      has_data <- colSums(!is.na(m[sp_rows, q$col, drop = FALSE])) > 0
      for (r in which(q$header_status == "blank" & has_data)) {
        qn <- suppressWarnings(as.integer(q$qnum[r]))
        if (is.na(qn)) { add_issue(block_id, "blank_header_with_data", paste("col", q$col[r])); next }
        nxt <- which(q$header_status == "ok" & q$col > q$col[r])[1]
        prv <- utils::tail(which(q$header_status == "ok" & q$col < q$col[r]), 1)
        t_inf <- if (!is.na(nxt) && q$quadrat[nxt] > qn) q$transect[nxt] else
                 if (length(prv) == 1) q$transect[prv] else NA_integer_
        q$transect[r] <- t_inf; q$quadrat[r] <- qn; q$header_status[r] <- "inferred"
        add_issue(block_id, "header_inferred",
                  sprintf("col %d -> transect %s quadrat %d", q$col[r], t_inf, qn))
      }
    }
    q_use <- dplyr::filter(q, header_status %in% c("ok", "inferred"))

    # ---- counts: matrix -> long --------------------------------------------
    if (n_sp > 0 && nrow(q_use) > 0) {
      vals <- m[sp_rows, q_use$col, drop = FALSE]
      long <- tibble::tibble(
        block_id  = block_id,
        col       = rep(q_use$col, each = length(sp_rows)),
        excel_row = rep(sp_rows, times = ncol(vals)),
        taxon_raw = rep(m[sp_rows, j0], times = ncol(vals)),
        value_raw = as.vector(vals)
      ) |>
        dplyr::filter(!is.na(value_raw)) |>
        dplyr::mutate(count = suppressWarnings(as.numeric(value_raw)))
      bad <- dplyr::filter(long, is.na(count))
      for (r in seq_len(nrow(bad))) {
        add_issue(block_id, "non_numeric_count",
                  sprintf("row %d col %d: %s = '%s'", bad$excel_row[r], bad$col[r], bad$taxon_raw[r], bad$value_raw[r]))
      }
      counts[[k]] <- long |>
        dplyr::filter(!is.na(count)) |>
        dplyr::left_join(dplyr::select(q_use, col, transect, quadrat, point_m), by = "col")
    }

    # ---- the sheet's own TOTALS section (verification only) -----------------
    tot_row <- which(stringr::str_detect(m[, j0], "^(?i)totals"))[1]
    if (!is.na(tot_row) && tot_row < nrow(m)) {
      below   <- m[(tot_row + 1):nrow(m), j0]
      n_tot   <- if (is.na(which(is.na(below))[1])) length(below) else which(is.na(below))[1] - 1
      if (n_tot > 0 && j0 < ncol(m)) {
        rr <- (tot_row + 1):(tot_row + n_tot)
        totals[[k]] <- tibble::tibble(block_id = block_id, taxon_raw = m[rr, j0],
                                      sheet_total = suppressWarnings(as.numeric(m[rr, j0 + 1])))
      }
    }

    quadrats[[k]] <- dplyr::transmute(q_use, block_id = block_id, transect, quadrat, point_m,
                                      header_status, source_col = col)
    surveys[[k]]  <- tibble::as_tibble(as.list(meta)) |>
      dplyr::mutate(block_id = block_id, source_file = fn, block = k, source_col = j0,
                    n_taxa_listed = n_sp, n_quadrat_cols = nrow(q_use), .before = 1)
    if (all(is.na(meta[c("date_raw", "site_raw", "school_raw")]))) {
      add_issue(block_id, "no_survey_metadata", "date, location and school all missing")
    }
  }

  list(
    surveys  = dplyr::bind_rows(surveys),
    quadrats = dplyr::bind_rows(quadrats),
    counts   = dplyr::bind_rows(counts),
    totals   = dplyr::bind_rows(totals),
    issues   = dplyr::bind_rows(issues)
  )
}

# -----------------------------------------------------------------------------
# import_legacy()
# Runs parse_legacy_file() on every legacy file in raw_dir, stacks results and
# converts metadata to proper types (raw text columns are kept alongside).
# -----------------------------------------------------------------------------
import_legacy <- function(raw_dir = here::here("Data", "Raw"),
                          pattern = "^OPIHI data (19|20)\\d{2}.*\\.xlsx$") {
  files <- list.files(raw_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No legacy files matched in ", raw_dir)
  parts <- purrr::map(files, parse_legacy_file)
  out   <- purrr::map(purrr::set_names(names(parts[[1]])),
                      function(nm) dplyr::bind_rows(purrr::map(parts, nm)))

  out$surveys <- out$surveys |>
    dplyr::mutate(
      date            = excel_to_date(date_raw),
      start_time      = excel_to_time(start_time_raw),
      stop_time       = excel_to_time(stop_time_raw),
      n_transects     = suppressWarnings(as.integer(n_transects_raw)),
      transect_length_m = suppressWarnings(as.numeric(stringr::str_extract(transect_length_raw, "[\\d.]+"))),
      point_spacing_m = suppressWarnings(as.numeric(stringr::str_extract(point_spacing_raw, "[\\d.]+"))),
      n_quad_per_t    = suppressWarnings(as.integer(n_quad_per_t_raw)),
      latitude        = suppressWarnings(as.numeric(latitude_raw)),
      longitude       = suppressWarnings(as.numeric(longitude_raw))
    )
  out
}

# -----------------------------------------------------------------------------
# check_legacy_totals()
# Verification: re-sum our parsed counts per block x taxon and compare with the
# TOTALS section each spreadsheet computed itself. 0 rows returned = perfect.
# -----------------------------------------------------------------------------
check_legacy_totals <- function(legacy) {
  ours <- legacy$counts |>
    dplyr::group_by(block_id, taxon_raw) |>
    dplyr::summarise(parsed_total = sum(count), .groups = "drop")
  legacy$totals |>
    dplyr::filter(!is.na(sheet_total)) |>
    dplyr::left_join(ours, by = c("block_id", "taxon_raw")) |>
    dplyr::mutate(parsed_total = dplyr::coalesce(parsed_total, 0)) |>
    dplyr::filter(abs(parsed_total - sheet_total) > 1e-9)
}
