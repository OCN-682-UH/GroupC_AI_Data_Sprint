# =============================================================================
# import_modern.R -- Family B importer: modern row-per-quadrat datasheets
#   OPIHI 2024 datasheet, OPIHI BIOL 499 2025 LT data,
#   OPIHI 2026 Longitudinal Monitoring Data, RAW Mauliola OPIHI Data 2026
#
# LAYOUT (tab "Quadrat PtCounts")
#   category row : labels above the taxon columns (Limu (Algae), Invertebrates,
#                  Substrate, Total). In 2024/2025 the algae block is unlabeled.
#   header row   : Date | Location | Data Entry By | Transect # |
#                  [Quadrat location (m)] | Position # on Transect (m) |
#                  <one column per taxon/substrate> | Total Points | notes...
#   data rows    : one row per quadrat, cell = points hit out of 25
#
# Other tabs in these files (Quadrat PercentCover, Grouped..., Sediment Data,
# Transect) are derived from this tab or use another method, so they are not
# imported here (team decision: quadrat point counts only).
#
# ROW RULES (each excluded row is logged with its reason)
#   example_row        : any cell contains "EXAMPLE" (student-template rows)
#   summary_row        : a "Column Totals" row
#   blank_template_row : no date AND every taxon cell blank or 0
#   no_taxon_data      : has date/site but every taxon cell blank (not entered)
#   invalid date       : if the rows before and after share the same date and
#                        site, that date is used and flagged "date_inferred";
#                        otherwise the date is NA and flagged "invalid_date"
#   blank taxon cell in a filled row = 0 (not detected)
# =============================================================================

# -----------------------------------------------------------------------------
# parse_modern_file()
# Input : path to one Family B .xlsx
# Output: list of tibbles
#   quadrats - one row per quadrat (sheet row) kept
#   counts   - one row per non-empty taxon cell: quadrat x taxon x count
#   metadata - the Metadata tab, one row per group/team, raw text kept
#   excluded - rows not imported, with the rule that excluded them
#   issues   - flags on kept rows (inferred dates, Total Points mismatches...)
# -----------------------------------------------------------------------------
parse_modern_file <- function(path, sheet = "Quadrat PtCounts") {
  fn <- basename(path)
  m  <- read_sheet_matrix(path, sheet)

  # ---- header + category rows ----------------------------------------------
  hi <- which(tolower(m[1:min(5, nrow(m)), 1]) == "date")[1]
  if (is.na(hi)) stop("No 'Date' header found in first 5 rows of ", fn)
  hdr <- m[hi, ]
  cat_row <- if (hi > 1) m[hi - 1, ] else rep(NA_character_, ncol(m))
  low <- tolower(dplyr::coalesce(hdr, ""))
  find_col <- function(re) which(stringr::str_detect(low, re))[1]

  j_date  <- find_col("^date")
  j_site  <- find_col("^location")
  j_entry <- find_col("^data entry")
  j_tran  <- find_col("^(# ?t$|transect)")
  j_qloc  <- find_col("^quadrat location")
  j_pos   <- find_col("^position")
  j_total <- find_col("^total points")
  if (is.na(j_total)) stop("No 'Total Points' column in ", fn)
  j_first <- max(c(j_tran, j_qloc, j_pos), na.rm = TRUE) + 1

  taxon_cols <- (j_first:(j_total - 1))[!is.na(hdr[j_first:(j_total - 1)])]

  # sheet category: carry each label right; unlabeled leading columns = algae
  cats <- cat_row[j_first:(j_total - 1)]
  cats <- tidyr::fill(tibble::tibble(x = cats), x)$x
  cats[is.na(cats)] <- "Limu (Algae)"
  names(cats) <- j_first:(j_total - 1)
  taxa <- tibble::tibble(col = taxon_cols, taxon_raw = hdr[taxon_cols],
                         sheet_category = unname(cats[as.character(taxon_cols)]))

  # ---- candidate rows (skip rows that are completely empty) ----------------
  rows <- (hi + 1):nrow(m)
  rows <- rows[rowSums(!is.na(m[rows, , drop = FALSE])) > 0]

  tv        <- m[rows, taxa$col, drop = FALSE]
  all_blank <- rowSums(!is.na(tv)) == 0
  blank_or0 <- rowSums(!is.na(tv) & tv != "0") == 0
  row_text  <- apply(m[rows, , drop = FALSE], 1, \(r) paste(r[!is.na(r)], collapse = " "))

  r <- tibble::tibble(
    excel_row = rows,
    date_raw  = m[rows, j_date],
    date      = excel_to_date(m[rows, j_date]),
    site_raw  = m[rows, j_site],
    entry_by_raw = if (!is.na(j_entry)) m[rows, j_entry] else NA_character_,
    transect  = suppressWarnings(as.integer(m[rows, j_tran])),
    point_m   = suppressWarnings(as.numeric(dplyr::coalesce(
                  if (!is.na(j_qloc)) m[rows, j_qloc] else NA_character_,
                  if (!is.na(j_pos))  m[rows, j_pos]  else NA_character_))),
    total_points_sheet = suppressWarnings(as.numeric(m[rows, j_total])),
    exclude = dplyr::case_when(
      stringr::str_detect(row_text, stringr::regex("example", ignore_case = TRUE))       ~ "example_row",
      stringr::str_detect(row_text, stringr::regex("column totals", ignore_case = TRUE)) ~ "summary_row",
      is.na(m[rows, j_date]) & blank_or0 ~ "blank_template_row",
      all_blank                          ~ "no_taxon_data",
      TRUE ~ NA_character_
    ),
    flag = NA_character_
  )

  # ---- repair invalid dates from identical neighbours ----------------------
  for (k in which(is.na(r$exclude) & is.na(r$date))) {
    prv <- utils::tail(which(!is.na(r$date) & seq_len(nrow(r)) < k), 1)
    nxt <- which(!is.na(r$date) & seq_len(nrow(r)) > k)[1]
    same <- length(prv) == 1 && !is.na(nxt) &&
      r$date[prv] == r$date[nxt] &&
      identical(r$site_raw[prv], r$site_raw[k]) && identical(r$site_raw[nxt], r$site_raw[k])
    if (same) { r$date[k] <- r$date[prv]; r$flag[k] <- "date_inferred" }
    else      { r$flag[k] <- "invalid_date" }
  }

  excluded <- r |>
    dplyr::filter(!is.na(exclude)) |>
    dplyr::transmute(source_file = fn, excel_row, reason = exclude, date_raw, site_raw)
  kept <- dplyr::filter(r, is.na(exclude))
  kept_idx <- match(kept$excel_row, rows)

  # ---- counts: wide -> long -----------------------------------------------
  vals <- tv[kept_idx, , drop = FALSE]
  counts <- tibble::tibble(
    excel_row = rep(kept$excel_row, times = ncol(vals)),
    col       = rep(taxa$col, each = nrow(vals)),
    value_raw = as.vector(vals)
  ) |>
    dplyr::filter(!is.na(value_raw)) |>
    dplyr::mutate(count = suppressWarnings(as.numeric(value_raw))) |>
    dplyr::left_join(taxa, by = "col")

  bad <- dplyr::filter(counts, is.na(count))
  counts <- dplyr::filter(counts, !is.na(count))

  pts <- counts |>
    dplyr::group_by(excel_row) |>
    dplyr::summarise(point_sum = sum(count), .groups = "drop")

  quadrats <- kept |>
    dplyr::left_join(pts, by = "excel_row") |>
    dplyr::mutate(point_sum = dplyr::coalesce(point_sum, 0),
                  quad_id = sprintf("%s|r%04d", fn, excel_row),
                  source_file = fn, .before = 1) |>
    dplyr::select(-exclude)

  issues <- dplyr::bind_rows(
    quadrats |> dplyr::filter(!is.na(flag)) |>
      dplyr::transmute(quad_id, issue = flag, detail = paste("date cell:", date_raw)),
    quadrats |> dplyr::filter(!is.na(total_points_sheet), abs(total_points_sheet - point_sum) > 1e-9) |>
      dplyr::transmute(quad_id, issue = "total_points_mismatch",
                       detail = sprintf("sheet says %g, cells sum to %g", total_points_sheet, point_sum)),
    bad |> dplyr::transmute(quad_id = sprintf("%s|r%04d", fn, excel_row), issue = "non_numeric_count",
                            detail = paste0(taxon_raw, " = '", value_raw, "'"))
  )

  counts <- counts |>
    dplyr::mutate(quad_id = sprintf("%s|r%04d", fn, excel_row)) |>
    dplyr::select(quad_id, excel_row, col, taxon_raw, sheet_category, count)

  # ---- Metadata tab (raw text; joined to surveys after site cleaning) ------
  metadata <- NULL
  if ("Metadata" %in% readxl::excel_sheets(path)) {
    mm <- read_sheet_matrix(path, "Metadata")
    nm <- stringr::str_replace_all(tolower(dplyr::coalesce(mm[1, ], "")), "[^a-z0-9#]+", "_")
    nm <- stringr::str_replace_all(nm, "#", "n")
    nm <- stringr::str_remove_all(nm, "^_+|_+$")
    nm[nm == ""] <- paste0("col", which(nm == ""))
    body <- mm[-1, , drop = FALSE]
    body <- body[rowSums(!is.na(body)) > 0, , drop = FALSE]
    metadata <- tibble::as_tibble(body, .name_repair = ~ make.unique(nm)) |>
      dplyr::mutate(source_file = fn, .before = 1)
  }

  list(quadrats = quadrats, counts = counts, metadata = metadata,
       excluded = excluded, issues = issues)
}

# -----------------------------------------------------------------------------
# import_modern()
# Runs parse_modern_file() on every Family B file and stacks the results.
# `files` = file names inside raw_dir.
# -----------------------------------------------------------------------------
import_modern <- function(files, raw_dir = here::here("Data", "Raw")) {
  parts <- purrr::map(file.path(raw_dir, files), parse_modern_file)
  purrr::map(purrr::set_names(names(parts[[1]])),
             function(nm) dplyr::bind_rows(purrr::map(parts, nm)))
}
