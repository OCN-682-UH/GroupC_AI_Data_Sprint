# =============================================================================
# import_classroom.R -- Family C importer: classroom "Field Study" workbooks
#   Copy of OPIHI Field Study 2023 - ALL Group Data P1_P2
#   Copy of OPIHI Field Study 2025 - ALL Group Data P1_P2 - use this one
#   OPIHI Field Study 2026 - Waipuilani Data
#
# LAYOUT (tabs "QUADRAT PT INTERCEPT ...")
#   a row of student-group labels (P1G1, P2G3, B4G1 ...), one column per group
#   an "Organism" label column, one row per taxon/substrate
#   a "TOTAL Points" row: points per group = (number of quadrats) x 25
#   Waipuilani only: a date written above each block of group columns
#
# RESOLUTION: each value is a GROUP TOTAL summed over that group's quadrats,
# so detection per quadrat cannot be recovered. Records are kept with
# resolution = "group_total" and n_quadrats = TOTAL Points / 25 (team decision).
#
# NOT IMPORTED: transect tabs (team decision: quadrat data only), tabs with
# "EXAMPLE" in the name, GROUP DATA / GROUPS template tabs, and the sheets'
# per-taxon TOTALS and percentage columns (calculated, sometimes broken).
# 2023 and 2025 contain no site or date anywhere: kept with year only, flagged.
# =============================================================================

CLASS_GROUP_RE <- "^[A-Z]\\d+G\\d+$"

# -----------------------------------------------------------------------------
# parse_classroom_tab()
# Input : character matrix of one tab, file name, tab name
# Output: list(groups, counts, issues)
# -----------------------------------------------------------------------------
parse_classroom_tab <- function(m, fn, tab) {
  # blank cells -> "" for pattern searches (dplyr::coalesce() is NOT used here:
  # on a matrix it treats each row as one record and leaves single NA cells)
  m0 <- m; m0[is.na(m0)] <- ""
  is_org <- matrix(stringr::str_detect(m0,
                                       stringr::regex("^organism", ignore_case = TRUE)), nrow(m))
  pos <- which(is_org, arr.ind = TRUE)
  if (nrow(pos) == 0) stop("No 'Organism' label in ", fn, " / ", tab)
  pos <- pos[order(pos[, 1], pos[, 2]), , drop = FALSE]
  hr <- pos[1, 1]; lc <- pos[1, 2]

  is_grp <- matrix(stringr::str_detect(m0, CLASS_GROUP_RE), nrow(m))
  gr    <- which.max(rowSums(is_grp))
  if (length(gr) == 0 || sum(is_grp[gr, ]) == 0) stop("No student-group labels (e.g. P1G1) in ", fn, " / ", tab)
  gcols <- which(is_grp[gr, ])

  # dates: each contiguous block of group columns takes the date written above it
  runs <- split(gcols, cumsum(c(1, diff(gcols) != 1)))
  grp_date <- rep(as.Date(NA), length(gcols))
  for (run in runs) {
    above <- if (gr > 1) as.vector(m[1:(gr - 1), min(run):max(run), drop = FALSE]) else character(0)
    d <- excel_to_date(above)
    d <- d[!is.na(d)]
    if (length(d) > 0) grp_date[gcols %in% run] <- d[1]
  }

  trow <- which(stringr::str_detect(dplyr::coalesce(m[, lc], ""),
                                    stringr::regex("^total points", ignore_case = TRUE)))
  trow <- trow[trow > hr][1]
  if (is.na(trow)) stop("No 'TOTAL Points' row in ", fn, " / ", tab)
  taxon_rows <- (hr + 1):(trow - 1)
  taxon_rows <- taxon_rows[!is.na(m[taxon_rows, lc])]

  groups <- tibble::tibble(
    group_uid    = sprintf("%s|%s|c%03d", fn, tab, gcols),
    source_file  = fn, source_tab = tab, source_col = gcols,
    date         = grp_date,
    group_label  = m[gr, gcols],
    total_points = suppressWarnings(as.numeric(m[trow, gcols])),
    n_quadrats   = total_points / 25
  )

  vals <- m[taxon_rows, gcols, drop = FALSE]
  counts <- tibble::tibble(
    group_uid = rep(groups$group_uid, each = length(taxon_rows)),
    excel_row = rep(taxon_rows, times = length(gcols)),
    taxon_raw = rep(m[taxon_rows, lc], times = length(gcols)),
    value_raw = as.vector(vals)
  ) |>
    dplyr::filter(!is.na(value_raw)) |>
    dplyr::mutate(count = suppressWarnings(as.numeric(value_raw)))

  # issues -----------------------------------------------------------------
  sums <- counts |>
    dplyr::filter(!is.na(count)) |>
    dplyr::group_by(group_uid) |>
    dplyr::summarise(point_sum = sum(count), .groups = "drop")
  chk <- dplyr::left_join(groups, sums, by = "group_uid") |>
    dplyr::mutate(point_sum = dplyr::coalesce(point_sum, 0))

  tab_period <- stringr::str_match(tab, "\\((P\\d)\\)")[, 2]
  issues <- dplyr::bind_rows(
    chk |> dplyr::filter(!is.na(total_points), abs(point_sum - total_points) > 1e-9) |>
      dplyr::transmute(group_uid, issue = "total_points_mismatch",
                       detail = sprintf("sheet says %g, cells sum to %g", total_points, point_sum)),
    chk |> dplyr::filter(is.na(total_points) | total_points %% 25 != 0) |>
      dplyr::transmute(group_uid, issue = "total_points_not_multiple_of_25",
                       detail = paste("TOTAL Points =", total_points)),
    if (!is.na(tab_period)) {
      groups |> dplyr::filter(!startsWith(group_label, tab_period)) |>
        dplyr::transmute(group_uid, issue = "group_label_mismatch",
                         detail = sprintf("label %s on tab for %s", group_label, tab_period))
    },
    counts |> dplyr::filter(is.na(count)) |>
      dplyr::transmute(group_uid, issue = "non_numeric_count",
                       detail = paste0(taxon_raw, " = '", value_raw, "'"))
  )

  list(groups = groups,
       counts = dplyr::filter(counts, !is.na(count)) |> dplyr::select(-value_raw),
       issues = issues)
}

# -----------------------------------------------------------------------------
# import_classroom()
# Input : file names inside raw_dir
# Output: list(groups, counts, issues). `year` comes from the file name;
#         `site_raw` is taken from the file name when it names a known site.
# -----------------------------------------------------------------------------
import_classroom <- function(files, raw_dir = here::here("Data", "Raw"),
                             known_sites = character(0)) {
  parts <- list()
  for (f in files) {
    path <- file.path(raw_dir, f)
    tabs <- readxl::excel_sheets(path)
    tabs <- tabs[stringr::str_detect(toupper(tabs), "^QUADRAT PT INTERCEPT") &
                 !stringr::str_detect(toupper(tabs), "EXAMPLE")]
    for (tab in tabs) {
      parts[[length(parts) + 1]] <- parse_classroom_tab(read_sheet_matrix(path, tab), f, tab)
    }
  }
  out <- purrr::map(c(groups = "groups", counts = "counts", issues = "issues"),
                    \(nm) dplyr::bind_rows(purrr::map(parts, nm)))

  file_site <- function(fn) {
    k   <- key_site(fn)
    hit <- known_sites[purrr::map_lgl(known_sites, \(s) stringr::str_detect(k, stringr::fixed(key_site(s))))]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  out$groups <- out$groups |>
    dplyr::mutate(year     = as.integer(stringr::str_extract(source_file, "20\\d{2}")),
                  site_raw = purrr::map_chr(source_file, file_site),
                  .after = source_tab)
  out
}
