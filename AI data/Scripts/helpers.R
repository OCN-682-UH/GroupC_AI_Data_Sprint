# =============================================================================
# helpers.R -- small utilities shared by every importer
# OPIHI AI-Assisted Community Data Sprint (OCN 682 / MBIO 612, Fall 2026)
# =============================================================================

# -----------------------------------------------------------------------------
# read_sheet_matrix()
# Reads ONE Excel tab as a plain character matrix, exactly as laid out in Excel.
#   - every cell read as text (no type guessing), so nothing is silently changed
#   - range anchored at A1 so row/column numbers match Excel (readxl would
#     otherwise skip leading blank rows/columns and shift everything)
#   - whitespace trimmed; empty strings become NA
# Input : path to .xlsx, sheet name or number
# Output: character matrix; m[i, j] == Excel row i, column j
# -----------------------------------------------------------------------------
read_sheet_matrix <- function(path, sheet = 1) {
  df <- readxl::read_excel(
    path, sheet = sheet,
    col_names = FALSE, col_types = "text",
    range = cellranger::cell_limits(c(1, 1), c(NA, NA)),
    .name_repair = "minimal"
  )
  m <- as.matrix(df)
  dimnames(m) <- NULL
  m[] <- stringr::str_squish(m)   # trims ends + collapses double spaces
  m[m == ""] <- NA
  m
}

# -----------------------------------------------------------------------------
# excel_to_date()
# Excel stores dates as serial numbers (days since 1899-12-30). Read as text
# they arrive as e.g. "42867". Some cells are typed text dates instead.
# Handles both. Anything unparseable -> NA (flagged later, never dropped).
# -----------------------------------------------------------------------------
excel_to_date <- function(x) {
  x   <- stringr::str_trim(as.character(x))
  out <- as.Date(rep(NA_character_, length(x)))
  num <- suppressWarnings(as.numeric(x))

  is_serial <- !is.na(num) & num > 20000 & num < 60000   # ~1954 to ~2064
  out[is_serial] <- as.Date(num[is_serial], origin = "1899-12-30")

  is_text <- !is_serial & !is.na(x) & is.na(num)
  if (any(is_text)) {
    out[is_text] <- as.Date(lubridate::parse_date_time(
      x[is_text], orders = c("Ymd HMS", "Ymd", "mdy", "mdY", "dmy"), quiet = TRUE
    ))
  }
  out
}

# -----------------------------------------------------------------------------
# excel_to_time()
# Returns 24-hour "HH:MM" text (team decision: always military time).
# Handles Excel time fractions (0.375 = 09:00), full date-time serials (only the
# fractional part is used), and typed text like "9:30", "9:30 AM", "13:00:00".
# -----------------------------------------------------------------------------
excel_to_time <- function(x) {
  x    <- stringr::str_trim(as.character(x))
  num  <- suppressWarnings(as.numeric(x))
  mins <- rep(NA_real_, length(x))

  is_num <- !is.na(num)
  mins[is_num] <- round((num[is_num] %% 1) * 24 * 60)

  is_text <- !is_num & !is.na(x)
  if (any(is_text)) {
    mt <- stringr::str_match(
      stringr::str_to_lower(x[is_text]),
      "^(\\d{1,2})[:.](\\d{2})(?::\\d{2})?\\s*(a\\.?m\\.?|p\\.?m\\.?)?"
    )
    h  <- as.numeric(mt[, 2])
    mm <- as.numeric(mt[, 3])
    ap <- mt[, 4]
    h  <- ifelse(!is.na(ap) & startsWith(ap, "p") & h < 12, h + 12, h)
    h  <- ifelse(!is.na(ap) & startsWith(ap, "a") & h == 12, 0, h)
    mins[is_text] <- h * 60 + mm
  }

  mins <- mins %% 1440
  out  <- sprintf("%02d:%02d", as.integer(mins %/% 60), as.integer(mins %% 60))
  out[is.na(mins)] <- NA_character_
  out
}

# -----------------------------------------------------------------------------
# export_raw_csv()
# Assignment step 2: convert proprietary Excel files to an open format.
# Writes every tab of every .xlsx in `raw_dir` to `csv_dir` as
#   <file name>__<tab name>.csv
# Raw Excel files are only READ, never modified.
# Tabs that cannot be read as cells (e.g. Excel "chart sheets", which hold only
# a graph) are skipped and logged with status "not_exported", never silently.
# Output: (invisibly) a tibble with one row per tab and its export status.
# -----------------------------------------------------------------------------
export_raw_csv <- function(raw_dir = here::here("Data", "Raw"),
                           csv_dir = here::here("Data", "Processed", "raw_csv")) {
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  files <- list_raw_xlsx(raw_dir)

  log <- purrr::map_dfr(files, function(f) {
    purrr::map_dfr(readxl::excel_sheets(f), function(sh) {
      out <- file.path(
        csv_dir,
        paste0(tools::file_path_sans_ext(basename(f)), "__",
               stringr::str_replace_all(sh, "[^A-Za-z0-9]+", "_"), ".csv")
      )
      res <- tryCatch({
        m <- read_sheet_matrix(f, sh)
        readr::write_csv(as.data.frame(m), out, na = "", col_names = FALSE)
        tibble::tibble(status = "exported", n_rows = nrow(m), n_cols = ncol(m),
                       csv = basename(out), note = NA_character_)
      }, error = function(e) {
        tibble::tibble(status = "not_exported", n_rows = NA_integer_, n_cols = NA_integer_,
                       csv = NA_character_,
                       note = paste("could not read as cells (likely a chart sheet):",
                                    conditionMessage(e)))
      })
      dplyr::bind_cols(tibble::tibble(file = basename(f), sheet = sh), res)
    })
  })
  invisible(log)
}

# -----------------------------------------------------------------------------
# list_raw_xlsx()
# Lists .xlsx files in a folder, ignoring the temporary "~$file.xlsx" lock files
# Excel creates while a workbook is open (they would otherwise break imports).
# -----------------------------------------------------------------------------
list_raw_xlsx <- function(raw_dir, full.names = TRUE) {
  list.files(raw_dir, pattern = "^[^~].*\\.xlsx$", full.names = full.names)
}
