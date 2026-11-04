# =============================================================================
# clean_build.R -- harmonise sites, taxa, teachers, times; deduplicate;
#                  cross-check compiled files; assemble the master tables.
# Every change is rule-based and logged; raw files are never modified.
# =============================================================================

# ---- 1. Name keys (must match the keys stored in Data/Lookups/*.csv) --------

# Taxon key: Unicode-normalised, lower case, single spaces, trailing "." dropped
key_taxon <- function(x) {
  x <- stringi::stri_trans_nfkc(as.character(x))
  x <- stringr::str_to_lower(stringr::str_squish(x))
  stringr::str_trim(stringr::str_remove(x, "\\.+$"))
}

# Site key: ʻokina / apostrophes / backticks removed, kahakō removed, lower case
key_site <- function(x) {
  x <- stringr::str_remove_all(as.character(x), "[ʻʼ'‘’`]")
  x <- stringi::stri_trans_general(x, "NFD; [:Nonspacing Mark:] Remove; NFC")
  stringr::str_to_lower(stringr::str_squish(x))
}

read_lookup <- function(path) readr::read_csv(path, col_types = readr::cols(.default = "c"), na = "")

# For names missing from a lookup: closest existing key by edit distance
suggest_match <- function(keys, lookup_keys) {
  if (length(keys) == 0) return(character(0))
  d <- utils::adist(keys, lookup_keys, ignore.case = TRUE)
  lookup_keys[apply(d, 1, which.min)]
}

# ---- 2. Site harmonisation ---------------------------------------------------
# Output: one row per distinct raw site string with its canonical site.
# Unmatched names are flagged with a suggested closest match (never guessed).
standardize_sites <- function(site_raw, site_xw) {
  tibble::tibble(site_raw = unique(site_raw)) |>
    dplyr::filter(!is.na(site_raw)) |>
    dplyr::mutate(site_key = key_site(site_raw)) |>
    dplyr::left_join(dplyr::select(site_xw, site_key, site_name, site_alias_english, island),
                     by = "site_key") |>
    dplyr::mutate(site_suggestion = dplyr::if_else(
      is.na(site_name), suggest_match(site_key, site_xw$site_key)[match(site_key, site_key)], NA_character_))
}

# ---- 3. Taxon harmonisation --------------------------------------------------
standardize_taxa <- function(taxon_raw, taxon_xw) {
  tibble::tibble(taxon_raw = unique(taxon_raw)) |>
    dplyr::filter(!is.na(taxon_raw)) |>
    dplyr::mutate(taxon_key = key_taxon(taxon_raw)) |>
    dplyr::left_join(dplyr::select(taxon_xw, taxon_key, taxon_clean, analysis_taxon, category),
                     by = "taxon_key") |>
    # rule: a taxon "name" that is only a number (e.g. "3", "3.0") is a data-entry error
    dplyr::mutate(
      is_number      = stringr::str_detect(taxon_key, "^[0-9]+(\\.[0-9]+)?$"),
      taxon_clean    = dplyr::if_else(is_number, "(invalid entry)", taxon_clean),
      analysis_taxon = dplyr::if_else(is_number, NA_character_, analysis_taxon),
      category       = dplyr::if_else(is_number, "Invalid", category)) |>
    dplyr::select(-is_number) |>
    dplyr::mutate(taxon_suggestion = dplyr::if_else(
      is.na(taxon_clean), suggest_match(taxon_key, taxon_xw$taxon_key)[match(taxon_key, taxon_key)], NA_character_))
}

# ---- 4. People: split, tidy, and assign anonymous IDs ------------------------
# split_names(): "Cam, Linnea" / "Patrick and Ashlan" / "A & B" -> tidy names
split_names <- function(x) {
  parts <- stringr::str_split(dplyr::coalesce(x, ""), ",|;|&|/|\\band\\b")
  purrr::map(parts, function(p) {
    p <- stringr::str_to_title(stringr::str_squish(p))
    unique(p[p != "" & !toupper(p) %in% c("NA", "N/A", "NONE")])
  })
}

# tidy_org(): legacy School names are NOT split ("Arts and Science" is one school)
tidy_org <- function(x) {
  x <- stringr::str_squish(x)
  dplyr::if_else(is.na(x) | x == "", NA_character_, x)
}

# ---- 5. Field times (team rule: 24-hour HH:MM; 01-05 h assumed PM) ----------
fix_field_time <- function(hhmm) {
  h <- suppressWarnings(as.integer(substr(hhmm, 1, 2)))
  flag <- dplyr::if_else(!is.na(h) & h >= 1 & h <= 5, "time_assumed_pm", NA_character_)
  new  <- dplyr::if_else(!is.na(flag), sprintf("%02d:%s", h + 12L, substr(hhmm, 4, 5)), hhmm)
  tibble::tibble(time = new, flag = flag)
}

# ---- 6. Distances -------------------------------------------------------------
haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(a))
}

# ---- 7. Compiled (Family D) trip list ----------------------------------------
# Reads Date / Location / School / coordinates (and whether any taxon count > 0)
# from the trip-total tab of a compiled workbook. Used ONLY for cross-checks.
read_compiled_trips <- function(path, sheet = "Totals_NO_Excl") {
  m  <- read_sheet_matrix(path, sheet)
  hi <- which(tolower(m[, 1]) == "date")[1]
  h  <- tolower(dplyr::coalesce(m[hi, ], ""))
  col <- function(re) which(stringr::str_detect(h, re))[1]
  j_notes <- col("^notes")
  body <- m[(hi + 1):nrow(m), , drop = FALSE]
  taxa <- body[, (j_notes + 1):ncol(m), drop = FALSE]
  tibble::tibble(
    excel_row = (hi + 1):nrow(m),
    date      = excel_to_date(body[, 1]),
    site_raw  = body[, col("^location")],
    school_raw = body[, col("^school")],
    latitude  = suppressWarnings(as.numeric(body[, col("^latitude")])),
    longitude = suppressWarnings(as.numeric(body[, col("^longitude")])),
    has_species_data = rowSums(matrix(suppressWarnings(as.numeric(taxa)), nrow(taxa)) > 0, na.rm = TRUE) > 0
  ) |>
    dplyr::filter(!is.na(date), !is.na(site_raw)) |>
    dplyr::mutate(source_file = basename(path), .before = 1)
}

# =============================================================================
# 8. build_opihi(): everything above, applied to the three imports
# Inputs : legacy / modern / classroom import lists, lookup tables, compiled trips
# Output : list of tibbles (surveys, samples, observations, master, exclusions,
#          qc_flags, crosscheck, unmatched, teacher_key)
# =============================================================================
build_opihi <- function(legacy, modern, classroom, site_xw, taxon_xw, compiled) {

  # ---------- sites (all raw strings from every source) ----------
  sites <- standardize_sites(c(legacy$surveys$site_raw, modern$quadrats$site_raw,
                               modern$metadata$location, classroom$groups$site_raw,
                               compiled$site_raw), site_xw)
  add_site <- function(df) dplyr::left_join(df, dplyr::select(sites, site_raw, site_name), by = "site_raw")

  exclusions <- list(); flags <- list()

  # ---------- A: legacy surveys ----------
  a_bad <- legacy$issues |> dplyr::filter(issue == "no_survey_metadata")
  exclusions$a_meta <- a_bad |>
    dplyr::transmute(level = "survey", source_file = stringr::str_remove(block_id, "\\|.*$"),
                     source_ref = block_id, reason = "no_survey_metadata", detail)

  st <- fix_field_time(legacy$surveys$start_time); et <- fix_field_time(legacy$surveys$stop_time)
  a_sv <- legacy$surveys |>
    dplyr::filter(!block_id %in% a_bad$block_id) |>
    add_site() |>
    dplyr::mutate(
      start_time = st$time[match(block_id, legacy$surveys$block_id)],
      stop_time  = et$time[match(block_id, legacy$surveys$block_id)],
      placeholder = start_time %in% "00:00" & (stop_time %in% c("01:00", "13:00")),
      time_flag = dplyr::case_when(placeholder ~ "placeholder_time_removed",
                                   !is.na(st$flag[match(block_id, legacy$surveys$block_id)]) |
                                   !is.na(et$flag[match(block_id, legacy$surveys$block_id)]) ~ "time_assumed_pm"),
      start_time = dplyr::if_else(placeholder, NA_character_, start_time),
      stop_time  = dplyr::if_else(placeholder, NA_character_, stop_time),
      teacher_names = purrr::map(tidy_org(school_raw), \(s) if (is.na(s)) character(0) else s)
    )

  # coordinate check: labelled site vs coordinates (reference = median per site
  # from legacy + compiled files). Relabel only on an exact-coordinate match.
  ref <- dplyr::bind_rows(
      dplyr::select(a_sv, site_name, latitude, longitude),
      add_site(compiled) |> dplyr::select(site_name, latitude, longitude)) |>
    dplyr::filter(!is.na(site_name), !is.na(latitude), !is.na(longitude)) |>
    dplyr::group_by(site_name) |>
    dplyr::summarise(ref_lat = stats::median(latitude), ref_lon = stats::median(longitude), .groups = "drop")
  nearest <- purrr::map2(a_sv$latitude, a_sv$longitude, function(la, lo) {
    if (is.na(la) || is.na(lo)) return(list(site = NA_character_, km = NA_real_))
    d <- haversine_km(la, lo, ref$ref_lat, ref$ref_lon)
    list(site = ref$site_name[which.min(d)], km = min(d))
  })
  a_sv <- a_sv |>
    dplyr::left_join(ref, by = "site_name") |>
    dplyr::mutate(
      km_to_label  = haversine_km(latitude, longitude, ref_lat, ref_lon),
      nearest_site = purrr::map_chr(nearest, "site"),
      km_nearest   = purrr::map_dbl(nearest, "km"),
      site_flag = dplyr::case_when(
        !is.na(km_to_label) & km_to_label > 5 & km_nearest < 0.2 ~ "site_corrected_by_coordinates",
        !is.na(km_to_label) & km_to_label > 5                     ~ "coordinates_conflict_with_site"),
      site_label_raw = site_name,
      site_name = dplyr::if_else(site_flag %in% "site_corrected_by_coordinates", nearest_site, site_name)
    )
  flags$a_site <- a_sv |> dplyr::filter(!is.na(site_flag)) |>
    dplyr::transmute(level = "survey", ref = block_id, flag = site_flag,
                     detail = sprintf("labelled %s, coordinates %.1f km from it; nearest site %s (%.2f km)",
                                      site_label_raw, km_to_label, nearest_site, km_nearest))
  flags$a_time <- a_sv |> dplyr::filter(!is.na(time_flag)) |>
    dplyr::transmute(level = "survey", ref = block_id, flag = time_flag,
                     detail = paste("start", start_time_raw, "stop", stop_time_raw))

  a_surveys <- a_sv |>
    dplyr::transmute(src_key = block_id, source_family = "A", source_file, resolution = "quadrat",
                     date, site_name, latitude, longitude, start_time, stop_time,
                     low_tide_time = NA_character_, teacher_names,
                     n_transects_reported = n_transects, n_quad_per_t_reported = n_quad_per_t,
                     transect_length_m, point_spacing_m)

  a_samples <- legacy$quadrats |>
    dplyr::filter(block_id %in% a_sv$block_id) |>
    dplyr::transmute(src_key = block_id, sample_key = sprintf("%s|c%04d", block_id, source_col),
                     sample_type = "quadrat", transect, quadrat, point_m, n_quadrats = 1,
                     sample_flag = dplyr::if_else(header_status == "inferred", "header_inferred", NA_character_))
  a_obs <- legacy$counts |>
    dplyr::filter(block_id %in% a_sv$block_id) |>
    dplyr::transmute(sample_key = sprintf("%s|c%04d", block_id, col), taxon_raw, count,
                     source_cell = sprintf("%s!R%dC%d", stringr::str_remove(block_id, "\\|.*$"), excel_row, col))

  # ---------- B: modern surveys ----------
  md <- modern$metadata |>
    dplyr::mutate(date = excel_to_date(date), site_raw = location) |>
    add_site() |>
    dplyr::mutate(low_tide_time = excel_to_time(low_tide_time),
                  start_time = excel_to_time(start_time), end_time = excel_to_time(end_time))
  b_q <- modern$quadrats |> add_site()
  b_sv <- b_q |>
    dplyr::distinct(source_file, date, site_name) |>
    dplyr::mutate(src_key = sprintf("%s|%s|%s", source_file, as.character(date), site_name))
  # metadata match: same file + site + date, else same file + site within 7 days
  b_meta <- purrr::pmap(b_sv, function(source_file, date, site_name, src_key) {
    cand <- md[md$source_file == source_file & md$site_name %in% site_name, ]
    exact <- cand[cand$date %in% date, ]
    if (nrow(exact) > 0) return(list(rows = exact, flag = NA_character_))
    near <- cand[!is.na(cand$date) & abs(as.numeric(cand$date - date)) <= 7, ]
    if (nrow(near) > 0) return(list(rows = near, flag = "metadata_date_mismatch"))
    list(rows = cand[0, ], flag = "no_metadata_match")
  })
  first_time <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA_character_ }
  b_sv <- b_sv |>
    dplyr::mutate(
      teacher_names = purrr::map(b_meta, \(b) as.character(sort(unique(unlist(split_names(b$rows$team_leader)))))),
      low_tide_raw  = purrr::map_chr(b_meta, \(b) first_time(b$rows$low_tide_time)),
      start_raw     = purrr::map_chr(b_meta, \(b) first_time(b$rows$start_time)),
      end_raw       = purrr::map_chr(b_meta, \(b) first_time(b$rows$end_time)),
      meta_flag     = purrr::map_chr(b_meta, "flag"))
  lt <- fix_field_time(b_sv$low_tide_raw); s2 <- fix_field_time(b_sv$start_raw); e2 <- fix_field_time(b_sv$end_raw)
  b_sv <- b_sv |>
    dplyr::mutate(low_tide_time = lt$time, start_time = s2$time, stop_time = e2$time,
                  time_flag = dplyr::if_else(!is.na(lt$flag) | !is.na(s2$flag) | !is.na(e2$flag),
                                             "time_assumed_pm", NA_character_))
  flags$b_meta <- b_sv |> dplyr::filter(!is.na(meta_flag)) |>
    dplyr::transmute(level = "survey", ref = src_key, flag = meta_flag,
                     detail = "metadata row matched on file + site within 7 days (data-row date kept)")
  flags$b_time <- b_sv |> dplyr::filter(!is.na(time_flag)) |>
    dplyr::transmute(level = "survey", ref = src_key, flag = time_flag,
                     detail = sprintf("low tide %s, start %s, end %s", low_tide_raw, start_raw, end_raw))

  b_surveys <- b_sv |>
    dplyr::transmute(src_key, source_family = "B", source_file, resolution = "quadrat", date, site_name,
                     latitude = NA_real_, longitude = NA_real_, start_time, stop_time, low_tide_time,
                     teacher_names, n_transects_reported = NA_integer_, n_quad_per_t_reported = NA_integer_,
                     transect_length_m = NA_real_, point_spacing_m = NA_real_)
  b_flag_by_q <- modern$issues |> dplyr::group_by(quad_id) |>
    dplyr::summarise(f = paste(unique(issue), collapse = ";"), .groups = "drop")
  b_samples <- b_q |>
    dplyr::transmute(src_key = sprintf("%s|%s|%s", source_file, as.character(date), site_name), sample_key = quad_id,
                     sample_type = "quadrat", transect, quadrat = NA_integer_, point_m, n_quadrats = 1,
                     sample_flag = b_flag_by_q$f[match(quad_id, b_flag_by_q$quad_id)])
  b_obs <- modern$counts |>
    dplyr::transmute(sample_key = quad_id, taxon_raw, count,
                     source_cell = sprintf("%s!R%dC%d", stringr::str_remove(quad_id, "\\|.*$"), excel_row, col))
  exclusions$b_rows <- modern$excluded |>
    dplyr::transmute(level = "row", source_file, source_ref = sprintf("%s|r%04d", source_file, excel_row),
                     reason, detail = paste("date:", date_raw, "site:", site_raw))

  # ---------- C: classroom group totals ----------
  c_g <- classroom$groups |> add_site()
  c_sv <- c_g |>
    dplyr::distinct(source_file, year, date, site_name) |>
    dplyr::mutate(src_key = sprintf("%s|%s", source_file, dplyr::coalesce(as.character(date), "no-date")))
  flags$c_nosite <- c_sv |> dplyr::filter(is.na(site_name) | is.na(date)) |>
    dplyr::transmute(level = "survey", ref = src_key, flag = "no_site_or_date_in_file",
                     detail = "classroom file names no site/date; year from file name only")
  c_surveys <- c_sv |>
    dplyr::transmute(src_key, source_family = "C", source_file, resolution = "group_total", date, site_name,
                     year_from_file = year, latitude = NA_real_, longitude = NA_real_,
                     start_time = NA_character_, stop_time = NA_character_, low_tide_time = NA_character_,
                     teacher_names = list(character(0)), n_transects_reported = NA_integer_,
                     n_quad_per_t_reported = NA_integer_, transect_length_m = NA_real_, point_spacing_m = NA_real_)
  c_flag_by_g <- classroom$issues |> dplyr::group_by(group_uid) |>
    dplyr::summarise(f = paste(unique(issue), collapse = ";"), .groups = "drop")
  c_samples <- c_g |>
    dplyr::transmute(src_key = sprintf("%s|%s", source_file, dplyr::coalesce(as.character(date), "no-date")),
                     sample_key = group_uid, sample_type = "student_group", group_label,
                     transect = NA_integer_, quadrat = NA_integer_, point_m = NA_real_, n_quadrats,
                     sample_flag = c_flag_by_g$f[match(group_uid, c_flag_by_g$group_uid)])
  c_obs <- classroom$counts |>
    dplyr::transmute(sample_key = group_uid, taxon_raw, count,
                     source_cell = sprintf("%s!R%d", stringr::str_remove(group_uid, "\\|c\\d+$"), excel_row))

  # ---------- stack ----------
  surveys <- dplyr::bind_rows(a_surveys, b_surveys, c_surveys) |>
    dplyr::mutate(year = dplyr::coalesce(as.integer(format(date, "%Y")), year_from_file)) |>
    dplyr::select(-year_from_file)
  samples <- dplyr::bind_rows(a_samples, b_samples, c_samples)
  obs     <- dplyr::bind_rows(a_obs, b_obs, c_obs)

  # ---------- taxa ----------
  taxa <- standardize_taxa(obs$taxon_raw, taxon_xw)
  obs <- obs |>
    dplyr::left_join(dplyr::select(taxa, taxon_raw, taxon_clean, analysis_taxon, category), by = "taxon_raw") |>
    dplyr::left_join(dplyr::select(samples, sample_key, src_key), by = "sample_key")
  exclusions$invalid_taxa <- obs |> dplyr::filter(category %in% "Invalid") |>
    dplyr::transmute(level = "observation", source_file = stringr::str_remove(source_cell, "!.*$"),
                     source_ref = source_cell, reason = "invalid_taxon_name",
                     detail = sprintf("taxon '%s' count %g", taxon_raw, count))
  obs <- obs |> dplyr::filter(!category %in% "Invalid")

  # ---------- duplicates ----------
  # (a) identical survey content in two places: fingerprint = sorted
  #     sample position + taxon + count. Keep the first (A, then B, then C).
  fp <- obs |>
    dplyr::left_join(dplyr::select(samples, sample_key, transect, quadrat, point_m, group_label), by = "sample_key") |>
    dplyr::mutate(item = paste(transect, quadrat, point_m, group_label, taxon_clean, count, sep = "|")) |>
    dplyr::group_by(src_key) |>
    dplyr::summarise(fingerprint = paste(sort(item), collapse = "#"), .groups = "drop")
  surveys <- surveys |>
    dplyr::left_join(fp, by = "src_key") |>
    dplyr::arrange(source_family, source_file, date) |>
    dplyr::mutate(dup_of = ifelse(!is.na(fingerprint) & duplicated(fingerprint),
                                  src_key[match(fingerprint, fingerprint)], NA_character_))
  exclusions$dup_surveys <- surveys |> dplyr::filter(!is.na(dup_of)) |>
    dplyr::transmute(level = "survey", source_file, source_ref = src_key, reason = "duplicate_survey",
                     detail = paste("identical content to", dup_of))
  # (b) duplicate quadrats inside one survey (same position + identical counts)
  sfp <- obs |>
    dplyr::group_by(sample_key) |>
    dplyr::summarise(content = paste(sort(paste(taxon_clean, count)), collapse = "#"), .groups = "drop")
  samples <- samples |>
    dplyr::left_join(sfp, by = "sample_key") |>
    dplyr::group_by(src_key, transect, quadrat, point_m, group_label, content) |>
    dplyr::mutate(dup_sample = sample_type == "quadrat" & !is.na(content) & dplyr::row_number() > 1) |>
    dplyr::ungroup()
  exclusions$dup_samples <- samples |> dplyr::filter(dup_sample) |>
    dplyr::transmute(level = "sample", source_file = stringr::str_remove(sample_key, "\\|.*$"),
                     source_ref = sample_key, reason = "duplicate_sample", detail = "identical position and counts")
  drop_src <- surveys$src_key[!is.na(surveys$dup_of)]
  surveys <- dplyr::filter(surveys, is.na(dup_of)) |> dplyr::select(-fingerprint, -dup_of)
  samples <- dplyr::filter(samples, !dup_sample, !src_key %in% drop_src) |> dplyr::select(-content, -dup_sample)
  obs     <- dplyr::filter(obs, sample_key %in% samples$sample_key)
  # (c) same site + date in more than one source file: flagged for review only
  flags$same_day <- surveys |> dplyr::filter(!is.na(date), !is.na(site_name)) |>
    dplyr::group_by(site_name, date) |> dplyr::filter(dplyr::n_distinct(source_file) > 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(level = "survey", ref = src_key, flag = "same_site_date_other_file",
                     detail = "another source file has a survey at this site on this date (different content)")

  # ---------- IDs (deterministic: ordered by date, site, source) ----------
  surveys <- surveys |>
    dplyr::arrange(date, site_name, source_family, source_file, src_key) |>
    dplyr::mutate(survey_id = sprintf("S%04d", dplyr::row_number()), .before = 1)
  samples <- samples |>
    dplyr::left_join(dplyr::select(surveys, src_key, survey_id), by = "src_key") |>
    dplyr::arrange(survey_id, transect, quadrat, point_m, group_label, sample_key) |>
    dplyr::group_by(survey_id) |>
    dplyr::mutate(sample_id = sprintf("%s-%03d", survey_id, dplyr::row_number())) |>
    dplyr::ungroup()

  # teacher IDs: order of first appearance in time; several names -> "T001;T004"
  people <- surveys |> dplyr::select(survey_id, teacher_names) |>
    tidyr::unnest_longer(teacher_names, values_to = "teacher_name") |>
    dplyr::filter(!is.na(teacher_name))
  teacher_key <- people |> dplyr::distinct(teacher_name) |>
    dplyr::mutate(teacher_id = sprintf("T%03d", dplyr::row_number()))
  tid <- people |> dplyr::left_join(teacher_key, by = "teacher_name") |>
    dplyr::group_by(survey_id) |>
    dplyr::summarise(teacher_id = paste(sort(unique(teacher_id)), collapse = ";"), .groups = "drop")
  surveys <- surveys |>
    dplyr::left_join(tid, by = "survey_id") |>
    dplyr::select(-teacher_names)

  # ---------- effort per sample ----------
  pts <- obs |> dplyr::group_by(sample_key) |>
    dplyr::summarise(point_sum = sum(count),
                     n_taxa = dplyr::n_distinct(taxon_clean[!category %in% "Substrate" & count > 0]),
                     .groups = "drop")
  samples <- samples |>
    dplyr::left_join(pts, by = "sample_key") |>
    dplyr::mutate(point_sum = dplyr::coalesce(point_sum, 0), n_taxa = dplyr::coalesce(n_taxa, 0L),
                  points_possible = 25 * n_quadrats)
  # add "points_over_25" to any existing flags for that quadrat
  samples$sample_flag <- purrr::map2_chr(samples$sample_flag,
    samples$sample_type == "quadrat" & samples$point_sum > samples$points_possible,
    \(f, over) { v <- c(if (!is.na(f)) strsplit(f, ";")[[1]], if (over) "points_over_25")
                 if (length(v)) paste(unique(v), collapse = ";") else NA_character_ })

  surveys <- surveys |>
    dplyr::left_join(samples |> dplyr::group_by(survey_id) |>
                       dplyr::summarise(n_samples = dplyr::n(), n_quadrats = sum(n_quadrats), .groups = "drop"),
                     by = "survey_id") |>
    dplyr::left_join(dplyr::distinct(site_xw, site_name, site_alias_english, island), by = "site_name")

  # ---------- cross-check against compiled (Family D) files ----------
  comp <- add_site(compiled) |> dplyr::distinct(source_file, date, site_name, has_species_data, .keep_all = TRUE)
  crosscheck <- comp |>
    dplyr::left_join(dplyr::distinct(surveys, date, site_name) |> dplyr::mutate(in_master = TRUE),
                     by = c("date", "site_name")) |>
    dplyr::mutate(in_master = dplyr::coalesce(in_master, FALSE),
                  status = dplyr::case_when(
                    in_master ~ "in master",
                    !has_species_data ~ "compiled only - trip has no species data (known gap)",
                    TRUE ~ "compiled only - HAS species data, check raw files")) |>
    dplyr::select(source_file, excel_row, date, site_raw, site_name, has_species_data, status)

  # ---------- master: one row per detection (count > 0) ----------
  obs <- obs |>
    dplyr::left_join(dplyr::select(samples, sample_key, sample_id, survey_id), by = "sample_key")
  master <- obs |>
    dplyr::filter(count > 0) |>
    dplyr::left_join(dplyr::select(surveys, survey_id, source_family, resolution, date, year, site_name,
                                   site_alias_english, island, latitude, longitude, teacher_id),
                     by = "survey_id") |>
    dplyr::left_join(dplyr::select(samples, sample_id, sample_type, group_label, transect, quadrat,
                                   point_m, n_quadrats, sample_flag), by = "sample_id") |>
    dplyr::transmute(survey_id, sample_id, date, year, site_name, site_alias_english, island,
                     latitude, longitude, teacher_id, resolution, sample_type, group_label,
                     transect, quadrat, point_m, n_quadrats, taxon_raw, taxon_clean, analysis_taxon,
                     category, count, points_possible = 25 * n_quadrats, sample_flag,
                     source_family, source_cell) |>
    dplyr::arrange(survey_id, sample_id, taxon_clean)

  qc_flags <- dplyr::bind_rows(
    dplyr::bind_rows(flags),
    samples |> dplyr::filter(!is.na(sample_flag)) |>
      dplyr::transmute(level = "sample", ref = sample_id, flag = sample_flag, detail = sample_key))

  unmatched <- dplyr::bind_rows(
    sites |> dplyr::filter(is.na(site_name)) |>
      dplyr::transmute(type = "site", raw = site_raw, suggestion = site_suggestion),
    taxa  |> dplyr::filter(is.na(taxon_clean)) |>
      dplyr::transmute(type = "taxon", raw = taxon_raw, suggestion = taxon_suggestion))

  list(
    surveys = surveys |> dplyr::select(survey_id, source_family, source_file, resolution, date, year,
                                       site_name, site_alias_english, island, latitude, longitude,
                                       teacher_id, start_time, stop_time, low_tide_time,
                                       n_transects_reported, n_quad_per_t_reported, transect_length_m,
                                       point_spacing_m, n_samples, n_quadrats, src_key),
    samples = samples |> dplyr::select(sample_id, survey_id, sample_type, group_label, transect, quadrat,
                                       point_m, n_quadrats, points_possible, point_sum, n_taxa,
                                       sample_flag, sample_key),
    master = master,
    exclusions = dplyr::bind_rows(exclusions),
    qc_flags = qc_flags,
    crosscheck = crosscheck,
    unmatched = unmatched,
    teacher_key = teacher_key
  )
}
