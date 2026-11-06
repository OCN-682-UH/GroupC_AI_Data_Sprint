# =============================================================================
# 1) Algae-only data sheet
# 2) Interactive Oahu map: native vs non-native algae at each survey location
# 3) Time series (native vs non-native) for each location
# =============================================================================
# Needs in the working directory: opihi_master.csv, algae_origin_lookup.csv
# Packages (install once): tidyverse, leaflet, htmlwidgets

library(tidyverse)
library(leaflet)
library(htmlwidgets)

# ---- 1. ALGAE-ONLY SHEET -------------------------------------------------------
raw    <- read_csv("opihi_master.csv", show_col_types = FALSE)
lookup <- read_csv("algae_origin_lookup.csv", show_col_types = FALSE) |>
  select(taxon_clean, origin_status)

algae <- raw |>
  filter(category == "Algae") |>
  left_join(lookup, by = "taxon_clean") |>
  mutate(origin_group = case_when(
    origin_status %in% c("Native", "Native (likely)")             ~ "Native",
    origin_status %in% c("Non-native", "Non-native (invasive)")   ~ "Non-native",
    TRUE ~ "Unresolved"   # cryptogenic, genus-level, turf, unknown groups
  ))

write_csv(algae, "algae_only.csv")

# Companion "effort" sheet: every quadrat that was surveyed, whether or not it
# had algae. Without this, quadrats with no algae vanish from an algae-only file
# and average cover would be overestimated.
samples <- raw |>
  filter(sample_type == "quadrat") |>          # drops student-group totals
  distinct(sample_id, survey_id, site_name, island, year,
           lat = latitude, lon = longitude, points_possible)

write_csv(samples, "quadrat_samples.csv")

# ---- 2. COVER PER QUADRAT, NATIVE vs NON-NATIVE -------------------------------
# Cover = points hit / points possible x 100. Quadrats with none of a group get 0.
oahu_samples <- samples |> filter(island == "Oahu")

group_points <- algae |>
  filter(sample_type == "quadrat", origin_group %in% c("Native", "Non-native")) |>
  group_by(sample_id, origin_group) |>
  summarise(points = sum(count, na.rm = TRUE), .groups = "drop")

quad_long <- oahu_samples |>
  crossing(origin_group = c("Native", "Non-native")) |>
  left_join(group_points, by = c("sample_id", "origin_group")) |>
  mutate(cover = pmin(100 * replace_na(points, 0) / points_possible, 100))

# ---- 3. INTERACTIVE MAP --------------------------------------------------------
# One circle per group per site, sized by mean % cover over all surveys there.
site_map <- quad_long |>
  group_by(site_name, lat, lon, origin_group) |>
  summarise(mean_cover = mean(cover),
            n_quadrats = n_distinct(sample_id),
            n_years    = n_distinct(year),
            first_year = min(year), last_year = max(year),
            .groups = "drop")

popup_tbl <- site_map |>
  select(site_name, origin_group, mean_cover) |>
  pivot_wider(names_from = origin_group, values_from = mean_cover,
              names_prefix = "cover_") |>
  rename(cover_nn = `cover_Non-native`, cover_nat = cover_Native)

site_map <- site_map |>
  left_join(popup_tbl, by = "site_name") |>
  mutate(
    radius = 4 + 18 * sqrt(mean_cover / max(mean_cover)),   # same scale for both groups
    popup  = paste0("<b>", site_name, "</b><br>",
                    "Native: ", round(cover_nat, 1), "% cover<br>",
                    "Non-native: ", round(cover_nn, 1), "% cover<br>",
                    "Quadrats: ", n_quadrats, "<br>",
                    "Survey years: ", n_years, " (", first_year, "-", last_year, ")")
  )

# Offset the two circles slightly east/west so both stay visible (~400 m).
# Set to 0 to stack them exactly on the site.
offset_deg <- 0.0035

native_pts <- filter(site_map, origin_group == "Native")
nonnat_pts <- filter(site_map, origin_group == "Non-native")

map <- leaflet() |>
  addProviderTiles(providers$Esri.WorldImagery) |>
  addCircleMarkers(data = native_pts, lng = ~lon - offset_deg, lat = ~lat,
                   radius = ~radius, color = "white", weight = 1,
                   fillColor = "#2ca25f", fillOpacity = 0.85,
                   popup = ~popup, label = ~site_name, group = "Native") |>
  addCircleMarkers(data = nonnat_pts, lng = ~lon + offset_deg, lat = ~lat,
                   radius = ~radius, color = "white", weight = 1,
                   fillColor = "#e34a33", fillOpacity = 0.85,
                   popup = ~popup, label = ~site_name, group = "Non-native") |>
  addLayersControl(overlayGroups = c("Native", "Non-native"),
                   options = layersControlOptions(collapsed = FALSE)) |>
  addLegend("bottomright",
            colors = c("#2ca25f", "#e34a33"),
            labels = c("Native (likely)", "Non-native (invasive)"),
            title = "Mean % cover", opacity = 1) |>
  fitBounds(-158.30, 21.25, -157.62, 21.72)    # whole island

map
saveWidget(map, "oahu_native_vs_nonnative_map.html", selfcontained = TRUE)

# ---- 4. TIME SERIES: NATIVE vs NON-NATIVE, ONE GRAPH PER SITE ------------------
site_year <- quad_long |>
  group_by(site_name, year, origin_group) |>
  summarise(n_quadrats = n_distinct(sample_id),
            mean_cover = mean(cover),
            se = sd(cover) / sqrt(n()),
            .groups = "drop")

# Only sites with enough survey years for a line to mean anything
min_years <- 3
sites_to_plot <- site_year |>
  group_by(site_name) |>
  summarise(n_years = n_distinct(year), .groups = "drop") |>
  filter(n_years >= min_years) |>
  pull(site_name)
sites_to_plot

cols <- c("Native" = "#2ca25f", "Non-native" = "#e34a33")

plot_site <- function(site) {
  site_year |>
    filter(site_name == site) |>
    ggplot(aes(year, mean_cover, color = origin_group)) +
    geom_linerange(aes(ymin = pmax(mean_cover - se, 0), ymax = mean_cover + se),
                   alpha = 0.4) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n_quadrats)) +
    scale_color_manual(values = cols, name = NULL) +
    scale_size_continuous(name = "Quadrats", range = c(1.5, 5)) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    labs(title = site, x = "Year", y = "Mean % cover per quadrat",
         caption = "Points = survey years; bars = +/- 1 SE; unresolved taxa excluded") +
    theme_minimal(base_size = 12)
}

# View a single site (swap the name to change locations)
plot_site("Leahi")

# Save one PNG per site
dir.create("site_timeseries", showWarnings = FALSE)
walk(sites_to_plot, \(s)
  ggsave(file.path("site_timeseries",
                   paste0(str_replace_all(tolower(s), "[^a-z0-9]+", "_"), ".png")),
         plot_site(s), width = 7, height = 4, dpi = 300))

# One-page overview of all selected sites
p_all <- site_year |>
  filter(site_name %in% sites_to_plot) |>
  ggplot(aes(year, mean_cover, color = origin_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = n_quadrats)) +
  facet_wrap(~ site_name, scales = "free_y") +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "Year", y = "Mean % cover per quadrat", size = "Quadrats") +
  theme_minimal(base_size = 11)
ggsave("site_timeseries/all_sites_overview.png", p_all,
       width = 11, height = 8, dpi = 300)
