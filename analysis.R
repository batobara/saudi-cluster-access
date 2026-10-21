# ==========================================================================
# After-hours primary care and hospital referral in Saudi Arabia's health
# clusters: a national cross-sectional study of access, placement, and
# opening hours. Batobara A, Alfaraidhy M, Lim C. Submitted to BMJ Global
# Health.
#
# Single-file analysis pipeline. Writes every main-text and supplementary
# table to output/tables (CSV) and every figure to output/figures (PNG);
# README.md maps each file to the manuscript.
#
# Manuscript questions and the sections that answer them:
#   Q1  Access to an open PHC through the week; night alternatives: D1, D2, D2b
#   Q2  Where overnight PHCs should be; conversions and physicians: D2c, D3,
#       D5, D9
#   Q3  Extending hours to midnight versus opening 24/7: D4
#   Q4  Referral within the cluster: 6, 7, 7b, 7c, 7d, 8, D6, D7
#   Tables 1-2 and Figure 3: "Manuscript exhibits for the decision analyses"
#   Sensitivity analyses: 7b, 11, D9
#
# Sections:
#   0. Setup and parameters
#   1. Facilities
#   2. Governorates and census population
#   3. Gridded population (WorldPop scaled to GASTAT 2022 citizen totals)
#   4. Travel times (dodgr street-network routing, no external server)
#   5-8. Descriptive access and referral
#  10. Tables, maps, cluster scorecard, per-cluster distributions, and
#      facility validation against the MoH Statistical Yearbook 2024
#  11. Sensitivity and supplement tables
#  12. Decision analyses (D1-D9) and the manuscript exhibits
#
# In the code, "tertiary" denotes the referral level: the medical cities and
# specialized hospitals (section 1). Some output file names keep the
# numbering of an earlier draft.
#
# Expensive steps cache to data/processed/ and beside the road extract as
# .rds; delete a cache file to force recomputation. Section 4 needs the OSM
# road extract (see its header); everything else is plain CRAN R.
# ==========================================================================

# ---- 0. Setup and parameters ---------------------------------------------

pkgs <- c("sf", "terra", "dplyr", "tidyr", "readr", "stringr", "dodgr",
          "osmextract", "units", "RANN", "curl", "ggplot2",
          "viridis", "ggrepel", "patchwork", "readxl")  # after 12: `::` only
new <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs[1:12], library, character.only = TRUE))

root  <- "."  # project root (the folder containing data/raw); auto-detected below
              # when the script is source()'d or Rscript'ed from elsewhere
if (!dir.exists(file.path(root, "data/raw"))) {
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)  # source() path
  if (is.null(f) || !nzchar(f)) {
    a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(a)) f <- sub("^--file=", "", a[1])              # Rscript path
  }
  cands <- character()
  if (!is.null(f) && nzchar(f)) {
    d <- dirname(normalizePath(f))
    cands <- c(d, dirname(d))        # script dir (deposit) or its parent (R/)
  }
  for (cand in cands)
    if (dir.exists(file.path(cand, "data/raw"))) { root <- cand; break }
  if (!dir.exists(file.path(root, "data/raw")))
    stop("data/raw not found. Either setwd() to the project root (the folder ",
         "containing data/raw) or run this file via source(\"<full path to ",
         "analysis.R>\") and it will locate the root itself.")
}
raw   <- file.path(root, "data/raw")
prc   <- file.path(root, "data/processed")
figs  <- file.path(root, "output/figures")
tabs  <- file.path(root, "output/tables")

PAR <- list(
  crs_m          = 32638,                 # UTM 38N (metric, covers KSA core)
  thresholds_phc = c(20, 30, 45, 60),     # minutes; 30 is primary
  thresholds_hosp= c(30, 60, 90, 120),    # minutes; 60 is primary
  knn_candidates = 10,                    # Euclidean candidates per origin
  pbf_path       = "~/osrm_gcc/gcc-states-latest.osm.pbf",  # Geofabrik Gulf
  hw_keep        = c("motorway", "motorway_link", "trunk", "trunk_link",
                     "primary", "primary_link", "secondary", "secondary_link",
                     "tertiary", "tertiary_link", "unclassified",
                     "residential", "track"),
  track_speed    = 40,   # km/h profile max_speed for unpaved tracks;
                         # dodgr's default car profile treats track as
                         # impassable, so this is set explicitly below.
                         # For motorcar, dodgr 0.4.3 gives roads without an
                         # OSM maxspeed tag the median tagged speed of their
                         # highway class and uses this profile value only
                         # when no road of the class is tagged, so the graph
                         # is the same at 30 or 40 (the cached graph, built
                         # 2026-09-02, was built with 30). The built graph
                         # runs tracks at 40 km/h (2,559,877 of 2,560,807
                         # track edges, 99.9% of track-km; checked 2026-10-09
                         # and 2026-10-20), the Saudi track speed in Weiss et
                         # al. 2020 Nat Med, Supplementary Table 1
  worldpop_url   = paste0("https://data.worldpop.org/GIS/Population/",
                          "Global_2000_2020_Constrained/2020/BSGM/SAU/",
                          "sau_ppp_2020_constrained.tif"),
  # Terrain sensitivity (section 11, run with SLOPE_SENS=1): linear speed
  # penalty from slope_thresh to slope_full grade, up to slope_maxpen
  # reduction. Base model already captures terrain via road geometry
  # (switchbacks) and class speeds; 3.1% of edge-km exceed 6% grade (SRTM).
  slope_thresh   = 0.04,
  slope_full     = 0.12,
  slope_maxpen   = 0.5
)

cache <- function(name, expr) {
  path <- file.path(prc, paste0(name, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  val <- force(expr)
  saveRDS(val, path)
  val
}

# ---- 1. Facilities (supply) ----------------------------------------------

# PHC opening hours as listed in the cluster directories (phc_opening_hours.csv;
# centers without a listing work regular hours, Sunday-Thursday 08:00-16:00).
# Open overnight on weekdays = "h24"; open beyond 16:00 or at weekends =
# "extended". These listed hours supersede the census hours_level column.
opening <- read_csv(file.path(raw, "phc_opening_hours.csv"), show_col_types = FALSE) |>
  select(facility_id, weekday_open, weekday_close, weekend_open, weekend_close)
# Referral level (revised 2026-10-20). No official list of tertiary hospitals
# exists. MoH referral guidance puts medical cities at the top level of care
# (level 4), and the 2014 regulation of medical cities and specialized
# hospitals designates specialized hospitals with them. The referral level
# ("tertiary" in this script) is therefore the cluster medical cities listed by
# the Health Holding Company plus King Fahad Specialist Hospital Dammam
# (referral_services.csv, referral_level = 1). The other directory entries the
# census codes 3ry (cardiac centers, specialist hospitals, centers of
# excellence) are specialized hospitals at the secondary level; their
# co-located entries still count as one site per campus (campus_entry).
# Advanced services (open-heart surgery, radiotherapy, comprehensive oncology,
# transplantation) were audited hospital by hospital in the same file
# (section 7d).
ref_svc <- read_csv(file.path(raw, "referral_services.csv"), show_col_types = FALSE)
referral_ids <- ref_svc$facility_id[ref_svc$referral_level == 1]
fac <- read_csv(file.path(raw, "list_of_healthcare_providers.csv"),
                show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon)) |>
  left_join(opening, by = "facility_id") |>
  mutate(
    tier = case_when(type == "PHC" ~ "primary",
                     facility_id %in% referral_ids ~ "tertiary",
                     TRUE ~ "secondary"),
    campus_entry = type == "3ry",
    hours = case_when(tier != "primary" ~ NA_character_,
                      weekday_open == 0 & weekday_close == 24 ~ "h24",
                      weekday_close > 16 | !is.na(weekend_open) ~ "extended",
                      TRUE ~ "regular"),
    # Service capability flags (scope strings verified against source sites)
    is_general    = tier == "secondary" &
                    str_detect(scope, "General|Emergency"),
    # Departments come from the hospital service audit recorded in the obgyn,
    # peds and psych columns (1 = department present, 0 = absent, blank = no
    # record either way): every general, specialist, and medical-city hospital was
    # checked against official and news sources, because hospital type alone
    # misses general hospitals with these departments. Obstetric and
    # pediatric access is measured to SECONDARY hospitals with the department;
    # the medical cities and specialized hospitals (the referral level) hold
    # these departments for advanced cases and enter as referral destinations
    # (section 7c). Section 7b adds hospitals with no record of the
    # department and tertiary centers as sensitivity analyses.
    serves_obgyn  = tier == "secondary" & obgyn %in% 1,
    serves_peds   = tier == "secondary" & peds %in% 1,
    serves_psych  = tier != "primary" & psych %in% 1,
    tert_obgyn    = tier == "tertiary" & obgyn %in% 1,
    tert_peds     = tier == "tertiary" & peds %in% 1
  )

fac_sf  <- st_as_sf(fac, coords = c("lon", "lat"), crs = 4326)
fac_m   <- st_transform(fac_sf, PAR$crs_m)

phc      <- fac_m |> filter(tier == "primary")
phc_ext  <- phc   |> filter(hours %in% c("extended", "h24"))
phc_24   <- phc   |> filter(hours == "h24")
hosp_gen <- fac_m |> filter(is_general)
hosp_ob  <- fac_m |> filter(serves_obgyn)
hosp_pd  <- fac_m |> filter(serves_peds)
hosp_psy <- fac_m |> filter(serves_psych)
tert     <- fac_m |> filter(tier == "tertiary")

# Co-located directory entries share a campus: single-linkage grouping of
# facilities within thresh_m. Entries the census codes as specialized (3ry)
# are counted as sites (28 entries on 16 campuses); section 10b checks this count against the
# site table behind the manuscript's facility counts.
campus_id <- function(lat, lon, thresh_m = 500) {
  n <- length(lat)
  if (n == 1) return(1L)
  p <- pi / 180
  dm <- outer(seq_len(n), seq_len(n), function(i, j) {
    a <- sin((lat[j] - lat[i]) * p / 2)^2 +
      cos(lat[i] * p) * cos(lat[j] * p) * sin((lon[j] - lon[i]) * p / 2)^2
    2 * 6371000 * asin(pmin(1, sqrt(a)))
  })
  cutree(hclust(as.dist(dm), method = "single"), h = thresh_m)
}
n_tert_sites <- fac |> filter(tier == "tertiary") |>
  distinct(cluster_en, lat, lon) |>
  group_by(cluster_en) |> mutate(site = campus_id(lat, lon)) |> ungroup() |>
  distinct(cluster_en, site) |> nrow()

message(sprintf(
  "Facilities: %d PHC (%d extended, %d 24h) | %d general | %d OB | %d peds | %d psych | %d tertiary",
  nrow(phc), nrow(phc_ext) - nrow(phc_24), nrow(phc_24), nrow(hosp_gen),
  nrow(hosp_ob), nrow(hosp_pd), nrow(hosp_psy), nrow(tert)))

# ---- 2. Governorates and census population (demand) ----------------------

gov <- st_read(file.path(raw, "governorate/Governorate.gpkg"), quiet = TRUE) |>
  st_transform(PAR$crs_m)
reg <- st_read(file.path(raw, "Regions/Regions.shp"), quiet = TRUE) |>
  st_transform(PAR$crs_m)

# Arabic name normalization for GASTAT <-> shapefile matching
norm_ar <- function(x) {
  x |>
    str_remove_all("\\(.*?\\)") |>            # drop parentheticals
    str_replace_all("[أإآ]", "ا") |>  # hamza forms -> alef
    str_replace_all("ة", "ه") |>    # taa marbuta -> haa
    str_replace_all("ى", "ي") |>    # alef maqsura -> yaa
    str_remove_all("[ً-ْ]") |>      # diacritics
    str_remove_all("\\s")           # census has stray spaces (e.g. Al Amwah)
}
gov$gov_key <- norm_ar(gov$Gov_AR)

# Census 2022: population by nationality x gender per city, summed to gov.
pop_city <- read_csv(
  file.path(raw, "PopulationbyNationalitybyRegionGovernorateCityandNationalityARCSV.csv"),
  show_col_types = FALSE)
names(pop_city) <- c("nationality", "governorate", "city", "region",
                     "gender", "pop")

pop_gov <- pop_city |>
  mutate(gov_key = norm_ar(governorate),
         saudi   = nationality == "سعودي") |>
  group_by(gov_key) |>
  summarise(pop_total = sum(pop),
            pop_saudi = sum(pop[saudi]), .groups = "drop")

# Census 2022: single-year age x nationality x gender per governorate.
pop_age <- read_csv(
  file.path(raw, "PopulationbydetailedAgebyRegionGovernorateNationalityandGenderARCSV.csv"),
  show_col_types = FALSE)
names(pop_age) <- c("nationality", "age", "age5", "governorate", "region",
                    "gender", "pop")

pop_dem <- pop_age |>
  mutate(gov_key = norm_ar(governorate),
         saudi   = nationality == "سعودي",
         female  = gender == "أنثى",
         age_n   = suppressWarnings(as.integer(age))) |>
  group_by(gov_key) |>
  summarise(
    women_1549_saudi = sum(pop[saudi & female & age_n >= 15 & age_n <= 49],
                           na.rm = TRUE),
    child_u15_saudi  = sum(pop[saudi & age_n < 15], na.rm = TRUE),
    adult_18p_saudi  = sum(pop[saudi & age_n >= 18], na.rm = TRUE),
    .groups = "drop")

gov_pop <- gov |>
  left_join(pop_gov, by = "gov_key") |>
  left_join(pop_dem, by = "gov_key")

unmatched <- gov_pop |> st_drop_geometry() |>
  filter(is.na(pop_total)) |> select(Gov_AR, Gov_EN, gov_key)
if (nrow(unmatched) > 0) {
  write_csv(unmatched, file.path(prc, "unmatched_governorates.csv"))
  warning(sprintf(
    "%d governorates unmatched to census names -> data/processed/unmatched_governorates.csv (fix names there, then rerun)",
    nrow(unmatched)))
}

# ---- 3. Gridded population -----------------------------------------------
# WorldPop 2020 constrained (100 m; population only in settled areas),
# aggregated to 1 km, supplies the within-governorate spatial distribution;
# cells are rescaled so each governorate sums to its GASTAT 2022 totals
# (total, Saudi, women 15-49 Saudi, children <15 Saudi).

wp_path <- file.path(raw, "sau_ppp_2020_constrained.tif")
if (!file.exists(wp_path)) {
  message("Downloading WorldPop constrained raster (~60 MB)...")
  curl::curl_download(PAR$worldpop_url, wp_path, quiet = TRUE)
}

grid <- cache("grid_population", {
  r  <- terra::aggregate(terra::rast(wp_path), fact = 10, fun = "sum",
                         na.rm = TRUE)   # 100 m -> 1 km
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "wp"
  df <- df[df$wp > 0, ]
  pts <- st_as_sf(df, coords = c("x", "y"), crs = 4326) |>
    st_transform(PAR$crs_m) |>
    st_join(gov_pop[, c("Gov_ID", "Gov_EN", "gov_key", "pop_total",
                        "pop_saudi", "women_1549_saudi", "child_u15_saudi")],
            join = st_within) |>
    filter(!is.na(Gov_ID))
  pts |>
    group_by(Gov_ID) |>
    mutate(w = wp / sum(wp)) |>
    ungroup() |>
    mutate(cell_total = w * pop_total,
           cell_saudi = w * pop_saudi,
           cell_women = w * women_1549_saudi,
           cell_child = w * child_u15_saudi,
           cell_id    = row_number()) |>
    select(cell_id, Gov_ID, Gov_EN, starts_with("cell_"))
})
# Saudi adults (18+) for the psychiatric-care denominator. Derived outside
# the cache: within a governorate the grid shares age structure, so
# cell_adult = cell_saudi x the governorate's Saudi adult share.
grid <- grid |>
  left_join(gov_pop |> st_drop_geometry() |>
              transmute(Gov_ID, adult_share = adult_18p_saudi / pop_saudi),
            by = "Gov_ID") |>
  mutate(cell_adult = cell_saudi * adult_share) |>
  select(-adult_share)

message(sprintf("Population grid: %s cells | Saudi total %.0f",
                format(nrow(grid), big.mark = ","),
                sum(grid$cell_saudi, na.rm = TRUE)))

# ---- 4. Travel times (dodgr street-network routing) ----------------------
# One-time input: download the Geofabrik Gulf extract to PAR$pbf_path
#   https://download.geofabrik.de/asia/gcc-states-latest.osm.pbf  (~250 MB)
# Routing runs entirely in R (dodgr): car travel times from OSM road
# geometry, the 'motorcar' weighting profile, and maxspeed/oneway tags.
# No server, no Docker, no Java. The road graph and the travel-time matrices
# are cached beside the pbf.

if (!file.exists(path.expand(PAR$pbf_path)))
  stop("Road extract not found at ", PAR$pbf_path,
       " - download it (see Section 4 header), then rerun.")

# Road ways are needed only to (re)build the routing graph, so they are read
# lazily inside the graph step below.
load_roads <- function() {
  roads <- cache("roads_saudi", {
    q <- sprintf(   # SELECT * keeps the geometry column (named columns drop it)
      "SELECT * FROM lines WHERE highway IN (%s)",
      paste(sprintf("'%s'", PAR$hw_keep), collapse = ", "))
    rd <- osmextract::oe_read(path.expand(PAR$pbf_path), layer = "lines",
                              extra_tags = c("oneway", "maxspeed", "lanes"),
                              query = q, quiet = TRUE)
    ksa <- st_union(gov) |> st_buffer(20000) |> st_transform(4326)
    st_filter(rd, ksa)              # keep Saudi roads (Gulf extract is wider)
  })
  message(sprintf("Road network: %s ways", format(nrow(roads), big.mark = ",")))
  roads
}

# Graph cache lives beside the pbf, outside the project folder: the file is
# ~1 GB. Keep every road component
# with >= 100 edges, not just the largest: inhabited islands with their own
# road networks (in practice one: the Farasan Islands in Jazan; the other
# retained components are disconnected mainland road pockets) would otherwise
# have their cells and facilities snapped across the sea. Cross-component
# routes return NA (unreachable).
graph_path <- file.path(dirname(path.expand(PAR$pbf_path)), "dodgr_graph.rds")
graph <- if (file.exists(graph_path)) readRDS(graph_path) else {
  # Motorcar profile with unpaved tracks passable (speed actually used: see
  # the note at PAR$track_speed), weighted like residential streets rather
  # than excluded. Must go through
  # a profile FILE: passing a profile data.frame makes dodgr skip travel
  # time computation ("graph has no time column").
  wpf <- file.path(prc, "wt_profiles.json")
  if (!file.exists(wpf)) {
    dodgr::write_dodgr_wt_profile(file.path(prc, "wt_profiles"))
    j <- jsonlite::fromJSON(wpf)
    i <- j$weighting_profiles$name == "motorcar" &
         j$weighting_profiles$way  == "track"
    j$weighting_profiles$value[i]     <- 0.5
    j$weighting_profiles$max_speed[i] <- PAR$track_speed
    jsonlite::write_json(j, wpf, pretty = TRUE)
  }
  g <- dodgr::weight_streetnet(load_roads(), wt_profile = "motorcar",
                               wt_profile_file = wpf)
  keep <- names(which(table(g$component) >= 100))
  g <- g[g$component %in% as.integer(keep), ]
  saveRDS(g, graph_path)
  g
}
verts   <- dodgr::dodgr_vertices(graph)
vert_xy <- as.matrix(verts[, c("x", "y")])

# Points attach only to vertices that can both reach and be reached from the
# main road network, or to vertices of the separate components kept above
# (such as the Farasan Islands). Otherwise the nearest vertex can be a one-way
# dead end in the OSM data, which left 75,000 citizens unroutable before
# 2026-10-21. Reachability is tested from one vertex in central Riyadh, on the
# graph and on its reverse.
mn_path <- file.path(dirname(graph_path), "main_network_vertices.rds")
v_ok <- if (file.exists(mn_path)) readRDS(mn_path) else NULL
if (is.null(v_ok) || length(v_ok) != nrow(verts)) {
  v_ok <- local({
    v_comp <- graph$component[match(verts$id, graph$from_id)]
    miss <- is.na(v_comp)
    v_comp[miss] <- graph$component[match(verts$id[miss], graph$to_id)]
    in_main <- v_comp %in% as.integer(names(which.max(table(graph$component))))
    hub <- verts$id[in_main][RANN::nn2(vert_xy[in_main, ], cbind(46.71, 24.69), k = 1)$nn.idx[1]]
    dodgr::dodgr_cache_off()
    swap <- function(d, a, b) { tmp <- d[[a]]; d[[a]] <- d[[b]]; d[[b]] <- tmp; d }
    rev_g <- swap(swap(swap(graph, "from_id", "to_id"), "from_lon", "to_lon"), "from_lat", "to_lat")
    ids <- verts$id[in_main]
    t_in  <- dodgr::dodgr_times(graph, from = hub, to = ids)     # hub to vertex
    t_out <- dodgr::dodgr_times(rev_g, from = hub, to = ids)     # vertex to hub
    ok <- !in_main
    ok[in_main] <- is.finite(t_in[1, match(ids, colnames(t_in))]) &
                   is.finite(t_out[1, match(ids, colnames(t_out))])
    message(sprintf("  main network: %.2f%% of main-component vertices reach and are reached",
                    100 * mean(ok[in_main])))
    ok
  })
  saveRDS(v_ok, mn_path)
}
ok_idx <- which(v_ok)

# Snap sf points to the nearest usable road-graph vertex; report worst snap distances.
snap_verts <- function(x, label) {
  ll <- st_coordinates(st_transform(x, 4326))
  nn <- RANN::nn2(vert_xy[ok_idx, , drop = FALSE], ll, k = 1)
  d_km <- nn$nn.dists[, 1] * 111    # degrees -> km (approximate)
  message(sprintf("  snap %-8s: median %.2f km, max %.1f km, >5 km: %d of %d",
                  label, median(d_km), max(d_km), sum(d_km > 5), nrow(x)))
  verts$id[ok_idx[nn$nn.idx[, 1]]]
}

# Travel time from each origin to its nearest destination (minutes).
# Euclidean k-NN pruning, then chunked dodgr many-to-many time matrices.
tt_nearest <- function(origins, dests, k = PAR$knn_candidates,
                       chunk = 1000, label = "dests") {
  k   <- min(k, nrow(dests))
  nn  <- RANN::nn2(st_coordinates(dests), st_coordinates(origins), k = k)
  o_v <- snap_verts(origins, "origins")
  d_v <- snap_verts(dests, label)
  out_t <- rep(NA_real_, nrow(origins))
  out_i <- rep(NA_integer_, nrow(origins))
  starts <- seq(1, nrow(origins), by = chunk)
  for (s in starts) {
    idx  <- s:min(s + chunk - 1, nrow(origins))
    cand <- sort(unique(as.vector(nn$nn.idx[idx, , drop = FALSE])))
    tmat <- dodgr::dodgr_times(graph, from = unique(o_v[idx]),
                               to = unique(d_v[cand])) / 60
    for (j in seq_along(idx)) {
      # index by vertex id: co-snapped points share matrix rows/columns
      tj  <- tmat[match(o_v[idx[j]], rownames(tmat)),
                  match(d_v[nn$nn.idx[idx[j], ]], colnames(tmat))]
      tj[!is.finite(tj)] <- NA
      if (all(is.na(tj))) next
      b   <- which.min(tj)
      out_t[idx[j]] <- tj[b]
      out_i[idx[j]] <- nn$nn.idx[idx[j], b]
    }
    message(sprintf("  ...%d / %d origins", max(idx), nrow(origins)))
  }
  list(minutes = out_t, dest_row = out_i)
}

grid_pts <- grid  # sf point layer

# Hospital-to-cell travel times for every hospital, computed once on the
# reversed road graph (one shortest-path tree per hospital gives each cell's
# time TO that hospital). The nearest hospital of any service subset is then
# a column minimum, so specialty services and their sensitivity sets need no
# further routing. Agrees with tt_nearest() for general hospitals in all but
# 0.01% of cells (differences only where k-NN pruning missed the road-nearest).
hosp_all <- fac_m |> filter(tier != "primary")
# ~150 MB, so cached beside the road graph
hc_path <- file.path(dirname(graph_path), "hospital_cell_times.rds")
hc <- if (file.exists(hc_path)) readRDS(hc_path) else NULL
if (is.null(hc) || !identical(rownames(hc), hosp_all$facility_id) ||
    ncol(hc) != nrow(grid_pts)) {
  hc <- local({
    dodgr::dodgr_cache_off()      # reversed copy must not reuse cached graph
    swap <- function(d, a, b) { tmp <- d[[a]]; d[[a]] <- d[[b]]; d[[b]] <- tmp; d }
    rev_g <- swap(swap(swap(graph, "from_id", "to_id"),
                       "from_lon", "to_lon"), "from_lat", "to_lat")
    cv <- snap_verts(grid_pts, "cells")
    hv <- snap_verts(hosp_all, "hospitals")
    tm <- dodgr::dodgr_times(rev_g, from = unique(hv), to = unique(cv)) / 60
    m  <- tm[match(hv, rownames(tm)), match(cv, colnames(tm)), drop = FALSE]
    m[!is.finite(m)] <- NA
    rownames(m) <- hosp_all$facility_id
    m
  })
  saveRDS(hc, hc_path)
}
# PHC-to-cell travel times, built the same way (~15 min, cached beside the
# road graph): any set of PHCs open at a given hour is evaluated from it.
pc_path <- file.path(dirname(graph_path), "phc_cell_times.rds")
pc <- if (file.exists(pc_path)) readRDS(pc_path) else NULL
# PHC x cell travel times on graph g (also used on the terrain-penalised
# graph in section 11).
phc_cell_matrix <- function(g) local({
    dodgr::dodgr_cache_off()
    swap <- function(d, a, b) { tmp <- d[[a]]; d[[a]] <- d[[b]]; d[[b]] <- tmp; d }
    rev_g <- swap(swap(swap(g, "from_id", "to_id"),
                       "from_lon", "to_lon"), "from_lat", "to_lat")
    cv  <- snap_verts(grid_pts, "cells")
    pv  <- snap_verts(phc, "phc")
    ucv <- unique(cv)
    m <- matrix(NA_real_, nrow(phc), nrow(grid_pts),
                dimnames = list(phc$facility_id, NULL))
    for (ix in split(seq_len(nrow(phc)), ceiling(seq_len(nrow(phc)) / 250))) {
      tm <- dodgr::dodgr_times(rev_g, from = unique(pv[ix]), to = ucv) / 60
      mm <- tm[match(pv[ix], rownames(tm)), match(cv, colnames(tm)), drop = FALSE]
      mm[!is.finite(mm)] <- NA
      m[ix, ] <- mm
    }
    m
  })
if (is.null(pc) || !identical(rownames(pc), phc$facility_id) ||
    ncol(pc) != nrow(grid_pts)) {
  pc <- phc_cell_matrix(graph)
  saveRDS(pc, pc_path)
}
colmin <- function(M) {
  if (!nrow(M)) return(rep(NA_real_, ncol(M)))
  x <- suppressWarnings(apply(M, 2, min, na.rm = TRUE))
  x[!is.finite(x)] <- NA
  x
}
t_ext_pc <- colmin(pc[phc$hours %in% c("extended", "h24"), , drop = FALSE])
t_24_pc  <- colmin(pc[phc$hours == "h24", , drop = FALSE])

nearest_of <- function(ids) {
  s <- hc[intersect(ids, rownames(hc)), , drop = FALSE]
  apply(s, 2, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
}

# Specialty routes record the facility set they were computed for, so a
# revised hospital service audit re-routes only the affected service.
tt_spec <- function(d, label) list(minutes  = nearest_of(d$facility_id),
                                   dest_ids = sort(d$facility_id))
spec_sets <- list(obgyn = hosp_ob, peds = hosp_pd, psych = hosp_psy)

tt <- cache("travel_times_grid", {
  list(
    phc_all = tt_nearest(grid_pts, phc),
    phc_ext = tt_nearest(grid_pts, phc_ext),
    phc_24  = tt_nearest(grid_pts, phc_24),
    gen     = tt_nearest(grid_pts, hosp_gen),
    obgyn   = tt_spec(hosp_ob,  "obgyn"),
    peds    = tt_spec(hosp_pd,  "peds"),
    psych   = tt_spec(hosp_psy, "psych"),
    tert    = tt_nearest(grid_pts, tert)
  )
})
for (nm in names(spec_sets)) {
  if (!identical(tt[[nm]]$dest_ids, sort(spec_sets[[nm]]$facility_id))) {
    message(sprintf("Re-routing %s: destination set changed (%d hospitals)",
                    nm, nrow(spec_sets[[nm]])))
    tt[[nm]] <- tt_spec(spec_sets[[nm]], nm)
    saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
  }
}
# The referral level records its facility set too (it changed on 2026-10-20,
# from 28 tertiary-coded facilities to the medical cities and specialized
# hospitals), so a revised set is re-routed rather than read from the cache.
if (!identical(tt$tert$dest_ids, sort(tert$facility_id))) {
  message(sprintf("Re-routing referral level: destination set changed (%d hospitals)", nrow(tert)))
  tt$tert <- tt_nearest(grid_pts, tert, label = "referral")
  tt$tert$dest_ids <- sort(tert$facility_id)
  saveRDS(tt, file.path(prc, "travel_times_grid.rds"))
}

acc <- grid |>
  mutate(t_phc   = tt$phc_all$minutes,  phc_row  = tt$phc_all$dest_row,
         t_ext   = t_ext_pc,                  # listed hours (see section 1)
         t_24    = t_24_pc,
         t_gen   = tt$gen$minutes,      gen_row  = tt$gen$dest_row,
         t_obgyn = tt$obgyn$minutes,
         t_peds  = tt$peds$minutes,
         t_psych = tt$psych$minutes,
         t_tert  = tt$tert$minutes,     tert_row = tt$tert$dest_row,
         cluster_phc  = phc$cluster_en[phc_row],       # registration proxy
         cluster_gen  = hosp_gen$cluster_en[gen_row],
         cluster_tert = tert$cluster_en[tert_row])

# Population-weighted summary helpers
wq <- function(x, w, p) {                       # weighted quantile
  ok <- !is.na(x) & !is.na(w); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]
  x[which.max(cumsum(w) / sum(w) >= p)]
}
cov_pct <- function(x, w, thr) sum(w[!is.na(x) & x <= thr]) / sum(w) * 100

# ---- 5. PHC access by hours category -------------------------------------

q1_national <- bind_rows(lapply(
  list(all = "t_phc", extended_or_24h = "t_ext", h24_only = "t_24"),
  function(v) tibble(
    median_min = wq(acc[[v]], acc$cell_saudi, .5),
    p90_min    = wq(acc[[v]], acc$cell_saudi, .9),
    !!!setNames(
      lapply(PAR$thresholds_phc,
             function(th) cov_pct(acc[[v]], acc$cell_saudi, th)),
      paste0("pct_within_", PAR$thresholds_phc, "min")))
), .id = "phc_tier")

q1_cluster <- acc |>
  st_drop_geometry() |>
  group_by(cluster = cluster_phc) |>
  summarise(saudi_pop     = sum(cell_saudi),
            median_min    = wq(t_phc, cell_saudi, .5),
            pct_within_30 = cov_pct(t_phc, cell_saudi, 30),
            pct_within_30_ext = cov_pct(t_ext, cell_saudi, 30),
            pct_within_30_24h = cov_pct(t_24, cell_saudi, 30),
            .groups = "drop")


# ---- 6. Q4: referral from PHC to general hospital ------------------------
# Nearest general hospital in the SAME cluster (stated proxy for the
# referral destination; referral directories are not published).

q2 <- cache("q2_referral", {
  res <- lapply(split(seq_len(nrow(phc)), phc$cluster_en), function(ix) {
    dst <- hosp_gen |> filter(cluster_en == phc$cluster_en[ix[1]])
    if (nrow(dst) == 0) return(tibble(phc_row = ix, t_ref = NA_real_))
    r <- tt_nearest(phc[ix, ], dst, k = min(5, nrow(dst)))
    tibble(phc_row = ix, t_ref = r$minutes)
  })
  bind_rows(res) |> arrange(phc_row)
})
phc$t_referral <- q2$t_ref

q2_cluster <- phc |> st_drop_geometry() |>
  group_by(cluster_en) |>
  summarise(n_phc = n(),
            median_ref_min = median(t_referral, na.rm = TRUE),
            p90_ref_min    = quantile(t_referral, .9, na.rm = TRUE),
            n_over_60      = sum(t_referral > 60, na.rm = TRUE),
            .groups = "drop")
# People exposed: each citizen takes the referral journey of their nearest PHC
# (the registration proxy used for cluster populations).
acc$t_referral <- phc$t_referral[acc$phc_row]
q2_pop <- acc |> st_drop_geometry() |> filter(!is.na(cluster_phc)) |>
  group_by(cluster_en = cluster_phc) |>
  summarise(citizens_ref_over_60 = sum(cell_saudi[t_referral > 60 & !is.na(t_referral)],
                                       na.rm = TRUE),
            pct_citizens_ref_over_60 = 100 * citizens_ref_over_60 / sum(cell_saudi, na.rm = TRUE),
            pop_median_ref_min = wq(t_referral, cell_saudi, .5),
            .groups = "drop")
q2_cluster <- q2_cluster |> left_join(q2_pop, by = "cluster_en")
message(sprintf(
  "  Q2: %d of %d PHCs refer >60 min; they serve %s citizens (%.2f%%)",
  sum(phc$t_referral > 60, na.rm = TRUE), nrow(phc),
  format(round(sum(q2_pop$citizens_ref_over_60)), big.mark = ","),
  100 * sum(q2_pop$citizens_ref_over_60) / sum(acc$cell_saudi, na.rm = TRUE)))

# ---- 7. Obstetric, pediatric, and psychiatric access (demand-matched) ----

q3 <- tibble(
  service    = c("general", "obstetric", "pediatric", "psychiatric"),
  denom      = c("cell_saudi", "cell_women", "cell_child", "cell_adult"),
  timevar    = c("t_gen", "t_obgyn", "t_peds", "t_psych")) |>
  rowwise() |>
  mutate(
    median_min    = wq(acc[[timevar]], acc[[denom]], .5),
    p90_min       = wq(acc[[timevar]], acc[[denom]], .9),
    pct_within_60 = cov_pct(acc[[timevar]], acc[[denom]], 60),
    extra_vs_general_median =
      median_min - wq(acc$t_gen, acc[[denom]], .5)) |>
  ungroup()

# ---- 7b. Sensitivity: hospital service audit -----------------------------
# Main results: secondary hospitals with a documented department. Alternative
# destination sets: hospitals with no record either way counted as having it
# (upper bound on access); tertiary centers added (their
# obstetric and pediatric departments serve advanced cases); and designated
# specialty hospitals only (maternity and children, children's, psychiatric).
fac_d <- st_drop_geometry(hosp_all)
sec_d <- fac_d$tier == "secondary"
svc_sets <- list(
  obstetric = list(denom = "cell_women", sets = list(
    documented         = hosp_ob$facility_id,
    documented_or_no_record = fac_d$facility_id[sec_d & (fac_d$obgyn %in% 1 | is.na(fac_d$obgyn))],
    including_tertiary = fac_d$facility_id[fac_d$obgyn %in% 1],
    designated_only    = fac_d$facility_id[sec_d & str_detect(fac_d$scope, "Maternity|Obstetrics")])),
  pediatric = list(denom = "cell_child", sets = list(
    documented         = hosp_pd$facility_id,
    documented_or_no_record = fac_d$facility_id[sec_d & (fac_d$peds %in% 1 | is.na(fac_d$peds))],
    including_tertiary = fac_d$facility_id[fac_d$peds %in% 1],
    designated_only    = fac_d$facility_id[sec_d & str_detect(fac_d$scope, "Children|Neonatal")])),
  psychiatric = list(denom = "cell_adult", sets = list(
    documented         = hosp_psy$facility_id,
    documented_or_no_record = fac_d$facility_id[fac_d$psych %in% 1 | is.na(fac_d$psych)],
    designated_only    = fac_d$facility_id[str_detect(fac_d$scope, "Psychiatry") |
                                             (fac_d$scope == "Medical City" & fac_d$psych %in% 1)])))
q3_sens <- bind_rows(lapply(names(svc_sets), function(s) {
  w <- acc[[svc_sets[[s]]$denom]]
  bind_rows(lapply(names(svc_sets[[s]]$sets), function(k) {
    ids <- svc_sets[[s]]$sets[[k]]
    t   <- nearest_of(ids)
    tibble(service = s, destination_set = k, n_hospitals = length(ids),
           median_min = wq(t, w, .5), p90_min = wq(t, w, .9),
           pct_within_30 = cov_pct(t, w, 30), pct_within_60 = cov_pct(t, w, 60),
           pct_within_120 = cov_pct(t, w, 120),
           extra_vs_general_median = wq(t, w, .5) - wq(acc$t_gen, w, .5))
  }))
}))
write_csv(q3_sens |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "q3_service_audit_sensitivity.csv"))

# ---- 7c. Q4: referral from general hospital to the referral level --------
# Referrals stay inside the accountable cluster, as in the PHC -> hospital
# step: each resident's chain runs through the nearest secondary hospital of
# their own cluster (cluster_phc) to the nearest tertiary center of that same
# cluster. General hospitals refer to any tertiary center (variant: medical
# cities and specialist hospitals only, since the other tertiary centers are
# single-specialty); obstetric and pediatric cases go to tertiary centers
# holding the department. Where the cluster has no such center, referral
# leaves the cluster network (for example to King Faisal Specialist Hospital
# and Research Centre, outside the facility census): those residents are
# counted as referred outside, not given a travel time, and the travel-time
# summaries cover residents referred within their cluster.
tert_ob <- fac_m |> filter(tert_obgyn)
tert_pd <- fac_m |> filter(tert_peds)
tert_ms <- tert  |> filter(scope %in% c("Medical City", "Specialist"))
ref_mat <- function(from_sf, to_sf) {                  # hospital x tertiary minutes
  fv <- snap_verts(from_sf, "referring")
  tv <- snap_verts(to_sf, "tertiary")
  tm <- dodgr::dodgr_times(graph, from = unique(fv), to = unique(tv)) / 60
  tm <- tm[match(fv, rownames(tm)), match(tv, colnames(tm)), drop = FALSE]
  tm[!is.finite(tm)] <- NA
  tm
}
pct_over <- function(x, w, thr) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  sum(w[ok & x > thr]) / sum(w[ok]) * 100
}
ref_within <- function(from_sf, to_sf) {
  tm  <- ref_mat(from_sf, to_sf)
  leg <- vapply(seq_len(nrow(tm)), function(i) {
    x <- tm[i, to_sf$cluster_en == from_sf$cluster_en[i]]
    if (!length(x) || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  }, numeric(1))
  tcell <- rep(NA_real_, nrow(acc))
  outside <- rep(FALSE, nrow(acc))
  for (k in unique(na.omit(acc$cluster_phc))) {
    ci <- which(acc$cluster_phc == k)
    hi <- which(from_sf$cluster_en == k)
    if (!any(to_sf$cluster_en == k)) { outside[ci] <- TRUE; next }
    if (!length(hi)) next
    s <- hc[from_sf$facility_id[hi], ci, drop = FALSE]   # nearest same-cluster hospital
    b <- apply(s, 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
    tcell[ci] <- leg[hi][b]
  }
  list(leg = leg, tcell = tcell, outside = outside,
       hosp_outside = !from_sf$cluster_en %in% to_sf$cluster_en)
}
ref_specs <- list(
  list(service = "general", from = hosp_gen, to = tert, denom = "cell_saudi"),
  list(service = "general_to_medcity_specialist", from = hosp_gen, to = tert_ms,
       denom = "cell_saudi"),
  list(service = "obstetric", from = hosp_ob, to = tert_ob, denom = "cell_women"),
  list(service = "pediatric", from = hosp_pd, to = tert_pd, denom = "cell_child"))
ref_cells <- list()
q3_ref <- bind_rows(lapply(ref_specs, function(s) {
  r <- ref_within(s$from, s$to)
  ref_cells[[s$service]] <<- r$tcell
  w  <- acc[[s$denom]]
  w0 <- ifelse(is.na(w) | is.na(acc$cluster_phc), 0, w)
  tibble(service = s$service, n_secondary = nrow(s$from),
         n_tertiary = nrow(s$to),
         n_clusters_with_tertiary = n_distinct(s$to$cluster_en),
         hosp_n_referred_outside  = sum(r$hosp_outside),
         hosp_median_min  = median(r$leg, na.rm = TRUE),
         hosp_p90_min     = unname(quantile(r$leg, .9, na.rm = TRUE)),
         hosp_n_over_120  = sum(r$leg > 120, na.rm = TRUE),
         pop_pct_referred_outside = 100 * sum(w0[r$outside]) / sum(w0),
         pop_median_min   = wq(r$tcell, w, .5),
         pop_p90_min      = wq(r$tcell, w, .9),
         pop_pct_over_60  = pct_over(r$tcell, w, 60),
         pop_pct_over_120 = pct_over(r$tcell, w, 120))
}))
write_csv(q3_ref |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "q3_referral_to_tertiary.csv"))
q3_ref_cluster <- acc |> st_drop_geometry() |>
  mutate(ref_gen = ref_cells$general,
         ref_gen_ms = ref_cells$general_to_medcity_specialist,
         ref_ob = ref_cells$obstetric, ref_pd = ref_cells$pediatric) |>
  filter(!is.na(cluster_phc)) |>
  group_by(cluster = cluster_phc) |>
  summarise(ref_gen_median = wq(ref_gen, cell_saudi, .5),
            ref_gen_pct_over_120 = pct_over(ref_gen, cell_saudi, 120),
            ref_gen_ms_median = wq(ref_gen_ms, cell_saudi, .5),
            ref_gen_ms_pct_over_120 = pct_over(ref_gen_ms, cell_saudi, 120),
            ref_ob_median = wq(ref_ob, cell_women, .5),
            ref_ob_pct_over_120 = pct_over(ref_ob, cell_women, 120),
            ref_pd_median = wq(ref_pd, cell_child, .5),
            ref_pd_pct_over_120 = pct_over(ref_pd, cell_child, 120),
            .groups = "drop") |>
  mutate(gen_tertiary_in_cluster = cluster %in% tert$cluster_en,
         ms_tertiary_in_cluster  = cluster %in% tert_ms$cluster_en,
         ob_tertiary_in_cluster  = cluster %in% tert_ob$cluster_en,
         pd_tertiary_in_cluster  = cluster %in% tert_pd$cluster_en)
write_csv(q3_ref_cluster |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "q3_referral_to_tertiary_by_cluster.csv"))

# ---- 7d. Q4: referral by advanced service ---------------------------------
# Medical cities are not interchangeable for every service, and some advanced
# services sit in other hospitals. Each advanced service was audited hospital by
# hospital from public sources (referral_services.csv): open-heart surgery,
# radiotherapy (linear accelerator on site), comprehensive oncology (inpatient
# medical oncology with surgical oncology or a named cancer center), and organ
# or stem-cell transplantation. As in 7c, each resident's nearest general
# hospital of their own cluster refers to the nearest hospital of the same
# cluster providing the service; residents of clusters without it are referred
# outside the cluster network. Direct travel from home to the nearest provider
# in any cluster is reported alongside.
svc_names <- c(open_heart_surgery = "Open-heart surgery", radiotherapy = "Radiotherapy",
               comprehensive_oncology = "Comprehensive oncology",
               transplantation = "Transplantation")
w_rt <- ifelse(is.na(acc$cell_saudi) | is.na(acc$cluster_phc), 0, acc$cell_saudi)
svc_runs <- lapply(names(svc_names), function(k) {
  to <- hosp_all |> filter(facility_id %in% ref_svc$facility_id[ref_svc[[k]] == 1])
  stopifnot(nrow(to) == sum(ref_svc[[k]] == 1))
  r  <- ref_within(hosp_gen, to)
  home <- nearest_of(to$facility_id)
  list(
    nat = tibble(service = svc_names[[k]], n_hospitals = nrow(to),
                 n_clusters_with = n_distinct(to$cluster_en),
                 pct_citizens_cluster_without = 100 * sum(w_rt[r$outside]) / sum(w_rt),
                 ref_median_min = wq(r$tcell, acc$cell_saudi, .5),
                 ref_pct_over_120 = pct_over(r$tcell, acc$cell_saudi, 120),
                 home_median_min = wq(home, acc$cell_saudi, .5),
                 home_pct_within_120 = cov_pct(home, acc$cell_saudi, 120)),
    cl = acc |> st_drop_geometry() |> mutate(tc = r$tcell) |>
      filter(!is.na(cluster_phc)) |>
      group_by(cluster = cluster_phc) |>
      summarise(median_min = wq(tc, cell_saudi, .5), .groups = "drop") |>
      mutate(service = svc_names[[k]], in_cluster = cluster %in% to$cluster_en))
})
g_ref <- q3_ref |> filter(service == "general")
ref_by_service <- bind_rows(
  tibble(service = "Medical city or specialized hospital", n_hospitals = nrow(tert),
         n_clusters_with = n_distinct(tert$cluster_en),
         pct_citizens_cluster_without = g_ref$pop_pct_referred_outside,
         ref_median_min = g_ref$pop_median_min, ref_pct_over_120 = g_ref$pop_pct_over_120,
         home_median_min = wq(acc$t_tert, acc$cell_saudi, .5),
         home_pct_within_120 = cov_pct(acc$t_tert, acc$cell_saudi, 120)),
  bind_rows(lapply(svc_runs, `[[`, "nat")))
ref_by_service_cl <- bind_rows(lapply(svc_runs, `[[`, "cl"))
write_csv(ref_by_service |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "referral_by_service.csv"))
write_csv(ref_by_service_cl |> mutate(median_min = round(median_min, 1)),
          file.path(tabs, "referral_by_service_cluster.csv"))
message(sprintf("  7d: clusters with each service: %s; citizens in clusters without: %s",
                paste(ref_by_service$n_clusters_with, collapse = "/"),
                paste(round(ref_by_service$pct_citizens_cluster_without, 1), collapse = "/")))

# ---- 8. Q4: referral level by cluster and cross-cluster nearest facility --

q4_tert_access <- acc |> st_drop_geometry() |>
  group_by(cluster = cluster_phc) |>
  summarise(median_tert_min = wq(t_tert, cell_saudi, .5),
            pct_within_120  = cov_pct(t_tert, cell_saudi, 120),
            .groups = "drop")

# Boundary alignment: assigned cluster proxied by nearest-PHC cluster
# (the model of care registers every person with a PHC); misalignment = nearest
# general hospital or tertiary center belongs to a different cluster.
q4_alignment <- acc |> st_drop_geometry() |>
  filter(!is.na(cluster_phc)) |>                  # citizens who can be assigned a cluster
  summarise(
    saudi_pop            = sum(cell_saudi),
    pct_gen_other_cluster  = 100 * sum(cell_saudi[cluster_gen  != cluster_phc],
                                       na.rm = TRUE) / saudi_pop,
    pct_tert_other_cluster = 100 * sum(cell_saudi[cluster_tert != cluster_phc],
                                       na.rm = TRUE) / saudi_pop)

# ---- 10. Tables and figures ----------------------------------------------

write_csv(q1_national, file.path(tabs, "q1_national_phc_by_hours.csv"))
write_csv(q1_cluster,  file.path(tabs, "q1_cluster_phc.csv"))
write_csv(q2_cluster,  file.path(tabs, "q2_referral_by_cluster.csv"))
write_csv(q3,          file.path(tabs, "q3_specialized_access.csv"))
write_csv(q4_tert_access |> filter(!is.na(cluster)), file.path(tabs, "q4_tertiary_access.csv"))
write_csv(q4_alignment,       file.path(tabs, "q4_boundary_alignment.csv"))

map_theme <- theme_void(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11,
                                  margin = margin(b = 4)),
        legend.title = element_text(size = 10),
        legend.key.height = grid::unit(4, "mm"))

# Binned travel-time scale: population cells drawn over a light backdrop.
t_labels <- c("<15", "15–30", "30–45", "45–60",
              "60–120", ">120")
bin_time <- function(t) {
  b <- cut(t, breaks = c(0, 15, 30, 45, 60, 120, Inf),
           labels = t_labels, include.lowest = TRUE)
  factor(ifelse(is.na(t), "No route", as.character(b)),
         levels = c(t_labels, "No route"))
}
pal_time <- c(setNames(viridis::viridis(length(t_labels), direction = -1),
                       t_labels), "No route" = "grey55")
basemap <- list(
  geom_sf(data = gov_pop, fill = "grey96", color = NA))        # land
# Borders go on top of the travel-time cells, so governorate lines stay
# visible in densely populated areas: governorates thin, regions thicker.
borders_top <- list(
  geom_sf(data = gov_pop, fill = NA, color = "grey45", linewidth = .1),
  geom_sf(data = reg, fill = NA, color = "grey15", linewidth = .3))
time_scale <- scale_color_manual(
  values = pal_time, drop = FALSE,
  guide = guide_legend(nrow = 1, override.aes = list(size = 3)))

# Fig 1 -- Q1: PHC travel time by working-hours tier, populated cells
# Panels are time-of-day scenarios (which PHCs are open), not facility
# categories: evening = extended-hours AND 24h PHCs both open.
tier_labs <- c(sprintf("Daytime: all PHCs (n = %s)", format(nrow(phc), big.mark = ",")),
               sprintf("Evening: open after 16:00 (n = %d)", nrow(phc_ext)),
               sprintf("Night: open overnight (n = %d)", nrow(phc_24)))
acc_tier <- bind_rows(
  acc |> mutate(panel = tier_labs[1], t = t_phc),
  acc |> mutate(panel = tier_labs[2], t = t_ext),
  acc |> mutate(panel = tier_labs[3], t = t_24)) |>
  mutate(panel = factor(panel, levels = tier_labs),
         t_bin = bin_time(t))
p1 <- ggplot() + basemap +
  geom_sf(data = acc_tier, aes(color = t_bin), shape = 15, size = .12) +
  borders_top +
  time_scale + labs(color = "Minutes to nearest PHC") +
  facet_wrap(~panel, nrow = 1) + map_theme
ggsave(file.path(figs, "map_phc_access.png"), p1, width = 13.5,
       height = 6, dpi = 300, bg = "white")

# Fig 2 -- Q3: hospital travel time by service, facilities overlaid
# Specialty panels name the general hospitals they include (as in the
# cluster travel-time figure).
svc_levels <- c(sprintf("General hospitals (n = %d)", nrow(hosp_gen)),
                sprintf(paste0("Hospitals with obstetric services (n = %d)\n",
                               "incl. general hospitals with obstetrics"),
                        nrow(hosp_ob)),
                sprintf(paste0("Hospitals with pediatric services (n = %d)\n",
                               "incl. general hospitals with pediatrics"),
                        nrow(hosp_pd)),
                sprintf(paste0("Hospitals with inpatient psychiatric ",
                               "services (n = %d)\nincl. general hospitals ",
                               "with psychiatry"), nrow(hosp_psy)))
acc_svc <- bind_rows(
  acc |> mutate(panel = svc_levels[1], t = t_gen),
  acc |> mutate(panel = svc_levels[2], t = t_obgyn),
  acc |> mutate(panel = svc_levels[3], t = t_peds),
  acc |> mutate(panel = svc_levels[4], t = t_psych)) |>
  mutate(panel = factor(panel, levels = svc_levels), t_bin = bin_time(t))
fac_svc <- bind_rows(
  hosp_gen |> mutate(panel = svc_levels[1]),
  hosp_ob  |> mutate(panel = svc_levels[2]),
  hosp_pd  |> mutate(panel = svc_levels[3]),
  hosp_psy |> mutate(panel = svc_levels[4])) |>
  mutate(panel = factor(panel, levels = svc_levels))
p2 <- ggplot() + basemap +
  geom_sf(data = acc_svc, aes(color = t_bin), shape = 15, size = .12) +
  borders_top +
  geom_sf(data = fac_svc, shape = 24, size = 1.2, fill = "white",
          color = "black", stroke = .4) +
  time_scale +
  labs(color = "Minutes to nearest hospital of that type",
       caption = "Triangles: the hospitals counted in each panel.") +
  facet_wrap(~panel, nrow = 2) + map_theme +
  theme(plot.caption = element_text(size = 10, hjust = .5))
ggsave(file.path(figs, "map_hospital_access.png"), p2, width = 10,
       height = 11, dpi = 300, bg = "white")

# Fig 3 -- Q4: tertiary travel time, tertiary centers overlaid
p3 <- ggplot() + basemap +
  geom_sf(data = acc |> mutate(t_bin = bin_time(t_tert)),
          aes(color = t_bin), shape = 15, size = .15) +
  borders_top +
  geom_sf(data = tert, shape = 23, size = 2, fill = "white",
          color = "black", stroke = .5) +
  time_scale +
  labs(title = sprintf("Medical cities and specialized hospitals (n = %d sites)", n_tert_sites),
       subtitle = "The four cluster medical cities and King Fahad Specialist Hospital Dammam",
       color = "Minutes to nearest medical city or specialized hospital",
       caption = "Diamonds: medical city and specialized hospital sites.") +
  map_theme +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = .5),
        plot.subtitle = element_text(size = 10.5, hjust = .5,
                                     margin = margin(t = 3, b = 6)),
        plot.caption = element_text(size = 10, hjust = .5),
        # the one-row legend is wider than this single map: title above the keys
        legend.title.position = "top",
        legend.title = element_text(size = 10, hjust = .5))
ggsave(file.path(figs, "map_tertiary_access.png"), p3, width = 8,
       height = 8.8, dpi = 300, bg = "white")


# ---- 10b. Cluster scorecard (the cluster as unit of comparison) ----------
# One row per cluster spanning all four sub-questions, plus cross-cluster
# dependence. Descriptive comparison only; the
# spatial support remains the 1-km grid, and cluster membership is applied
# as an attribute after routing, never as a constraint.

# Count facility SITES, not directory rows. Entries the census codes as
# specialized (3ry) list hospital-campus centers separately: Al Baha's ten centers of
# excellence share the King Fahad Hospital - Al Bahah campus (identical
# coordinates); Hail's cardiac and diabetes-and-endocrine centers sit on the
# King Salman Specialist Hospital campus (120-300 m from the hospital pin);
# Qassim's specialist medical center sits on the King Fahd Specialist
# Hospital campus (~240 m). Entries of the same cluster within 500 m of one
# another are merged into one site; the smallest distance between genuinely
# distinct specialized sites is 2.8 km (Eastern: King Fahad Specialist Hospital
# Dammam vs Saud Al Babtain Cardiac Center), so the threshold is unambiguous. The merge applies to
# these specialized entries ONLY: adjacent PHCs
# and secondary hospitals (e.g. a maternity hospital beside a general
# hospital) are distinct facilities and keep exact-coordinate deduplication.
# One named exception on the secondary tier: the Al Makhwah endocrinology
# and diabetes center (ALB-010) is an outpatient center of excellence hosted
# INSIDE Al Makhwah General Hospital (ALB-014, 4 m away) and is counted with
# its host as one secondary site. Coordinates are kept as geocoded, so
# routing and travel times are unaffected.
#
# The referral level ("tertiary" in this script) is the designated medical
# cities and specialized hospitals set in section 1 (Supplementary Methods);
# the four advanced services are recorded hospital by hospital in
# referral_services.csv. The other specialized entries are secondary-level
# hospital sites.
# campus_id() is defined in section 1.

fac_sites <- fac |>
  filter(facility_id != "ALB-010") |>  # hosted inside Al Makhwah GH; one site
  distinct(cluster_en, tier, campus_entry, lat, lon) |>
  group_by(cluster_en, tier, campus_entry) |>
  mutate(site = if (first(campus_entry)) campus_id(lat, lon)
                else seq_len(n())) |>
  ungroup() |>
  distinct(cluster_en, tier, campus_entry, site)

fac_n <- fac_sites |>
  count(cluster_en, tier) |>
  pivot_wider(names_from = tier, values_from = n, values_fill = 0,
              names_prefix = "n_") |>
  left_join(fac |> filter(tier == "primary", hours == "h24") |>
              count(cluster_en, name = "n_phc_24h"), by = "cluster_en") |>
  mutate(n_phc_24h = coalesce(n_phc_24h, 0L))
stopifnot(n_tert_sites == sum(fac_n$n_tertiary))   # section 1 count
stopifnot(sum(fac_n$n_primary)   == 2169,
          sum(fac_n$n_secondary) == 293,
          sum(fac_n$n_tertiary)  == 5,
          all(fac_n$n_tertiary[fac_n$cluster_en %in%
                c("Riyadh First", "Riyadh Second", "Makkah", "Madinah", "Eastern")] == 1))

scorecard <- acc |> st_drop_geometry() |>
  filter(!is.na(cluster_phc)) |>
  group_by(cluster = cluster_phc) |>
  summarise(
    saudi_pop         = sum(cell_saudi),
    med_phc_min       = wq(t_phc, cell_saudi, .5),
    pct30_phc         = cov_pct(t_phc, cell_saudi, 30),
    pct30_phc_24h     = cov_pct(t_24,  cell_saudi, 30),
    med_gen_min       = wq(t_gen, cell_saudi, .5),
    pct60_gen         = cov_pct(t_gen, cell_saudi, 60),
    pct60_obgyn_women = cov_pct(t_obgyn, cell_women, 60),
    pct60_peds_child  = cov_pct(t_peds, cell_child, 60),
    pct60_psych_adult = cov_pct(t_psych, cell_adult, 60),
    med_tert_min      = wq(t_tert, cell_saudi, .5),
    pct120_tert       = cov_pct(t_tert, cell_saudi, 120),
    pct_gen_other     = 100 * sum(cell_saudi[cluster_gen != cluster_phc],
                                  na.rm = TRUE) / saudi_pop,
    pct_tert_other    = 100 * sum(cell_saudi[cluster_tert != cluster_phc],
                                  na.rm = TRUE) / saudi_pop,
    .groups = "drop") |>
  left_join(q2_cluster |>
              select(cluster_en, median_ref_min, ref_over_60 = n_over_60),
            by = c("cluster" = "cluster_en")) |>
  left_join(fac_n, by = c("cluster" = "cluster_en")) |>
  mutate(saudi_per_phc = saudi_pop / n_primary)   # vs typical care-team panel sizes

write_csv(scorecard, file.path(tabs, "cluster_scorecard.csv"))

# Cluster comparison figure: five coverage indicators as a table-style
# heatmap in one sequential hue (darker = more of the population covered),
# with each value printed in its cell. Tertiary access is left out: six
# clusters have no tertiary center, so it is compared as referral outside
# the cluster instead (section 7c).
sc_ind <- c(
  "PHC within 30 min\n(daytime)"                         = "pct30_phc",
  "PHC open at night\nwithin 30 min"                     = "pct30_phc_24h",
  "Hospital with obstetric\nservices within 60 min\n(women 15–49)" =
    "pct60_obgyn_women",
  "Hospital with pediatric\nservices within 60 min\n(children under 15)" =
    "pct60_peds_child",
  "Hospital with inpatient\npsychiatric services\nwithin 60 min (adults 18+)" =
    "pct60_psych_adult")
sc_long <- scorecard |>
  select(cluster, all_of(sc_ind)) |>
  pivot_longer(-cluster, names_to = "indicator", values_to = "pct") |>
  mutate(cluster = factor(cluster, levels = scorecard |>
                            arrange(pct30_phc_24h) |> pull(cluster)),
         indicator = factor(indicator, levels = names(sc_ind)),
         ink = ifelse(pct >= 55, "white", "grey10"))   # text contrast >= ~4:1
sc_ramp <- c("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b")
p6 <- ggplot(sc_long, aes(indicator, cluster, fill = pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f", pct), color = ink), size = 3.4) +
  scale_fill_gradientn(colours = sc_ramp, limits = c(0, 100),
                       breaks = seq(0, 100, 25),
                       labels = function(x) paste0(x, "%"),
                       name = "Share of the population covered") +
  scale_color_identity() +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       caption = paste(strwrap(paste(
         "Clusters ordered by coverage of a PHC open at night. PHC columns:",
         "all Saudi citizens; hospital columns: the Saudi population that",
         "needs the service."), width = 110), collapse = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 9, lineheight = .95),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom", legend.title.position = "top",
        legend.key.width = unit(2.4, "cm"), legend.key.height = unit(.35, "cm"),
        plot.caption = element_text(size = 8.5, hjust = 0, color = "grey30"),
        plot.caption.position = "plot")
ggsave(file.path(figs, "fig_cluster_scorecard.png"), p6, width = 8.5,
       height = 9.5, dpi = 300, bg = "white")

# ---- 10c. Reference maps: administrative geography and cluster coverage --

pal20 <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
           "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f",
           "#cab2d6", "#8dd3c7", "#bebada", "#fa8072", "#80b1d3",
           "#fdb462", "#b3de69", "#fccde5", "#bc80bd", "#ffed6f")

# Map A: the 13 administrative regions and their 150 governorates
reg_xy <- st_coordinates(st_point_on_surface(st_geometry(reg)))
pA <- ggplot() +
  geom_sf(data = gov, aes(fill = Region_EN), color = "white",
          linewidth = .15, alpha = .85, show.legend = FALSE) +
  geom_sf(data = reg, fill = NA, color = "grey25", linewidth = .45) +
  ggrepel::geom_text_repel(
    data = reg |> st_drop_geometry() |>
      mutate(x = reg_xy[, 1], y = reg_xy[, 2]),
    aes(x, y, label = Region_EN), size = 3.1, fontface = "bold",
    color = "grey10", bg.color = "white", bg.r = .12, seed = 1) +
  scale_fill_manual(values = pal20[1:13]) +
  map_theme
ggsave(file.path(figs, "map_regions_governorates.png"), pA, width = 8,
       height = 8.5, dpi = 300, bg = "white")

# Map B: the 20 health clusters and the governorates they cover.
# Cluster membership is facility-level; each governorate is filled with
# the cluster holding the majority of its facilities, and governorates
# whose facilities span more than one cluster are starred.
gov_clu <- fac |>
  count(governorate_en, cluster_en) |>
  group_by(governorate_en) |>
  mutate(n_clu = n_distinct(cluster_en)) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(governorate_en, cluster_en, n_clu)
gov_cl_sf <- gov |> left_join(gov_clu, by = c("Gov_EN" = "governorate_en"))
clu_diss <- gov_cl_sf |> filter(!is.na(cluster_en)) |>
  group_by(cluster_en) |> summarise(.groups = "drop")
# Drop interior rings from the dissolved cluster outlines: the union of the
# ArcGIS-edited governorate borders leaves sliver interior rings (785 rings,
# all < 0.007 km2; verified 2026-09-12 that no genuine holes or fully
# enclosed cross-cluster enclaves exist), which would otherwise be stroked
# dark inside the clusters. Exact — exterior boundaries are untouched.
drop_holes <- function(geom) {
  crs <- sf::st_crs(geom)
  sf::st_sfc(lapply(geom, function(g) {
    polys <- sf::st_cast(sf::st_sfc(g, crs = crs), "POLYGON")
    sf::st_multipolygon(lapply(polys, function(p) list(p[[1]])))
  }), crs = crs)
}
sf::st_geometry(clu_diss) <- drop_holes(sf::st_geometry(clu_diss))
stopifnot(sum(vapply(st_geometry(clu_diss), function(g) {
  polys <- st_cast(st_sfc(g, crs = st_crs(clu_diss)), "POLYGON")
  sum(vapply(polys, function(p) length(p) - 1L, 1L))
}, 1L)) == 0)  # zero interior rings remain
# Label anchors: the most interior point of each sizeable part of a cluster
# (both parts of Makkah get a label; small islands do not). The four small
# west-coast clusters are labelled in the Red Sea with leader lines.
clu_parts <- suppressWarnings(st_cast(clu_diss, "POLYGON"))
clu_parts$area <- as.numeric(st_area(clu_parts))
clu_parts <- clu_parts |>
  group_by(cluster_en) |>
  filter(area >= .25 * max(area), area >= 2e9) |>
  ungroup()
# nQuadSegs = 0 returns one 2-point line per part: circle center, then a
# boundary point
ic <- st_coordinates(st_inscribed_circle(st_geometry(clu_parts),
                                         dTolerance = 1000, nQuadSegs = 0))
lab_xy <- ic[!duplicated(ic[, "L1"]), c("X", "Y"), drop = FALSE]
lab_df <- tibble(cluster_en = clu_parts$cluster_en,
                 x = lab_xy[, 1], y = lab_xy[, 2])
west_small <- c("Jeddah First", "Jeddah Second", "Makkah", "Al Baha")
north_out  <- "Northern Borders"      # narrow strip: label outside, with a line
lab_in  <- lab_df |> filter(!cluster_en %in% c(west_small, north_out)) |>
  mutate(y = y + ifelse(cluster_en == "Hafar Al Batin", 45e3, 0))  # clear *
lab_out <- lab_df |> filter(cluster_en %in% west_small)
lab_north <- lab_df |> filter(cluster_en %in% north_out)
x_sea   <- min(lab_out$x) - 230e3              # label column in the Red Sea
# neighbors get clearly different fills (Makkah vs Jeddah Second, Jazan vs
# Najran were near-identical in pal20)
cl_names <- sort(unique(clu_diss$cluster_en))
clu_cols <- setNames(pal20[seq_along(cl_names)], cl_names)
clu_cols[c("Makkah", "Jazan")] <- c("#e7298a", "#ffff99")
split_gov <- gov_cl_sf |> filter(n_clu > 1)
pB <- ggplot() +
  geom_sf(data = gov_cl_sf, aes(fill = cluster_en), color = "white",
          linewidth = .12, alpha = .9, show.legend = FALSE) +
  geom_sf(data = clu_diss, fill = NA, color = "grey15", linewidth = .45) +
  geom_sf(data = st_point_on_surface(st_geometry(split_gov)),
          shape = 8, size = 1.8, color = "black", stroke = .6) +
  geom_label(data = lab_in, aes(x, y, label = cluster_en), size = 3,
             fontface = "bold", color = "grey10",
             fill = scales::alpha("white", .8), linewidth = 0,
             label.padding = unit(.12, "lines")) +
  ggrepel::geom_label_repel(
    data = lab_out, aes(x, y, label = cluster_en), size = 3,
    fontface = "bold", color = "grey10", fill = "white",
    nudge_x = x_sea - lab_out$x, xlim = c(NA, x_sea), direction = "y",
    hjust = 1, min.segment.length = 0, segment.color = "grey15",
    segment.size = .35, box.padding = .35, seed = 1) +
  ggrepel::geom_label_repel(
    data = lab_north, aes(x, y, label = cluster_en), size = 3,
    fontface = "bold", color = "grey10", fill = "white",
    nudge_x = 90e3, nudge_y = 170e3, min.segment.length = 0,
    segment.color = "grey15", segment.size = .35, seed = 1) +
  scale_fill_manual(values = clu_cols, na.value = "grey92") +
  labs(caption = paste0(
    "* Governorate hosts facilities of more than one cluster; fill shows ",
    "the majority cluster.\nEastern and Makkah clusters each cover two ",
    "separate areas, labelled in both. Lines point to clusters too small ",
    "or narrow to label inside.")) +
  map_theme +
  theme(plot.caption = element_text(size = 8, hjust = 0))
ggsave(file.path(figs, "map_clusters_governorates.png"), pB, width = 8,
       height = 8.5, dpi = 300, bg = "white")

# ---- 10d. Per-cluster travel-time distributions --------------------------
# Population-weighted travel-time quantiles per cluster for every
# destination, with demand-matched weights (women 15-49 for maternity,
# children under 15 for pediatric, adults 18+ for psychiatric,
# Saudi citizens otherwise).

dest_set <- tribble(
  ~var,      ~w_var,       ~destination,
  "t_phc",   "cell_saudi", "PHC, daytime (all)",
  "t_ext",   "cell_saudi", "PHC, evening (open after 16:00)",
  "t_24",    "cell_saudi", "PHC, night (open overnight)",
  "t_gen",   "cell_saudi", "General hospital",
  "t_obgyn", "cell_women", "Hospital with obstetric services",
  "t_peds",  "cell_child", "Hospital with pediatric services",
  "t_psych", "cell_adult", "Hospital with inpatient psychiatric services",
  "t_tert",  "cell_saudi", "Medical city or specialized hospital")
# Figure panel titles only (the CSV keeps the short names): facility counts,
# that the specialty sets include general hospitals with the department, and
# whose travel time is shown when it is not all citizens. Tertiary centers
# are counted as sites (n_tert_sites, section 1).
dest_set$panel <- c(
  sprintf("PHC, daytime (all)\nn = %s", format(nrow(phc), big.mark = ",")),
  sprintf("PHC, evening (open after 16:00)\nn = %d", nrow(phc_ext)),
  sprintf("PHC, night (open overnight)\nn = %d", nrow(phc_24)),
  sprintf("General hospital\nn = %d", nrow(hosp_gen)),
  sprintf(paste0("Hospital with obstetric services\nn = %d, incl. general ",
                 "hospitals with obstetrics\nfor women 15–49"), nrow(hosp_ob)),
  sprintf(paste0("Hospital with pediatric services\nn = %d, incl. general ",
                 "hospitals with pediatrics\nfor children under 15"),
          nrow(hosp_pd)),
  sprintf(paste0("Hospital with inpatient psychiatric services\nn = %d, ",
                 "incl. general hospitals with psychiatry\nfor adults 18+"),
          nrow(hosp_psy)),
  sprintf("Medical city or\nspecialized hospital\nn = %d sites", n_tert_sites))

acc_clu <- acc |> st_drop_geometry() |> filter(!is.na(cluster_phc))
clu_tt <- bind_rows(lapply(seq_len(nrow(dest_set)), function(i) {
  v <- dest_set$var[i]; wv <- dest_set$w_var[i]
  acc_clu |>
    group_by(cluster = cluster_phc) |>
    summarise(destination  = dest_set$destination[i],
              p25_min      = wq(.data[[v]], .data[[wv]], .25),
              median_min   = wq(.data[[v]], .data[[wv]], .50),
              p75_min      = wq(.data[[v]], .data[[wv]], .75),
              p90_min      = wq(.data[[v]], .data[[wv]], .90),
              pct_no_route = 100 * sum(.data[[wv]][is.na(.data[[v]])],
                                       na.rm = TRUE) /
                sum(.data[[wv]], na.rm = TRUE),
              .groups = "drop")
})) |>
  mutate(across(where(is.numeric), ~ round(.x, 1)))
write_csv(clu_tt |> arrange(destination, cluster),
          file.path(tabs, "cluster_travel_times.csv"))

# Interval figure: median dot with p25-p90 bar, clusters ordered by their
# night-PHC median (the indicator with the widest between-cluster spread).
# Tertiary centers are left out of the figure (six clusters have none; the
# CSV keeps them).
ord <- clu_tt |> filter(destination == "PHC, night (open overnight)") |>
  arrange(median_min) |> pull(cluster)
p7 <- clu_tt |>
  filter(destination != "Medical city or specialized hospital") |>
  mutate(cluster = factor(cluster, levels = ord),
         panel = factor(dest_set$panel[match(destination,
                                             dest_set$destination)],
                        levels = dest_set$panel[dest_set$var != "t_tert"])) |>
  ggplot(aes(median_min, cluster)) +
  geom_linerange(aes(xmin = p25_min, xmax = p90_min),
                 color = "grey60", linewidth = .8) +
  geom_point(size = 1.9, color = "#1b7837") +
  facet_wrap(~panel, nrow = 2, scales = "free_x") +
  labs(x = "Minutes (dot = population-weighted median; bar = p25–p90)",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(figs, "fig_cluster_travel_times.png"), p7, width = 12.5,
       height = 8.5, dpi = 300, bg = "white")

# Small-multiple maps: one panel per cluster, cells colored by binned
# travel time, cluster boundary from the facility-majority dissolve.
# Every panel is a framed square of the same size, so each cluster name
# sits directly above its own map; neighboring clusters are drawn in light
# grey for context (this also shows clusters in two separate parts).
clu_map <- function(v, legend_lab, file) {
  cls <- sort(clu_diss$cluster_en)
  plots <- lapply(cls, function(cl) {
    cells <- acc |> filter(cluster_phc == cl) |>
      mutate(t_bin = bin_time(.data[[v]]))
    bound <- clu_diss |> filter(cluster_en == cl)
    b1 <- st_bbox(bound); b2 <- st_bbox(cells)   # cells may lie outside
    xr <- range(b1[c(1, 3)], b2[c(1, 3)])
    yr <- range(b1[c(2, 4)], b2[c(2, 4)])
    half <- max(diff(xr), diff(yr)) / 2 * 1.06   # square extent, 6% margin
    ggplot() +
      geom_sf(data = clu_diss, fill = "grey98", color = "grey78",
              linewidth = .15) +
      geom_sf(data = bound, fill = "grey90", color = "grey25",
              linewidth = .35) +
      # legend from the first panel only; otherwise patchwork shows
      # near-duplicate legends side by side when any panel differs
      geom_sf(data = cells, aes(color = t_bin), shape = 15, size = .25,
              show.legend = identical(cl, cls[1])) +
      scale_color_manual(values = pal_time, drop = FALSE,
                         name = legend_lab,
                         guide = guide_legend(nrow = 1,
                           override.aes = list(size = 3))) +
      coord_sf(xlim = mean(xr) + c(-half, half),
               ylim = mean(yr) + c(-half, half), expand = FALSE) +
      labs(title = cl) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(size = 9, face = "bold",
                                      hjust = .5, margin = margin(b = 3)),
            plot.background = element_rect(fill = "white",
                                           color = "grey60",
                                           linewidth = .4),
            plot.margin = margin(4, 4, 4, 4))
  })
  # legend in its own row below the grid, so no panel frame runs into it
  pw <- (patchwork::wrap_plots(plots, ncol = 5) / patchwork::guide_area()) +
    patchwork::plot_layout(guides = "collect", heights = c(1, .05)) +
    patchwork::plot_annotation(
      caption = paste(strwrap(paste(
        "Dark outline: governorates where the cluster runs most of its",
        "facilities; other clusters in light grey. Colored cells: 1-km cells",
        "whose nearest PHC belongs to the cluster, so in governorates shared",
        "by two clusters (for example Jeddah and Riyadh) some cells fall",
        "outside the outline."), width = 150), collapse = "\n"),
      theme = theme(plot.caption = element_text(size = 9, hjust = 0)))
  ggsave(file.path(figs, file), pw, width = 12, height = 11.5,
         dpi = 300, bg = "white")
}
clu_map("t_24",  "Minutes to nearest PHC open overnight",
        "map_cluster_phc24.png")
clu_map("t_gen", "Minutes to nearest general hospital",
        "map_cluster_gen.png")

# ---- 10e. Physicians and facility validation (MoH Yearbook 2024) ---------
# MoH Statistical Yearbook 2024 (data/raw): PHC physicians (table 2-45,
# national total used to cost the conversions), physicians by nationality
# and grade (2-20, with the GP total checked against 2-29), Saudization of
# PHC physicians (2-43), and facility counts by cluster (2-4), which
# validate the facility file. Parsed values are checked against the
# yearbook's own totals (stopifnot checks below).

yb <- file.path(raw, "Statistical-Yearbook-2024 1.xlsx")
num <- function(x) suppressWarnings(as.numeric(x))

# Branch column order in table 2-45 (5 stacked blocks x 4 branches,
# S/NS/Total triples in cols 3-14; category total row = header + offset).
br_2445 <- c("Riyadh", "The Holy Capital", "Jeddah", "Taif",
             "Madinah", "Qassim", "Eastern", "Al Ahsa",
             "Hafr Al Batin", "Asir", "Bishah", "Tabuk",
             "Hail", "Northern Borders", "Jazan", "Najran",
             "Al Baha", "Al Jouf", "Qurayyat", "Qunfudah")
tot_cols <- c(5, 8, 11, 14)

d45 <- suppressMessages(readxl::read_excel(yb, "2-45", col_names = FALSE))
grab45 <- function(off) setNames(unlist(lapply(c(4, 30, 56, 82, 108),
  function(h) num(unlist(d45[h + off, tot_cols])))), br_2445)
dent_br   <- grab45(9)                      # PHC dentists (Yearbook lists them among PHC "physicians"; cf. 2-43)
phcdoc_br <- grab45(24) - dent_br           # PHC physicians excluding dentists
stopifnot(sum(dent_br) == 2408)
# Yearbook grand total (row 158) counts dentists among PHC physicians
stopifnot(sum(phcdoc_br) + sum(dent_br) == num(d45[158, 3]) + num(d45[158, 7]))

# 2-29: physicians by specialty (3 stacked blocks x 7 branches, Total col
# every 4th from col 5; specialty matched by Arabic label within block).
d29 <- suppressMessages(readxl::read_excel(yb, "2-29", col_names = FALSE))
br_29 <- list(c("Riyadh", "The Holy Capital", "Jeddah", "Taif",
                "Madinah", "Qassim", "Eastern"),
              c("Al Ahsa", "Hafr Al Batin", "Asir", "Bishah",
                "Tabuk", "Hail", "Northern Borders"),
              c("Jazan", "Najran", "Al Baha", "Al Jouf",
                "Qurayyat", "Qunfudah", "KSA_TOTAL"))
grab29 <- function(label) {
  out <- unlist(lapply(seq_along(c(4, 47, 90)), function(b) {
    h <- c(4, 47, 90)[b]
    rows <- (h + 4):min(h + 43, nrow(d29))
    hit <- rows[which(trimws(unlist(d29[rows, 1])) == label)]
    stopifnot(length(hit) == 1)
    setNames(num(unlist(d29[hit, seq(5, 29, 4)])), br_29[[b]])
  }))
  stopifnot(sum(out[names(out) != "KSA_TOTAL"]) == out[["KSA_TOTAL"]])
  out[names(out) != "KSA_TOTAL"]
}
gp_br    <- grab29("عام")   # general practitioners; Yearbook records all at resident grade

# National specialty-by-nationality table 2-20 (Saudi/non-Saudi x grade), cited in
# Methods and Discussion. Cols: 2-6 resident (SM,SF,NSM,NSF,Total), 7-11 registrar,
# 12-16 consultant, 17 grand total. Row 7 = GP, row 32 = family medicine.
d20 <- suppressMessages(readxl::read_excel(yb, "2-20", col_names = FALSE))
n20 <- function(r, cols) sum(sapply(cols, function(cc) num(d20[r, cc])))
gp_saudi_20  <- n20(7, c(2, 3))                       # 11,453
gp_total_20  <- num(d20[7, 6])                        # 20,606
fm_reg_s_20  <- n20(32, c(7, 8));  fm_reg_t_20 <- num(d20[32, 11])
fm_con_s_20  <- n20(32, c(12, 13)); fm_con_t_20 <- num(d20[32, 16])
stopifnot(gp_total_20 == sum(gp_br),                  # 2-20 matches 2-29 GP total
          gp_saudi_20 == 11453, fm_reg_t_20 == 3133, fm_con_t_20 == 1702)
cat(sprintf(paste0("Table 2-20 nationality splits: GP %.1f%% Saudi; ",
                   "FM registrars %.1f%% Saudi; FM consultants %.1f%% Saudi; ",
                   "qualified FM (reg+cons) %.1f%% Saudi\n"),
    100 * gp_saudi_20 / gp_total_20,
    100 * fm_reg_s_20 / fm_reg_t_20,
    100 * fm_con_s_20 / fm_con_t_20,
    100 * (fm_reg_s_20 + fm_con_s_20) / (fm_reg_t_20 + fm_con_t_20)))

# Five-year Saudization trend in PHCs, national (table 2-43; dentistry excluded).
# Rows: 9/11 dentistry Saudi/Total, 18/20 family medicine Saudi/Total,
# 24/26 all-specialty Saudi/Total; cols 3:7 = years 2020..2024.
d43 <- suppressMessages(readxl::read_excel(yb, "2-43", col_names = FALSE))
n43 <- function(r) num(unlist(d43[r, 3:7]))
phc_saudi_43 <- n43(24) - n43(9)
phc_total_43 <- n43(26) - n43(11)
fm_saudi_43  <- n43(18); fm_total_43 <- n43(20)
stopifnot(phc_total_43[5] == 9693, phc_saudi_43[5] == 5255,
          phc_total_43[5] == sum(phcdoc_br))          # 2-43 matches 2-45
cat(sprintf(paste0("Table 2-43 Saudization 2020-2024: PHC physicians (excl. dentists) ",
                   "%.1f%% -> %.1f%% Saudi; PHC family medicine %.1f%% -> %.1f%% Saudi\n"),
    100 * phc_saudi_43[1] / phc_total_43[1], 100 * phc_saudi_43[5] / phc_total_43[5],
    100 * fm_saudi_43[1] / fm_total_43[1],  100 * fm_saudi_43[5] / fm_total_43[5]))

# Facility validation at true cluster level: yearbook table 2-4 vs ours
d4 <- suppressMessages(readxl::read_excel(yb, "2-4", col_names = FALSE))
map4 <- c("تجمع الرياض الصحي الأول"   = "Riyadh First",
          "تجمع الرياض الصحي الثاني"  = "Riyadh Second",
          "تجمع الرياض الصحي الثالث"  = "Riyadh Third",
          "إجمالي تجمع مكة المكرمة الصحي" = "Makkah",
          "تجمع جدة الصحي الأول"      = "Jeddah First",
          "تجمع جدة الصحي الثاني"     = "Jeddah Second",
          "تجمع الطائف الصحي"         = "Al Taif",
          "تجمع المدينة المنورة الصحي" = "Madinah",
          "تجمع القصيم الصحي"         = "Qassim",
          "تجمع الشرقية الصحي"        = "Eastern",
          "تجمع الأحساء الصحي"        = "Al Ahsa",
          "تجمع حفر الباطن الصحي"     = "Hafar Al Batin",
          "إجمالي تجمع عسير الصحي"    = "Aseer",
          "تجمع تبوك الصحي"           = "Tabuk",
          "تجمع حائل الصحي"           = "Hail",
          "تجمع الحدود الشمالية الصحي" = "Northern Borders",
          "تجمع جازان الصحي"          = "Jazan",
          "تجمع نجران الصحي"          = "Najran",
          "تجمع الباحة الصحي"         = "Al Baha",
          "إجمالي تجمع الجوف الصحي"   = "Al Jouf")
idx4 <- match(names(map4), trimws(unlist(d4[, 1])))
stopifnot(!anyNA(idx4))
phc_pop <- tibble(cluster = unname(map4),
                  phc_yearbook  = num(unlist(d4[idx4, 2])),
                  hosp_yearbook = num(unlist(d4[idx4, 3]))) |>
  left_join(fac |> filter(tier == "primary") |>
              count(cluster_en, name = "n_phc_ours"),
            by = c("cluster" = "cluster_en")) |>
  left_join(fac_sites |> filter(tier != "primary") |>   # campus-merged sites
              count(cluster_en, name = "n_hosp_ours"),
            by = c("cluster" = "cluster_en")) |>
  left_join(acc |> st_drop_geometry() |> filter(!is.na(cluster_phc)) |>
              group_by(cluster = cluster_phc) |>
              summarise(saudi_pop = sum(cell_saudi), .groups = "drop"),
            by = "cluster") |>
  mutate(phc_per_10k_saudi = round(1e4 * n_phc_ours / saudi_pop, 2),
         saudis_per_phc    = round(saudi_pop / n_phc_ours))
write_csv(phc_pop, file.path(tabs, "cluster_phc_population.csv"))


# ---- 11. Sensitivity: terrain-gradient speed penalty ---------------------
# Run with:  SLOPE_SENS=1 Rscript R/analysis.R
# The base model reflects topography through the road network itself:
# mountain roads carry their true switchback geometry (longer paths) and
# lower highway classes (lower profile speeds), and roadless terrain is
# unreachable. What the base model omits is gradient-induced slowdown on
# steep segments. This sensitivity extracts edge gradients from the SRTM
# 30-arcsec DEM (geodata) and applies a linear speed penalty between
# PAR$slope_thresh and PAR$slope_full grade (max PAR$slope_maxpen
# reduction, i.e. speed halved at >=12% grade), then recomputes travel
# times to the three key destination sets.

if (nzchar(Sys.getenv("SLOPE_SENS"))) {
  message("Section 11: terrain-gradient sensitivity")
  dem <- geodata::elevation_30s(country = "SAU",
                                path = dirname(path.expand(PAR$pbf_path)))
  el1 <- terra::extract(dem, cbind(graph$from_lon, graph$from_lat))[, 1]
  el2 <- terra::extract(dem, cbind(graph$to_lon,   graph$to_lat))[, 1]
  grade <- abs(el2 - el1) / pmax(graph$d, 1)
  pen <- 1 - PAR$slope_maxpen *
    pmin(pmax(grade - PAR$slope_thresh, 0) /
           (PAR$slope_full - PAR$slope_thresh), 1)
  pen[!is.finite(pen)] <- 1          # edges outside DEM: no penalty
  time0 <- graph$time; time_w0 <- graph$time_weighted   # restored below
  message(sprintf("  edge-km with grade > %.0f%%: %.2f%% (penalised)",
                  100 * PAR$slope_thresh,
                  100 * sum(graph$d[pen < 1], na.rm = TRUE) / sum(graph$d)))
  graph$time          <- graph$time / pen
  graph$time_weighted <- graph$time_weighted / pen

  tt_s <- cache("travel_times_slope", {
    list(phc  = tt_nearest(grid_pts, phc),
         gen  = tt_nearest(grid_pts, hosp_gen),
         tert = tt_nearest(grid_pts, tert))
  })
  if (!identical(tt_s$tert$dest_ids, sort(tert$facility_id))) {   # see section 4
    tt_s$tert <- tt_nearest(grid_pts, tert, label = "referral")
    tt_s$tert$dest_ids <- sort(tert$facility_id)
    saveRDS(tt_s, file.path(prc, "travel_times_slope.rds"))
  }

  sens <- bind_rows(lapply(list(
    phc  = list(b = acc$t_phc,  s = tt_s$phc$minutes,  thr = 30),
    gen  = list(b = acc$t_gen,  s = tt_s$gen$minutes,  thr = 60),
    tert = list(b = acc$t_tert, s = tt_s$tert$minutes, thr = 120)),
    function(z) tibble(
      threshold_min    = z$thr,
      median_base      = wq(z$b, acc$cell_saudi, .5),
      median_slope     = wq(z$s, acc$cell_saudi, .5),
      p90_base         = wq(z$b, acc$cell_saudi, .9),
      p90_slope        = wq(z$s, acc$cell_saudi, .9),
      pct_within_base  = cov_pct(z$b, acc$cell_saudi, z$thr),
      pct_within_slope = cov_pct(z$s, acc$cell_saudi, z$thr))),
    .id = "destination")
  write_csv(sens, file.path(tabs, "sensitivity_slope_penalty.csv"))

  # PHC x cell times on the penalised graph, for the decision outputs (D9)
  pcs_path <- file.path(dirname(graph_path), "phc_cell_times_slope.rds")
  if (!file.exists(pcs_path)) saveRDS(phc_cell_matrix(graph), pcs_path)
  # Cluster-level shift for PHC access (mountain clusters move most), from the
  # same PHC x cell matrices as the decision analyses (D3, D9)
  sens_cl <- acc |> st_drop_geometry() |>
    mutate(t_base = colmin(pc), t_slope = colmin(readRDS(pcs_path))) |>
    group_by(cluster = cluster_phc) |>
    summarise(pct30_base  = cov_pct(t_base,  cell_saudi, 30),
              pct30_slope = cov_pct(t_slope, cell_saudi, 30),
              .groups = "drop") |>
    mutate(delta_pp = pct30_slope - pct30_base)
  write_csv(sens_cl, file.path(tabs, "sensitivity_slope_by_cluster.csv"))
  graph$time <- time0; graph$time_weighted <- time_w0   # later sections route on the base graph
  message("  slope sensitivity written")
}

# ---- 11. Supplement sensitivity tables (thresholds, residents, assignment)
# S1: threshold coverage grid, primary (Saudi / demand-matched) denominators.
# S2: the same grid under the total-resident denominator (the cluster model's
#     mandate covers all residents).
# S11: cluster-assignment sensitivity — the study assigns each grid cell to
#     the cluster of its nearest PHC (registration proxy); the alternative
#     assigns it to the cluster operating the majority of facilities in its
#     governorate. GASTAT publishes no population by health cluster, so this
#     internal check replaces external validation of cluster totals.

dest_spec <- tribble(
  ~col,      ~destination,                 ~wcol_primary, ~denom_primary,
  "t_phc",   "PHC, any",                   "cell_saudi",  "Saudi citizens",
  "t_ext",   "PHC, open after 16:00",      "cell_saudi",  "Saudi citizens",
  "t_24",    "PHC, open overnight",        "cell_saudi",  "Saudi citizens",
  "t_gen",   "General hospital",           "cell_saudi",  "Saudi citizens",
  "t_obgyn", "Hospital with obstetric services", "cell_women",  "Saudi women 15-49",
  "t_peds",  "Hospital with pediatric services", "cell_child",  "Saudi children <15",
  "t_psych", "Hospital with inpatient psychiatric services", "cell_adult",  "Saudi adults 18+",
  "t_tert",  "Medical city or specialized hospital",            "cell_saudi",  "Saudi citizens")

thresh_grid <- function(wcol_field, denom_field) {
  acc_df <- acc |> st_drop_geometry()
  bind_rows(lapply(seq_len(nrow(dest_spec)), function(i) {
    x <- acc_df[[dest_spec$col[i]]]
    w <- acc_df[[if (wcol_field == "total") "cell_total" else dest_spec$wcol_primary[i]]]
    tibble(destination = dest_spec$destination[i],
           denominator = if (denom_field == "total") "Total residents"
                         else dest_spec$denom_primary[i],
           median = round(wq(x, w, .5), 1), p90 = round(wq(x, w, .9), 1),
           pct20  = round(cov_pct(x, w, 20), 1),
           pct30  = round(cov_pct(x, w, 30), 1),
           pct45  = round(cov_pct(x, w, 45), 1),
           pct60  = round(cov_pct(x, w, 60), 1),
           pct120 = round(cov_pct(x, w, 120), 1))
  }))
}
write_csv(thresh_grid("primary", "primary"),
          file.path(tabs, "table_S1_threshold_coverage_saudi.csv"))
write_csv(thresh_grid("total", "total"),
          file.path(tabs, "table_S2_threshold_coverage_total_residents.csv"))

gov_major <- fac |>
  count(governorate_en, cluster_en) |>
  group_by(governorate_en) |>
  slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |>
  select(governorate_en, cl_gov = cluster_en)
assign_cmp <- acc |> st_drop_geometry() |>
  left_join(gov_major, by = c("Gov_EN" = "governorate_en")) |>
  filter(!is.na(cluster_phc), !is.na(cl_gov))     # routable cells, both methods
s11 <- full_join(
  assign_cmp |> group_by(cluster = cluster_phc) |>
    summarise(saudi_nearest_phc = round(sum(cell_saudi))),
  assign_cmp |> group_by(cluster = cl_gov) |>
    summarise(saudi_gov_majority = round(sum(cell_saudi))),
  by = "cluster") |>
  mutate(diff = saudi_gov_majority - saudi_nearest_phc,
         diff_pct = round(100 * diff / saudi_nearest_phc, 1)) |>
  arrange(cluster)
write_csv(s11, file.path(tabs, "table_S11_cluster_assignment_sensitivity.csv"))
message(sprintf(
  "  S11: %.1f%% of citizens assigned to the same cluster by both methods",
  100 * sum(assign_cmp$cell_saudi[assign_cmp$cluster_phc == assign_cmp$cl_gov]) /
        sum(assign_cmp$cell_saudi)))


# S12: region-level census validation. GASTAT publishes population by
# governorate and by region, but not by health cluster, so cluster totals
# have no direct official comparator. Aggregating the nearest-PHC cluster
# populations to the 13 administrative regions gives the coarsest geography
# at which the assignment can be checked against census totals: any
# divergence measures citizens pulled across a region border by the
# nearest-PHC rule. Governorate totals match GASTAT by construction
# (section 3) and are not a check.
cl_region <- c(
  "Riyadh First" = "Riyadh", "Riyadh Second" = "Riyadh",
  "Riyadh Third" = "Riyadh",
  "Makkah" = "Makkah Al Mukarramah", "Jeddah First" = "Makkah Al Mukarramah",
  "Jeddah Second" = "Makkah Al Mukarramah", "Al Taif" = "Makkah Al Mukarramah",
  "Eastern" = "Eastern", "Al Ahsa" = "Eastern", "Hafar Al Batin" = "Eastern",
  "Madinah" = "Al Madinah Al Munawarah", "Qassim" = "Al Qassim",
  "Hail" = "Hail", "Tabuk" = "Tabuk", "Al Jouf" = "Al Jouf",
  "Northern Borders" = "Northern Borders", "Aseer" = "Aseer",
  "Al Baha" = "Al Baha", "Jazan" = "Jazan", "Najran" = "Najran")
reg_cells <- acc |> st_drop_geometry() |>
  left_join(gov |> st_drop_geometry() |> distinct(Gov_ID, Region_EN),
            by = "Gov_ID")
s12 <- full_join(
  reg_cells |> group_by(region = Region_EN) |>
    summarise(census_saudi   = round(sum(cell_saudi)),
              routable_saudi = round(sum(cell_saudi[!is.na(cluster_phc)]))),
  reg_cells |> filter(!is.na(cluster_phc)) |>
    group_by(region = cl_region[cluster_phc]) |>
    summarise(assigned_saudi = round(sum(cell_saudi)), .groups = "drop"),
  by = "region") |>
  mutate(diff = assigned_saudi - routable_saudi,
         diff_pct = round(100 * diff / routable_saudi, 2)) |>
  arrange(region)
write_csv(s12, file.path(tabs, "table_S12_region_census_validation.csv"))
message(sprintf(
  "  S12: %.2f%% of routable citizens assigned across a region border",
  100 * sum(reg_cells$cell_saudi[!is.na(reg_cells$cluster_phc) &
        cl_region[reg_cells$cluster_phc] != reg_cells$Region_EN]) /
        sum(reg_cells$cell_saudi[!is.na(reg_cells$cluster_phc)])))


# ---- Descriptive summary tables (earlier draft) --------------------------
# Facility and access summaries from an earlier draft, kept because they pin
# the facility counts; the current Tables 1-2 are built in "Manuscript
# exhibits for the decision analyses". Site counts
# fold co-located directory entries into one site per campus: identical
# coordinates collapse the ten Al Baha centers of excellence, and four
# named satellite entries are folded into their host hospitals (two centers
# at King Salman Specialist Hospital in Hail, one at King Fahd Specialist
# Hospital in Buraidah, one outpatient center at Al Makhwah General
# Hospital). The stopifnot() calls pin the documented totals (2,169 PHCs,
# 282 secondary sites, 16 tertiary sites); a future data edit that changes
# them fails loudly here.
satellites <- c("Diabetes and Endocrine Center",                    # Hail
                "Heart Center in King Salman Specialist Hospital",  # Hail
                "Specialist Medical Centre King Fahad Hospital",    # Qassim
                "Endocrinology and Diabetes Center - Al Makhwah")   # Al Baha
stopifnot(sum(fac$name_en %in% satellites) == 4)
fac_sites <- fac |>
  filter(!name_en %in% satellites) |>
  distinct(tier, campus_entry, lat, lon, .keep_all = TRUE)
stopifnot(sum(fac_sites$tier == "primary")   == 2169,
          sum(fac_sites$tier == "secondary") == 293,
          sum(fac_sites$tier == "tertiary")  == 5)

site_counts <- fac_sites |>
  group_by(cluster = cluster_en) |>
  summarise(n_phc  = sum(tier == "primary"),
            n_ext  = sum(tier == "primary" & hours == "extended"),
            n_24h  = sum(tier == "primary" & hours == "h24"),
            n_hosp = sum(tier == "secondary"),
            n_tert = sum(tier == "tertiary"), .groups = "drop")
t1 <- acc |> st_drop_geometry() |>
  filter(!is.na(cluster_phc)) |>
  group_by(cluster = cluster_phc) |>
  summarise(saudi_pop = round(sum(cell_saudi)), .groups = "drop") |>
  left_join(site_counts, by = "cluster") |>
  mutate(per_phc = round(saudi_pop / n_phc)) |>
  arrange(cluster)
t1 <- bind_rows(t1, t1 |> summarise(cluster = "National",
  across(c(saudi_pop, n_phc, n_ext, n_24h, n_hosp, n_tert), sum)) |>
  mutate(per_phc = round(saudi_pop / n_phc)))
write_csv(t1, file.path(tabs, "table1_clusters.csv"))

t2_dest <- c("PHC, any"              = "PHC, any (daytime)",
             "PHC, open after 16:00" = "PHC, open after 16:00 (evening)",
             "PHC, open overnight"   = "PHC, open overnight (night)",
             "General hospital"      = "General hospital",
             "Hospital with obstetric services" = "Hospital with obstetric services",
             "Hospital with pediatric services" = "Hospital with pediatric services",
             "Hospital with inpatient psychiatric services" = "Hospital with inpatient psychiatric services",
             "Medical city or specialized hospital"       = "Medical city or specialized hospital")
t2_den  <- c("Saudi citizens"     = "Saudi citizens",
             "Saudi women 15-49"  = "Saudi women 15\u201349",
             "Saudi children <15" = "Saudi children <15",
             "Saudi adults 18+"   = "Saudi adults \u226518")
t2 <- thresh_grid("primary", "primary") |>
  transmute(destination = unname(t2_dest[destination]),
            denom = unname(t2_den[denominator]),
            median, p90, pct30, pct60, pct120)
stopifnot(!anyNA(t2$destination), !anyNA(t2$denom))
write_csv(t2, file.path(tabs, "table2_national_access.csv"))
message("  Manuscript tables 1-2 written.")

# ---- 12. Decision analyses ---------------------------------------------------
# Uses the objects built above (acc, phc, hosp_*, tert*, hc, pc, ref_mat, wq,
# cov_pct, ...).
#
# D1  Coverage through the week, from each PHC's listed opening hours
# D2  Night-time gap types: open PHC / MoH general hospital ED only / neither
# D2b Night care for children and for women in labor
# D2c Night load: citizens per overnight PHC
# D3  Overnight PHC placement: re-siting the night network (nationally and
#     within each cluster, with and without protection); conversions needed
#     for national and cluster targets, and for citizens with no night option
# D4  Extending hours versus opening 24/7, per added opening hour
# D5  Physicians needed for the conversions
# D6  Cluster borders: nearest facility versus nearest same-cluster facility
# D7  Hospital departments: additions that bring the most people within 60 min
# D8  In-cluster referral hub where a cluster lacks the service (not reported)
# D9  Exact integer-programming solutions; decision outputs at 20 and 45 min
#     and with the terrain penalty

message("---- Decision analyses ----")
suppressPackageStartupMessages(library(Matrix))
thr_phc  <- PAR$thresholds_phc[2]      # 30 min
thr_hosp <- PAR$thresholds_hosp[2]     # 60 min

# PHC x cell travel times (pc) and colmin() come from analysis.R section 4.
cells_of <- lapply(seq_len(nrow(pc)), function(j) which(pc[j, ] <= thr_phc))
A <- sparseMatrix(i = rep(seq_along(cells_of), lengths(cells_of)),
                  j = unlist(cells_of), x = 1, dims = dim(pc))
wS <- ifelse(is.na(acc$cell_saudi), 0, acc$cell_saudi)
wT <- ifelse(is.na(acc$cell_total), 0, acc$cell_total)
covered <- function(idx) { v <- logical(ncol(pc)); for (j in idx) v[cells_of[[j]]] <- TRUE; v }
pct <- function(v, w = wS) 100 * sum(w[v]) / sum(w)
stopifnot(abs(pct(covered(which(phc$hours == "h24"))) -
              cov_pct(acc$t_24, acc$cell_saudi, thr_phc)) < 0.2)
er <- !is.na(acc$t_gen) & acc$t_gen <= thr_phc          # hospital ER, open 24/7

# ---- D1. Coverage through the week ---------------------------------------
# Opening hours from the cluster directories (data/raw/phc_opening_hours.csv):
# centers without a listing keep regular hours (Sunday-Thursday 08:00-16:00,
# closed Friday-Saturday; MoH working-hours schedule, 2016). The week splits
# into ten blocks whose boundaries are every listed opening and closing time,
# so coverage is exact per block.
oh <- read_csv(file.path(raw, "phc_opening_hours.csv"), show_col_types = FALSE)
oh <- oh[match(phc$facility_id, oh$facility_id), ]
stopifnot(!anyNA(oh$weekday_open))
blocks <- tibble(day  = rep(c("Weekday", "Weekend"), each = 5),
                 from = rep(c(0, 8, 16, 21, 23), 2),
                 to   = rep(c(8, 16, 21, 23, 24), 2)) |>
  mutate(n_days = ifelse(day == "Weekday", 5, 2),
         hours_week = (to - from) * n_days,
         label = sprintf("%s %02d:00-%02d:00", day, from, to))
open_blocks <- function(wo, wc, eo, ec) sapply(seq_len(nrow(blocks)), function(b) {
  o <- if (blocks$day[b] == "Weekday") wo else eo
  c <- if (blocks$day[b] == "Weekday") wc else ec
  !is.na(o) & o <= blocks$from[b] & c >= blocks$to[b]
})
O <- open_blocks(oh$weekday_open, oh$weekday_close, oh$weekend_open, oh$weekend_close)
cov_b <- lapply(seq_len(nrow(blocks)), function(b) covered(which(O[, b])))
d1 <- blocks |>
  mutate(n_phc_open    = colSums(O),
         pct_phc       = sapply(cov_b, pct),
         pct_phc_or_er = sapply(cov_b, function(v) pct(v | er)),
         pct_er_only   = sapply(cov_b, function(v) pct(!v & er)),
         pct_neither   = sapply(cov_b, function(v) pct(!v & !er)),
         pct_phc_all_residents = sapply(cov_b, function(v) pct(v, wT)))
write_csv(d1 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d1_coverage_by_time.csv"))
# Sensitivity: the 22 PHCs whose census lists them as extended or 24-hour but
# which have no directory listing are set to regular hours (lower bound).
unl  <- grepl("^default for census", oh$hours_source)
O_lo <- open_blocks(ifelse(unl, 8, oh$weekday_open), ifelse(unl, 16, oh$weekday_close),
                    ifelse(unl, NA_real_, oh$weekend_open), ifelse(unl, NA_real_, oh$weekend_close))
sens_unlisted <- blocks |>
  transmute(label, hours_week, n_phc_open = colSums(O), pct_phc = d1$pct_phc,
            n_phc_open_lower = colSums(O_lo),
            pct_phc_lower = sapply(seq_len(nrow(blocks)), function(b) pct(covered(which(O_lo[, b]))))) |>
  mutate(change = pct_phc_lower - pct_phc)
write_csv(sens_unlisted |> mutate(across(where(is.double), ~ round(.x, 2))),
          file.path(tabs, "sensitivity_unlisted_hours.csv"))
message(sprintf("  D1 sensitivity: %d unlisted PHCs set to regular hours; largest coverage change %.2f points",
                sum(unl), min(sens_unlisted$change)))
day_blk <- which(blocks$day == "Weekday" & blocks$from == 8)
wk_avg  <- sum(d1$hours_week * d1$pct_phc) / 168
ah_avg  <- sum(d1$hours_week[-day_blk] * d1$pct_phc[-day_blk]) / sum(d1$hours_week[-day_blk])
message(sprintf("  D1: 30-min PHC coverage averaged over the week %.1f%%; outside weekday daytime %.1f%%",
                wk_avg, ah_avg))

hourly <- expand_grid(day = c("Weekday", "Weekend"), hour = 0:24) |>
  mutate(b = mapply(function(d, h) which(blocks$day == d & blocks$from <= min(h, 23) &
                                           blocks$to > min(h, 23)), day, hour))
hourly <- bind_rows(
  hourly |> mutate(measure = "Open PHC", pct = d1$pct_phc[b]),
  hourly |> mutate(measure = "Open PHC or MoH general hospital ED", pct = d1$pct_phc_or_er[b])) |>
  mutate(day = factor(recode(day, Weekday = "Sunday to Thursday", Weekend = "Friday and Saturday"),
                      levels = c("Sunday to Thursday", "Friday and Saturday")))
pD1 <- ggplot(hourly, aes(hour, pct, colour = measure)) +
  geom_step(linewidth = 0.9) +
  facet_wrap(~day) +
  scale_x_continuous(breaks = seq(0, 24, 4),
                     labels = c("00", "04", "08", "12", "16", "20", "24")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = "Hour of day", y = "Saudi citizens within 30 min (%)", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(figs, "d1_coverage_by_hour.png"), pD1, width = 9, height = 4.5,
       dpi = 300, bg = "white")

# ---- D2. Night-time gap types (Sunday-Thursday, 00:00-08:00) -------------
b_night <- which(blocks$day == "Weekday" & blocks$from == 0)
b_wkend_night <- which(blocks$day == "Weekend" & blocks$from == 0)
gap_levels <- c("MoH PHC open at night within 30 min", "MoH general hospital ED only", "Neither MoH option within 30 min")
gap_of <- function(v) factor(ifelse(v, gap_levels[1], ifelse(er, gap_levels[2], gap_levels[3])),
                             levels = gap_levels)
acc$night_gap <- gap_of(cov_b[[b_night]])
night_open_d2 <- which(O[, b_night])
acc$wkend_night_gap <- gap_of(cov_b[[b_wkend_night]])
share <- function(g, lev, w) 100 * sum(w[g == lev], na.rm = TRUE) / sum(w, na.rm = TRUE)
d2 <- acc |> st_drop_geometry() |> filter(!is.na(cluster_phc)) |>
  group_by(cluster = cluster_phc) |>
  summarise(saudi = sum(cell_saudi, na.rm = TRUE),
            n_phc_open_night = sum(phc$cluster_en[night_open_d2] == cluster_phc[1]),
            pct_24h_phc   = share(night_gap, gap_levels[1], cell_saudi),
            pct_24h_phc_all = share(night_gap, gap_levels[1], cell_total),  # all residents
            pct_er_only   = share(night_gap, gap_levels[2], cell_saudi),
            pct_neither   = share(night_gap, gap_levels[3], cell_saudi),
            weekend_pct_24h_phc = share(wkend_night_gap, gap_levels[1], cell_saudi),
            .groups = "drop")
d2 <- bind_rows(d2, acc |> st_drop_geometry() |>
  summarise(cluster = "National", saudi = sum(cell_saudi, na.rm = TRUE),
            n_phc_open_night = length(night_open_d2),
            pct_24h_phc = share(night_gap, gap_levels[1], cell_saudi),
            pct_24h_phc_all = share(night_gap, gap_levels[1], cell_total),
            pct_er_only = share(night_gap, gap_levels[2], cell_saudi),
            pct_neither = share(night_gap, gap_levels[3], cell_saudi),
            weekend_pct_24h_phc = share(wkend_night_gap, gap_levels[1], cell_saudi)))
write_csv(d2 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d2_night_gap_by_cluster.csv"))
message(sprintf("  D2: at night %.1f%% 24-hour PHC, %.1f%% ER only, %.1f%% neither (%s citizens)",
                tail(d2$pct_24h_phc, 1), tail(d2$pct_er_only, 1), tail(d2$pct_neither, 1),
                format(round(sum(wS[acc$night_gap == gap_levels[3]])), big.mark = ",")))
pD2 <- ggplot() + basemap +
  geom_sf(data = acc |> filter(cell_saudi > 0), aes(color = night_gap), shape = 15, size = .12) +
  borders_top +
  scale_color_manual(values = c("#1B9E77", "#E6AB02", "#D95F02"), drop = FALSE,
                     guide = guide_legend(override.aes = list(size = 3))) +
  labs(color = "At night (00:00-08:00, Sunday-Thursday)") + map_theme
ggsave(file.path(figs, "d2_night_gap_map.png"), pD2, width = 8, height = 7.5,
       dpi = 300, bg = "white")

# ---- D2b. Night care for children and for women in labor -----------------
# Hospital pediatric and obstetric departments work around the clock, so at
# night a sick child can be seen at a 24-hour PHC or at a hospital with a
# pediatric department, and a woman in labor needs a hospital with obstetric
# services (secondary or tertiary).
t_pd_any <- nearest_of(c(hosp_pd$facility_id, tert_pd$facility_id))
t_ob_any <- nearest_of(c(hosp_ob$facility_id, tert_ob$facility_id))
wC <- ifelse(is.na(acc$cell_child), 0, acc$cell_child)
wW <- ifelse(is.na(acc$cell_women), 0, acc$cell_women)
phc_n <- cov_b[[b_night]]
pd30  <- !is.na(t_pd_any) & t_pd_any <= thr_phc
d2b <- tibble(
  group = c(rep("Saudi children under 15", 4), rep("Saudi women 15-49", 2)),
  measure = c("24-hour PHC within 30 min at night",
              "Pediatric hospital only within 30 min at night",
              "Neither within 30 min at night",
              "Pediatric hospital within 60 min (any hour)",
              "Obstetric hospital within 30 min (any hour)",
              "Obstetric hospital within 60 min (any hour)"),
  pct = c(pct(phc_n, wC), pct(!phc_n & pd30, wC), pct(!phc_n & !pd30, wC),
          pct(!is.na(t_pd_any) & t_pd_any <= thr_hosp, wC),
          pct(!is.na(t_ob_any) & t_ob_any <= thr_phc, wW),
          pct(!is.na(t_ob_any) & t_ob_any <= thr_hosp, wW)))
write_csv(d2b |> mutate(pct = round(pct, 1)), file.path(tabs, "d2b_night_children_women.csv"))
message(sprintf("  D2b: children at night: %.1f%% 24-hour PHC, %.1f%% pediatric hospital only, %.1f%% neither",
                d2b$pct[1], d2b$pct[2], d2b$pct[3]))

# ---- D2c. Night load on overnight PHCs -------------------------------------
# Coverage counts whether an open PHC is within 30 min, not how many people
# share it. Each citizen within 30 min of an open PHC at night is assigned to
# the nearest one, in any cluster; a center's load is the citizens assigned.
pn <- pc[night_open_d2, , drop = FALSE]
near_n <- apply(pn, 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
in30 <- cov_b[[b_night]]
stopifnot(!anyNA(near_n[in30]))
nl_cluster <- phc$cluster_en[night_open_d2]
night_load <- phc[night_open_d2, ] |> st_drop_geometry() |>
  select(facility_id, name_en, cluster = cluster_en)
night_load$load <- vapply(seq_along(night_open_d2), function(k)
  sum(wS[in30 & near_n %in% k]), 0)
night_load$load_other_clusters <- vapply(seq_along(night_open_d2), function(k)
  sum(wS[in30 & near_n %in% k & !(acc$cluster_phc %in% nl_cluster[k])]), 0)
stopifnot(abs(sum(night_load$load) - sum(wS[in30])) < 1)
d2c <- night_load |> group_by(cluster) |>
  summarise(phcs_open_night = n(), citizens_per_night_phc = mean(load),
            max_load = max(load), pct_load_other_clusters = 100 * sum(load_other_clusters) / max(sum(load), 1),
            .groups = "drop")
write_csv(night_load |> arrange(desc(load)) |> mutate(across(where(is.double), round)),
          file.path(tabs, "d2c_night_load_by_phc.csv"))
nl_origin <- bind_rows(lapply(seq_along(night_open_d2), function(k) {
  ci <- which(in30 & near_n %in% k)
  if (!length(ci)) return(NULL)
  tibble(facility_id = phc$facility_id[night_open_d2[k]], cluster = nl_cluster[k],
         residents_of = acc$cluster_phc[ci], w = wS[ci]) |>
    group_by(facility_id, cluster, residents_of) |> summarise(citizens = sum(w), .groups = "drop")
}))
write_csv(nl_origin |> arrange(facility_id, desc(citizens)) |> mutate(citizens = round(citizens)),
          file.path(tabs, "d2c_night_load_origin.csv"))
write_csv(d2c |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d2c_night_load_by_cluster.csv"))
message(sprintf("  D2c: citizens per overnight PHC: national %s; by cluster %s (%s) to %s (%s)",
                format(round(mean(night_load$load)), big.mark = ","),
                format(round(min(d2c$citizens_per_night_phc)), big.mark = ","),
                d2c$cluster[which.min(d2c$citizens_per_night_phc)],
                format(round(max(d2c$citizens_per_night_phc)), big.mark = ","),
                d2c$cluster[which.max(d2c$citizens_per_night_phc)]))

# ---- D3. Overnight PHC placement (weekday night) --------------------------
# Greedy maximal covering: each step adds the PHC that brings the most
# uncovered citizens within 30 min. `base` marks cells already served (for
# the ER-aware objective, cells within 30 min of a hospital ER).
greedy <- function(start = integer(0), n_max = Inf, target = Inf,
                   pool = seq_len(nrow(A)), w = wS, base = logical(ncol(A))) {
  cvd <- covered(start) | base
  sel <- start
  blocked <- rep(TRUE, nrow(A)); blocked[pool] <- FALSE; blocked[start] <- TRUE
  steps <- list(tibble(n_added = 0L, facility_id = NA_character_, pct = pct(cvd, w)))
  while (length(sel) - length(start) < n_max && tail(steps, 1)[[1]]$pct < target) {
    g <- as.numeric(A %*% (w * !cvd)); g[blocked] <- -1
    j <- which.max(g); if (g[j] <= 0) break
    sel <- c(sel, j); blocked[j] <- TRUE; cvd[cells_of[[j]]] <- TRUE
    steps[[length(steps) + 1]] <- tibble(n_added = length(sel) - length(start),
                                         facility_id = phc$facility_id[j], pct = pct(cvd, w))
  }
  list(sel = sel, path = bind_rows(steps))
}
night_open <- which(O[, b_night])
n_needed <- function(path, tgt) path$n_added[which(path$pct >= tgt)[1]]
g_resite <- greedy(n_max = length(night_open))
g_add    <- greedy(start = night_open, target = 90)
g_er     <- greedy(start = night_open, base = er, target = 97)
t_resite <- colmin(pc[g_resite$sel, , drop = FALSE])
t_now    <- colmin(pc[night_open, , drop = FALSE])
d3 <- tibble(
  scenario = c("Current night network", "Same number, best locations",
               "Current + conversions to 80%", "Current + conversions to 85%",
               "Current + conversions to 90%", "All PHCs open at night (ceiling)",
               "Current, PHC or hospital ED", "PHC or ED + conversions to 96%",
               "PHC or ED + conversions to 97%"),
  n_24h_phc = c(length(night_open), length(g_resite$sel),
                length(night_open) + c(n_needed(g_add$path, 80), n_needed(g_add$path, 85),
                                       n_needed(g_add$path, 90)), nrow(phc),
                length(night_open), length(night_open) + n_needed(g_er$path, 96),
                length(night_open) + n_needed(g_er$path, 97)),
  pct_within_30 = c(pct(cov_b[[b_night]]), tail(g_resite$path$pct, 1),
                    sapply(c(80, 85, 90), function(t) g_add$path$pct[n_needed(g_add$path, t) + 1]),
                    pct(covered(seq_len(nrow(pc)))), pct(cov_b[[b_night]] | er),
                    sapply(c(96, 97), function(t) g_er$path$pct[n_needed(g_er$path, t) + 1])),
  median_min = c(wq(t_now, acc$cell_saudi, .5), wq(t_resite, acc$cell_saudi, .5),
                 rep(NA, 7)),
  p90_min = c(wq(t_now, acc$cell_saudi, .9), wq(t_resite, acc$cell_saudi, .9), rep(NA, 7)))
write_csv(d3 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d3_night_siting.csv"))
added <- setdiff(g_add$sel, night_open)
write_csv(g_add$path |> left_join(phc |> st_drop_geometry() |>
                                    select(facility_id, name_en, governorate_en, cluster_en, hours),
                                  by = "facility_id") |>
            mutate(pct = round(pct, 2)),
          file.path(tabs, "d3_night_conversions_order.csv"))
# The from-scratch placement itself (Supplementary Table S16, part B): the 90
# PHCs in the order chosen, their current hours, citizens each one adds, and
# cumulative coverage.
write_csv(g_resite$path[-1, ] |>
            mutate(citizens_added = round(diff(g_resite$path$pct) / 100 * sum(wS))) |>
            left_join(phc |> st_drop_geometry() |>
                        select(facility_id, name_en, governorate_en, cluster_en, hours),
                      by = "facility_id") |>
            mutate(pct = round(pct, 2)),
          file.path(tabs, "d3_night_resite_from_scratch.csv"))
message(sprintf("  D3: night 30-min coverage %.1f%% now; same %d centres best placed %.1f%%; 80/85/90%% need %d/%d/%d conversions",
                pct(cov_b[[b_night]]), length(night_open), tail(g_resite$path$pct, 1),
                n_needed(g_add$path, 80), n_needed(g_add$path, 85), n_needed(g_add$path, 90)))

# Cluster targets: conversions each cluster needs, from its own PHCs, for 80%
# of its citizens within 30 min of a night PHC (any cluster's open PHCs count).
cl_runs <- lapply(sort(unique(na.omit(acc$cluster_phc))), function(k) {
  wk <- ifelse(acc$cluster_phc %in% k, wS, 0)
  pool <- which(phc$cluster_en == k)
  now  <- pct(cov_b[[b_night]], wk)
  ceil <- pct(covered(union(night_open, pool)), wk)
  g <- if (now >= 80) NULL else greedy(start = night_open, pool = pool, w = wk, target = 80)
  list(row = tibble(cluster = k, pct_now = now, conversions_for_80 = if (is.null(g)) 0L else
                      if (tail(g$path$pct, 1) >= 80) tail(g$path$n_added, 1) else NA_integer_,
                    pct_ceiling_own_phcs = ceil),
       added = if (is.null(g)) integer(0) else setdiff(g$sel, night_open))
})
d3c <- bind_rows(lapply(cl_runs, `[[`, "row"))
cl_added <- unlist(lapply(cl_runs, `[[`, "added"))
pct_all80 <- pct(covered(c(night_open, cl_added)))
kept_resite <- sum(g_resite$sel %in% night_open)
write_csv(d3c |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d3_night_cluster_targets.csv"))
message(sprintf("  D3: every cluster to 80%%: %d conversions in total (national %.1f%%); %d clusters cannot reach 80%% with their own PHCs; best placement keeps %d of %d current centres",
                sum(d3c$conversions_for_80, na.rm = TRUE), pct_all80, sum(is.na(d3c$conversions_for_80)),
                kept_resite, length(night_open)))

# Moves within clusters, at no cost: each cluster keeps its number of overnight
# PHCs, and a move shifts one night shift from a center to another PHC of the
# same cluster. Moves are made one at a time, each the one that raises national
# coverage most, until no move gains a citizen (an exchange heuristic for the
# maximal covering model). cnt = open PHCs within 30 min of each cell.
cnt <- as.numeric(colSums(A[night_open, , drop = FALSE]))
cur <- night_open
mv <- list(tibble(move = 0L, cluster = NA_character_, from_id = NA_character_,
                  to_id = NA_character_, pct = pct(cnt > 0)))
repeat {
  g_now <- as.numeric(A %*% (wS * (cnt == 0)))
  best <- 0; bi <- NA_integer_; bj <- NA_integer_
  for (i in cur) {
    cand <- setdiff(which(phc$cluster_en == phc$cluster_en[i]), cur)
    if (!length(cand)) next
    u <- cells_of[[i]][cnt[cells_of[[i]]] == 1]            # cells only i covers
    h <- if (length(u)) as.numeric(A[cand, u, drop = FALSE] %*% wS[u]) else 0
    net <- g_now[cand] + h - sum(wS[u])
    k <- which.max(net)
    if (net[k] > best) { best <- net[k]; bi <- i; bj <- cand[k] }
  }
  if (is.na(bi) || best < 1) break
  cnt[cells_of[[bi]]] <- cnt[cells_of[[bi]]] - 1
  cnt[cells_of[[bj]]] <- cnt[cells_of[[bj]]] + 1
  cur <- c(setdiff(cur, bi), bj)
  mv[[length(mv) + 1]] <- tibble(move = length(mv), cluster = phc$cluster_en[bi],
                                 from_id = phc$facility_id[bi], to_id = phc$facility_id[bj],
                                 pct = pct(cnt > 0))
}
d3m <- bind_rows(mv)
stopifnot(all(diff(d3m$pct) > 0), abs(tail(d3m$pct, 1) - pct(covered(cur))) < 1e-9,
          identical(table(phc$cluster_en[cur]), table(phc$cluster_en[night_open])))
nm_of  <- function(id) phc$name_en[match(id, phc$facility_id)]
gov_of <- function(id) phc$governorate_en[match(id, phc$facility_id)]
write_csv(d3m |> mutate(from_name = nm_of(from_id), from_governorate = gov_of(from_id),
                        to_name = nm_of(to_id), to_governorate = gov_of(to_id), pct = round(pct, 2)),
          file.path(tabs, "d3_night_moves_within_clusters.csv"))
cov_cl <- function(v) sapply(sort(unique(na.omit(acc$cluster_phc))), function(k)
  pct(v, ifelse(acc$cluster_phc %in% k, wS, 0)))
d3m_cl <- tibble(cluster = sort(unique(na.omit(acc$cluster_phc))),
                 phcs_open_night = as.integer(table(factor(phc$cluster_en[night_open], levels = cluster))),
                 moves = as.integer(table(factor(d3m$cluster[-1], levels = cluster))),
                 pct_now = cov_cl(cov_b[[b_night]]), pct_after_moves = cov_cl(cnt > 0))
write_csv(d3m_cl |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d3_night_moves_by_cluster.csv"))
mv_at <- function(n) d3m$pct[min(n, nrow(d3m) - 1) + 1]
# Citizens gaining and losing a night PHC within 30 min after the first six
# moves and after all moves, and the load on the busiest center afterwards.
set_after <- function(n) {
  s <- night_open
  for (m in seq_len(n)) s <- c(setdiff(s, match(d3m$from_id[m + 1], phc$facility_id)),
                               match(d3m$to_id[m + 1], phc$facility_id))
  s
}
d3m_eff <- bind_rows(lapply(c(6L, nrow(d3m) - 1L), function(n) {
  s <- set_after(n); aft <- covered(s); now <- cov_b[[b_night]]
  nn <- apply(pc[s, , drop = FALSE], 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
  ld <- vapply(seq_along(s), function(k) sum(wS[aft & nn %in% k]), 0)
  tibble(moves = n, pct = pct(aft), citizens_gaining = sum(wS[!now & aft]),
         citizens_losing = sum(wS[now & !aft]), losing_with_er_30 = sum(wS[now & !aft & er]),
         max_load = max(ld), max_load_phc = phc$name_en[s[which.max(ld)]],
         max_load_now = max(night_load$load))
}))
stopifnot(abs(d3m_eff$pct - c(mv_at(6), tail(d3m$pct, 1))) < 1e-9)
write_csv(d3m_eff |> mutate(pct = round(pct, 1), across(c(citizens_gaining, citizens_losing,
                                                     losing_with_er_30, max_load, max_load_now), round)),
          file.path(tabs, "d3_night_moves_effects.csv"))
message(sprintf("  D3: moves within clusters (no added hours): %.1f%% after %d moves (6 moves %.1f%%, 10 moves %.1f%%); %d of %d current centres kept",
                tail(d3m$pct, 1), nrow(d3m) - 1, mv_at(6), mv_at(10),
                sum(cur %in% night_open), length(night_open)))
message(sprintf("      first 6 moves: +%s / -%s citizens (%.1f%% of losers keep an ER within 30 min); busiest centre after all moves %s (now %s)",
                format(round(d3m_eff$citizens_gaining[1]), big.mark = ","),
                format(round(d3m_eff$citizens_losing[1]), big.mark = ","),
                100 * d3m_eff$losing_with_er_30[1] / d3m_eff$citizens_losing[1],
                format(round(d3m_eff$max_load[2]), big.mark = ","),
                format(round(d3m_eff$max_load_now[2]), big.mark = ",")))

# Protected moves within clusters: as above, but no move may leave a citizen
# with neither a night PHC nor a hospital ER within 30 min. Protected cells are
# those now within 30 min of a night PHC and not within 30 min of an ER.
prot  <- cov_b[[b_night]] & !er
cnt_p <- as.numeric(colSums(A[night_open, , drop = FALSE]))
cur_p <- night_open
mvp <- list(tibble(move = 0L, cluster = NA_character_, from_id = NA_character_,
                   to_id = NA_character_, pct = pct(cnt_p > 0)))
repeat {
  g_now <- as.numeric(A %*% (wS * (cnt_p == 0)))
  best <- 0; bi <- NA_integer_; bj <- NA_integer_
  for (i in cur_p) {
    cand <- setdiff(which(phc$cluster_en == phc$cluster_en[i]), cur_p)
    if (!length(cand)) next
    u  <- cells_of[[i]][cnt_p[cells_of[[i]]] == 1]          # cells only i covers
    uP <- u[prot[u]]                                        # ... that are protected
    h  <- if (length(u)) as.numeric(A[cand, u, drop = FALSE] %*% wS[u]) else 0
    net <- g_now[cand] + h - sum(wS[u])
    if (length(uP))                                         # j must cover them all
      net[as.numeric(A[cand, uP, drop = FALSE] %*% rep(1, length(uP))) < length(uP)] <- -Inf
    k <- which.max(net)
    if (net[k] > best) { best <- net[k]; bi <- i; bj <- cand[k] }
  }
  if (is.na(bi) || best < 1) break
  cnt_p[cells_of[[bi]]] <- cnt_p[cells_of[[bi]]] - 1
  cnt_p[cells_of[[bj]]] <- cnt_p[cells_of[[bj]]] + 1
  cur_p <- c(setdiff(cur_p, bi), bj)
  mvp[[length(mvp) + 1]] <- tibble(move = length(mvp), cluster = phc$cluster_en[bi],
                                   from_id = phc$facility_id[bi], to_id = phc$facility_id[bj],
                                   pct = pct(cnt_p > 0))
}
d3mp <- bind_rows(mvp)
stopifnot(all(cnt_p[prot] > 0), abs(tail(d3mp$pct, 1) - pct(covered(cur_p))) < 1e-9,
          identical(table(phc$cluster_en[cur_p]), table(phc$cluster_en[night_open])))
write_csv(d3mp |> mutate(from_name = nm_of(from_id), from_governorate = gov_of(from_id),
                         to_name = nm_of(to_id), to_governorate = gov_of(to_id), pct = round(pct, 2)),
          file.path(tabs, "d3_night_moves_protected.csv"))

# Gains and losses of each relocation scenario: citizens who gain or lose a
# night PHC within 30 min, losers split by whether a hospital ER is within 30
# min, and the busiest center (citizens within 30 min for whom it is nearest).
reloc_eff <- function(s, label) {
  aft <- covered(s); now <- cov_b[[b_night]]
  nn <- apply(pc[s, , drop = FALSE], 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
  ld <- vapply(seq_along(s), function(k) sum(wS[aft & nn %in% k]), 0)
  tibble(scenario = label, shifts_in_new_phc = sum(!s %in% night_open), pct = pct(aft),
         citizens_gaining = sum(wS[!now & aft]), citizens_losing = sum(wS[now & !aft]),
         losing_keep_er = sum(wS[now & !aft & er]), losing_neither = sum(wS[now & !aft & !er]),
         max_load = max(ld), max_load_phc = phc$name_en[s[which.max(ld)]])
}
d3_reloc <- bind_rows(
  reloc_eff(g_resite$sel, "Same 90 PHCs, best placement found nationally"),
  reloc_eff(set_after(6), "Within clusters, first 6 moves"),
  reloc_eff(cur, "Within clusters, all moves"),
  reloc_eff(cur_p, "Within clusters, all moves, no one left without a night option"))
stopifnot(abs(d3_reloc$pct - c(tail(g_resite$path$pct, 1), mv_at(6), tail(d3m$pct, 1),
                               tail(d3mp$pct, 1))) < 1e-9, d3_reloc$losing_neither[4] == 0)
write_csv(d3_reloc |> mutate(pct = round(pct, 1), across(c(citizens_gaining:max_load), round)),
          file.path(tabs, "d3_night_relocation_effects.csv"))
message(sprintf("  D3: best national placement: +%s / -%s citizens (%s with neither); protected moves within clusters: %.1f%% after %d moves (+%s / -%s, all keep an ER)",
                format(round(d3_reloc$citizens_gaining[1]), big.mark = ","),
                format(round(d3_reloc$citizens_losing[1]), big.mark = ","),
                format(round(d3_reloc$losing_neither[1]), big.mark = ","),
                tail(d3mp$pct, 1), nrow(d3mp) - 1,
                format(round(d3_reloc$citizens_gaining[4]), big.mark = ","),
                format(round(d3_reloc$citizens_losing[4]), big.mark = ",")))

curve <- bind_rows(
  g_resite$path |> transmute(n = n_added, pct, path = "Best locations chosen from scratch"),
  g_add$path |> transmute(n = length(night_open) + n_added, pct,
                          path = "Current centers plus the best conversions"))
pD3 <- ggplot(curve, aes(n, pct, colour = path)) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.9) +
  annotate("point", x = length(night_open), y = pct(cov_b[[b_night]]), size = 2.5) +
  annotate("text", x = length(night_open) + 3, y = pct(cov_b[[b_night]]) - 3, hjust = 0,
           size = 3.2, label = sprintf("Current: %d centers, %.1f%%", length(night_open),
                                        pct(cov_b[[b_night]]))) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "PHCs open at night", y = "Saudi citizens within 30 min at night (%)", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")

# ---- D4. Extending hours versus opening 24/7 ------------------------------
# Four actions, each adding opening hours to one PHC. Gains are citizen-hours
# per week within 30 min of an open PHC, over all ten blocks, so evening and
# night hours count equally. Each strategy adds its best PHCs one at a time.
wk_hours <- as.numeric(O %*% blocks$hours_week)
is_regular <- oh$weekday_open == 8 & oh$weekday_close == 16 & is.na(oh$weekend_open)
is_ext     <- !is_regular & !(oh$weekday_open == 0 & oh$weekday_close == 24)
is_wd_ext  <- is_ext & is.na(oh$weekend_open)
tgt <- list(
  midnight  = open_blocks(8, 24, NA, NA),             # Sun-Thu 08:00-24:00
  all_night = open_blocks(0, 24, 0, 24),              # every day, all day
  weekend   = open_blocks(8, 23, 8, 23))              # extended hours on Fri-Sat too
strategies <- list(
  list(name = "Extend regular PHCs to midnight (Sun-Thu)", cand = which(is_regular), to = "midnight"),
  list(name = "Open regular PHCs 24/7", cand = which(is_regular), to = "all_night"),
  list(name = "Open extended-hours PHCs 24/7", cand = which(is_ext), to = "all_night"),
  list(name = "Open weekday-only extended PHCs on Fri-Sat", cand = which(is_wd_ext), to = "weekend"))
run_strategy <- function(s, budget = 6000, bw = rep(1, nrow(blocks)),
                         Am = A, cb = cov_b, co = cells_of) {
  hw <- blocks$hours_week * bw
  newb <- t(sapply(s$cand, function(j) tgt[[s$to]] & !O[j, ]))
  addh <- as.numeric(newb %*% blocks$hours_week)
  Ac   <- Am[s$cand, , drop = FALSE]
  cvd  <- cb
  G    <- sapply(seq_len(nrow(blocks)), function(b) as.numeric(Ac %*% (wS * !cvd[[b]])))
  used <- rep(FALSE, length(s$cand)); hrs <- 0; gain_tot <- 0
  steps <- list(tibble(strategy = s$name, n = 0L, added_hours = 0, citizen_hours = 0))
  while (hrs < budget) {
    gain <- as.numeric((G * newb) %*% hw); gain[used | addh <= 0] <- -1
    i <- which.max(gain); if (gain[i] <= 0) break
    used[i] <- TRUE; hrs <- hrs + addh[i]; gain_tot <- gain_tot + gain[i]
    j <- s$cand[i]
    for (b in which(newb[i, ])) {
      cvd[[b]][co[[j]]] <- TRUE
      G[, b] <- as.numeric(Ac %*% (wS * !cvd[[b]]))
    }
    steps[[length(steps) + 1]] <- tibble(strategy = s$name, n = sum(used), added_hours = hrs,
                                         citizen_hours = gain_tot)
  }
  bind_rows(steps)
}
d4_path <- bind_rows(lapply(strategies, run_strategy))
# Stress test: night hours (00:00-08:00) count 40% as much as other hours,
# since fewer people seek care then (no hourly demand data are published).
night_w <- ifelse(blocks$from == 0, 0.4, 1)
d4_path_w <- bind_rows(lapply(strategies, run_strategy, bw = night_w))
# Lower bound: night hours count for nothing, so only evening and weekend
# hours gained by the same picks count.
d4_path_w0 <- bind_rows(lapply(strategies, run_strategy, bw = ifelse(blocks$from == 0, 0, 1)))
at_budget <- function(p, h) { i <- max(which(p$added_hours <= h)); p$citizen_hours[i] / max(p$added_hours[i], 1) }
d4 <- bind_rows(lapply(split(d4_path, d4_path$strategy), function(p) tibble(
  strategy = p$strategy[1],
  citizen_hours_per_hour_at_1000 = at_budget(p, 1000),
  citizen_hours_per_hour_at_3000 = at_budget(p, 3000),
  centres_at_3000 = p$n[max(which(p$added_hours <= 3000))]))) |>
  left_join(bind_rows(lapply(split(d4_path_w, d4_path_w$strategy), function(p) tibble(
    strategy = p$strategy[1],
    night_weighted_0.4_per_hour_at_1000 = at_budget(p, 1000)))), by = "strategy") |>
  left_join(bind_rows(lapply(split(d4_path_w0, d4_path_w0$strategy), function(p) tibble(
    strategy = p$strategy[1],
    night_weighted_0_per_hour_at_1000 = at_budget(p, 1000)))), by = "strategy") |>
  arrange(desc(citizen_hours_per_hour_at_1000))
write_csv(d4 |> mutate(across(where(is.double), ~ round(.x, 0))),
          file.path(tabs, "d4_hours_strategies.csv"))
write_csv(d4_path |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d4_hours_strategies_path.csv"))
message("  D4: citizen-hours gained per added opening hour (first 1,000 hours/week):")
for (i in seq_len(nrow(d4)))
  message(sprintf("      %-48s %7.0f  (night weighted 0.4: %7.0f)", d4$strategy[i],
                  d4$citizen_hours_per_hour_at_1000[i], d4$night_weighted_0.4_per_hour_at_1000[i]))
pD4 <- ggplot(d4_path, aes(added_hours, citizen_hours / 1e6, colour = strategy)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#009E73", "#CC79A7", "#0072B2", "#D55E00")) +
  labs(x = "Opening hours added per week", y = "Citizen-hours of 30-min access gained per week (millions)",
       colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")

# ---- D5. Physicians needed -------------------------------------------------
# Physicians on duty during every added opening hour: one (minimum) or two (the
# male-and-female pairing of Makkah's on-duty roster, SPA 2017); 40-hour weeks.
d5 <- bind_rows(lapply(c(80, 85, 90), function(tg) {
  k <- n_needed(g_add$path, tg)
  ids <- match(g_add$path$facility_id[seq_len(k) + 1], phc$facility_id)
  h <- sum(168 - wk_hours[ids])
  tibble(target_pct = tg, conversions = k, added_hours_per_week = h,
         pct_reached = g_add$path$pct[k + 1],
         physicians_on_duty = c(1, 2), physician_fte = h * c(1, 2) / 40,
         physicians_minimum = ceiling(h * c(1, 2) / 40),    # minimums round up
         pct_of_phc_physicians = 100 * ceiling(h * c(1, 2) / 40) / sum(phcdoc_br))
}))
write_csv(d5 |> mutate(across(where(is.double), ~ round(.x, 2))),
          file.path(tabs, "d5_physicians_needed.csv"))
d5_90 <- d5 |> filter(target_pct == 90)
message(sprintf("  D5: %d conversions (90%%) add %.0f opening hours/week = %.0f-%.0f physician FTE (%.1f-%.1f%% of %s PHC physicians)",
                d5_90$conversions[1], d5_90$added_hours_per_week[1], d5_90$physician_fte[1],
                d5_90$physician_fte[2], d5_90$pct_of_phc_physicians[1], d5_90$pct_of_phc_physicians[2],
                format(sum(phcdoc_br), big.mark = ",")))

# Conversions aimed at citizens with no night option (g_er, D3): citizens
# within 30 min of a hospital ER already count as served, so conversions go
# where neither a PHC nor an ER is within 30 min.
k_er  <- n_needed(g_er$path, 96)
ids_er <- match(g_er$path$facility_id[seq_len(k_er) + 1], phc$facility_id)
ids_90 <- match(g_add$path$facility_id[seq_len(h90_k <- n_needed(g_add$path, 90)) + 1], phc$facility_id)
h_er  <- sum(168 - wk_hours[ids_er])
er_now <- pct(cov_b[[b_night]] | er)
er_96  <- g_er$path$pct[k_er + 1]
write_csv(g_er$path |> left_join(phc |> st_drop_geometry() |>
                                   select(facility_id, name_en, governorate_en, cluster_en, hours), by = "facility_id") |>
            mutate(pct = round(pct, 2)),
          file.path(tabs, "d3_night_conversions_no_option.csv"))
message(sprintf("  D5b: no night option: PHC or ER within 30 min %.1f%% now; %d conversions reach %.1f%% (+%s citizens), %.0f opening hours/week = %.0f-%.0f FTE; %d of them among the %d conversions for 90%% PHC coverage",
                er_now, k_er, er_96, format(round((er_96 - er_now) / 100 * sum(wS)), big.mark = ","),
                h_er, h_er / 40, h_er * 2 / 40, length(intersect(ids_er, ids_90)), h90_k))

# ---- D6. Cluster borders -----------------------------------------------------
# Extra travel when residents must use a facility of their own cluster rather
# than the nearest one, and the neighbouring clusters involved.
gov_major <- fac |> count(governorate_en, cluster_en) |> group_by(governorate_en) |>
  slice_max(n, n = 1, with_ties = FALSE) |> ungroup() |> select(governorate_en, cl_gov = cluster_en)
cl_gov <- gov_major$cl_gov[match(acc$Gov_EN, gov_major$governorate_en)]
border <- function(M, fac_cluster, cell_cluster, label, mask = rep(TRUE, ncol(M))) {
  near <- apply(M, 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
  t_any <- colmin(M)
  t_own <- rep(NA_real_, ncol(M))
  for (k in unique(na.omit(cell_cluster))) {
    ci <- which(cell_cluster == k); ri <- which(fac_cluster == k)
    if (length(ri)) t_own[ci] <- colmin(M[ri, ci, drop = FALSE])
  }
  ok <- !is.na(cell_cluster) & !is.na(t_any)
  other <- ok & (fac_cluster[near] != cell_cluster) %in% TRUE
  none  <- ok & is.na(t_own)
  extra <- t_own - t_any
  pairs <- tibble(own = cell_cluster, nearest = fac_cluster[near], w = wS, extra, other) |>
    filter(other, !is.na(extra)) |>
    group_by(own, nearest) |> summarise(citizens = sum(w), median_extra_min = wq(extra, w, .5),
                                        .groups = "drop") |>
    mutate(facility = label)
  # Shares that depend on the resident's cluster use citizens who can be routed
  # and assigned a cluster; the coverage loss uses all citizens, like coverage.
  # `mask` restricts the summary to some cells (sensitivity below).
  pct_rt <- function(v) 100 * sum(wS[v & mask]) / sum(wS[ok & mask])
  pct_m  <- function(v) 100 * sum(wS[v & mask]) / sum(wS[mask])
  sel <- other & !none & mask
  list(summary = tibble(facility = label,
                        pct_nearest_in_other_cluster = pct_rt(other),
                        median_extra_min = wq(extra[sel], wS[sel], .5),
                        pct_own_over_30_other_within_30 = pct_m(other & (is.na(t_own) | t_own > thr_phc) %in% TRUE &
                                                                  (t_any <= thr_phc) %in% TRUE),
                        pct_no_facility_in_own_cluster = pct_rt(none)),
       pairs = pairs)
}
bd <- list(
  border(pc[night_open, , drop = FALSE], phc$cluster_en[night_open], acc$cluster_phc,
         "Night PHC (00:00-08:00)"),
  border(pc, phc$cluster_en, cl_gov, "Daytime PHC, governorate-based clusters"),
  border(hc[hosp_gen$facility_id, , drop = FALSE], hosp_gen$cluster_en, acc$cluster_phc,
         "General hospital ED"),
  border(hc[tert$facility_id, , drop = FALSE], tert$cluster_en, acc$cluster_phc,
         "Medical city or specialized hospital"))
d6 <- bind_rows(lapply(bd, `[[`, "summary"))
d6_pairs <- bind_rows(lapply(bd, `[[`, "pairs")) |> arrange(facility, desc(citizens))
write_csv(d6 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d6_border_costs.csv"))
write_csv(d6_pairs |> mutate(citizens = round(citizens), median_extra_min = round(median_extra_min, 1)),
          file.path(tabs, "d6_cross_border_pairs.csv"))
# Sensitivity: the same border costs for citizens outside the Riyadh and Jeddah
# governorates. Both are split between clusters, so membership there rests
# entirely on the nearest-PHC rule and cannot be checked against census totals.
stopifnot(all(c("Riyadh", "Jeddah") %in% acc$Gov_EN))
out_rj <- !(acc$Gov_EN %in% c("Riyadh", "Jeddah"))
bd_rj <- list(
  border(pc[night_open, , drop = FALSE], phc$cluster_en[night_open], acc$cluster_phc,
         "Night PHC (00:00-08:00)", mask = out_rj),
  border(pc, phc$cluster_en, cl_gov, "Daytime PHC, governorate-based clusters", mask = out_rj),
  border(hc[hosp_gen$facility_id, , drop = FALSE], hosp_gen$cluster_en, acc$cluster_phc,
         "General hospital ED", mask = out_rj),
  border(hc[tert$facility_id, , drop = FALSE], tert$cluster_en, acc$cluster_phc,
         "Medical city or specialized hospital", mask = out_rj))
night_now_rj <- pct(cov_b[[b_night]], wS * out_rj)
d6_rj <- bind_rows(lapply(bd_rj, `[[`, "summary")) |>
  mutate(night_coverage_now = c(night_now_rj, NA, NA, NA),
         night_coverage_own_cluster_only = night_coverage_now - pct_own_over_30_other_within_30,
         citizens_in_scope = sum(wS[out_rj]), pct_of_all_citizens = 100 * sum(wS[out_rj]) / sum(wS))
write_csv(d6_rj |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d6_border_costs_excl_riyadh_jeddah.csv"))
message(sprintf("  D6 excluding Riyadh and Jeddah governorates (%.1f%% of citizens): only night PHC within 30 min across a border %.1f%% (all citizens %.1f%%); night coverage %.1f%% -> %.1f%% own cluster only",
                d6_rj$pct_of_all_citizens[1], d6_rj$pct_own_over_30_other_within_30[1],
                d6$pct_own_over_30_other_within_30[1], night_now_rj, d6_rj$night_coverage_own_cluster_only[1]))
message(sprintf("  D6: nearest night PHC in another cluster for %.1f%% (median %.1f min closer)",
                d6$pct_nearest_in_other_cluster[1], d6$median_extra_min[1]))

# ---- D7. Hospital departments ------------------------------------------------
# General hospitals without the department, added one at a time by how many of
# the demand population each brings within 60 min.
dept_add <- function(have_ids, cand_ids, w, label, targets) {
  cells_h <- lapply(cand_ids, function(id) which(hc[id, ] <= thr_hosp))
  w0 <- ifelse(is.na(w), 0, w)
  cvd <- !is.na(nearest_of(have_ids)) & nearest_of(have_ids) <= thr_hosp
  now <- pct(cvd, w0)
  Ah <- sparseMatrix(i = rep(seq_along(cells_h), lengths(cells_h)), j = unlist(cells_h),
                     x = 1, dims = c(length(cand_ids), ncol(hc)))
  used <- rep(FALSE, length(cand_ids)); path <- tibble(n = 0L, hospital = NA_character_, pct = now)
  while (tail(path$pct, 1) < max(targets)) {
    g <- as.numeric(Ah %*% (w0 * !cvd)); g[used] <- -1
    i <- which.max(g); if (g[i] <= 0) break
    used[i] <- TRUE; cvd[cells_h[[i]]] <- TRUE
    path <- bind_rows(path, tibble(n = sum(used), hospital = cand_ids[i], pct = pct(cvd, w0)))
  }
  names_of <- function(ids) hosp_all$name_en[match(ids, hosp_all$facility_id)]
  tibble(service = label, pct_now = now,
         target = targets,
         hospitals_needed = sapply(targets, function(t) path$n[which(path$pct >= t)[1]]),
         pct_reached = sapply(targets, function(t) path$pct[which(path$pct >= t)[1]]),
         first_five = paste(names_of(na.omit(path$hospital[2:6])), collapse = "; "))
}
cand_of <- function(flag) hosp_gen$facility_id[!(hosp_gen[[flag]] %in% 1)]
d7 <- bind_rows(
  dept_add(hosp_psy$facility_id, cand_of("psych"), acc$cell_adult, "Psychiatric inpatient", c(90, 95)),
  dept_add(hosp_ob$facility_id,  cand_of("obgyn"), acc$cell_women, "Obstetric", c(98, 99)),
  dept_add(hosp_pd$facility_id,  cand_of("peds"),  acc$cell_child, "Pediatric", c(98, 99)))
write_csv(d7 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d7_department_additions.csv"))
message(sprintf("  D7: psychiatric %.1f%% within 60 min; +%s hospitals for 90%%, +%s for 95%%",
                d7$pct_now[1], d7$hospitals_needed[1], d7$hospitals_needed[2]))

# ---- D8. In-cluster referral hub (not reported in the manuscript) ---------
# In clusters without a tertiary center for the service, the natural upgrade is
# the cluster's designated maternity or children's hospital (largest by beds;
# otherwise the largest hospital with the department). Reported: the referral
# journey from each resident's nearest same-cluster hospital to that hub.
hub_runs <- function(sec, tert_set, w, label, desig) {
  lapply(setdiff(sort(unique(sec$cluster_en)), tert_set$cluster_en), function(k) {
    h <- sec |> filter(cluster_en == k)
    d <- which(str_detect(h$scope, desig))
    pool <- if (length(d)) d else seq_len(nrow(h))
    beds <- suppressWarnings(as.numeric(h$bed_capacity[pool]))
    c  <- pool[which.max(replace(beds, is.na(beds), -1))]
    tm <- ref_mat(h, h[c, ])
    ci <- which(acc$cluster_phc == k)
    s  <- hc[h$facility_id, ci, drop = FALSE]
    b  <- apply(s, 2, function(x) if (all(is.na(x))) NA_integer_ else which.min(x))
    x  <- tm[b, 1]
    ww <- ifelse(is.na(w[ci]), 0, w[ci])
    list(row = tibble(service = label, cluster = k, hub = h$name_en[c],
                      hub_type = if (length(d)) "designated specialty hospital" else "largest hospital with the department",
                      pop_median_min = wq(x, ww, .5),
                      pop_pct_over_120 = 100 * sum(ww[!is.na(x) & x > 120]) / sum(ww[!is.na(x)]),
                      pop_pct_unreachable = 100 * sum(ww[is.na(x)]) / sum(ww)),
         ci = ci, x = x)
  })
}
hr_ob <- hub_runs(hosp_ob, tert_ob, acc$cell_women, "Obstetric", "Maternity|Obstetrics")
hr_pd <- hub_runs(hosp_pd, tert_pd, acc$cell_child, "Pediatric", "Children|Neonatal")
d8 <- bind_rows(lapply(c(hr_ob, hr_pd), `[[`, "row"))
within2h <- function(times, w) {
  ok <- !is.na(acc$cluster_phc)
  100 * sum(w[ok & !is.na(times) & times <= 120], na.rm = TRUE) / sum(w[ok], na.rm = TRUE)
}
with_hubs <- function(base, runs) { x <- base; for (r in runs) x[r$ci] <- r$x; x }
d8_nat <- tibble(
  service = c("General", "Obstetric", "Pediatric"),
  pct_in_cluster_within_120_now = c(within2h(ref_cells$general, acc$cell_saudi),
                                    within2h(ref_cells$obstetric, acc$cell_women),
                                    within2h(ref_cells$pediatric, acc$cell_child)),
  pct_in_cluster_within_120_with_hubs = c(NA,
                                          within2h(with_hubs(ref_cells$obstetric, hr_ob), acc$cell_women),
                                          within2h(with_hubs(ref_cells$pediatric, hr_pd), acc$cell_child)))
write_csv(d8_nat |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d8_in_cluster_referral_national.csv"))
write_csv(d8 |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "d8_tertiary_hub_candidates.csv"))
message(sprintf("  D8: in-cluster tertiary referral within 2 h: obstetric %.1f%% now, %.1f%% with regional hubs; pediatric %.1f%% now, %.1f%% with hubs",
                d8_nat$pct_in_cluster_within_120_now[2], d8_nat$pct_in_cluster_within_120_with_hubs[2],
                d8_nat$pct_in_cluster_within_120_now[3], d8_nat$pct_in_cluster_within_120_with_hubs[3]))

# ---- Manuscript exhibits for the decision analyses ---------------------------
# Table 1: night-time access by cluster; Table 2 source: decision options with
# their effect and cost; Figure 3: placement curve (A) and hours strategies (B).
t1 <- d2 |> filter(cluster != "National") |>
  left_join(phc |> st_drop_geometry() |> count(cluster = cluster_en, name = "phcs"), by = "cluster") |>
  left_join(d3c |> select(cluster, conversions_for_80), by = "cluster") |>
  left_join(q3_ref_cluster |> select(cluster, gen_tertiary_in_cluster), by = "cluster") |>
  left_join(d2c |> select(cluster, citizens_per_night_phc), by = "cluster") |>
  transmute(cluster, saudi = round(saudi), phcs, phcs_open_night = n_phc_open_night,
            citizens_per_night_phc = round(citizens_per_night_phc),
            pct_open_phc_30_night = round(pct_24h_phc, 1),
            pct_open_phc_30_night_all = round(pct_24h_phc_all, 1), pct_er_only = round(pct_er_only, 1),
            pct_neither = round(pct_neither, 1), conversions_for_80,
            tertiary_in_cluster = ifelse(gen_tertiary_in_cluster, "Yes", "No"))
t1 <- bind_rows(t1, tibble(cluster = "National", saudi = sum(t1$saudi), phcs = nrow(phc),
                           phcs_open_night = length(night_open),
                           citizens_per_night_phc = round(mean(night_load$load)),
                           pct_open_phc_30_night = round(tail(d2$pct_24h_phc, 1), 1),
                           pct_open_phc_30_night_all = round(tail(d2$pct_24h_phc_all, 1), 1),
                           pct_er_only = round(tail(d2$pct_er_only, 1), 1),
                           pct_neither = round(tail(d2$pct_neither, 1), 1),
                           conversions_for_80 = sum(d3c$conversions_for_80, na.rm = TRUE),
                           tertiary_in_cluster = sprintf("%d of 20", sum(t1$tertiary_in_cluster != "No"))))
write_csv(t1, file.path(tabs, "table1_decision_clusters.csv"))

h80 <- d5 |> filter(target_pct == 80); h85 <- d5 |> filter(target_pct == 85); h90 <- d5 |> filter(target_pct == 90)
add_all80 <- sum(168 - wk_hours[cl_added])
own_only <- pct(cov_b[[b_night]]) - d6$pct_own_over_30_other_within_30[1]
# Tertiary referral: clusters without a tertiary center refer outside the
# cluster network (the D8 obstetric/pediatric hubs are not in the exhibits).
ref_gen <- q3_ref |> filter(service == "general")
opt <- tibble(
  option = c("current", "resite_same_number", "resite_within_clusters", "resite_within_clusters_first_6",
             "resite_within_clusters_protected", "convert_7", "convert_20", "convert_50", "every_cluster_80", "convert_no_night_option_96",
             "own_cluster_only_at_night", "extend_to_midnight",
             "open_all_night_regular", "open_all_night_extended", "weekend_opening",
             "referral_agreements_outside", "psych_9_hospitals"),
  centres = c(length(night_open), length(g_resite$sel), sum(!cur %in% night_open), 6L,
              sum(!cur_p %in% night_open), h80$conversions[1],
              h85$conversions[1], h90$conversions[1], length(cl_added), k_er, NA, NA, NA, NA, NA,
              20 - ref_gen$n_clusters_with_tertiary, d7$hospitals_needed[1]),
  added_hours_week = c(0, 0, 0, 0, 0, h80$added_hours_per_week[1], h85$added_hours_per_week[1],
                       h90$added_hours_per_week[1], add_all80, h_er, 0, NA, NA, NA, NA, NA, NA),
  physician_fte_1 = c(0, 0, 0, 0, 0, h80$physician_fte[1], h85$physician_fte[1], h90$physician_fte[1],
                      add_all80 / 40, h_er / 40, 0, NA, NA, NA, NA, NA, NA),
  physician_fte_2 = c(0, 0, 0, 0, 0, h80$physician_fte[2], h85$physician_fte[2], h90$physician_fte[2],
                      add_all80 * 2 / 40, h_er * 2 / 40, 0, NA, NA, NA, NA, NA, NA),
  metric = c(rep("pct citizens within 30 min of an open PHC at night", 9),
             "pct citizens within 30 min of an open PHC or hospital ED at night",
             "pct citizens within 30 min of an open PHC at night",
             rep("citizen-hours gained per added opening hour (first 1,000 h/week)", 4),
             "pct citizens referred outside the cluster network (clusters without a medical city or specialized hospital; centres = clusters)",
             "pct adults within 60 min of inpatient psychiatric care"),
  before = c(NA, rep(pct(cov_b[[b_night]]), 8), er_now, pct(cov_b[[b_night]]), NA, NA, NA, NA,
             ref_gen$pop_pct_referred_outside, d7$pct_now[1]),
  after = c(pct(cov_b[[b_night]]), tail(g_resite$path$pct, 1), tail(d3m$pct, 1), mv_at(6),
            tail(d3mp$pct, 1), h80$pct_reached[1], h85$pct_reached[1], h90$pct_reached[1],
            pct_all80, er_96, own_only,
            d4$citizen_hours_per_hour_at_1000[match(sapply(strategies, `[[`, "name"), d4$strategy)][c(1, 2, 3, 4)],
            NA, d7$pct_reached[1]))
lose <- d3_reloc[c(1, 3, 2, 4), ]                    # relocation rows 2-5 of opt
opt <- opt |> mutate(physicians_min_1 = ceiling(added_hours_week / 40),   # minimums round up
                     physicians_min_2 = ceiling(added_hours_week * 2 / 40))
opt <- opt |> mutate(citizens_losing_night_phc = c(NA, lose$citizens_losing, rep(NA, n() - 5)),
                     losing_with_neither = c(NA, lose$losing_neither, rep(NA, n() - 5)))
write_csv(opt |> mutate(across(where(is.double), ~ round(.x, 1))),
          file.path(tabs, "table2_decision_options.csv"))

fig3 <- patchwork::wrap_plots(
  pD3 + labs(title = "A  Placement of overnight PHCs"),
  pD4 + labs(title = "B  Extending hours versus opening 24/7"), ncol = 1)
ggsave(file.path(figs, "fig3_night_decisions.png"), fig3, width = 8.5, height = 12,
       dpi = 300, bg = "white")
message("  Decision exhibits written (table1_decision_clusters, table2_decision_options, fig3_night_decisions).")

# ---- D9. Exact solutions and robustness of the decision outputs ------------
# Greedy selection (D3) reaches at least 63% of the optimum; the one-at-a-time
# moves have no such guarantee and can stall under the protection rule. D9
# solves the same problems exactly as integer programs (maximal covering,
# HiGHS): the best placement of the current overnight PHCs nationally and
# within clusters (each cluster keeps its count), the within-cluster placement
# with protection (no citizen now within 30 min of a night PHC and not of a
# hospital ER may lose it), the fewest conversions for 80%, 85% and 90%
# coverage, for 80% in every cluster, and for citizens with no night option. It then recomputes the decision outputs at 20 and 45 min and on
# the terrain-penalised network (section 11; run once with SLOPE_SENS=1).
suppressPackageStartupMessages(library(highs))
cover_of <- function(M, thr) {
  co <- lapply(seq_len(nrow(M)), function(j) which(M[j, ] <= thr))
  list(cells = co, A = sparseMatrix(i = rep(seq_along(co), lengths(co)), j = unlist(co),
                                    x = 1, dims = dim(M)))
}
# Cells reached by the same set of PHCs are merged into one pattern.
patterns <- function(Am, w, prot) {
  p <- Am@p
  key <- vapply(seq_len(ncol(Am)), function(c)
    if (p[c + 1] > p[c]) paste(Am@i[(p[c] + 1):p[c + 1]], collapse = ",") else "", "")
  ok <- nzchar(key)
  f  <- factor(key[ok])
  rep_col <- which(ok)[match(levels(f), key[ok])]
  list(B = t(Am[, rep_col, drop = FALSE]), W = as.numeric(tapply(w[ok], f, sum)),
       prot = as.logical(tapply(prot[ok], f, any)))
}
# Maximal covering MIP. mode "cover": maximise covered citizens; mode "count":
# fewest PHCs open at night (current ones fixed open) reaching `target`.
mclp <- function(pt, rows = NULL, rlo = NULL, rhi = NULL, xlo = 0, xhi = 1, protect = FALSE,
                 mode = "cover", target = NA, time = 900) {
  B <- pt$B; n <- ncol(B); P <- nrow(B)
  Z <- function(k) Matrix(0, k, P, sparse = TRUE)
  Acon <- cbind(-B, Diagonal(P)); lo <- rep(-Inf, P); hi <- rep(0, P)
  if (protect && any(pt$prot)) {
    Acon <- rbind(Acon, cbind(B[pt$prot, , drop = FALSE], Z(sum(pt$prot))))
    lo <- c(lo, rep(1, sum(pt$prot))); hi <- c(hi, rep(Inf, sum(pt$prot)))
  }
  if (!is.null(rows)) {
    Acon <- rbind(Acon, cbind(rows, Z(nrow(rows)))); lo <- c(lo, rlo); hi <- c(hi, rhi)
  }
  if (mode == "count") {
    Acon <- rbind(Acon, sparseMatrix(i = rep(1, P), j = n + seq_len(P), x = pt$W, dims = c(1, n + P)))
    lo <- c(lo, target); hi <- c(hi, Inf)
  }
  s <- highs_solve(L = if (mode == "cover") c(rep(0, n), pt$W) else c(rep(1, n), rep(0, P)),
                   lower = c(rep_len(xlo, n), rep(0, P)), upper = c(rep_len(xhi, n), rep(1, P)),
                   A = Acon, lhs = lo, rhs = hi, types = c(rep("I", n), rep("C", P)),
                   maximum = mode == "cover",
                   control = highs_control(time_limit = time, log_to_console = FALSE))
  x <- which(s$primal_solution[seq_len(n)] > 0.5)
  list(x = x, status = s$status_message, gap = s$info$mip_gap, bound = s$info$mip_dual_bound,
       covered = sum(pt$W[as.numeric(B[, x, drop = FALSE] %*% rep(1, length(x))) > 0]))
}
cl_rows <- function(cl) {
  u <- sort(unique(cl))
  sparseMatrix(i = match(cl, u), j = seq_along(cl), x = 1, dims = c(length(u), length(cl)))
}
n_night_cl <- as.numeric(table(factor(phc$cluster_en[night_open], levels = sort(unique(phc$cluster_en)))))
W_all <- sum(wS)
d9_picks <- list()
# Decision outputs for one travel-time matrix and threshold.
decision_outputs <- function(M, thr, t_er, label, exact = TRUE) {
  cv <- cover_of(M, thr)
  cov_alt <- function(idx) { v <- logical(ncol(M)); for (j in idx) v[cv$cells[[j]]] <- TRUE; v }
  er_alt <- !is.na(t_er) & t_er <= thr
  now <- cov_alt(night_open)
  # greedy conversions from the current network, as in D3
  cvd <- now; picks <- integer(0); path <- pct(cvd)
  while (tail(path, 1) < 90) {
    g <- as.numeric(cv$A %*% (wS * !cvd)); g[c(night_open, picks)] <- -1
    j <- which.max(g); if (g[j] <= 0) break
    picks <- c(picks, j); cvd[cv$cells[[j]]] <- TRUE; path <- c(path, pct(cvd))
  }
  need <- function(tgt) if (any(path >= tgt)) which(path >= tgt)[1] - 1L else NA_integer_
  base7 <- match(g_add$path$facility_id[2:8], phc$facility_id)
  d9_picks[[label]] <<- tibble(scenario = label, order = seq_len(min(7, length(picks))),
                               facility_id = phc$facility_id[picks[seq_len(min(7, length(picks)))]]) |>
    mutate(name_en = nm_of(facility_id), governorate_en = gov_of(facility_id),
           cluster_en = phc$cluster_en[match(facility_id, phc$facility_id)])
  # hours strategies (Q3): extend to midnight versus 24/7, first 1,000 hours
  cb_alt <- lapply(seq_len(nrow(blocks)), function(b) cov_alt(which(O[, b])))
  q3 <- sapply(strategies, function(s)
    at_budget(run_strategy(s, budget = 1200, Am = cv$A, cb = cb_alt, co = cv$cells), 1000))
  first7 <- picks[seq_len(min(7, length(picks)))]
  out <- tibble(scenario = label, threshold_min = thr,
                night_coverage = pct(now), night_er_only = pct(!now & er_alt),
                night_neither = pct(!now & !er_alt),
                conversions_for_80 = need(80), conversions_for_90 = need(90),
                first7_same_as_base = sum(picks[1:7] %in% base7),
                first7_cities_shared = length(intersect(gov_of(phc$facility_id[first7]),
                                                        gov_of(phc$facility_id[base7]))),
                q3_midnight = q3[1], q3_24_7 = q3[2], q3_24_7_extended = q3[3], q3_weekend = q3[4],
                q3_ratio = q3[2] / q3[1])
  if (exact) {
    pt <- patterns(cv$A, wS, now & !er_alt)
    pr <- mclp(pt, rows = cl_rows(phc$cluster_en), rlo = n_night_cl, rhi = n_night_cl, protect = TRUE)
    m80 <- mclp(pt, xlo = as.numeric(seq_len(nrow(M)) %in% night_open), mode = "count",
                target = 0.80 * W_all)
    out <- out |> mutate(protected_relocation_exact = 100 * pr$covered / W_all,
                         protected_status = pr$status, protected_gap = pr$gap,
                         min_conversions_80_exact = length(m80$x) - length(night_open),
                         min80_status = m80$status, min80_gap = m80$gap)
  }
  out
}

# Exact solutions at the primary threshold
pt30 <- patterns(A, wS, cov_b[[b_night]] & !er)
ex <- list(
  national   = mclp(pt30, rows = Matrix(1, 1, nrow(A), sparse = TRUE),
                    rlo = length(night_open), rhi = length(night_open)),
  within     = mclp(pt30, rows = cl_rows(phc$cluster_en), rlo = n_night_cl, rhi = n_night_cl),
  protected  = mclp(pt30, rows = cl_rows(phc$cluster_en), rlo = n_night_cl, rhi = n_night_cl,
                    protect = TRUE))
for (tg in c(80, 85, 90))
  ex[[paste0("min_", tg)]] <- mclp(pt30, xlo = as.numeric(seq_len(nrow(A)) %in% night_open),
                                   mode = "count", target = tg / 100 * W_all)
d9_exact <- tibble(
  problem = c("Same number, best placement nationally", "Same number, within clusters",
              "Same number, within clusters, protected",
              "Fewest conversions for 80%", "Fewest conversions for 85%", "Fewest conversions for 90%"),
  heuristic = c(tail(g_resite$path$pct, 1), tail(d3m$pct, 1), tail(d3mp$pct, 1),
                n_needed(g_add$path, 80), n_needed(g_add$path, 85), n_needed(g_add$path, 90)),
  exact = c(sapply(ex[1:3], function(e) 100 * e$covered / W_all),
            sapply(ex[4:6], function(e) length(e$x) - length(night_open))),
  status = sapply(ex, `[[`, "status"), mip_gap = sapply(ex, `[[`, "gap"),
  shifts_moved = c(sapply(ex[1:3], function(e) sum(!e$x %in% night_open)), NA, NA, NA))
ex_prot <- ex$protected$x
# Each cluster to 80% from its own PHCs (D3 cluster targets) and the fewest
# conversions for 96% within 30 min of an open PHC or hospital ER (D5b).
xlo_now <- as.numeric(seq_len(nrow(A)) %in% night_open)
drop0 <- function(pt) { k <- pt$W > 0; list(B = pt$B[k, , drop = FALSE], W = pt$W[k], prot = pt$prot[k]) }
cl_ex <- lapply(sort(unique(na.omit(acc$cluster_phc))), function(k) {
  wk <- ifelse(acc$cluster_phc %in% k, wS, 0)
  if (pct(cov_b[[b_night]], wk) >= 80) return(list(n = 0L, status = "Optimal", gap = 0))
  e <- mclp(drop0(patterns(A, wk, logical(ncol(A)))), xlo = xlo_now,
            xhi = as.numeric(seq_len(nrow(A)) %in% union(night_open, which(phc$cluster_en == k))),
            mode = "count", target = 0.80 * sum(wk))
  list(n = length(e$x) - length(night_open), status = e$status, gap = e$gap)
})
e_er <- mclp(drop0(patterns(A, ifelse(er, 0, wS), logical(ncol(A)))), xlo = xlo_now,
             mode = "count", target = 0.96 * W_all - sum(wS[er]))
cl_status <- unique(sapply(cl_ex, `[[`, "status"))
d9_exact <- d9_exact |> mutate(
  citizens_losing = c(sapply(ex[1:3], function(e) sum(wS[cov_b[[b_night]] & !covered(e$x)])), NA, NA, NA),
  losing_neither  = c(sapply(ex[1:3], function(e) sum(wS[cov_b[[b_night]] & !covered(e$x) & !er])), NA, NA, NA))
d9_exact <- bind_rows(d9_exact, tibble(
  problem = c("Fewest conversions for 80% in every cluster (total)",
              "Fewest conversions for 96% within 30 min of an open PHC or ED"),
  heuristic = c(sum(d3c$conversions_for_80, na.rm = TRUE), k_er),
  exact = c(sum(sapply(cl_ex, `[[`, "n")), length(e_er$x) - length(night_open)),
  status = c(paste(cl_status, collapse = "; "), e_er$status),
  mip_gap = c(max(sapply(cl_ex, `[[`, "gap")), e_er$gap)))
write_csv(d9_exact |> mutate(across(where(is.double), ~ round(.x, 2))),
          file.path(tabs, "d9_exact_solutions.csv"))
message("  D9 exact: ", paste(sprintf("%s %.2f vs %.2f (%s)", d9_exact$problem, d9_exact$exact,
                                     d9_exact$heuristic, d9_exact$status), collapse = "; "))

# Decision outputs at 20 and 45 min and on the terrain-penalised network
d9_sens <- bind_rows(
  decision_outputs(pc, 30, acc$t_gen, "Primary (30 min)", exact = FALSE) |>
    mutate(protected_relocation_exact = 100 * ex$protected$covered / W_all,
           protected_status = ex$protected$status, protected_gap = ex$protected$gap,
           min_conversions_80_exact = length(ex$min_80$x) - length(night_open),
           min80_status = ex$min_80$status, min80_gap = ex$min_80$gap),
  decision_outputs(pc, 20, acc$t_gen, "20 min"),
  decision_outputs(pc, 45, acc$t_gen, "45 min"))
pcs_path <- file.path(dirname(graph_path), "phc_cell_times_slope.rds")
tts_path <- file.path(prc, "travel_times_slope.rds")
if (file.exists(pcs_path) && file.exists(tts_path)) {
  pc_s <- readRDS(pcs_path); tt_s9 <- readRDS(tts_path)
  stopifnot(identical(dim(pc_s), dim(pc)))
  d9_sens <- bind_rows(d9_sens, decision_outputs(pc_s, 30, tt_s9$gen$minutes, "Terrain penalty (30 min)"))
  rm(pc_s); invisible(gc())
}
write_csv(d9_sens |> mutate(across(where(is.double), ~ round(.x, 2))),
          file.path(tabs, "d9_decision_sensitivity.csv"))
write_csv(bind_rows(d9_picks), file.path(tabs, "d9_first_conversions.csv"))
message("  D9 sensitivity written: ", paste(d9_sens$scenario, collapse = ", "))

# Submission copies of the main figures, named by number: BMJ takes TIFF,
# EPS, JPEG or PDF (not PNG), at least 600 dpi for line art, one file per
# figure. Supplementary figures stay inside the supplementary file.
sub_dir <- file.path(root, "submission")
dir.create(sub_dir, showWarnings = FALSE)
for (f in list(list("Figure_1.tiff", pD1, 9, 4.5), list("Figure_2.tiff", pD2, 8, 7.5),
               list("Figure_3.tiff", fig3, 8.5, 12)))
  ggsave(file.path(sub_dir, f[[1]]), f[[2]], width = f[[3]], height = f[[4]],
         dpi = 600, compression = "lzw", bg = "white")
message("  Submission figures written (submission/Figure_1-3.tiff, 600 dpi).")

message("Done. Tables in output/tables, figures in output/figures.")
