# ==========================================================================
# Geographic accessibility to the three tiers of care under Saudi Arabia's
# health cluster (Accountable Care Organization) model
#
# Single-file analysis pipeline. Sections:
#   0. Setup and parameters
#   1. Facilities (supply)
#   2. Governorates and census population (demand)
#   3. Gridded population (WorldPop scaled to GASTAT 2022 totals)
#   4. Travel times (dodgr street-network routing, no external server)
#   5. Q1 PHC access by working-hours tier
#   6. Q2 PHC -> general hospital referral burden
#   7. Q3 Specialized hospital access (demand-matched denominators)
#   8. Q4 Tertiary completeness and cluster boundary alignment
#   9. Inequality (Gini, Theil decomposition) and spatial clustering (LISA)
#  9b. E2SFCA accessibility index (sensitivity to nearest-facility metric)
#  10. Tables and figures
#
# Expensive steps cache to data/processed/ as .rds; delete a cache file to
# force recomputation. Section 4 needs the OSM road extract (see its header);
# everything else is plain CRAN R.
# ==========================================================================

# ---- 0. Setup and parameters ---------------------------------------------

pkgs <- c("sf", "terra", "dplyr", "tidyr", "readr", "stringr", "dodgr",
          "osmextract", "units", "ineq", "spdep", "RANN", "curl", "ggplot2",
          "viridis", "ggrepel", "patchwork", "readxl")  # after 14: `::` only
new <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs[1:14], library, character.only = TRUE))

root  <- "."  # project root (cluster_access_paper/); auto-detected below
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
  # Last-resort fallback (author's machine) so the script also runs when the
  # code is pasted into a console; inert elsewhere -- edit to your own path
  # or simply run from the project root.
  cands <- c(cands, path.expand(paste0("~/Library/CloudStorage/",
    "OneDrive-Personal/Documents/1 Code/cluster_access_paper")))
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
  track_speed    = 30,   # km/h on unpaved tracks (AccessMod convention);
                         # dodgr's default car profile treats track as
                         # impassable, so this is set explicitly below
  worldpop_url   = paste0("https://data.worldpop.org/GIS/Population/",
                          "Global_2000_2020_Constrained/2020/BSGM/SAU/",
                          "sau_ppp_2020_constrained.tif"),
  e2sfca_breaks  = c(0, 15, 30, 60),      # minutes (Luo & Qi 2009 zones)
  e2sfca_weights = c(1, 0.68, 0.22),
  empanel_ratio  = 2500,                  # care-team panel: conservative upper bound of typical primary care panel sizes (Raffoul et al. 2016, JABFM)
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

fac <- read_csv(file.path(raw, "list_of_healthcare_providers.csv"),
                show_col_types = FALSE) |>
  filter(!is.na(lat), !is.na(lon)) |>
  mutate(
    tier = case_when(type == "PHC" ~ "primary",
                     type == "2ry" ~ "secondary",
                     type == "3ry" ~ "tertiary"),
    hours = case_when(hours_level == 0 ~ "regular",
                      hours_level == 1 ~ "extended",
                      hours_level == 2 ~ "h24"),
    # Service capability flags (scope strings verified against source sites)
    is_general    = tier == "secondary" &
                    str_detect(scope, "General|Emergency"),
    serves_obgyn  = tier == "secondary" &
                    str_detect(scope, "Maternity|Obstetrics"),
    serves_peds   = tier == "secondary" &
                    str_detect(scope, "Children|Neonatal"),
    serves_psych  = tier == "secondary" & str_detect(scope, "Psychiatry")
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
# No server, no Docker, no Java. The pbf sits outside the OneDrive-synced
# tree; the derived road graph is cached in data/processed/.

if (!file.exists(path.expand(PAR$pbf_path)))
  stop("Road extract not found at ", PAR$pbf_path,
       " - download it (see Section 4 header), then rerun.")

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

# Graph cache lives beside the pbf, outside the OneDrive-synced tree: the
# file is ~1 GB and cloud sync stalls the write. Keep every road component
# with >= 100 edges, not just the largest: inhabited islands with their own
# road networks (in practice one: the Farasan Islands in Jazan; the other
# retained components are disconnected mainland road pockets) would otherwise
# have their cells and facilities snapped across the sea. Cross-component
# routes return NA (unreachable).
graph_path <- file.path(dirname(path.expand(PAR$pbf_path)), "dodgr_graph.rds")
graph <- if (file.exists(graph_path)) readRDS(graph_path) else {
  # Motorcar profile with unpaved tracks passable at PAR$track_speed,
  # weighted like residential streets rather than excluded. Must go through
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
  g <- dodgr::weight_streetnet(roads, wt_profile = "motorcar",
                               wt_profile_file = wpf)
  keep <- names(which(table(g$component) >= 100))
  g <- g[g$component %in% as.integer(keep), ]
  saveRDS(g, graph_path)
  g
}
verts   <- dodgr::dodgr_vertices(graph)
vert_xy <- as.matrix(verts[, c("x", "y")])

# Snap sf points to nearest road-graph vertex; report worst snap distances.
snap_verts <- function(x, label) {
  ll <- st_coordinates(st_transform(x, 4326))
  nn <- RANN::nn2(vert_xy, ll, k = 1)
  d_km <- nn$nn.dists[, 1] * 111    # degrees -> km (approximate)
  message(sprintf("  snap %-8s: median %.2f km, max %.1f km, >5 km: %d of %d",
                  label, median(d_km), max(d_km), sum(d_km > 5), nrow(x)))
  verts$id[nn$nn.idx[, 1]]
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

tt <- cache("travel_times_grid", {
  list(
    phc_all = tt_nearest(grid_pts, phc),
    phc_ext = tt_nearest(grid_pts, phc_ext),
    phc_24  = tt_nearest(grid_pts, phc_24),
    gen     = tt_nearest(grid_pts, hosp_gen),
    obgyn   = tt_nearest(grid_pts, hosp_ob),
    peds    = tt_nearest(grid_pts, hosp_pd),
    psych   = tt_nearest(grid_pts, hosp_psy),
    tert    = tt_nearest(grid_pts, tert)
  )
})

acc <- grid |>
  mutate(t_phc   = tt$phc_all$minutes,  phc_row  = tt$phc_all$dest_row,
         t_ext   = tt$phc_ext$minutes,
         t_24    = tt$phc_24$minutes,
         t_gen   = tt$gen$minutes,      gen_row  = tt$gen$dest_row,
         t_obgyn = tt$obgyn$minutes,
         t_peds  = tt$peds$minutes,
         t_psych = tt$psych$minutes,
         t_tert  = tt$tert$minutes,     tert_row = tt$tert$dest_row,
         cluster_phc  = phc$cluster_en[phc_row],       # empanelment proxy
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

# ---- 5. Q1: PHC access by working-hours tier -----------------------------

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

# Empanelment benchmark: Saudi population per PHC vs a 2,500-patient panel
q1_empanel <- acc |>
  st_drop_geometry() |>
  group_by(Gov_EN) |>
  summarise(saudi_pop = sum(cell_saudi), .groups = "drop") |>
  left_join(fac |> filter(tier == "primary") |>
              count(governorate_en, name = "n_phc"),
            by = c("Gov_EN" = "governorate_en")) |>
  mutate(pop_per_phc   = saudi_pop / n_phc,
         teams_needed  = ceiling(saudi_pop / PAR$empanel_ratio))

# ---- 6. Q2: PHC -> general hospital referral burden ----------------------
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

# ---- 7. Q3: specialized access, demand-matched denominators --------------

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

# ---- 8. Q4: tertiary completeness and boundary alignment -----------------

q4_tert_access <- acc |> st_drop_geometry() |>
  group_by(cluster = cluster_phc) |>
  summarise(median_tert_min = wq(t_tert, cell_saudi, .5),
            pct_within_120  = cov_pct(t_tert, cell_saudi, 120),
            .groups = "drop")

# Boundary alignment: assigned cluster proxied by nearest-PHC cluster
# (empanelment assigns every person to a PHC); misalignment = nearest
# general hospital or tertiary center belongs to a different cluster.
q4_alignment <- acc |> st_drop_geometry() |>
  summarise(
    saudi_pop            = sum(cell_saudi),
    pct_gen_other_cluster  = 100 * sum(cell_saudi[cluster_gen  != cluster_phc],
                                       na.rm = TRUE) / saudi_pop,
    pct_tert_other_cluster = 100 * sum(cell_saudi[cluster_tert != cluster_phc],
                                       na.rm = TRUE) / saudi_pop)

# ---- 9. Inequality and spatial clustering --------------------------------

gov_acc <- acc |> st_drop_geometry() |>
  group_by(Gov_ID, Gov_EN) |>
  summarise(saudi_pop  = sum(cell_saudi),
            med_phc    = wq(t_phc, cell_saudi, .5),
            med_gen    = wq(t_gen, cell_saudi, .5),
            cluster    = { tb <- table(cluster_phc)
                           if (length(tb)) names(which.max(tb))
                           else NA_character_ },
            .groups = "drop")

gini_phc <- ineq::Gini(gov_acc$med_phc)
gini_gen <- ineq::Gini(gov_acc$med_gen)

# Theil decomposition (between vs within clusters), population-weighted
theil_decomp <- function(x, w, g) {
  ok <- !is.na(x) & x > 0 & !is.na(g); x <- x[ok]; w <- w[ok]; g <- g[ok]
  mu <- weighted.mean(x, w); s <- w / sum(w)
  Tt <- sum(s * (x / mu) * log(x / mu))
  bg <- tapply(seq_along(x), g, function(i) {
    sg <- sum(s[i]); mg <- weighted.mean(x[i], w[i])
    sg * (mg / mu) * log(mg / mu)
  })
  c(total = Tt, between = sum(bg), within = Tt - sum(bg))
}
theil_phc <- theil_decomp(gov_acc$med_phc, gov_acc$saudi_pop, gov_acc$cluster)

# Moran's I and LISA on governorate median PHC time
gov_sf <- gov_pop |> left_join(gov_acc |> select(Gov_ID, med_phc, med_gen),
                               by = "Gov_ID") |>
  filter(!is.na(med_phc))
nb  <- spdep::poly2nb(gov_sf, queen = TRUE)
lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
moran_phc <- spdep::moran.test(gov_sf$med_phc, lw, zero.policy = TRUE)
lisa      <- spdep::localmoran(gov_sf$med_phc, lw, zero.policy = TRUE)
gov_sf$lisa_I <- lisa[, "Ii"]
gov_sf$lisa_p <- lisa[, "Pr(z != E(Ii))"]

# ---- 9b. E2SFCA (sensitivity to the nearest-facility metric) -------------
# Enhanced two-step floating catchment area (Luo & Qi 2009 zones/weights in
# PAR): step 1 gives each facility a supply-to-weighted-demand ratio within
# its zoned catchment; step 2 sums those ratios over facilities reachable
# from each cell. Supply: 1 per PHC (no staffing data), beds for hospitals.

e2sfca_zone_w <- function(t) {
  w  <- rep(0, length(t))
  br <- PAR$e2sfca_breaks
  for (z in seq_along(PAR$e2sfca_weights))
    w[!is.na(t) & t >= br[z] & t <= br[z + 1]] <- PAR$e2sfca_weights[z]
  w
}

e2sfca <- function(fac_sf, supply, pop, label) {
  supply <- ifelse(is.na(supply) | supply <= 0,
                   median(supply[supply > 0], na.rm = TRUE), supply)
  f_v <- snap_verts(fac_sf, label)
  o_v <- snap_verts(grid_pts, "origins")
  uo  <- unique(o_v)
  trip <- list()
  for (ix in split(seq_len(nrow(fac_sf)),
                   ceiling(seq_len(nrow(fac_sf)) / 200))) {
    tm <- dodgr::dodgr_times(graph, from = unique(f_v[ix]), to = uo) / 60
    for (j in ix) {
      tj <- tm[match(f_v[j], rownames(tm)), match(o_v, colnames(tm))]
      w  <- e2sfca_zone_w(tj)
      nz <- which(w > 0)
      if (length(nz))
        trip[[length(trip) + 1]] <- data.frame(fac = j, cell = nz, w = w[nz])
    }
    message(sprintf("  E2SFCA %s: %d / %d facilities", label,
                    max(ix), nrow(fac_sf)))
  }
  tr  <- bind_rows(trip)
  dem <- tr |> mutate(wp = w * pop[cell]) |>
    group_by(fac) |> summarise(D = sum(wp), .groups = "drop")
  Rj  <- setNames(rep(0, nrow(fac_sf)), seq_len(nrow(fac_sf)))
  Rj[as.character(dem$fac)] <- supply[dem$fac] / pmax(dem$D, 1)
  ai  <- tr |> mutate(a = w * Rj[as.character(fac)]) |>
    group_by(cell) |> summarise(A = sum(a), .groups = "drop")
  out <- rep(NA_real_, nrow(grid_pts))    # NA: no facility within 60 min
  out[ai$cell] <- ai$A
  out
}

spai <- cache("e2sfca_spai", {
  list(phc = e2sfca(phc, rep(1, nrow(phc)), grid$cell_saudi, "phc"),
       gen = e2sfca(hosp_gen, hosp_gen$bed_capacity, grid$cell_saudi, "gen"))
})

gov_spai <- acc |> st_drop_geometry() |>
  mutate(spai_phc = spai$phc, spai_gen = spai$gen) |>
  group_by(Gov_ID, Gov_EN) |>
  summarise(spai_phc = weighted.mean(spai_phc, cell_saudi, na.rm = TRUE),
            spai_gen = weighted.mean(spai_gen, cell_saudi, na.rm = TRUE),
            .groups = "drop")

# ---- 10. Tables and figures ----------------------------------------------

write_csv(q1_national, file.path(tabs, "q1_national_phc_by_hours.csv"))
write_csv(q1_cluster,  file.path(tabs, "q1_cluster_phc.csv"))
write_csv(q1_empanel,  file.path(tabs, "q1_empanelment_benchmark.csv"))
write_csv(q2_cluster,  file.path(tabs, "q2_referral_by_cluster.csv"))
write_csv(q3,          file.path(tabs, "q3_specialized_access.csv"))
write_csv(q4_tert_access |> filter(!is.na(cluster)), file.path(tabs, "q4_tertiary_access.csv"))
write_csv(q4_alignment,       file.path(tabs, "q4_boundary_alignment.csv"))
write_csv(gov_spai,           file.path(tabs, "e2sfca_spai_by_governorate.csv"))
write_csv(tibble(metric  = c("gini_phc", "gini_gen",
                             names(theil_phc), "moran_I", "moran_p"),
                 value   = c(gini_phc, gini_gen, theil_phc,
                             moran_phc$estimate[1], moran_phc$p.value)),
          file.path(tabs, "q_equity_metrics.csv"))

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
  geom_sf(data = gov_pop, fill = "grey96", color = "grey88",
          linewidth = .05),
  geom_sf(data = reg, fill = NA, color = "grey60", linewidth = .2))
time_scale <- scale_color_manual(
  values = pal_time, drop = FALSE,
  guide = guide_legend(nrow = 1, override.aes = list(size = 3)))

# Fig 1 -- Q1: PHC travel time by working-hours tier, populated cells
# Panels are time-of-day scenarios (which PHCs are open), not facility
# categories: evening = extended-hours AND 24h PHCs both open.
tier_labs <- c("Daytime — all PHCs (n = 2,169)",
               "Evening — extended-hours and 24-hour open (n = 392)",
               "Night — 24-hour only (n = 87)")
acc_tier <- bind_rows(
  acc |> mutate(panel = tier_labs[1], t = t_phc),
  acc |> mutate(panel = tier_labs[2], t = t_ext),
  acc |> mutate(panel = tier_labs[3], t = t_24)) |>
  mutate(panel = factor(panel, levels = tier_labs),
         t_bin = bin_time(t))
p1 <- ggplot() + basemap +
  geom_sf(data = acc_tier, aes(color = t_bin), shape = 15, size = .12) +
  time_scale + labs(color = "Minutes to nearest PHC") +
  facet_wrap(~panel, nrow = 1) + map_theme
ggsave(file.path(figs, "map_phc_access.png"), p1, width = 13.5,
       height = 6, dpi = 300, bg = "white")

# Fig 2 -- Q3: hospital travel time by service, facilities overlaid
svc_levels <- c("General (n = 222)", "Maternity & children (n = 18)",
                "Pediatric (n = 18)", "Psychiatric (n = 19)")
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
  geom_sf(data = fac_svc, shape = 24, size = 1.2, fill = "white",
          color = "black", stroke = .4) +
  time_scale + labs(color = "Minutes to nearest hospital") +
  facet_wrap(~panel, nrow = 2) + map_theme
ggsave(file.path(figs, "map_hospital_access.png"), p2, width = 10,
       height = 11, dpi = 300, bg = "white")

# Fig 3 -- Q4: tertiary travel time, tertiary centers overlaid
p3 <- ggplot() + basemap +
  geom_sf(data = acc |> mutate(t_bin = bin_time(t_tert)),
          aes(color = t_bin), shape = 15, size = .15) +
  geom_sf(data = tert, shape = 23, size = 2, fill = "white",
          color = "black", stroke = .5) +
  time_scale + labs(color = "Minutes to nearest tertiary center") +
  map_theme
ggsave(file.path(figs, "map_tertiary_access.png"), p3, width = 8,
       height = 8.5, dpi = 300, bg = "white")

# Fig 4 -- E2SFCA SPAI quintiles by governorate (PHC and general hospital)
q5cut <- function(x) cut(x, quantile(x, seq(0, 1, .2), na.rm = TRUE),
                         include.lowest = TRUE,
                         labels = c("Q1 (lowest)", "Q2", "Q3", "Q4",
                                    "Q5 (highest)"))
gov_spai_j <- gov_pop |>
  left_join(gov_spai |> select(-Gov_EN), by = "Gov_ID")
gov_spai_sf <- bind_rows(
  gov_spai_j |> mutate(panel = "PHC", quint = q5cut(spai_phc)),
  gov_spai_j |> mutate(panel = "General hospital",
                       quint = q5cut(spai_gen))) |>
  mutate(panel = factor(panel, levels = c("PHC", "General hospital")))
p4 <- ggplot(gov_spai_sf) +
  geom_sf(aes(fill = quint), color = "grey70", linewidth = .08) +
  geom_sf(data = reg, fill = NA, color = "grey45", linewidth = .25) +
  scale_fill_viridis_d(name = "E2SFCA accessibility (quintile)",
                       option = "mako", direction = -1, begin = 0.25,
                       end = 0.95, na.value = "grey85",
                       guide = guide_legend(nrow = 1)) +
  facet_wrap(~panel, nrow = 1) + map_theme
ggsave(file.path(figs, "map_spai_governorate.png"), p4, width = 12,
       height = 6.5, dpi = 300, bg = "white")

# Fig 5 -- Governorate median PHC travel time (choropleth). The LISA test
# above yields no significant local cluster of long times (no governorate
# with p < 0.05 and Ii > 0; global Moran's I reported in q_equity_metrics):
# long median times are dispersed pockets, which this map displays; the
# LISA statistics are reported in the figure caption.
med_labels <- c("<3", "3–5", "5–10", "10–15",
                "15–30", ">30")
gov_sf$med_bin <- cut(gov_sf$med_phc, breaks = c(0, 3, 5, 10, 15, 30, Inf),
                      labels = med_labels, include.lowest = TRUE)
p5 <- ggplot(gov_sf) +
  geom_sf(aes(fill = med_bin), color = "grey70", linewidth = .08) +
  geom_sf(data = reg, fill = NA, color = "grey45", linewidth = .25) +
  scale_fill_manual(name = "Median PHC travel time (minutes)",
                    values = setNames(viridis::viridis(length(med_labels),
                                                       direction = -1),
                                      med_labels),
                    drop = FALSE, guide = guide_legend(nrow = 1)) +
  map_theme
ggsave(file.path(figs, "map_gov_med_phc.png"), p5, width = 8, height = 8.5,
       dpi = 300, bg = "white")

# ---- 10b. Cluster scorecard (the cluster as unit of comparison) ----------
# One row per cluster spanning all four sub-questions, plus within-cluster
# inequality and cross-cluster dependence. Descriptive comparison only; the
# spatial support remains the 1-km grid, and cluster membership is applied
# as an attribute after routing, never as a constraint.

gini_w <- function(x, w) {                 # population-weighted Gini
  ok <- is.finite(x) & is.finite(w) & w > 0; x <- x[ok]; w <- w[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]
  p <- cumsum(w) / sum(w); L <- cumsum(w * x) / sum(w * x)
  1 - sum(diff(c(0, p)) * (L + c(0, head(L, -1))))
}

# Count facility SITES, not directory rows. Tertiary directories list
# hospital-campus centres as separate entries: Al Baha's ten centres of
# excellence share the King Fahad Hospital - Al Bahah campus (identical
# coordinates); Hail's cardiac and diabetes-and-endocrine centres sit on the
# King Salman Specialist Hospital campus (120-300 m from the hospital pin);
# Qassim's specialist medical centre sits on the King Fahd Specialist
# Hospital campus (~240 m). Entries of the same cluster within 500 m of one
# another are merged into one site; the smallest distance between genuinely
# distinct tertiary sites is 2.8 km (Eastern: King Fahad Specialist Hospital
# Dammam vs Saud Al Babtain Cardiac Center), so the threshold is unambiguous. The merge applies to the tertiary tier ONLY: adjacent PHCs
# and secondary hospitals (e.g. a maternity hospital beside a general
# hospital) are distinct facilities and keep exact-coordinate deduplication.
# One named exception on the secondary tier: the Al Makhwah endocrinology
# and diabetes centre (ALB-010) is an outpatient centre of excellence hosted
# INSIDE Al Makhwah General Hospital (ALB-014, 4 m away) and is counted with
# its host as one secondary site. Coordinates are kept as geocoded, so
# routing and travel times are unaffected.
#
# Tier classification (tertiary vs secondary) is analyst-assigned; the
# directories do not label tiers. Two official anchors ground the scheme:
# (1) the Health Holding Company describes each cluster as comprising
# primary care centres, hospitals, and medical cities / specialised
# hospitals (health.sa/en/clusters), and (2) the reform's model of care
# organises cluster services as integrated levels of care, with primary
# care as first contact and referral upward to secondary and tertiary
# services (see manuscript references) - i.e. tier is a level of care,
# defined by function, not by premises. Medical cities and specialised
# hospitals are a distinct MoH category under the 2014 Law of Medical
# Cities and Specialized Hospitals. Tertiary therefore = the medical-city / specialised-hospital
# stratum acting as the cluster referral terminus with advanced specialised
# services (open-heart surgery, transplantation, radiation/comprehensive
# oncology), verified facility by facility. Centres of excellence
# delivering such services from a host-hospital campus are tertiary
# (function over building); specialty hospitals and centres without such
# services are secondary even when named "specialist" or "centre of
# excellence" (e.g. the Al Makhwah outpatient endocrine centre, the King
# Salman kidney centre in Riyadh, and the Qatif hereditary blood diseases
# hospital, which offers chronic hemoglobinopathy care but not bone-marrow
# transplantation).
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

fac_sites <- fac |>
  filter(facility_id != "ALB-010") |>  # hosted inside Al Makhwah GH; one site
  distinct(cluster_en, tier, lat, lon) |>
  group_by(cluster_en, tier) |>
  mutate(site = if (first(tier) == "tertiary") campus_id(lat, lon)
                else seq_len(n())) |>
  ungroup() |>
  distinct(cluster_en, tier, site)

fac_n <- fac_sites |>
  count(cluster_en, tier) |>
  pivot_wider(names_from = tier, values_from = n, values_fill = 0,
              names_prefix = "n_") |>
  left_join(fac |> filter(tier == "primary", hours == "h24") |>
              count(cluster_en, name = "n_phc_24h"), by = "cluster_en") |>
  mutate(n_phc_24h = coalesce(n_phc_24h, 0L))
stopifnot(sum(fac_n$n_primary)   == 2169,
          sum(fac_n$n_secondary) == 282,
          sum(fac_n$n_tertiary)  == 16,
          fac_n$n_tertiary[fac_n$cluster_en == "Hail"] == 1,
          fac_n$n_tertiary[fac_n$cluster_en == "Qassim"] == 1,
          fac_n$n_tertiary[fac_n$cluster_en == "Al Baha"] == 1,
          fac_n$n_tertiary[fac_n$cluster_en == "Eastern"] == 2,
          fac_n$n_tertiary[fac_n$cluster_en == "Riyadh Second"] == 1)

scorecard <- acc |> st_drop_geometry() |>
  mutate(spai_phc_cell = spai$phc) |>
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
    med_tert_min      = wq(t_tert, cell_saudi, .5),
    pct120_tert       = cov_pct(t_tert, cell_saudi, 120),
    pct_gen_other     = 100 * sum(cell_saudi[cluster_gen != cluster_phc],
                                  na.rm = TRUE) / saudi_pop,
    pct_tert_other    = 100 * sum(cell_saudi[cluster_tert != cluster_phc],
                                  na.rm = TRUE) / saudi_pop,
    gini_phc_within   = gini_w(t_phc, cell_saudi),
    spai_phc          = weighted.mean(spai_phc_cell, cell_saudi,
                                      na.rm = TRUE),
    .groups = "drop") |>
  left_join(q2_cluster |>
              select(cluster_en, median_ref_min, ref_over_60 = n_over_60),
            by = c("cluster" = "cluster_en")) |>
  left_join(fac_n, by = c("cluster" = "cluster_en")) |>
  mutate(saudi_per_phc = saudi_pop / n_primary)   # vs typical care-team panel sizes

write_csv(scorecard, file.path(tabs, "cluster_scorecard.csv"))

# Cluster comparison figure: the four coverage indicators, one dot plot
sc_long <- scorecard |>
  select(cluster,
         `PHC within 30 min`                       = pct30_phc,
         `24-hour PHC within 30 min`               = pct30_phc_24h,
         `Maternity hospital within 60 min (women 15–49)` =
           pct60_obgyn_women,
         `Tertiary center within 120 min`          = pct120_tert) |>
  pivot_longer(-cluster, names_to = "indicator", values_to = "pct") |>
  mutate(cluster = factor(cluster, levels = scorecard |>
                            arrange(pct30_phc_24h) |> pull(cluster)),
         indicator = factor(indicator, levels = c(
           "PHC within 30 min", "24-hour PHC within 30 min",
           "Maternity hospital within 60 min (women 15–49)",
           "Tertiary center within 120 min")))
p6 <- ggplot(sc_long, aes(pct, cluster, color = indicator)) +
  geom_line(aes(group = cluster), color = "grey80", linewidth = .4) +
  geom_point(size = 2.4) +
  scale_color_manual(values = c("#1b7837", "#d95f02", "#7570b3",
                                "#e7298a"), name = NULL) +
  scale_x_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Saudi-citizen population covered", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2))
ggsave(file.path(figs, "fig_cluster_scorecard.png"), p6, width = 8.5,
       height = 7, dpi = 300, bg = "white")

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
clu_xy <- st_coordinates(st_point_on_surface(st_geometry(clu_diss)))
split_gov <- gov_cl_sf |> filter(n_clu > 1)
pB <- ggplot() +
  geom_sf(data = gov_cl_sf, aes(fill = cluster_en), color = "white",
          linewidth = .12, alpha = .9, show.legend = FALSE) +
  geom_sf(data = clu_diss, fill = NA, color = "grey15", linewidth = .45) +
  geom_sf(data = st_point_on_surface(st_geometry(split_gov)),
          shape = 8, size = 1.8, color = "black", stroke = .6) +
  ggrepel::geom_text_repel(
    data = clu_diss |> st_drop_geometry() |>
      mutate(x = clu_xy[, 1], y = clu_xy[, 2]),
    aes(x, y, label = cluster_en), size = 3, fontface = "bold",
    color = "grey10", bg.color = "white", bg.r = .12, seed = 1) +
  scale_fill_manual(values = pal20, na.value = "grey92") +
  labs(caption = paste("* Governorate hosts facilities of more than one",
                       "cluster; fill shows the majority cluster.")) +
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
  "t_phc",   "cell_saudi", "PHC — daytime (all)",
  "t_ext",   "cell_saudi", "PHC — evening (extended + 24h)",
  "t_24",    "cell_saudi", "PHC — night (24h only)",
  "t_gen",   "cell_saudi", "General hospital",
  "t_obgyn", "cell_women", "Maternity & children hospital",
  "t_peds",  "cell_child", "Pediatric hospital",
  "t_psych", "cell_adult", "Psychiatric hospital",
  "t_tert",  "cell_saudi", "Tertiary center")

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
ord <- clu_tt |> filter(destination == "PHC — night (24h only)") |>
  arrange(median_min) |> pull(cluster)
p7 <- clu_tt |>
  mutate(cluster = factor(cluster, levels = ord),
         destination = factor(destination,
                              levels = dest_set$destination)) |>
  ggplot(aes(median_min, cluster)) +
  geom_linerange(aes(xmin = p25_min, xmax = p90_min),
                 color = "grey60", linewidth = .8) +
  geom_point(size = 1.9, color = "#1b7837") +
  facet_wrap(~destination, nrow = 2, scales = "free_x") +
  labs(x = "Minutes (dot = population-weighted median; bar = p25–p90)",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(figs, "fig_cluster_travel_times.png"), p7, width = 11,
       height = 7.5, dpi = 300, bg = "white")

# Small-multiple maps: one panel per cluster, cells colored by binned
# travel time, cluster boundary from the facility-majority dissolve.
clu_map <- function(v, legend_lab, file) {
  cls <- sort(clu_diss$cluster_en)
  plots <- lapply(cls, function(cl) {
    cells <- acc |> filter(cluster_phc == cl) |>
      mutate(t_bin = bin_time(.data[[v]]))
    bound <- clu_diss |> filter(cluster_en == cl)
    bb <- st_bbox(bound)
    ggplot() +
      geom_sf(data = bound, fill = "grey96", color = "grey40",
              linewidth = .3) +
      # legend from the first panel only; otherwise patchwork shows
      # near-duplicate legends side by side when any panel differs
      geom_sf(data = cells, aes(color = t_bin), shape = 15, size = .25,
              show.legend = identical(cl, cls[1])) +
      scale_color_manual(values = pal_time, drop = FALSE,
                         name = legend_lab,
                         guide = guide_legend(nrow = 1,
                           override.aes = list(size = 3))) +
      coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)]) +
      labs(title = cl) +
      theme_void(base_size = 9) +
      theme(plot.title = element_text(size = 8, face = "bold",
                                      hjust = .5))
  })
  pw <- patchwork::wrap_plots(plots, ncol = 5, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(figs, file), pw, width = 12, height = 10.5,
         dpi = 300, bg = "white")
}
clu_map("t_24",  "Minutes to nearest 24-hour PHC",
        "map_cluster_phc24.png")
clu_map("t_gen", "Minutes to nearest general hospital",
        "map_cluster_gen.png")

# ---- 10e. Workforce and capacity linkage (MoH Yearbook 2024) -------------
# Links spatial access to staffing: MoH Statistical Yearbook 2024
# (data/raw). Manpower tables (2-24, 2-29, 2-45) are published by the 20
# former MoH branches/offices; the yearbook crosswalk (sheet "regions")
# maps branches onto clusters, but Riyadh's three clusters share one
# branch and Jeddah's two share another, so workforce indicators use 17
# cluster groups. Facility counts (table 2-4) ARE at true cluster level
# and validate the facility file. Parsed values verified against the
# yearbook's own national-total blocks (stopifnot checks below).

yb <- file.path(raw, "Statistical-Yearbook-2024 1.xlsx")
num <- function(x) suppressWarnings(as.numeric(x))

# Branch column order in tables 2-24 / 2-45 (5 stacked blocks x 4 branches,
# S/NS/Total triples in cols 3-14; category total row = header + offset).
br_2445 <- c("Riyadh", "The Holy Capital", "Jeddah", "Taif",
             "Madinah", "Qassim", "Eastern", "Al Ahsa",
             "Hafr Al Batin", "Asir", "Bishah", "Tabuk",
             "Hail", "Northern Borders", "Jazan", "Najran",
             "Al Baha", "Al Jouf", "Qurayyat", "Qunfudah")
tot_cols <- c(5, 8, 11, 14)

d24 <- suppressMessages(readxl::read_excel(yb, "2-24", col_names = FALSE))
grab24 <- function(off) setNames(unlist(lapply(c(4, 33, 62, 91, 120),
  function(h) num(unlist(d24[h + off, tot_cols])))), br_2445)
phys_br  <- grab24(6)                       # physicians (total row)
stopifnot(sum(phys_br)  == num(d24[155, 3]) + num(d24[155, 7]))

d45 <- suppressMessages(readxl::read_excel(yb, "2-45", col_names = FALSE))
grab45 <- function(off) setNames(unlist(lapply(c(4, 30, 56, 82, 108),
  function(h) num(unlist(d45[h + off, tot_cols])))), br_2445)
dent_br   <- grab45(9)                      # PHC dentists (Yearbook lists them among PHC "physicians"; cf. 2-43)
phcdoc_br <- grab45(24) - dent_br           # PHC physicians excluding dentists
stopifnot(sum(dent_br) == 2408)
phcgpfm_br <- grab45(6) + grab45(18)        # GP + family medicine
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
obgyn_br <- grab29("نساء وولادة")
peds_br  <- grab29("طب أطفال")
psych_br <- grab29("أمراض نفسية")
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

# Branch -> cluster-group crosswalk (yearbook sheet: former directorates)
br2grp <- c(Riyadh = "Riyadh (3 clusters)", `The Holy Capital` = "Makkah",
            Qunfudah = "Makkah", Jeddah = "Jeddah (2 clusters)",
            Taif = "Al Taif", Madinah = "Madinah", Qassim = "Qassim",
            Eastern = "Eastern", `Al Ahsa` = "Al Ahsa",
            `Hafr Al Batin` = "Hafar Al Batin", Asir = "Aseer",
            Bishah = "Aseer", Tabuk = "Tabuk", Hail = "Hail",
            `Northern Borders` = "Northern Borders", Jazan = "Jazan",
            Najran = "Najran", `Al Baha` = "Al Baha",
            `Al Jouf` = "Al Jouf", Qurayyat = "Al Jouf")
clu2grp <- function(cl) case_when(
  grepl("^Riyadh", cl) ~ "Riyadh (3 clusters)",
  grepl("^Jeddah", cl) ~ "Jeddah (2 clusters)",
  TRUE ~ cl)
agg_br <- function(v) tapply(v, br2grp[names(v)], sum)

wf_br <- tibble(grp = names(agg_br(phys_br)),
                physicians = as.numeric(agg_br(phys_br)),
                phc_physicians = as.numeric(agg_br(phcdoc_br)),
                phc_gp_fm      = as.numeric(agg_br(phcgpfm_br)),
                obgyn = as.numeric(agg_br(obgyn_br)),
                peds  = as.numeric(agg_br(peds_br)),
                psych = as.numeric(agg_br(psych_br)),
                gp    = as.numeric(agg_br(gp_br)))

# Demand denominators and access medians per cluster group (our pipeline)
grp_acc <- acc |> st_drop_geometry() |> filter(!is.na(cluster_phc)) |>
  group_by(grp = clu2grp(cluster_phc)) |>
  summarise(saudi_pop = sum(cell_saudi),
            women_pop = sum(cell_women),
            child_pop = sum(cell_child),
            med_phc_min   = wq(t_phc,   cell_saudi, .5),
            med_24h_min   = wq(t_24,    cell_saudi, .5),
            med_obgyn_min = wq(t_obgyn, cell_women, .5),
            med_peds_min  = wq(t_peds,  cell_child, .5),
            med_psych_min = wq(t_psych, cell_adult, .5),
            .groups = "drop")
grp_phc <- fac |> filter(tier == "primary") |>
  count(grp = clu2grp(cluster_en), name = "n_phc")

workforce <- grp_acc |>
  left_join(grp_phc, by = "grp") |>
  left_join(wf_br,   by = "grp") |>
  mutate(phc_per_10k_saudi   = 1e4 * n_phc / saudi_pop,
         saudis_per_phc      = saudi_pop / n_phc,
         phys_per_10k_saudi  = 1e4 * physicians / saudi_pop,
         phc_phys_per_phc    = phc_physicians / n_phc,
         phc_phys_per_10k    = 1e4 * phc_physicians / saudi_pop,
         obgyn_per_10k_women = 1e4 * obgyn / women_pop,
         peds_per_10k_child  = 1e4 * peds / child_pop,
         psych_per_100k      = 1e5 * psych / saudi_pop,
         gp_share_pct        = 100 * gp / physicians) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))
write_csv(workforce, file.path(tabs, "workforce_by_cluster_group.csv"))

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

# Figure: workforce supply vs spatial access, one panel per service.
# x = population-weighted median travel time, y = staffing rate.
sup_panels <- bind_rows(
  workforce |> transmute(grp, x = med_phc_min, y = phc_phys_per_10k,
    panel = "PHC physicians per 10,000 vs minutes to PHC"),
  workforce |> transmute(grp, x = med_obgyn_min, y = obgyn_per_10k_women,
    panel = "ObGyn per 10,000 women 15–49 vs minutes to maternity hospital"),
  workforce |> transmute(grp, x = med_peds_min, y = peds_per_10k_child,
    panel = "Pediatricians per 10,000 children vs minutes to pediatric hospital"),
  workforce |> transmute(grp, x = med_psych_min, y = psych_per_100k,
    panel = "Psychiatrists per 100,000 vs minutes to psychiatric hospital")) |>
  mutate(panel = factor(panel, levels = unique(panel)))
med_ref <- sup_panels |> group_by(panel) |>
  summarise(mx = median(x, na.rm = TRUE), my = median(y, na.rm = TRUE),
            .groups = "drop")
p8 <- ggplot(sup_panels, aes(x, y)) +
  geom_vline(data = med_ref, aes(xintercept = mx),
             linetype = 2, color = "grey65") +
  geom_hline(data = med_ref, aes(yintercept = my),
             linetype = 2, color = "grey65") +
  geom_point(color = "#1b7837", size = 2) +
  ggrepel::geom_text_repel(aes(label = grp), size = 2.6, seed = 1,
                           max.overlaps = 20) +
  facet_wrap(~panel, scales = "free") +
  labs(x = "Population-weighted median travel time (minutes)",
       y = "Staffing rate (MoH Statistical Yearbook 2024)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5))
ggsave(file.path(figs, "fig_supply_vs_access.png"), p8, width = 10.5,
       height = 8, dpi = 300, bg = "white")

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

  # Cluster-level shift for PHC access (mountain clusters move most)
  sens_cl <- acc |> st_drop_geometry() |>
    mutate(t_slope = tt_s$phc$minutes) |>
    group_by(cluster = cluster_phc) |>
    summarise(pct30_base  = cov_pct(t_phc,   cell_saudi, 30),
              pct30_slope = cov_pct(t_slope, cell_saudi, 30),
              .groups = "drop") |>
    mutate(delta_pp = pct30_slope - pct30_base)
  write_csv(sens_cl, file.path(tabs, "sensitivity_slope_by_cluster.csv"))
  message("  slope sensitivity written")
}

# ---- 11. Supplement sensitivity tables (S1, S2, S11) ---------------------
# S1: threshold coverage grid, primary (Saudi / demand-matched) denominators.
# S2: the same grid under the total-resident denominator (manual defines the
#     eligible population as all residents).
# S11: cluster-assignment sensitivity — the study assigns each grid cell to
#     the cluster of its nearest PHC (empanelment proxy); the alternative
#     assigns it to the cluster operating the majority of facilities in its
#     governorate. GASTAT publishes no population by health cluster, so this
#     internal check replaces external validation of cluster totals.

dest_spec <- tribble(
  ~col,      ~destination,                 ~wcol_primary, ~denom_primary,
  "t_phc",   "PHC, any",                   "cell_saudi",  "Saudi citizens",
  "t_ext",   "PHC, extended or 24-hour",   "cell_saudi",  "Saudi citizens",
  "t_24",    "PHC, 24-hour",               "cell_saudi",  "Saudi citizens",
  "t_gen",   "General hospital",           "cell_saudi",  "Saudi citizens",
  "t_obgyn", "Obstetric hospital",         "cell_women",  "Saudi women 15-49",
  "t_peds",  "Paediatric hospital",        "cell_child",  "Saudi children <15",
  "t_psych", "Psychiatric hospital",       "cell_adult",  "Saudi adults 18+",
  "t_tert",  "Tertiary centre",            "cell_saudi",  "Saudi citizens")

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


# ---- Manuscript main-text tables ----------------------------------------
# Table 1 and Table 2 exactly as they appear in the manuscript, so the
# pipeline reproduces the main text as well as the supplement. Site counts
# fold co-located directory entries into one site per campus: identical
# coordinates collapse the ten Al Baha centres of excellence, and four
# named satellite entries are folded into their host hospitals (two centres
# at King Salman Specialist Hospital in Hail, one at King Fahd Specialist
# Hospital in Buraidah, one outpatient centre at Al Makhwah General
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
  distinct(tier, lat, lon, .keep_all = TRUE)
stopifnot(sum(fac_sites$tier == "primary")   == 2169,
          sum(fac_sites$tier == "secondary") == 282,
          sum(fac_sites$tier == "tertiary")  == 16)

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

t2_dest <- c("PHC, any"                 = "PHC, any (daytime)",
             "PHC, extended or 24-hour" = "PHC, extended or 24-hour (evening)",
             "PHC, 24-hour"             = "PHC, 24-hour (night)",
             "General hospital"         = "General hospital",
             "Obstetric hospital"    = "Hospital with obstetric services",
             "Paediatric hospital"   = "Hospital with paediatric services",
             "Psychiatric hospital"  = "Hospital with psychiatric services",
             "Tertiary centre"          = "Tertiary centre")
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

message("Done. Tables in output/tables, figures in output/figures.")
