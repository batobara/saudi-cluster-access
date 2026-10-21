# After-hours primary care and hospital referral in Saudi Arabia's health clusters

Code and derived tables for the manuscript:

> *After-hours primary care and hospital referral in Saudi Arabia's health clusters: a national cross-sectional
> study of access, placement, and opening hours.* Batobara A, Alfaraidhy M, Lim C. Submitted to BMJ Global Health.

**Version 2.0 (October 2026).** This version reproduces the submitted manuscript: PHC opening hours from the
cluster directories (accessed 28 September 2026; 90 PHCs open on weekday nights), coverage for every hour of the
week, overnight PHC placement, hours strategies, cluster borders, and referral. Version 1 (5 September 2026)
carried an earlier analysis of the same facility census and an earlier manuscript title.

## Files

| File | Description | GitHub | Zenodo |
|---|---|:---:|:---:|
| `README.md` | This file. | ✓ | ✓ |
| `analysis.R` | The complete analysis pipeline in one R script. Writes every main-text and supplementary table (CSV) and figure (PNG). | ✓ | ✓ |
| `renv.lock` | R 4.6.0 and the 93 CRAN packages (with dependencies) at the versions used, including `highs` for the exact integer programs. `renv::restore()` rebuilds the environment. | ✓ | ✓ |
| `LICENSE` | MIT license for the code (on Zenodo, the license is set in the record). | ✓ | |
| `referral_services.csv` | Input: the referral level (medical cities and specialized hospitals) and the hospitals with a documented advanced service, with public sources. | | ✓ |
| `derived_tables.zip` | All 53 derived tables (CSV) written by `analysis.R`, mapped to the manuscript below. | | ✓ |

The figures are written by `analysis.R` and appear in the article and its supplement; they are not deposited
separately.

## Where each table and figure comes from

Section names refer to the `# ----` headers in `analysis.R`. Some CSV names keep the numbering of an earlier
draft; this map gives the current numbers.

| Manuscript item | Output file(s) | Section of `analysis.R` |
|---|---|---|
| Table 1 | `table1_decision_clusters.csv` | Manuscript exhibits for the decision analyses |
| Table 2 | `table2_decision_options.csv` | Manuscript exhibits for the decision analyses |
| Figure 1 | `d1_coverage_by_hour.png` | D1. Coverage through the week |
| Figure 2 | `d2_night_gap_map.png` | D2. Night-time gap types |
| Figure 3 | `fig3_night_decisions.png` | Manuscript exhibits for the decision analyses |
| Table S1 | PHC opening hours (restricted input `phc_opening_hours.csv`) | 1. Facilities; D1 |
| Table S2 | Hospital service audit (restricted input `hospital_services_evidence.csv`) | 1. Facilities |
| Table S3 | `cluster_phc_population.csv` | 10e. Physicians and facility validation |
| Table S4 | `table_S11_cluster_assignment_sensitivity.csv` | 11. Supplement sensitivity tables |
| Table S5 | `table_S12_region_census_validation.csv` | 11. Supplement sensitivity tables |
| Table S6 | `table_S1_threshold_coverage_saudi.csv` | 11. Supplement sensitivity tables |
| Table S7 | `d1_coverage_by_time.csv`, `sensitivity_unlisted_hours.csv`, `d2_night_gap_by_cluster.csv`, `d2b_night_children_women.csv` | D1, D2, D2b |
| Table S8 | `d3_night_siting.csv`, `d3_night_resite_from_scratch.csv`, `d3_night_conversions_order.csv`, `d3_night_moves_within_clusters.csv`, `d3_night_cluster_targets.csv`, `d3_night_conversions_no_option.csv`, `d9_exact_solutions.csv`, `d9_decision_sensitivity.csv`, `d9_first_conversions.csv`, `d3_night_relocation_effects.csv` (part G; protected moves in `d3_night_moves_protected.csv`) | D3. Overnight PHC placement; D5; D9 |
| Table S9 | `d4_hours_strategies.csv`, `d5_physicians_needed.csv` | D4, D5 |
| Table S10 | `q2_referral_by_cluster.csv` | 6. Q4: referral from PHC to general hospital |
| Table S11 | `q3_referral_to_tertiary.csv`, `q3_referral_to_tertiary_by_cluster.csv`, `referral_by_service.csv`, `referral_by_service_cluster.csv` | 7c. Q4: referral from general hospital to the referral level; 7d. Q4: referral by advanced service |
| Table S12 | `d6_border_costs.csv`, `d6_cross_border_pairs.csv`, `d6_border_costs_excl_riyadh_jeddah.csv` | D6. Cluster borders |
| Table S13 | `d7_department_additions.csv`, `q3_service_audit_sensitivity.csv` | D7. Hospital departments; 7b. Sensitivity: hospital service audit |
| Table S14 | `sensitivity_slope_penalty.csv` | 11. Sensitivity: terrain-gradient speed penalty |
| Table S15 | `sensitivity_slope_by_cluster.csv` | 11. Sensitivity: terrain-gradient speed penalty |
| Table S16 | `table_S2_threshold_coverage_total_residents.csv` | 11. Supplement sensitivity tables |
| Table S17 | `cluster_travel_times.csv` | 10d. Per-cluster travel-time distributions |
| Figure S1 | `map_regions_governorates.png` | 10c. Reference maps |
| Figure S2 | `map_clusters_governorates.png` | 10c. Reference maps |
| Figure S3 | `fig_cluster_scorecard.png` | 10b. Cluster scorecard |
| Figure S4 | `fig_cluster_travel_times.png` | 10d. Per-cluster travel-time distributions |
| Figure S5 | `map_phc_access.png` | 10. Tables and figures |
| Figure S6 | `map_cluster_phc24.png` | 10d. Per-cluster travel-time distributions |
| Figure S7 | `map_cluster_gen.png` | 10d. Per-cluster travel-time distributions |
| Figure S8 | `map_hospital_access.png` | 10. Tables and figures |
| Figure S9 | `map_tertiary_access.png` | 10. Tables and figures |

In the code, "tertiary" denotes the referral level: the medical cities and specialized hospitals.

## Restricted inputs

Three files compiled for this study are not openly distributed: the facility census
(`list_of_healthcare_providers.csv`: 2,169 PHCs and 311 hospital records on 298 sites, with coordinates,
cluster, scope, and the hospital service flags), the PHC opening hours (`phc_opening_hours.csv`, 2,169 rows,
with the directory text for each PHC open beyond regular hours), and the hospital service audit
(`hospital_services_evidence.csv`, 311 records, with sources and quoted evidence). They are available from the
corresponding author on reasonable request for non-commercial academic use; the restricted-access record is
<https://doi.org/10.5281/zenodo.22358580>.

## Public inputs (download from source)

| Input | Source | Used for |
|---|---|---|
| WorldPop constrained population 2020, Saudi Arabia (`sau_ppp_2020_constrained.tif`) | <https://www.worldpop.org> | Population grid, aggregated to ~1 km and scaled to GASTAT governorate totals of Saudi citizens |
| GASTAT 2022 census: population by nationality, and by detailed age, gender, and nationality, by region and governorate | <https://portal.saudicensus.sa> | Citizen denominators; women 15 to 49, children under 15, adults 18 and over; all-resident sensitivity |
| Governorate and region boundaries (`governorate/Governorate.gpkg`, `Regions/Regions.shp`) | GASTAT | Reporting geography |
| OpenStreetMap, GCC states extract (2 September 2026) | <https://download.geofabrik.de/asia/gcc-states-latest.osm.pbf> | Routing graph |
| MoH Statistical Yearbook 2024 (`Statistical-Yearbook-2024 1.xlsx`) | <https://www.moh.gov.sa> | Facility validation (table 2-4); physicians posted to PHCs |
| SRTM 30 arc-second elevation (via the `geodata` package) | downloaded by the script | Terrain sensitivity only |

## How to reproduce

1. Make a project folder with this layout: `R/analysis.R`, `data/raw/`, and `data/processed/`. Put
   `referral_services.csv`, the public inputs, and the restricted inputs in `data/raw/`.
2. Download the OpenStreetMap extract and set `PAR$pbf_path` in section 0 of `analysis.R` to its location. The
   routing graph and the travel-time matrices are cached as `.rds` files in the same folder.
3. Restore the R environment with `renv::restore()` (or install the packages in `renv.lock`), then run
   `Rscript R/analysis.R` from the project folder. Outputs go to `output/tables/` and `output/figures/`. With the
   caches in place the run takes about 12 minutes; from scratch, several hours, dominated by routing.
4. Run `SLOPE_SENS=1 Rscript R/analysis.R` once for the terrain sensitivity (Supplementary Tables S14, S15, and
   S8I and S8J); it routes a second travel-time matrix on the slope-penalized network and caches it.

Routing runs entirely in R (dodgr); no routing server is needed. Travel times are computed on the full national
road network, and cluster membership is applied afterward as an attribute, never as a routing constraint.

**Unpaved tracks.** The routing graph runs tracks at 40 km/h, the Saudi track speed in Weiss et al. 2020
(Nature Medicine, Supplementary Table 1). `PAR$track_speed` makes tracks passable in the dodgr profile; for a car
profile, dodgr 0.4.3 gives roads without a speed tag the median tagged speed of their road class and uses the
profile value only when no road of the class is tagged, so the graph is the same whether the value is 30 (as when
the cached graph was built on 2 September 2026) or 40 (as now). In the built graph, 99.9% of the 336,680 km of
tracks run at 40 km/h.

## License and citation

- Code (`analysis.R`): MIT license.
- Derived tables and `referral_services.csv` (Zenodo deposit): CC-BY 4.0; please cite the manuscript when reusing.
- Restricted inputs: available on request, as above.
- OpenStreetMap data © OpenStreetMap contributors, ODbL. WorldPop data CC-BY 4.0. GASTAT and MoH statistics remain
  subject to their publishers' terms.

## Contact

Corresponding author: Ahmed Batobara, Mel and Enid Zuckerman College of Public Health, University of Arizona.
Email: abatobara@arizona.edu · ORCID: [0009-0005-7607-1213](https://orcid.org/0009-0005-7607-1213)
