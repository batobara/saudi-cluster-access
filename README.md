# Travel-time access to three tiers of care under Saudi Arabia's health-cluster reform

Reproducibility deposit for the manuscript:

> *Travel-time access to three tiers of care under Saudi Arabia's health-cluster
> reform: a national ecological study.* Submitted to BMJ Global Health.

## Contents

| Path | Description |
|---|---|
| `analysis.R` | The complete analysis pipeline in a single R script. Reproduces every main-text and supplementary table and figure from the inputs below. |
| `table_S1_…` – `table_S12_….csv` | Supplementary tables as submitted (CSV). |
| `figure_S1_…` – `figure_S9_….png` | Supplementary figures as submitted (PNG). |
| `supplementary_material.md` / `.docx` | Captions and notes for all supplementary tables and figures. |

The GitHub repository carries `analysis.R` and this README; the supplementary
tables, figures, and caption files listed above are archived together with the
code in the project's Zenodo deposit ([OSF/Zenodo DOI TBD]).

The facility census is **not distributed** in either location (see Facility
dataset availability below).

## Facility dataset availability

The facility census (`list_of_healthcare_providers.csv`) is not distributed in
the GitHub repository or the Zenodo deposit. It comprises 2,480 rows covering 2,467 health-cluster facility
sites (2,169 primary health care centers, 282 secondary hospitals, 16 tertiary
sites; directory entries sharing one hospital campus are counted as one site
each: ten centers of excellence at King Fahad Hospital in Al Baha, two centers
at King Salman Specialist Hospital in Hail, one center at King Fahd Specialist
Hospital in Buraidah, and one outpatient center at Al Makhwah General Hospital)
with coordinates, cluster assignment, tier, specialty scope, working-hours
classification, and the source URL for each geocode. It was validated against
MoH Statistical Yearbook 2024 table 2-4 (2,169 vs 2,172 PHCs, per-cluster
differences of six or fewer; 298 vs 290 hospital sites, per-cluster differences
of two or fewer; supplementary Table S8).

The dataset was compiled and verified facility by facility by the authors over
several months from public MoH directories. It is **available from the
corresponding author on reasonable request** for non-commercial academic use,
with citation of the manuscript. All derived tables (including every per-cluster
and per-governorate result) are openly included in the Zenodo deposit.

## How to reproduce

1. Request the facility census `list_of_healthcare_providers.csv` from the
   corresponding author (see Contact) and place it under `data/raw/` in the
   project root (the script expects `data/raw/` and `data/processed/`
   directories).
2. Download the public inputs listed below into `data/raw/`.
3. Download the Geofabrik Gulf road extract
   (<https://download.geofabrik.de/asia/gcc-states-latest.osm.pbf>, ~250 MB) and set
   `PAR$pbf_path` in section 0 of `analysis.R` to its location.
4. Run `analysis.R` from the project root. Expensive steps (population grid,
   routing graph, travel-time matrices) cache to `data/processed/` as `.rds`;
   delete a cache file to force recomputation. The full run from scratch takes
   several hours, dominated by the dodgr routing step.

Routing runs entirely in R (dodgr street-network routing); no external routing
server is required. Travel times are computed on the full national road network;
cluster membership is applied afterward as an attribute, never as a routing
constraint.

## Public inputs (not redistributed here; download from source)

| Input | Source | Used for |
|---|---|---|
| WorldPop constrained population 2020, Saudi Arabia, ~100 m (`sau_ppp_2020_constrained.tif`) | <https://www.worldpop.org> | Population grid, aggregated to ~1 km and scaled to GASTAT governorate totals |
| GASTAT 2022 census tables: population by nationality, gender, and detailed age, by region and governorate | <https://portal.saudicensus.sa> | Saudi-citizen denominators; demand-matched subgroups (women 15–49, children under 15, adults 18+); total-resident sensitivity denominator |
| OpenStreetMap road network, Gulf (GCC states) extract | <https://download.geofabrik.de/asia/gcc-states-latest.osm.pbf> (accessed 2 September 2026) | dodgr routing graph |
| MoH Statistical Yearbook 2024 | <https://www.moh.gov.sa> | External validation of the facility census (table 2-4); workforce indicators (tables 2-24, 2-45, 2-43, 2-29, 2-20) |
| Administrative boundaries (13 regions, 150 governorates) | GASTAT | Reporting geography |

## Software versions

Analysis run under **R 4.6.0** on macOS with the following packages:

sf 1.1-1, terra 1.9-27, dplyr 1.2.1, tidyr 1.3.2, readr 2.2.0, stringr 1.6.0,
dodgr 0.4.3, osmextract 0.6.0, units 1.0-1, ineq 0.2-13, spdep 1.4-2,
RANN 2.6.3, curl 7.1.0, ggplot2 4.0.3, viridis 0.6.5, ggrepel 0.9.8,
patchwork 1.3.2, readxl 1.5.0.

OSM extract date: 2 September 2026 (Geofabrik GCC states).

## Notes on the facility census

- `hours_level` codes the working-hours tier: 0 = regular hours (1,777 PHCs),
  1 = extended hours (305), 2 = 24-hour (87 PHCs; all hospitals operate 24/7).
- `type` codes the tier: PHC, 2ry (secondary), 3ry (tertiary). Tier is
  analyst-assigned: the source directories do not label tiers. The scheme
  follows two official anchors: the Health Holding Company describes each
  cluster as comprising primary care centers, hospitals, and medical
  cities / specialized hospitals (health.sa/en/clusters), and the reform's
  model of care organizes cluster services as integrated levels of care,
  with primary care as first contact and referral upward to secondary and
  tertiary services; tier is thus a level of care, defined by function
  rather than premises. Medical cities and specialized hospitals are a
  distinct MoH category under the 2014 Law of Medical Cities and
  Specialized Hospitals. Tertiary =
  the medical-city / specialized-hospital stratum acting as the cluster
  referral terminus with advanced specialized services (open-heart
  surgery, transplantation, radiation or comprehensive oncology),
  verified facility by facility; centers of excellence delivering such
  services from a host-hospital campus are coded 3ry, while specialty
  hospitals and centers without such services are coded 2ry even when
  named "specialist" or "center of excellence" (the `notes` column
  records the rationale for these entries).
- Coordinates were geocoded from MoH directories and verified facility by
  facility; `source_url` records the map link used for each.
- The study covers health-cluster facilities only. Other government
  providers (National Guard, military, university hospitals) serve specific
  populations outside the cluster system and define the scope of the study.

## License and citation

- Code (`analysis.R`): MIT license.
- Derived tables and figures (`table_S*.csv`, `figure_S*.png`, in the Zenodo
  deposit): CC-BY 4.0; please cite the manuscript when reusing.
- Facility dataset: not openly distributed; available from the
  corresponding author on reasonable request for non-commercial academic use,
  with citation of the manuscript (see Facility dataset availability above).
- OpenStreetMap data © OpenStreetMap contributors, ODbL. WorldPop data CC-BY 4.0.
  GASTAT and MoH statistics remain subject to their publishers' terms.

