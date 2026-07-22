# Sparidae length-based stock assessment pipeline

*Length-based stock status of exploited seabreams (Sparidae) in British Gibraltar Territorial Waters: a single-season assessment of catch composition, community diversity and spawning potential.*

The work presented here sets out to address the data-poverty of British Gibraltar Territorial Waters,
where the exploited seabreams (family Sparidae) have gone almost entirely unassessed, by building a
length-based spawning-potential-ratio (LB-SPR) assessment amended to represent spawning potential more
accurately in sequentially hermaphroditic fishes. Assembled for a single-season master's thesis in
Gibraltar, the pipeline is written to run on any comparable length dataset. Given a season of length
records (with weight, sex and location where those were recorded), it works through catch composition
and community diversity, length-based indicators, the **LBSPR** fit itself, the
**sequential-hermaphroditism amendment** to that ratio, a length-converted catch-curve cross-check on
mortality, diel (day/night) and length–weight sub-analyses, a tag–recapture summary, and a simulation
layer that tests the amended estimator against a known operating model and stress-tests it under four
kinds of misspecification.

The amendment is the central contribution, because the standard method quietly gets sparids wrong.
Standard LBSPR assumes a gonochoristic stock, one whose sex ratio holds fixed at about half female
across all lengths, so the female fraction is a constant that cancels out of the spawning-potential
calculation. That assumption fails for the protandrous and protogynous species that change sex, where
the proportion female shifts with length, so the amendment recomputes spawning potential with a
proportion-female-at-length weighting, applied **as a post-processing reweight of egg output**. The
underlying LBSPR fit (F/M, the selectivity lengths SL50/SL95, the package's own SPR) is **left alone**;
switching the sex layer off, or applying it to a gonochore, reproduces standard LBSPR exactly. The pipeline is modular and
threshold-gated, so it runs on sparse single-season data and grows as more accumulates, and it carries
over to other data-poor sparid or hermaphroditic fisheries given a life-history table and catch records
in the same format.

---

## 1. What the pipeline does

| Stage | What it computes |
|-------|------------------|
| Data assembly | Sorts out species names and the fork-length / total-length scale (FL = TL × 0.92), merges the season catch with the historical tagging series, and writes one tidy dataset. |
| Catch composition & diversity | Species composition, rank–abundance, Hill numbers, rarefaction / sample coverage (via **vegan** and **iNEXT**). |
| Length-based indicators | Standard length-based reference indicators per species (subject to a minimum-sample-size gate). |
| Mortality cross-check | Empirical natural-mortality envelope and a length-converted catch curve (both via **TropFishR**) as an independent check on the assessment. |
| LBSPR assessment | Spawning potential ratio from length composition via the **LBSPR** package (Hordyk), with life-history uncertainty propagated by a length bootstrap and a parameter Monte-Carlo draw. |
| Sequential-hermaphroditism amendment | A sex-structured re-weighting of egg output layered on top of the LBSPR fit. |
| Diel & length–weight | Day/night composition and length–weight relationships (each auto-activates only when the data clears its threshold). |
| Tagging | Capture / recapture summary (movement/growth analysis is present but intentionally commented out). |
| Operating-model validation | A recovery simulation: simulate a known stock, sample lengths, re-estimate through the **deployed** LBSPR + amendment path, and measure the estimator's error. |
| Misspecification stress battery | Four stressors, namely dome-shaped selectivity, compensatory (plastic) sex change, biased life-history and non-equilibrium recruitment, each of which breaks one LBSPR assumption on the *truth* side while the estimator is left naive. |

---

## 2. Methods provenance and package credits

The statistics are not reimplemented from scratch. The heavy lifting is done by established,
peer-reviewed packages, and the code here is a thin sex-structured layer sitting on top of them. Credit
where it is due:

- **LBSPR** (Hordyk et al.), the core length-based spawning-potential-ratio fit.
- **TropFishR**, length-frequency handling, the length-converted catch curve, and the empirical
  natural-mortality (M) envelope (`M_empirical`).
- **fishmethods**, an alternative empirical natural-mortality estimator.
- **rfishbase** and **FishLife** (Thorson), life-history cross-checks and the multivariate
  parameter-correlation structure used in the parameter Monte-Carlo.
- **vegan** and **iNEXT**, community diversity, Hill numbers, rarefaction and sample coverage.
- Supporting: **tidyverse**, **here**, **janitor**, **skimr**, **sf**, **leaflet**, **geosphere**,
  **ggrepel**, **scales**, **htmlwidgets**, **magick**, **jsonlite**.

The sequential-hermaphroditism amendment is a post-processing reweight and doesn't touch any
package's estimator.

---

## 3. Repository layout

```
sparidae-lbspr/
├── README.md                     this file
├── LICENSE                       MIT
├── CITATION.cff                  how to cite
├── .gitignore                    code-only: excludes data and generated outputs (see §4)
├── renv.lock                     pinned package versions (renv)
├── renv/ , .Rprofile             renv activation machinery
├── sparidae-lbspr.Rproj          RStudio project anchor (defines the here() root)
├── Sparid_LBSPR_LifeHistory_trinomial_Rready.csv   authoritative life-history inputs
├── R/
│   ├── 00_run_all.R              master runner (sources 01→04 in order)
│   ├── 01_combine_data.R         data assembly  -> combined_records/combined_dataset.csv
│   ├── 02_analysis.R             assessment engine -> results/ + graphs and maps/
│   ├── 03_stress_test.R          misspecification battery (+ figures); uses 02's functions
│   └── 04_catalog.R              photographic specimen catalogue -> catalog/
├── sex-structured-lbspr/  self-contained extract of the sex-change amendment (README + code + runnable demo) for reuse on its own
├── all_records/         raw season records (bream_*.csv) + photos/   [not shipped, supply your own]
├── historical_tagging/  historical tagging series (*.csv)            [not shipped, supply your own]
├── combined_records/    generated by R/01                            [generated by a run]
├── results/             generated tables (*.csv)                      [generated by a run]
├── graphs and maps/     generated figures (*.png) + leaflet map       [generated by a run]
└── catalog/             generated fish_catalog.html                   [generated by a run]
```

The layout is deliberately flatter than a proper R package. Every data path is resolved with `here()`
relative to the `.Rproj` at the repo root, so the scripts sit under `R/` without any edits, but the data
and output folders have to stay at the root where `here()` expects them.

---

## 4. Expected input data and the folder contract

The code expects this structure at the repo root:

| Path | Contents |
|------|----------|
| `all_records/bream_*.csv` | the season's catch records (one or more `bream_…` CSVs). |
| `all_records/photos/` | specimen photographs, consumed by `R/04_catalog.R`. |
| `historical_tagging/*.csv` | the historical tag/recapture series merged by `R/01_combine_data.R`. |
| `Sparid_LBSPR_LifeHistory_trinomial_Rready.csv` | one-row-per-species life-history table (growth, maturity, M/K, CVLinf, length–weight a/b, sex system, length-at-sex-change). Root level. |

---

## 5. Installation

- You'll want **R** ≥ 4.5 (developed on 4.5.2).
- Package versions are pinned with **renv**, so from the project root one command restores the exact
  library:

  ```r
  renv::restore()
  ```

- **FishLife** is GitHub-only (not on CRAN). `renv::restore()` installs it from the source recorded in
  `renv.lock`. If you're setting the environment up by hand instead, install it with:

  ```r
  install.packages("pak")
  pak::pkg_install("james-thorson/FishLife")
  ```

`renv::restore()` also builds a project-local library, and `renv/` and `.Rprofile` activate it
automatically when you open the project.

---

## 6. How to run

Run the whole script pipeline from the project root in one session:

```r
source(here::here("R", "00_run_all.R"))
```

(or run `R/01_…` through `R/04_…` one at a time, **in order**, `03` uses functions defined in `02`, so
it has to follow it in the same session). That writes the combined dataset, every assessment table and
figure, the photographic catalogue, and the **stress-battery caches** (`results/om_stress_*.csv`).

### The slow stages, and the caching that spares you
- The operating-model **recovery simulation** and the four-stressor **misspecification battery** (dome
  selectivity, compensatory sex change, life-history bias, recruitment variability) are the slow parts.
  They are `RUN_OM_*`-flag-gated and **cache their outputs to `results/`**, read back by the figure and
  table code, so once the caches exist you never pay the simulation cost again until you delete them.
- `FINAL_RUN` sets the Monte-Carlo depth. The slow step is the stress battery at `FINAL_RUN = TRUE`:
  1000 replicates per severity level, each fitting a full `LBSPRfit`. Give it time, it can run a while;
  drop `FINAL_RUN` to `FALSE` for fast iteration.

---

## 7. Outputs

Everything below is **regenerated by a run and git-ignored** (not tracked):

- `results/*.csv`, the assessment tables: the per-species assessment estimates (`assessment_results.csv`)
  and the pooled 2026+historical version, the LBSPR bootstrap and parameter-MC CIs, the sex-structured
  SPR, the LD50/LD95 sex-change-length sensitivity, length-based indicators, diversity / Hill numbers /
  sample coverage, the sourced life-history table and its FishBase audit, the catch-curve and mortality
  cross-checks, day/night composition, the tagging summary, and the operating-model / stress-battery
  caches.
- `graphs and maps/*.png`, the figures (length-frequency, SPR status, size-vs-reference, rank-abundance,
  day/night, tagged fish, fishing pressure, and the rest) plus a `leaflet.html` map.
- `catalog/fish_catalog.html`, the photographic specimen catalogue.

A few sections (diel, length–weight, tag movement) only produce output once the data clears their
thresholds or when recaptures exist, so if they're missing on a given dataset that's expected, not a bug.

---

## 8. Citation

Cite this pipeline using the metadata in [`CITATION.cff`](CITATION.cff).

---

## 9. License

Released under the **MIT License** (see [`LICENSE`](LICENSE)).

---

### Reproducibility notes
- The Monte-Carlo stages are seeded, so given the same input data and package versions the seeded tables
  reproduce exactly. The one exception is the FishBase audit, which is fetched live over the network
  (via `rfishbase`) and can vary with the remote database.
