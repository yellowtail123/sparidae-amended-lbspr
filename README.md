# Sparidae length-based stock assessment pipeline

LB-SPR stock assessment for seabreams (Sparidae), amended to handle fish that change sex.

Standard LB-SPR assumes half the fish at any given length are female. Most sparids are protandrous
or protogynous, so the proportion female shifts with length and the standard calculation gets
spawning potential wrong. This adds a sex-at-length correction.

Built for a single-season assessment of British Gibraltar Territorial Waters. The code runs on any
comparable length data.

## Try the amendment

Self-contained in `sex-structured-lbspr/`, no data required:

```sh
cd sex-structured-lbspr
Rscript example.R
```

It runs on three real sparids, one per sex system, and needs only **tibble** and **dplyr**. Run it
from inside that folder so the project's renv setup stays out of the way. The mathematics is in
[`sex-structured-lbspr/README.md`](sex-structured-lbspr/README.md).

## How it works

Egg output is reweighted by a proportion-female-at-length ogive, applied after the fit rather than
inside it. The LB-SPR fit is untouched, so F/M, selectivity and the package's own SPR come through
unchanged. Apply it to a gonochore and you get standard LB-SPR back exactly.

For protogynous species it also computes a male capacity ratio and reports the lower of the two.
Fishing the large end of a protogynous stock removes males, which standard SPR cannot see.

## Using it on your own species

Three fields, in the life-history table rather than the catch records:

| Field | Meaning |
|---|---|
| `sex_system` | `protandry`, `protogyny`, `gonochore` or `rudimentary` |
| `LD50_sexchange_cm_TL` | length at 50% sex change |
| `LD95_sexchange_cm_TL` | length at 95%, defaults to 1.10 × LD50 |

**No per-fish sex data is needed.** Leave `LD50` blank and it falls back to standard LB-SPR instead
of failing. Anything outside those four `sex_system` values stops the run rather than guessing.

Report **`SPR_bind`**: the female egg-based SPR, or under protogyny the lower of that and the male
ratio.

## Running the pipeline

**1. Install.** R ≥ 4.5; `renv.lock` pins R 4.5.2 and 193 packages.

```r
renv::restore()
install.packages("patchwork")                       # see below
pak::pkg_install("james-thorson/FishLife")          # GitHub-only; renv will not fetch it
```

Two packages need a hand. **FishLife** is GitHub-only, and skipping it costs you only the correlated
parameter draws, which fall back to independent ones with a message. **patchwork** is listed as a
hard requirement by the preflight but is not in `renv.lock`, so `renv::restore()` alone leaves the
preflight failing on it. Installing it is the quick fix. It is only genuinely needed by `05`, which
warns and writes its panels separately without it, so the alternative is to move `"patchwork"` from
`hard` to `soft` in the package block of `00_run_all.R`.

If you install by hand rather than through renv, note that `ggrepel` is required by `02` and is not
part of the tidyverse metapackage, nor is it named in the preflight list. `renv::restore()` picks it
up; a manual install from the preflight's error message will not.

**2. Put your data in place.** Neither folder ships here and neither is optional.

```
all_records/bream_*.csv        catch records, one file per export; all are read and pooled
historical_tagging/*.csv       tagging series; every CSV in the folder is read
```

Lengths must be **fork length in cm** in `length_true_cm`, with `length_type` set to `FL`. Species
names must match the trinomials in `Sparid_LBSPR_LifeHistory_trinomial_Rready.csv`; a name that does
not match gets no life-history row and is quietly not assessed, so check the readiness table if a
species you expected is missing. A species needs **20 measured lengths** to be assessed, 10 for the
length-based indicators.

**3. Run.**

```r
source(here::here("R", "00_run_all.R"))
```

A preflight checks the folders, the life-history table and every required package before anything
runs, so it cannot stop forty minutes in on something it could have caught at the start.

**Order is not optional.** `02` reads `combined_records/combined_dataset.csv`, which `01` writes, and
uses it to gate the season analysis; run `02` on its own and the season *n* comes out one higher than
the reported figure, with a warning. `03` and `05` both call functions defined in `02` and must
follow it **in the same session**. `04` is standalone.

| Script | Does | Needs |
|---|---|---|
| `01_combine_data.R` | Pools the season exports and the tagging series into one length-reconciled file | both data folders |
| `02_analysis.R` | Catch composition, length indicators, mortality envelope, the LBSPR fit with bootstrap and parameter Monte-Carlo intervals, the amendment, catch-curve cross-check, diel, length–weight, tagging, operating model | `01` |
| `03_stress_test.R` | Misspecification battery: dome selectivity, compensatory sex change, biased life history, recruitment variability | `02`, same session |
| `04_catalog.R` | Photographic catalogue, one self-contained HTML file | magick, jsonlite |
| `05_manuscript_outputs.R` | Combined-dataset assessment, capture/recapture, closed-loop recovery, ogive sensitivity | `02`, same session |

**4. Collect.** `results/` holds the CSV tables, `graphs and maps/` the figures and the interactive
map, plus `combined_records/` and `catalog/`. All four are regenerated and git-ignored. Every run
writes `results/sessionInfo.txt`, so a figure can be traced back to the session that made it.

### Expect hours, not minutes

Two thousand bootstrap replicates and two thousand parameter Monte-Carlo draws per assessed species,
a thousand-replicate recovery simulation, and a thousand-replicate stress battery at every severity
level of four stressors. Turn it down while you are setting up:

| Switch | Where | Effect |
|---|---|---|
| `FINAL_RUN <- FALSE` | `02_analysis.R` | Halves the bootstrap and Monte-Carlo draws to 1000 |
| `RUN_OM_SIM <- FALSE` | `02_analysis.R` | Skips the recovery simulation; the figure rebuilds from cache |
| `RUN_OM_* <- FALSE` | `03_stress_test.R` | Skips a stressor; same caching |
| `RUN_CATALOG <- FALSE` | `00_run_all.R` | Skips `04`, the photo re-encoding |

The simulation stages cache to `results/om_*.csv` and the figures are rebuilt from whatever is on
disk, so **delete the stale CSV before re-running a stressor with changed settings**; the figures
will not tell you they came from the previous run. Keep `FINAL_RUN` the same in `02` and in
`analysisfinal.Rmd`: both write the same files into `results/`, so if they disagree, the intervals
you report depend on which ran last.

### When it stops

| Message | Cause |
|---|---|
| `Preflight failed. Fix these before running:` | Listed underneath: a missing folder, the life-history table, or a hard package. `patchwork` is the usual first one; see Install |
| `No files matching '^bream_.*\.csv$' found in .../all_records` | Empty or misnamed data folder; the same error names `historical_tagging` |
| `05_manuscript_outputs.R needs 02_analysis.R sourced first in this session` | `05` run on its own; source `02` first |
| `sex must be gonochore, rudimentary, protandry, or protogyny` | A `sex_system` value the ogive cannot read |
| `This script needs the 'magick' and 'jsonlite' packages` | `04` without ImageMagick behind magick |
| `[skip] TropFishR not installed` | Mortality envelope and catch curve skipped; `MK_cv` falls back to the global default |
| `[skip] rfishbase not installed` | FishBase cross-check skipped; nothing else changes |

Optional packages skip their own section and say so. Hard ones stop the preflight.

## Layout

```
sparidae-amended-lbspr/
├── README.md                                       this file
├── LICENSE                                         MIT
├── CITATION.cff                                    how to cite
├── .gitignore                                      code-only policy
├── .Rprofile                                       activates renv on open
├── renv.lock                                       pinned package versions (R 4.5.2)
├── renv/
│   ├── activate.R                                  renv bootstrap
│   └── .gitignore
├── sparidae-lbspr.Rproj                            project anchor, defines the here() root
├── Sparid_LBSPR_LifeHistory_trinomial_Rready.csv   life-history inputs
├── life histories/                                 reference copies, not read by the code
│   ├── Sparid_LBSPR_Params.xlsx                    superseded Excel workbook
│   └── biometrics.csv                              human-facing parameter sheet
├── R/
│   ├── 00_run_all.R                                master runner, with preflight
│   ├── 01_combine_data.R                           data assembly
│   ├── 02_analysis.R                               assessment engine
│   ├── 03_stress_test.R                            misspecification battery
│   ├── 04_catalog.R                                photographic catalogue
│   └── 05_manuscript_outputs.R                     combined assessment, recovery, ogive sensitivity
└── sex-structured-lbspr/                           the amendment, standalone
    ├── README.md                                   method write-up with the mathematics
    ├── sex_structured_lbspr.R                      the functions
    ├── example.R                                   runnable demonstration
    └── example_life_history.csv                    three species' parameters
```

No field data is committed here, and the `.gitignore` is written to keep it that way.

## Caveats

Standard LB-SPR assumptions carry over: equilibrium, asymptotic selectivity. Life-history values are
literature-derived with a single FL:TL ratio across species, so the intervals are wide. The male
floor is a precautionary bound, not an estimate of fertilisation success. Two intervals are reported
per species; the parameter Monte-Carlo is the wider and the more honest at these sample sizes, and
is the one to read as the working uncertainty.

## Credits

A thin sex-structured layer over established packages, not a reimplementation:

| Package | Author | Used for |
|---|---|---|
| [LBSPR](https://github.com/AdrianHordyk/LBSPR) | Hordyk | the LB-SPR fit itself |
| [TropFishR](https://github.com/tokami/TropFishR) | Mildenberger et al. | catch curve, natural-mortality envelope |
| [FishLife](https://github.com/James-Thorson-NOAA/FishLife) | Thorson | life-history priors |
| [rfishbase](https://github.com/ropensci/rfishbase) | Boettiger et al. | FishBase cross-checks |

Everything else is declared at the top of each script and pinned in `renv.lock`.

Cite this pipeline using [`CITATION.cff`](CITATION.cff). Released under the MIT License (see
[`LICENSE`](LICENSE)).

## References

- Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N. & Prince, J. (2015a) Some explorations of the
  life history ratios to describe length composition, spawning-per-recruit, and the spawning
  potential ratio. *ICES Journal of Marine Science* 72(1): 204–216.
  [doi](https://doi.org/10.1093/icesjms/fst235)
- Hordyk, A., Ono, K., Valencia, S., Loneragan, N. & Prince, J. (2015b) A novel length-based
  empirical estimation method of spawning potential ratio (SPR). *ICES Journal of Marine Science*
  72(1): 217–231. [doi](https://doi.org/10.1093/icesjms/fsu004)
- Hordyk, A.R., Ono, K., Prince, J.D. & Walters, C.J. (2016) A simple length-structured model based
  on life history ratios and incorporating size-dependent selectivity. *Canadian Journal of
  Fisheries and Aquatic Sciences* 73(12): 1787–1799.
  [doi](https://doi.org/10.1139/cjfas-2015-0422)
- Mildenberger, T.K., Taylor, M.H. & Wolff, M. (2017) TropFishR: an R package for fisheries analysis
  with length-frequency data. *Methods in Ecology and Evolution* 8(11): 1520–1527.
  [doi](https://doi.org/10.1111/2041-210X.12791)
- Thorson, J.T., Munch, S.B., Cope, J.M. & Gao, J. (2017) Predicting life history parameters for all
  fishes worldwide. *Ecological Applications* 27(8): 2262–2276.
  [doi](https://doi.org/10.1002/eap.1606)
- Boettiger, C., Lang, D.T. & Wainwright, P.C. (2012) rfishbase: exploring, manipulating and
  visualizing FishBase data from R. *Journal of Fish Biology* 81(6): 2030–2039.
  [doi](https://doi.org/10.1111/j.1095-8649.2012.03464.x)
- Goodyear, C.P. (1993) Spawning stock biomass per recruit in fisheries management: foundation and
  current use. *Canadian Special Publication of Fisheries and Aquatic Sciences* 120: 67–81.
