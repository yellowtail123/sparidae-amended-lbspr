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
| `sex_system` | `protandry`, `protogyny`, or anything else for a gonochore |
| `LD50_sexchange_cm_TL` | length at 50% sex change |
| `LD95_sexchange_cm_TL` | length at 95%, defaults to 1.10 × LD50 |

**No per-fish sex data is needed.** Leave `LD50` blank and it falls back to standard LB-SPR instead
of failing.

Report **`SPR_bind`**: the female egg-based SPR, or under protogyny the lower of that and the male
ratio.

## The full pipeline

`R/00_run_all.R` runs the whole analysis: catch composition, length indicators, mortality
cross-checks, the LBSPR fit with bootstrap intervals, the amendment, diel and length–weight,
tagging, and a photo catalogue. `R/03_stress_test.R` tests it against a known operating model under
four kinds of misspecification.

```r
renv::restore()
source(here::here("R", "00_run_all.R"))
```

R ≥ 4.5. You supply the catch records, as `all_records/bream_*.csv`; they do not ship here. A
historical tagging series in `historical_tagging/*.csv` is optional and that step is skipped if the
folder is empty. The life-history table does ship and shows the expected format. Lengths must be
**fork length in cm**, species names must match that table's trinomials, and a species needs 20
measured lengths before it is assessed.

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
│   ├── 00_run_all.R                                master runner
│   ├── 01_combine_data.R                           data assembly
│   ├── 02_analysis.R                               assessment engine
│   ├── 03_stress_test.R                            misspecification battery
│   └── 04_catalog.R                                photographic catalogue
└── sex-structured-lbspr/                           the amendment, standalone
    ├── README.md                                   method write-up with the mathematics
    ├── sex_structured_lbspr.R                      the functions
    ├── example.R                                   runnable demonstration
    └── example_life_history.csv                    three species' parameters
```

Output goes to `results/`, `graphs and maps/` and `catalog/`, all regenerated and git-ignored.

## Caveats

Standard LB-SPR assumptions carry over: equilibrium, asymptotic selectivity. Life-history values are
literature-derived with a single FL:TL ratio across species, so the intervals are wide. The male
floor is a precautionary bound, not an estimate of fertilisation success.

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
