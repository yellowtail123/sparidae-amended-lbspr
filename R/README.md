# R/ — the pipeline

Five scripts, run in order. Full usage, data requirements and troubleshooting are in the
[project README](../README.md); this file is the map.

| Script | Does | Needs |
|---|---|---|
| `00_run_all.R` | Master runner. Preflights folders, life-history table and packages, then sources 01–05 in order and writes `results/sessionInfo.txt` | — |
| `01_combine_data.R` | Pools the season exports and the tagging series into one length-reconciled file, collapsing records logged twice by different programmes | both data folders |
| `02_analysis.R` | The assessment engine: catch composition, length indicators, mortality envelope, the LB-SPR fit with bootstrap and parameter Monte-Carlo intervals, the sex-structured amendment, catch-curve cross-check, diel, length–weight, tagging, and the operating model | `01` |
| `03_stress_test.R` | Misspecification battery: dome selectivity, compensatory sex change, biased life history, recruitment variability | `02`, same session |
| `04_catalog.R` | Photographic catalogue as one self-contained HTML file | magick, jsonlite |
| `05_manuscript_outputs.R` | Combined-dataset assessment, capture/recapture, closed-loop recovery, ogive sensitivity, panel figures | `02`, same session |

**Order is not optional.** `02` reads the file `01` writes and uses it to gate the season analysis.
`03` and `05` call functions defined in `02` and must follow it in the same R session. `04` is
standalone and can be skipped with `RUN_CATALOG <- FALSE`.

Run-control switches sit in a block at the top of `02_analysis.R` and `03_stress_test.R`. Turn
`FINAL_RUN` off while setting up; a full run takes hours.

The sex-structured amendment itself is not here. It is standalone, with its own write-up and a
runnable example, in [`../sex-structured-lbspr/`](../sex-structured-lbspr/).
