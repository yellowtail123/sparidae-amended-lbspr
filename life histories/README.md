# life histories/ — reference copies

**Nothing in this folder is read by the pipeline.** The code reads
[`../Sparid_LBSPR_LifeHistory_trinomial_Rready.csv`](../Sparid_LBSPR_LifeHistory_trinomial_Rready.csv)
in the project root. These are the human-facing sources behind it, kept so the parameter choices can
be audited rather than taken on trust.

| File | What it is |
|---|---|
| `biometrics.csv` | The working parameter sheet. One row per species, carrying the values the pipeline uses alongside the columns it does not read: the source of each maturity ogive and each growth and mortality estimate, sex ratios, spawning season, fecundity, and a free-text confidence note |
| `Sparid_LBSPR_Params.xlsx` | The superseded Excel workbook the CSV was derived from. Retained for provenance only |

## How the values were chosen

Each parameter was resolved to the geographically nearest defensible published source, working
outward from the Strait of Gibraltar. Published total lengths were converted to fork length by a
single factor of 0.92 so that parameters and length observations share one basis.

Two flag columns matter when reading any result. `M_estimated` marks natural mortality values
inferred from the mortality-to-growth invariant rather than reported by a source; `K_estimated`
does the same for growth. Six of the nineteen species carry an inferred mortality and one an
inferred growth coefficient, so those rows are assumptions, not measurements.

The `confidence_note` column records what is doubtful about a row, including cases where the
adopted value sits at the edge of what the source supports. Read it before quoting a species.

## A caution on the sex-change ogive

`LD50_sexchange_cm_TL` is populated from the literature. `LD95_sexchange_cm_TL` is **not measured for
any species here**: it is set at 1.10 × LD50 throughout. The ogive width is therefore an assumption
applied uniformly across all twelve functional hermaphrodites, and the pipeline sweeps it as a
sensitivity axis for that reason.
