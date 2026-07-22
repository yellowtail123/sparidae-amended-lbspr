# Sex-structured LB-SPR, the sequential-hermaphroditism amendment

A small, self-contained implementation of a **sex-structured spawning potential ratio (SPR)**
for length-based stock assessment of **sequentially hermaphroditic** fishes (protandrous and
protogynous species, e.g. many seabreams, Sparidae).

Standard length-based SPR (LB-SPR) accumulates per-recruit egg production as
`maturity(L) × Lᵇ`, with no sex term: it assumes a gonochoristic stock whose sex ratio is fixed
(about half female) at every length, so the female fraction is a constant that cancels out. For
sequential hermaphrodites that assumption is wrong:

- **Protogyny** (female → male): the largest fish are *male*, so size-selective fishing removes
  much of the population's *fertilisation capacity* even while egg output looks adequate.
- **Protandry** (male → female): the largest fish *are* the egg-producers, so fishing hits egg
  output harder than a sex-blind model implies.

This amendment reweights egg output by a **proportion-female-at-length ogive** and adds a
**precautionary male-capacity floor**. It is applied as a **post-processing step**: the LB-SPR
fit itself (fishing mortality `F/M`, gear selectivity `SL50`/`SL95`, and the package SPR) is
**left completely unchanged**, because maturity and sex play no part in the LB-SPR likelihood.
With the sex layer disabled, or for a gonochoristic (separate-sex) species, it reduces
**exactly** to standard LB-SPR, so it is a *generalisation*, not a competing method.

> **Provenance.** These functions are copied **verbatim** from the thesis analysis pipeline
> (`../R/02_analysis.R`, section 10b), which remains the canonical source. Nothing in the
> thesis was moved or removed; this folder is a standalone extract so the method can be read and
> reused on its own. The typeset mathematics is given in §§1–2 below.

---

## 1. Standard LB-SPR (the baseline)

LB-SPR (Hordyk et al. 2015a, 2015b, 2016; R package `LBSPR`) exploits the result that, at
equilibrium and under von Bertalanffy growth, the standardised length composition of an
exploited stock is governed by two dimensionless ratios, the life-history ratio `M/K` and the
ratio of fishing to natural mortality `F/M`, together with the gear-selectivity ogive. Both
maturity and selectivity are logistic in length:

$$
S_L = \left[\,1 + \exp\!\left(-\ln(19)\,\frac{L - SL_{50}}{SL_{95} - SL_{50}}\right)\right]^{-1},
$$

with the maturity ogive $m(L)$ taking the same form with $L_{50}, L_{95}$. Fitting the expected
to the observed length composition yields `F/M` and the selectivity parameters, and the spawning
potential ratio is the ratio of expected lifetime egg production per recruit of the fished
stock to that of the unfished stock:

$$
\mathrm{SPR} = \frac{E_F}{E_0},
\qquad
E = \sum_{L} N(L)\,m(L)\,W(L),\quad W(L) = a\,L^{b},
$$

where $N(L)$ is survivorship-at-length (fished $N_F$ or unfished $N_0$) and $W(L)=aL^b$ is
mass-at-length. Status is read against the conventional reference points $\mathrm{SPR}=0.40$
(target) and $\mathrm{SPR}=0.20$ (limit) (Goodyear 1993), with $F/M = 1$ marking fishing
mortality equal to natural mortality.

## 2. The amendment (sex-structured SPR)

Write $\psi_f(L)$ for the **proportion female at length** $L$. Egg production is reweighted so
that only the female fraction contributes eggs:

$$
E_F = \sum_{L} N_F(L)\,m(L)\,\psi_f(L)\,W(L),
\qquad
E_0 = \sum_{L} N_0(L)\,m(L)\,\psi_f(L)\,W(L),
\qquad
\mathrm{SPR}_{f} = \frac{E_F}{E_0}.
$$

The proportion-female ogive is a logistic on the **functional sex-change** lengths
$L_{\Delta 50}, L_{\Delta 95}$ (`LD50`, `LD95`), on the same length scale as $L_\infty$:

$$
\psi_f(L)=
\begin{cases}
\tfrac{1}{2} & \text{gonochoristic or rudimentary,}\\[4pt]
\left[\,1+\exp\!\left(-\ln(19)\,\dfrac{L-L_{\Delta 50}}{L_{\Delta 95}-L_{\Delta 50}}\right)\right]^{-1} & \text{protandrous (ascending),}\\[10pt]
1-\left[\,1+\exp\!\left(-\ln(19)\,\dfrac{L-L_{\Delta 50}}{L_{\Delta 95}-L_{\Delta 50}}\right)\right]^{-1} & \text{protogynous (descending).}
\end{cases}
$$

- **Gonochores and rudimentary hermaphrodites** take a *constant* female fraction $\tfrac12$.
  Because a constant factor cancels from the ratio $E_F/E_0$, $\mathrm{SPR}_f$ **reduces exactly
  to the standard SPR**; the amendment changes nothing where there is no functional sex change.
- **Male-capacity floor.** In protogyny the large (male) fish are removed first, so egg-based
  SPR alone can be optimistic. A male analogue replaces the female fraction with its complement
  $1-\psi_f(L)$, giving a mature-male reproductive-capacity ratio $\mathrm{SPR}_{m}$. The status
  reported for protogynous species is the more precautionary of the two:

$$
\mathrm{SPR}_{\mathrm{bind}} =
\begin{cases}
\min\!\left(\mathrm{SPR}_{f},\,\mathrm{SPR}_{m}\right) & \text{protogyny,}\\[2pt]
\mathrm{SPR}_{f} & \text{otherwise.}
\end{cases}
$$

  The male floor is deliberately conservative: it *bounds the risk* of male limitation but does
  not estimate fertilisation success, because the mating function relating sex ratio to realised
  reproduction is unknown for these species (Heppell et al. 2006; Easter & White 2020).

- **Fit-invariance.** Maturity and sex never enter the LB-SPR likelihood, which depends only on
  the length data and the selectivity and mortality parameters. The amendment therefore takes
  the fitted $F/M$, $SL_{50}$, $SL_{95}$ **unchanged** and recomputes only the SPR, re-implementing
  the growth-type-group per-recruit calculation of the package for that purpose (following the
  established practice of building size-specific sex transformation into per-recruit analyses for
  hermaphroditic stocks, Shepherd 1993; Punt et al. 1993; Benvenuto et al. 2017).

## 3. Files

| File | What it is |
|---|---|
| `sex_structured_lbspr.R` | The method. `source()` it to get `spr_sex_structured()`, `fit_lbspr_one()` and `make_lf()`. Pure base R + `tibble`/`dplyr`; the core reweight needs no `LBSPR`. |
| `example.R` | A runnable, deterministic demonstration on three real sparids (one per sex system). |
| `example_life_history.csv` | The three species' life-history parameters (fork-length scale) so the folder runs on its own. |
| `README.md` | This file. |

## 4. Usage

```r
source("sex_structured_lbspr.R")

# Apply the amendment to a fitted stock. FM, SL50, SL95 come from a standard LBSPR fit.
sx <- spr_sex_structured(
  FM = 1.5, SL50 = 20.8, SL95 = 24.9,          # fitted fishing pressure + selectivity
  Linf = 47.4, MK = 1.5, L50 = 20.8, L95 = 24.9,
  sex_system = "protogyny", LD50 = 25.8, LD95 = 28.3,
  FecB = 3, CVLinf = 0.10, MaleExp = 3)
sx$SPR_gono_check   # standard LBSPR (replica of the package SPR)
sx$SPR_fem          # female (egg-based) SPR
sx$SPR_male         # male reproductive-capacity ratio (the floor)
sx$SPR_bind         # the status to report
```

`spr_sex_structured()` arguments: `FM`, `SL50`, `SL95` (from the fit); `Linf`, `MK`, `L50`,
`L95` (growth + maturity); `sex_system` (`"protandry"`, `"protogyny"`, or anything else →
gonochore); `LD50`, `LD95` (functional sex-change lengths; `LD95` defaults to `1.10 × LD50`);
`FecB` (fecundity–length exponent, default 3); `CVLinf` (default 0.10); `MaleExp` (mass exponent
for male capacity, default 3); optional `anchor_SPR` (the package SPR, to express the correction
as a ratio applied to it).

To run the whole **fit → reweight** path from measured lengths (requires the `LBSPR` package):

```r
res <- fit_lbspr_one(L, lh_row)   # res$SPR = standard LBSPR; res$SPR_bind = amended status
```

Run the demonstration:

```sh
Rscript example.R
```

## 5. References

- Hordyk, A., Ono, K., Valencia, S., Loneragan, N. & Prince, J. (2015a) *Some explorations of
  the life history ratios to describe length composition, spawning-per-recruit, and the
  spawning potential ratio.* ICES Journal of Marine Science.
- Hordyk, A., Ono, K., Sainsbury, K., Loneragan, N. & Prince, J. (2015b) *A novel length-based
  empirical estimation method of spawning potential ratio (SPR), and tests of its performance,
  for small-scale, data-poor fisheries.* ICES Journal of Marine Science.
- Hordyk, A.R., Ono, K., Prince, J.D. & Walters, C.J. (2016) *A simple length-structured model
  based on life history ratios and incorporating size-dependent selectivity: application to
  spawning potential ratios for data-poor stocks.* Canadian Journal of Fisheries and Aquatic
  Sciences.
- Hordyk, A. (`LBSPR`) *LBSPR: Length-Based Spawning Potential Ratio.* R package.
- Goodyear, C.P. (1993) *Spawning stock biomass per recruit in fisheries management: foundation
  and current use.*
- Buxton, C.D. & Garratt, P.A. (1990) *Alternative reproductive styles in seabreams (Pisces:
  Sparidae).* Environmental Biology of Fishes.
- Shepherd, G.R. (1993) *Length-based analyses of yield and spawning biomass per recruit for
  black sea bass, a protogynous hermaphrodite.*
- Punt, A.E., Garratt, P.A. & Govender, A. (1993) *On an approach for applying per-recruit
  methods to a protogynous hermaphrodite.* South African Journal of Marine Science.
- Provost, M.M. & Jensen, O.P. (2015) *The impacts of fishing on hermaphroditic species and
  treatment of sex change in stock assessments.* Fisheries.
- Heppell, S.S., Heppell, S.A., Coleman, F.C. & Koenig, C.C. (2006) *Models to compare
  management options for a protogynous fish.* Ecological Applications.
- Easter, E.E. & White, J.W. (2020) *Influence of protogynous sex change on recovery of fish
  populations within marine protected areas.* Ecological Applications.
- Benvenuto, C., Coscia, I., Chopelet, J., Sala-Bozano, M. & Mariani, S. (2017) *Ecological and
  evolutionary consequences of alternative sex-change pathways in fish.* Scientific Reports.

---

*Extracted from the BGTW Sparidae length-based stock assessment (LBSPR with a
sequential-hermaphroditism amendment). Licensed as the parent repository (see `../LICENSE`).*
