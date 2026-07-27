# =====================================================================
# Bream catch analysis  —  Sparidae field season 2026
# =====================================================================
# RE-RUN WORKFLOW (unchanged)
#   Drop any new CatchApp export named  bream_<date>_<time>.csv  into the
#   all_records/ folder and drop matching photos into  photos/ , then source
#   this file again. The loader globs every bream_*.csv automatically, so NO
#   code edits are needed as data arrives. Figures are written to
#   graphs and maps/  and result tables to  results/  with stable filenames,
#   so each re-run overwrites the previous outputs in place.
#
#   Every analytical section below is GUARDED: if a package is missing, a
#   parameter is absent, or a model fails to converge, that section is skipped
#   with a message and the rest of the script still runs to completion.
#
# CONTENTS  (top-to-bottom; ORDER MATTERS - later sections use objects built earlier)
#   packages        attach + optional-package guards (sets have_vegan / have_inext / have_lbspr ...)
#   run-control     ALL tunable knobs in one block, grouped (assessment / uncertainty /
#                   length-type / amendment / simulation / timezone)
#   visual style    one shared ggplot theme + colour palette for every figure
#   1   Load        auto-discovers every bream_*.csv in all_records/
#   2   Clean       types, length reconciliation, tag status
#   3   Length-estimate tables   (drive the grey angler-estimate overlays)
#   4   Plots       catch / length / date figures on the shared theme
#   5   Map         leaflet catch map
#   7   Life-history  parameter table + sex_info + per-species ASSESSMENT-READINESS report
#   8   Diversity   Shannon, Simpson, Hill numbers, coverage-based rarefaction
#   9   Size structure + length-based indicators (Froese / LBI)
#   10  LB-SPR      SPR + F/M with bootstrap AND parameter-MC 95% CIs (two status figures)
#   10b Sex-structured SPR   re-weights SPR by sex-at-length for hermaphrodites
#                   -- ALWAYS-ON production code; the amendment this script is built around
#   11  Diurnal vs nocturnal   (auto-activates once night records exist)
#   12  Length-weight relationship + condition (auto-activates once weights exist)
#   13  Tagging     live effort summary + movement analysis (COMMENTED until returns)
#   §Z  Operating-model simulation   GATED validation of the amendment (switch RUN_OM_SIM)
#                   -- defined unconditionally, runs only when on; NOTHING depends on it
#   14  Save figures -> graphs and maps/   (result tables are saved to results/ as they are made)
#
# PLACEMENT LOGIC: the amendment (10b) is production code that later sections consume, so it
# lives inside the assessment; the simulation (§Z) only validates it and nothing reads its
# output, so it sits at the foot of the file behind a switch.
#
# NOTE ON METHODS: this version runs LB-SPR only. LIME and LBB were removed -
# both are GitHub/JAGS/TMB-only with fragile installs, and for a SINGLE season
# they estimate essentially the same SPR/F-M as LBSPR (their advantage is
# multi-year data), so they added install pain without new information here.
# =====================================================================

# ---- packages -------------------------------------------------------
# install.packages(c("tidyverse","here","janitor","skimr","sf","leaflet",
#                     "vegan","iNEXT","LBSPR","ggrepel","htmlwidgets","geosphere"))
# install.packages("pak")
# pak::pkg_install("james-thorson/FishLife")
# install.packages(c("fishmethods", "TropFishR"))
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, ggplot2, readr, purrr, stringr, forcats, lubridate
  library(here)
  library(janitor)
  library(skimr)
  library(sf)
  library(leaflet)
  library(FishLife)
  library(fishmethods)
  library(TropFishR)
})

# Attach optional packages if present; otherwise mark the section to skip.
need <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (ok) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  else message(sprintf("[skip] package '%s' not installed - its section will be skipped.", pkg))
  ok
}
have_vegan <- need("vegan")
have_inext <- need("iNEXT")
have_lbspr <- need("LBSPR")
have_fishmethods <- need("fishmethods")   # empirical M estimators (Section 7b)
have_fishlife  <- requireNamespace("FishLife",  quietly = TRUE)   # multivariate LH covariance (GitHub: james-thorson/FishLife)
have_rfishbase <- requireNamespace("rfishbase", quietly = TRUE)   # FishBase audit (CRAN)
have_tropfishr <- requireNamespace("TropFishR", quietly = TRUE)   # length-converted catch curve (CRAN)
have_geo   <- requireNamespace("geosphere", quietly = TRUE)   # for commented movement block

select <- dplyr::select
filter <- dplyr::filter

# ---- FishLife Install ----
# Sanity check that FishLife is properly installed and is sourcing the correct data
# library(FishLife)
# Taxa <- Search_species(Genus = "Spondyliosoma", Species = "cantharus")$match_taxonomy
# P <- Plot_taxa(Taxa, mfrow = c(2, 2))
# P[[1]]$Cov_pred[c("Loo","K","M","Lm"), c("Loo","K","M","Lm")]

# ---- run-control knobs (edit these, not the code below) -------------
# - core assessment defaults -
BINWIDTH      <- 1      # length-bin width (cm) for histograms & assessment models
MK_DEFAULT    <- 1.5    # Beverton-Holt default M/K used when M not in the database
L95_FACTOR    <- 1.10   # FALLBACK ONLY: workbook ships L95 (=1.20*L50); used only if a cell is blank
# - minimum sample sizes that gate each model -
MIN_N_ASSESS  <- 20     # min MEASURED lengths before a species enters LBSPR
FM_CEIL       <- 4      # reliability ceiling: an estimated F/M above this flags an unstable data-poor fit
MIN_N_IND     <- 10     # min MEASURED lengths before a species gets length-indicators
# - uncertainty: length bootstrap + parameter Monte-Carlo (life-history) -
# Monte-Carlo reps. FINAL_RUN = TRUE for camera-ready (2000); FALSE for fast iteration (1000).
# 9999 was overkill: a percentile CI is already stable by ~2000 full GTG fits, at a fraction of the runtime.
FINAL_RUN     <- FALSE
N_BOOT        <- if (FINAL_RUN) 2000L else 1000L   # length-bootstrap reps for LBSPR
N_PARAM_MC    <- if (FINAL_RUN) 2000L else 1000L   # parameter Monte-Carlo draws (life-history)
PARAM_MC_NESTED <- FALSE # TRUE -> also resample lengths within each draw (parameter + sampling = total uncertainty)
# Multiplicative draw spread (lognormal sigma on the log scale, ~ a CV) for each Monte-Carlo input.
# Data-poor LBSPR is dominated by life-history uncertainty (Hordyk 2014; Medeiros-Leal 2023), which the
# length-only bootstrap cannot see; these knobs inject it. Widen/narrow per how well each input is known.
# MK now governs the M (numerator) draw; K is drawn separately and M/K is formed per draw
PARAM_MC_CV   <- c(MK = 0.15, K = 0.10, Linf = 0.10, CVLinf = 0.20, L50 = 0.10, LD50 = 0.10, LD95 = 0.10)
# - section activation thresholds (diel / weight) -
NIGHT_MIN_N   <- 8      # min night records before the diurnal/nocturnal section activates
WEIGHT_MIN_N  <- 10     # min weights before the length-weight section activates

# Length-type reconciliation. Your catch lengths are FORK LENGTH (length_type = "FL").
# The literature Linf and L50 below are almost all TOTAL LENGTH (TL). LBSPR compares
# your data to Linf/L50, so the two must be on the same scale. With CONVERT_TL_TO_FL
# = TRUE the tabulated TL parameters are multiplied by TL_FL_RATIO to put them on the
# FL scale that matches your data. 0.92 is a typical sparid FL/TL ratio; refine per
# species with your own FL-TL regressions when you have enough paired measurements.
CONVERT_TL_TO_FL <- TRUE
TL_FL_RATIO      <- 0.92   # FL ≈ TL * 0.92  for deep-bodied seabreams

# Sex-structured SPR (sequential hermaphroditism amendment to LB-SPR). Standard LBSPR
# treats every mature fish as an egg-producer; for protandrous/protogynous sparids that
# is wrong. With this ON, SPR is recomputed with a proportion-female-at-length weighting
# (see section 10b). The LBSPR FIT (F/M, SL50, SL95) is NOT affected - only the SPR.
# Sex system and length-at-sex-change live in the sex_info table in section 7.
SEX_STRUCTURED_SPR <- TRUE   # FALSE -> behave exactly like standard LBSPR
USE_FISHLIFE_CORR  <- TRUE   # TRUE -> FishLife correlation in the parameter MC; FALSE -> independent draws
RUN_OM_SIM <- FALSE  # operating-model validation (slow). FALSE for normal runs; TRUE only to re-validate §Z.
LDELTA_FACTOR      <- 1.10   # LD95 = LDELTA_FACTOR * LD50 when only LD50 is provided (ogive width)

# Local timezone for the diurnal/nocturnal split. Your timestamps are stored in UTC;
# Gibraltar is UTC+1 (winter) / UTC+2 (summer). Day/night must be judged in LOCAL time
# (see section 11), or an early-morning daytime catch gets mislabelled as night.
LOCAL_TZ <- "Europe/Gibraltar"

# - life-history CSV: authoritative parameter source, ingested in Section 7 -
# Flat one-row-per-species table. Already on trinomial names, ships CVLinf / FecB /
# length-weight a,b as columns, and carries M_estimated / K_estimated as logicals
# (no dagger parsing needed). Verified value-for-value against the former Excel
# workbook, so swapping the source is output-neutral.
# Authoritative parameter file. This MUST be the Rready export: it is the only one that
# ships the length-weight a / b / lw_src columns that Section 7 requires (coalesce(b, 3),
# coalesce(a, 0.01)) and it carries the accepted trinomials. The human-facing biometrics.csv
# is value-identical for columns 1-27 but LACKS a / b / lw_src, so pointing here at it would
# halt Section 7 with "object 'b' not found". If you keep the file under a different name in
# the project root, change only this string.
PARAM_CSV <- here("Sparid_LBSPR_LifeHistory_trinomial_Rready.csv")

dir.create(here("graphs and maps"), showWarnings = FALSE)
dir.create(here("results"),         showWarnings = FALSE)

# =====================================================================
# VISUAL STYLE  — one theme + palette shared by every figure
# =====================================================================
# A single warm, colour-blind-aware palette and a clean minimal theme are applied
# to every plot so the report reads as one coherent set, with plain-language titles
# and captions aimed at readers who do not know the underlying analyses.
BREAM <- list(
  gold      = "#E1A140",   # primary: measured fish
  gold_dark = "#9C6B1E",   # outlines / accents
  estimate  = "#9AA7B1",   # angler size-estimates (muted slate)
  ink       = "#33312E",   # text / points
  teal      = "#2E7D8A",   # secondary accent
  line      = "#B9C2C9"    # connector lines
)
REF_COLS  <- c(Lm = "#2E8B57", Lopt = "#1F6FB2", Lmega = "#7B4FA3", Lmean = "#D1495B")
ZONE_FILL <- c(Overfished = "#F4CFC6", Cautionary = "#FBE7B0", Healthy = "#CBE5D4")
DIEL_COLS <- c(day = "#F2A93B", night = "#3B4A6B")

theme_bream <- function(base_size = 13) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title            = element_text(face = "bold", size = rel(1.22), hjust = 0,
                                           margin = margin(b = 3), colour = BREAM$ink),
      plot.subtitle         = element_text(size = rel(0.92), hjust = 0, colour = "grey38",
                                           margin = margin(b = 12), lineheight = 1.15),
      plot.caption          = element_text(size = rel(0.78), hjust = 0, colour = "grey45",
                                           margin = margin(t = 12), lineheight = 1.2),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title            = element_text(colour = "grey30", size = rel(0.95)),
      axis.title.x          = element_text(margin = margin(t = 8)),
      axis.title.y          = element_text(margin = margin(r = 8), angle = 90),
      axis.text             = element_text(colour = "grey45"),
      panel.grid.minor      = element_blank(),
      panel.grid.major      = element_line(colour = "grey91", linewidth = 0.4),
      panel.spacing         = unit(1.1, "lines"),
      strip.background      = element_rect(fill = "grey95", colour = NA),
      strip.text            = element_text(face = "bold", size = rel(0.85), colour = "grey25",
                                           margin = margin(t = 5, b = 5)),
      legend.position       = "top",
      legend.title          = element_text(face = "bold", size = rel(0.85), colour = "grey30"),
      legend.text           = element_text(size = rel(0.85), colour = "grey35"),
      plot.margin           = margin(16, 20, 14, 16)
    )
}
theme_set(theme_bream())

# =====================================================================
# 1. LOAD  — auto-discovers every bream_*.csv in all_records/
# =====================================================================
raw <- list.files(here("all_records"), pattern = "^bream_.*\\.csv$", full.names = TRUE) |>
  set_names(basename) |>
  map(\(f) read_csv(f, na = c("", "NA", "NR"), show_col_types = FALSE)) |>
  list_rbind(names_to = "source_file") |>
  clean_names() |>
  distinct(fish_id, .keep_all = TRUE)

if (interactive()) { glimpse(raw); print(skim(raw)) }   # console-only; skipped under batch source()

# =====================================================================
# 2. CLEAN
# =====================================================================
cleaned <- raw |>
  mutate(
    had_existing_tag = as.logical(had_existing_tag),
    was_tagged       = as.logical(was_tagged),
    has_photo        = as.logical(has_photo),
    
    date          = ymd(date),
    timestamp_utc = as_datetime(timestamp_utc),
    
    length_source = case_when(
      !is.na(length_true_cm)  ~ "measured",
      !is.na(length_estimate) ~ "estimated",
      TRUE                    ~ "none"
    ),
    length    = coalesce(as.character(length_true_cm), length_estimate),
    length_cm = coalesce(
      length_true_cm,
      map_dbl(str_split(length_estimate, "-"), \(x) suppressWarnings(mean(as.numeric(x))))
    )
  )

cleaned <- cleaned |>
  mutate(
    common_name   = fct_infreq(common_name),
    release_fate  = fct_relevel(release_fate, "Kept", "Released Alive"),
    length_source = factor(length_source, levels = c("measured", "estimated", "none")),
    tag_status = case_when(
      was_tagged       ~ "newly tagged",
      had_existing_tag ~ "recapture (existing tag)",
      TRUE             ~ "untagged"
    )
  )

# NOTE: the original script dropped the recapture/movement/growth columns here.
# They are RETAINED now because section 13 (movement) needs them once returns
# arrive. Only genuinely empty admin columns are dropped, and any_of() means the
# drop never errors if a future export omits one of them.
cleaned <- cleaned |>
  select(-any_of(c(
    "notes", "tag_condition", "fishing_method", "bait_lure", "depth_m",
    "habitat", "tide_state", "session_id", "session_start_utc", "session_end_utc",
    "num_anglers"
  )))

cleaned_sf <- cleaned |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# species lookup: scientific <-> common name, used to label facets/axes with the
# friendly common names readers recognise while the data stays keyed by binomials.
sp_lookup <- cleaned |> distinct(scientific_name, common_name)
sci2com <- sp_lookup |>
  transmute(scientific_name = as.character(scientific_name),
            common_name     = as.character(common_name)) |>
  tibble::deframe()

# helper: finite MEASURED fork lengths for one species (used by the indicator and
# assessment loops, so the same extraction is not repeated/duplicated inline).
measured_lengths <- function(sp) {
  v <- cleaned |>
    filter(scientific_name == sp, length_source == "measured") |>
    pull(length_true_cm)
  v[is.finite(v)]
}

# Sanity checks (console-only; skipped in batch source())
if (interactive()) {
  print(nrow(cleaned)); print(n_distinct(cleaned$common_name))
  print(cleaned |> count(length_source))
  print(cleaned |> summarise(n_dates = n_distinct(date)))
}

# =====================================================================
# 3. LENGTH-ESTIMATE TABLES  (drive the grey estimate overlays)
# =====================================================================
estimates_by_sp <- cleaned |>
  filter(length_source == "estimated") |>
  count(common_name, length_estimate) |>
  mutate(
    lo = as.numeric(str_extract(length_estimate, "^\\d+")),
    hi = as.numeric(str_extract(length_estimate, "\\d+$")),
    lo = coalesce(lo, 0)
  ) |>
  filter(!is.na(hi))

estimates <- estimates_by_sp |>
  summarise(n = sum(n), .by = c(lo, hi))

# =====================================================================
# 4. PLOTS  (restyled with the shared theme + plain-language captions)
# =====================================================================
measured <- filter(cleaned, length_source == "measured")

histogramlength <- ggplot() +
  geom_histogram(data = measured, aes(length_true_cm), binwidth = BINWIDTH,
                 fill = BREAM$gold, colour = BREAM$gold_dark, linewidth = 0.25) +
  geom_rect(data = estimates, aes(xmin = lo, xmax = hi, ymin = 0, ymax = n),
            fill = BREAM$estimate, colour = "white", alpha = 0.55, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "How big are the seabreams being caught?",
       subtitle = "Length-frequency of the whole catch (all species pooled)",
       x = "Fork length (cm)", y = "Number of fish",
       caption = "Gold bars = individually measured fish. Grey blocks = fish the angler size-estimated, drawn across the\nlength range reported. Bar width = 1 cm.")

lengthfrequency <- ggplot() +
  geom_histogram(data = measured, aes(length_true_cm), binwidth = BINWIDTH,
                 fill = BREAM$gold, colour = BREAM$gold_dark, linewidth = 0.2) +
  geom_rect(data = estimates_by_sp, aes(xmin = lo, xmax = hi, ymin = 0, ymax = n),
            fill = BREAM$estimate, colour = "white", alpha = 0.55, linewidth = 0.25) +
  facet_wrap(~common_name, ncol = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Length-frequency by species",
       subtitle = "Same data as above, split per species (note the differing fish counts)",
       x = "Fork length (cm)", y = "Number of fish",
       caption = "Gold = measured fish; grey = angler size-estimates across their reported range.")

sp_counts <- cleaned |> count(common_name, scientific_name, sort = TRUE)

catchesspecies <- sp_counts |>
  ggplot(aes(n, fct_reorder(common_name, n))) +
  geom_col(fill = BREAM$gold, width = 0.72) +
  geom_text(aes(label = n), hjust = -0.35, size = 3.3, colour = "grey35") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Which seabreams dominate the catch?",
       subtitle = "Number of fish recorded per species",
       x = "Number of fish", y = NULL,
       caption = "Two tournaments plus recreational logging; black seabream and the two Diplodus dominate.")

sp_grouping <- cleaned |>
  group_by(common_name) |>
  summarise(
    n = n(),
    mean_len = mean(length_true_cm, na.rm = TRUE),
    sd_len   = sd(length_true_cm, na.rm = TRUE),
    min_len  = suppressWarnings(min(length_true_cm, na.rm = TRUE)),
    max_len  = suppressWarnings(max(length_true_cm, na.rm = TRUE)),
    n_kept     = sum(release_fate == "Kept"),
    n_released = sum(release_fate == "Released Alive"),
    .groups = "drop"
  ) |>
  arrange(common_name)

lengthspecies <- cleaned |>
  filter(length_source == "measured") |>
  ggplot(aes(length_true_cm,
             fct_reorder(common_name, length_true_cm, median, .na_rm = TRUE))) +
  geom_boxplot(fill = scales::alpha(BREAM$gold, 0.5), colour = BREAM$gold_dark,
               outlier.shape = NA, width = 0.6, linewidth = 0.4) +
  geom_jitter(height = 0.16, size = 1.9, alpha = 0.5, colour = BREAM$ink) +
  labs(title = "Size range of each species",
       subtitle = "Box = middle 50% of fish, line = median, dots = individual fish",
       x = "Fork length (cm)", y = NULL,
       caption = "Species ordered by median length. Only measured fish are shown.")

catchdate <- cleaned |> count(date) |>
  ggplot(aes(date, n)) +
  geom_col(fill = BREAM$teal, width = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "When were fish caught?",
       subtitle = "Records per day",
       x = NULL, y = "Number of fish",
       caption = "Spikes are the two tournament days; smaller bars are recreational sessions.")

lengthwater <- cleaned |>
  filter(is.finite(length_true_cm), is.finite(water_temp_c)) |>
  ggplot(aes(water_temp_c, length_true_cm)) +
  geom_point(size = 2.6, alpha = 0.55, colour = BREAM$teal) +
  labs(title = "Does fish size track water temperature?",
       subtitle = "Each point is one measured fish",
       x = "Water temperature (°C)", y = "Fork length (cm)",
       caption = "The two temperature clusters are the two sampling periods, not a biological gradient — interpret with care.")

# =====================================================================
# 5. MAP  (your original, unchanged)
# =====================================================================
pal <- leaflet::colorFactor("viridis", cleaned_sf$common_name)
catch_map <- leaflet(cleaned_sf) |>
  addProviderTiles("CartoDB.Positron") |>
  addCircleMarkers(radius = 6, fillColor = ~pal(common_name),
                   color = "#14100b", weight = 0.5, fillOpacity = 0.8,
                   popup = ~paste0(common_name, "<br>", length_true_cm, "cm")) |>
  addLegend(pal = pal, values = ~common_name, title = "Species")

# #####################################################################
# #####################  ADDED ANALYSIS SECTIONS  #####################
# #####################################################################

# =====================================================================
# 7. LIFE-HISTORY PARAMETERS  (ingested from the workbook)  +  READINESS REPORT
# =====================================================================
# The life-history CSV biometrics.csv is the single authoritative source for growth,
# mortality, maturity, the sex-change ogive, the per-species CVLinf and FecB, and the
# length-weight a/b. It is a flat one-row-per-species table on trinomial names. This
# block reads it, renames the descriptive column headers to the short names the rest
# of the pipeline expects, reconciles length type, and rebuilds every derived column
# the downstream sections consume. The object it produces, `life_history`, is a
# drop-in for the former workbook-plus-tribble path, so sections 8 onward are unchanged.
#
# WHAT CHANGED VS THE OLD WORKBOOK PATH (all output-neutral; values verified identical):
#   - source is a CSV, not Excel: no readxl, no sheet/skip, no dagger parsing
#     (M_estimated / K_estimated arrive as logicals straight from the file).
#   - names are already trinomials: no name_fix recode.
#   - a / b / lw_src ship in the file: no inline lw_supp tribble.
#   - CVLinf and FecB ship per species: fully externalised, no hard-coded constants.
#
# UNITS & CONVENTIONS (unchanged):
#   - Lengths cm. Weight W(g) = a*L^b; a cancels in the SPR ratio (default 0.01),
#     b drives LBSPR Wbeta and the protogyny male exponent, b = 3 isometric only
#     where no power-law fit exists.
#   - ALL CSV Linf/L50/LD50 are TOTAL LENGTH; CONVERT_TL_TO_FL / TL_FL_RATIO
#     (top of file) reconcile them with the FORK-LENGTH catch data.
#   - M_estimated TRUE where mortality is estimated, not directly reported.
#     MK (= M/K) is the LB-SPR input.

life_history <- readr::read_csv(PARAM_CSV, show_col_types = FALSE) |>
  filter(!is.na(scientific_name)) |>
  # rename the descriptive CSV headers to the short names the pipeline uses
  rename(
    Linf = Linf_cm_TL,
    K    = K_per_yr,
    M    = M_per_yr,
    MK   = M_over_K,
    L50  = L50_maturity_cm_TL,
    L95  = L95_maturity_cm_TL,
    LD50 = LD50_sexchange_cm_TL,
    LD95 = LD95_sexchange_cm_TL
  ) |>
  mutate(
    M_estimated = as.logical(M_estimated),
    K_estimated = as.logical(K_estimated),
    b  = coalesce(b, 3),                            # isometric fallback where unverified
    a  = coalesce(a, 0.01),                         # cancels in the SPR ratio
    M  = coalesce(M, MK_DEFAULT * K),               # file fills all 19; fallback only
    MK = coalesce(M / K, MK, MK_DEFAULT),           # recompute precisely; file MK as backup
    CVLinf = coalesce(CVLinf, 0.10),                # per-species from file; 0.10 fallback
    FecB   = coalesce(FecB, 3),                     # per-species from file; 3 (isometric) fallback
    # L50 fallback via Froese & Binohlan (2000); file fills all 19 -> no-op guard
    L50_src = if_else(is.na(L50) & is.finite(Linf), "est(Linf, Froese-Binohlan)", "file"),
    L50     = if_else(is.na(L50) & is.finite(Linf), 10^(0.8979 * log10(Linf) - 0.0782), L50),
    # file ships L95 and LD95 as VALUES; the *_FACTOR knobs are FALLBACK ONLY,
    # used solely if a cell is empty (LD95 is empty for the non-sex-changers, by design).
    L95  = coalesce(L95,  L50  * L95_FACTOR),
    LD95 = coalesce(LD95, LD50 * LDELTA_FACTOR),
    sex_system = coalesce(sex_system, "gonochore")  # unlisted -> ssr cancels to standard LBSPR
  ) |>
  mutate(
    # length-type reconciliation: file is all TL -> put every length on the FL catch scale
    across(c(Linf, L50, L95, LD50, LD95),
           \(x) if (CONVERT_TL_TO_FL) x * TL_FL_RATIO else x),
    Lopt         = Linf * 3 / (3 + MK),             # Froese optimum (post-conversion scale)
    Lmega        = 1.10 * Lopt,
    length_basis = if (CONVERT_TL_TO_FL) "FL (converted from TL)" else "TL (as published)"
  )

# join-key guard: a silent NA from a name mismatch would drop a species out of
# ready_lbspr with no error. Warn loudly (do not halt) if a caught species is unmatched.
catch_sp <- unique(as.character(cleaned$scientific_name))
orphans  <- setdiff(catch_sp, life_history$scientific_name)
if (length(orphans))
  warning("[life-history] no parameter row for caught species: ",
          paste(orphans, collapse = ", "),
          " - check trinomial spelling against biometrics.csv.", call. = FALSE)

# =====================================================================
# 7b. NATURAL MORTALITY ENVELOPE  (TropFishR::M_empirical, multi-estimator)
#     DROP-IN REPLACEMENT for the former single-Pauly (fishmethods) block.
#
#     Design: the adopted workbook M stays the HEADLINE input and is NOT
#     overwritten. This block only (i) builds a cross-check envelope from
#     several published empirical estimators and (ii) derives a per-species
#     draw CV that the parameter Monte-Carlo (Section 10) uses to widen the
#     M draw for species whose M was estimated rather than reported.
#
#     Package-first: every estimator is TropFishR::M_empirical (already
#     attached), not hand-coded. TropFishR bundles the updated Then et al.
#     (2015) reformulations alongside Pauly (1980), which is what turns the
#     old single value (CV undefined -> widening never fired) into a real
#     envelope with a real CV.
#
#     Fixes vs the previous version:
#       - the `within_2x` column is now actually created (old code created
#         `in_envelope` but select()-ed `within_2x`, which halted the run
#         whenever fishmethods was attached);
#       - estimated-M species with no computable CV now fall back to a WIDER
#         fixed CV rather than silently inheriting the reported-M default.
#
#     Caveat to verify once: the exact method strings and argument names in
#     TropFishR::M_empirical have drifted slightly across versions. Each call
#     is individually guarded, so a renamed method is SKIPPED rather than
#     fatal, and the surviving methods are listed in the `methods` column.
#     Confirm against ?TropFishR::M_empirical in your installed version and
#     adjust the strings below if the `methods` column comes back short.
# =====================================================================

T_ANNUAL_C <- 18   # regional mean annual habitat temperature (Strait / Alboran ~17-18 C)

# Proxy longevity from growth: the VBGF age at which L = 0.95 * Linf (Taylor 1958),
# with t0 assumed 0. Used ONLY to unlock the tmax-based estimators for species with
# no reported maximum age. These are K-correlated with the growth-based estimators,
# so they narrow rather than widen the true spread; they are included for coverage
# and flagged in `methods`, not treated as independent evidence.
tmax_from_K <- function(K) if (is.finite(K) && K > 0) log(1 - 0.95) / (-K) else NA_real_

# Empirical-M envelope for one species. Linf arrives on the FL (catch) scale; the
# empirical predictors are calibrated on published TL growth, so Linf is taken back
# to TL here (K is a rate and is scale-free).
m_envelope <- function(Linf_FL, K) {
  na <- tibble(est_min = NA_real_, est_med = NA_real_, est_max = NA_real_,
               est_cv = NA_real_, n_est = 0L, methods = NA_character_)
  if (!have_tropfishr || !is.finite(Linf_FL) || !is.finite(K)) return(na)
  Linf_TL <- if (CONVERT_TL_TO_FL) Linf_FL / TL_FL_RATIO else Linf_FL
  tmx     <- tmax_from_K(K)

  # one guarded call per estimator -> a wrong/renamed method string is skipped
  try_m <- function(method, ...) {
    v <- tryCatch(suppressWarnings(TropFishR::M_empirical(..., method = method)),
                  error = function(e) NULL)
    if (is.null(v)) return(NULL)
    v <- as.numeric(v)[1]
    if (is.finite(v) && v > 0) setNames(v, method) else NULL
  }

  ests <- c(
    # growth-only (need Linf + K): available for every species
    try_m("Pauly_Linf",  Linf = Linf_TL, K_l = K, temp = T_ANNUAL_C),
    try_m("Then_growth", Linf = Linf_TL, K_l = K),
    # longevity-based (need tmax): from the growth proxy, so K-correlated
    if (is.finite(tmx)) try_m("Then_tmax",      tmax = tmx)            else NULL,
    if (is.finite(tmx)) try_m("Hoenig",         tmax = tmx)            else NULL,
    if (is.finite(tmx)) try_m("AlversonCarney", K_l  = K, tmax = tmx)  else NULL
  )
  ests <- ests[is.finite(ests) & ests > 0]
  if (!length(ests)) return(na)

  tibble(est_min = min(ests), est_med = median(ests), est_max = max(ests),
         est_cv  = if (length(ests) > 1) stats::sd(ests) / mean(ests) else NA_real_,
         n_est   = length(ests),
         methods = paste(names(ests), collapse = ";"))
}

if (have_tropfishr) {
  M_cross_check <- life_history |>
    transmute(scientific_name, sex_system, M_workbook = M, M_estimated, K, Linf) |>
    mutate(env = purrr::pmap(list(Linf, K), m_envelope)) |>
    tidyr::unnest(env) |>
    mutate(within_2x    = is.finite(est_med) &
                          M_workbook >= 0.5 * est_med & M_workbook <= 2 * est_med,
           ratio_to_med = round(M_workbook / est_med, 2)) |>
    select(scientific_name, sex_system, M_workbook, M_estimated,
           est_min, est_med, est_max, est_cv, n_est, methods,
           within_2x, ratio_to_med)

  message("\n====== NATURAL MORTALITY CROSS-CHECK (workbook M vs TropFishR envelope) ======")
  print(M_cross_check, n = Inf)
  message("within_2x = TRUE means the adopted M sits within a factor of two of the estimator median.")
  message("`methods` lists which estimators actually returned a value for that species.")
  message("=============================================================================\n")
  write_csv(M_cross_check, here("results", "M_cross_check.csv"))

  # Per-species MK draw CV for the parameter Monte-Carlo. Reported M keeps the global
  # default. Estimated M is ALWAYS drawn at least as wide as EST_MK_CV_FLOOR, i.e. never
  # tighter than reported M, and wider still where the estimator envelope genuinely
  # warrants it, capped at 0.40 so one wild estimator cannot blow the interval open. The
  # empirical estimators are largely K-driven, so a small envelope CV is expected and the
  # floor, not the envelope, usually governs; the envelope earns its keep as the
  # face-validity check (within_2x) rather than as a precise uncertainty source.
  EST_MK_CV_FLOOR <- 0.20
  mk_cv_tbl <- M_cross_check |>
    transmute(scientific_name,
              MK_cv = dplyr::case_when(
                M_estimated & is.finite(est_cv) ~ pmin(pmax(est_cv, EST_MK_CV_FLOOR), 0.40),
                M_estimated                     ~ EST_MK_CV_FLOOR,
                TRUE                            ~ PARAM_MC_CV[["MK"]]))

  life_history <- life_history |>
    left_join(mk_cv_tbl, by = "scientific_name") |>
    mutate(MK_cv = coalesce(MK_cv, PARAM_MC_CV[["MK"]]))
} else {
  message("[skip] TropFishR not installed - M cross-check skipped; MK_cv set to the global default.")
  life_history <- life_history |> mutate(MK_cv = PARAM_MC_CV[["MK"]])
}

# Auditable appendix: one row per species, every value traceable, now carrying the
# workbook's growth_mortality_source / confidence_note, the retained lw_src, the
# M_estimated / K_estimated flags, and the per-species MK_cv used by the Monte-Carlo.
write_csv(life_history, here("results", "life_history_sourced.csv"))

# =====================================================================
# 7c. FISHLIFE LIFE-HISTORY CORRELATION  (covariance for the parameter MC)
# =====================================================================
# The parameter Monte-Carlo (section 10) formerly drew Linf, K, M and L50 as INDEPENDENT
# lognormal jitters, which manufactures biologically impossible combinations (high-Linf +
# high-K + low-M) because it ignores the invariants that make those traits co-vary. FishLife
# (Thorson et al. 2017) supplies the missing dependence: its multivariate model returns a
# predictive covariance among (Loo, K, M, Lm) for any taxon. We keep the CURATED workbook
# values as the draw CENTRES and PARAM_MC_CV as the marginal SPREADS, borrowing ONLY
# FishLife's correlation - curated means, FishLife dependence (Nadon et al. 2016).
#
# FishLife is a GitHub package: remotes::install_github("james-thorson/FishLife"). If it is
# absent, a species lookup fails, or USE_FISHLIFE_CORR is FALSE, draw_lh falls back to
# independent draws, so the pipeline is unchanged in that case. get_fl_corr() is the single
# place to adjust should a future FishLife API differ.

# correlated standard normals via Cholesky (base R; no extra dependency). A tiny ridge
# rescues a numerically non-positive-definite correlation block.
rmvn_corr <- function(R) {
  U <- tryCatch(chol(R), error = function(e) chol(R + diag(1e-6, ncol(R))))
  as.numeric(crossprod(U, stats::rnorm(ncol(R))))     # t(U) %*% z  ~  N(0, R)
}

fl_traits <- c(Linf = "Loo", K = "K", M = "M", L50 = "Lm")   # our name -> FishLife trait
get_fl_corr <- function(sciname) {
  if (!have_fishlife) return(NULL)
  parts <- strsplit(trimws(sciname), "\\s+")[[1]]
  if (length(parts) < 2) return(NULL)
  tryCatch({
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)   # swallow Plot_taxa's device
    taxa <- FishLife::Search_species(Genus = parts[1], Species = parts[2])$match_taxonomy
    pred <- FishLife::Plot_taxa(taxa)
    Cov  <- pred[[1]]$Cov_pred
    if (!all(fl_traits %in% rownames(Cov))) return(NULL)
    R <- stats::cov2cor(Cov[fl_traits, fl_traits, drop = FALSE])
    dimnames(R) <- list(names(fl_traits), names(fl_traits))          # -> Linf / K / M / L50
    R
  }, error = function(e) NULL)
}

fl_corr <- list()
if (isTRUE(USE_FISHLIFE_CORR) && have_fishlife) {
  fl_corr <- setNames(lapply(life_history$scientific_name, get_fl_corr),
                      life_history$scientific_name)
  n_ok <- sum(!vapply(fl_corr, is.null, logical(1)))
  message(sprintf("[FishLife] correlation obtained for %d of %d species; the rest use independent draws.",
                  n_ok, length(fl_corr)))
} else {
  message("[FishLife] correlated draws OFF (package absent or USE_FISHLIFE_CORR = FALSE); MC uses independent draws.")
}

# =====================================================================
# 7d. FISHBASE AUDIT  (rfishbase; non-blocking QC, never feeds estimation)
# =====================================================================
# Compares the adopted workbook values against FishBase central tendencies and writes
# results/fishbase_audit.csv. FishBase is pan-global, so deviation from a region-matched
# workbook value is EXPECTED and is a sanity flag, not an error. Comparison is on the
# published TL scale (workbook values taken back from FL); K and M are scale-free.
if (have_rfishbase) {
  fb_audit <- tryCatch({
    fb_key <- sub("^(\\S+\\s+\\S+).*$", "\\1", life_history$scientific_name)  # trinomial -> binomial
    pg <- rfishbase::popgrowth(fb_key)
    mt <- rfishbase::maturity(fb_key)
    pg_s <- pg |> dplyr::group_by(Species) |>
      dplyr::summarise(fb_Linf = stats::median(Loo, na.rm = TRUE),
                       fb_K    = stats::median(K,   na.rm = TRUE),
                       fb_M    = stats::median(M,   na.rm = TRUE), .groups = "drop")
    mt_s <- mt |> dplyr::group_by(Species) |>
      dplyr::summarise(fb_L50 = stats::median(Lm, na.rm = TRUE), .groups = "drop")
    inv <- if (CONVERT_TL_TO_FL) TL_FL_RATIO else 1
    tibble(scientific_name = life_history$scientific_name, Species = fb_key,
           wb_Linf = life_history$Linf / inv, wb_K = life_history$K,
           wb_M = life_history$M, wb_L50 = life_history$L50 / inv) |>
      left_join(pg_s, by = "Species") |>
      left_join(mt_s, by = "Species") |>
      mutate(Linf_ratio = round(wb_Linf / fb_Linf, 2),
             K_ratio    = round(wb_K    / fb_K,    2),
             M_ratio    = round(wb_M    / fb_M,    2),
             L50_ratio  = round(wb_L50  / fb_L50,  2),
             flag = dplyr::if_else(
               pmin(Linf_ratio, K_ratio, M_ratio, L50_ratio, na.rm = TRUE) < 0.7 |
               pmax(Linf_ratio, K_ratio, M_ratio, L50_ratio, na.rm = TRUE) > 1.4,
               "review", "ok"))
  }, error = function(e) { message("[FishBase audit] skipped: ", conditionMessage(e)); NULL })
  if (!is.null(fb_audit)) {
    write_csv(fb_audit, here("results", "fishbase_audit.csv"))
    message(sprintf("[FishBase audit] written for %d species; %d flagged 'review' (expected for region-matched values).",
                    nrow(fb_audit), sum(fb_audit$flag == "review", na.rm = TRUE)))
  }
} else {
  message("[skip] rfishbase not installed - FishBase audit skipped.")
}

# join measured-length counts and decide what is assessable
n_measured_sp <- cleaned |>
  filter(length_source == "measured") |>
  count(scientific_name, name = "n_measured")

readiness <- life_history |>
  left_join(n_measured_sp, by = "scientific_name") |>
  mutate(
    n_measured  = coalesce(n_measured, 0L),
    has_Linf    = is.finite(Linf),
    has_L50     = is.finite(L50),
    has_K       = is.finite(K),
    ready_lbspr = has_Linf & has_L50 & n_measured >= MIN_N_ASSESS,
    blocker = case_when(
      n_measured < MIN_N_ASSESS & (!has_Linf | !has_L50) ~ "too few fish AND missing Linf/L50",
      n_measured < MIN_N_ASSESS                          ~ sprintf("too few fish (%d < %d)", n_measured, MIN_N_ASSESS),
      !has_Linf & !has_L50                               ~ "missing Linf AND L50",
      !has_L50                                           ~ "missing L50",
      !has_Linf                                          ~ "missing Linf",
      TRUE                                               ~ "ready"
    )
  ) |>
  arrange(desc(n_measured))

message("\n================ ASSESSMENT-READINESS REPORT ================")
readiness |>
  select(scientific_name, n_measured, has_Linf, has_L50, ready_lbspr, blocker) |>
  print(n = Inf)
message("=============================================================\n")
write_csv(readiness, here("results", "assessment_readiness.csv"))

# =====================================================================
# 8. CATCH COMPOSITION  — rank-abundance + sample coverage
# =====================================================================
# Restrict to the BGTW Sparidae pool and drop unidentified records so the composition
# figure and the coverage statistic are robust to any non-sparid bycatch or NA species
# in a future export. Membership is taken from the life-history table (the single source
# of truth for the 19-species pool) and matched at GENUS level, which is immune to the
# binomial/trinomial subspecies naming difference (e.g. "Diplodus sargus sargus").
sparid_genera <- unique(stringr::word(life_history$scientific_name, 1))
abund <- cleaned |>
  dplyr::filter(!is.na(scientific_name),
                stringr::word(scientific_name, 1) %in% sparid_genera) |>
  count(scientific_name) |>
  arrange(desc(n))
abund_vec <- setNames(abund$n, abund$scientific_name)

# Shannon / Simpson / Hill indices are still computed and ARCHIVED to results/ for the
# record, but are no longer reported: a size- and target-selective rod-and-reel and
# citizen catch characterises composition, not standing-assemblage diversity, so no
# index table and no rarefaction/extrapolation figure are drawn from them.
diversity_summary <- NULL
if (have_vegan) {
  H  <- vegan::diversity(abund_vec, index = "shannon")
  D1 <- vegan::diversity(abund_vec, index = "simpson")
  D2 <- vegan::diversity(abund_vec, index = "invsimpson")
  J  <- H / log(length(abund_vec))
  diversity_summary <- tibble(
    metric = c("Species richness (q0)", "Shannon H'", "Hill q1 = exp(H')",
               "Gini-Simpson (1-D)", "Inverse Simpson (Hill q2)", "Pielou evenness J"),
    value  = c(length(abund_vec), H, exp(H), D1, D2, J)
  )
  write_csv(diversity_summary, here("results", "diversity_summary.csv"))
}

# Sample coverage (Chao & Jost, 2012) is kept as a single completeness statistic
# reported in prose alongside observed sparid richness; the accumulation CURVES
# (size- and coverage-based rarefaction/extrapolation) are no longer drawn.
inext_asy <- NULL
sample_coverage <- NA_real_
obs_richness    <- length(abund_vec)
if (have_inext) {
  ix <- iNEXT::iNEXT(list("BGTW Sparidae" = as.numeric(abund_vec)), q = c(0, 1, 2),
                     datatype = "abundance", endpoint = 2 * sum(abund_vec))
  inext_asy       <- ix$AsyEst
  sample_coverage <- ix$DataInfo$SC
  obs_richness    <- ix$DataInfo$S.obs
  write_csv(as_tibble(inext_asy),   here("results", "hill_numbers_asymptotic.csv"))
  write_csv(as_tibble(ix$DataInfo), here("results", "sample_coverage.csv"))
}

# Rank-abundance (Whittaker) plot - dependency-free, common-name labels.
# Pick the label layer explicitly (ggrepel if available, else plain text) rather than
# relying on operator precedence inside the ggplot chain.
rank_labels <- if (requireNamespace("ggrepel", quietly = TRUE)) {
  ggrepel::geom_text_repel(aes(label = common_name), size = 3, colour = "grey30",
                           max.overlaps = 20, seed = 1, min.segment.length = 0)
} else {
  geom_text(aes(label = common_name), size = 3, hjust = -0.05, colour = "grey30")
}
rankabund <- abund |>
  left_join(sp_lookup, by = "scientific_name") |>
  mutate(rank = row_number(), prop = n / sum(n)) |>
  ggplot(aes(rank, prop)) +
  geom_line(colour = BREAM$line, linewidth = 0.7) +
  geom_point(size = 2.8, shape = 21, fill = BREAM$gold, colour = BREAM$gold_dark, stroke = 0.7) +
  rank_labels +
  scale_y_log10(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "A few common species, a long tail of rare ones",
       subtitle = "Share of the catch by species, ranked most to least common (log scale)",
       x = "Species rank", y = "Share of catch",
       caption = "A steep drop after the top few species is typical of a diverse, lightly-structured assemblage.")

# =====================================================================
# 9. PER-SPECIES SIZE STRUCTURE  +  LENGTH-BASED INDICATORS (Froese / LBI)
# =====================================================================
# Indicators computed on MEASURED lengths only (you measure without size bias;
# angler estimate-ranges are too coarse for these ratios).
estimate_Lc <- function(x, binwidth = BINWIDTH) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NA_real_)
  br <- seq(0, ceiling(max(x)) + binwidth, by = binwidth)
  h <- hist(x, breaks = br, plot = FALSE)
  peak <- which.max(h$counts)
  asc <- h$counts[seq_len(peak)]
  idx <- which(asc >= 0.5 * max(h$counts))[1]      # ascending limb half-peak
  h$mids[idx]
}

ind_species <- readiness |>
  filter(n_measured >= MIN_N_IND, is.finite(Linf), is.finite(L50)) |>
  pull(scientific_name)

length_indicators <- map_dfr(ind_species, function(sp) {
  lh <- readiness |> filter(scientific_name == sp)
  L  <- measured_lengths(sp)
  Lc <- estimate_Lc(L)
  tibble(
    scientific_name = sp, n = length(L),
    Lmean = mean(L), Lmedian = median(L),
    Lm = lh$L50, Lopt = lh$Lopt, Lmega = lh$Lmega, Lc = Lc, Linf = lh$Linf,
    Pmat  = mean(L >= lh$L50),                                   # proportion mature
    Popt  = mean(L >= 0.9 * lh$Lopt & L <= 1.1 * lh$Lopt),       # within +/-10% of Lopt
    Pmega = mean(L > lh$Lmega),                                  # mega-spawners
    Lmean_over_Lopt = mean(L) / lh$Lopt
  ) |>
    mutate(Pobj = Pmat + Popt + Pmega)   # Cope & Punt composite (each ~1 is healthy)
})
if (nrow(length_indicators)) {
  print(length_indicators)
  write_csv(length_indicators, here("results", "length_based_indicators.csv"))
}

# annotated per-species length distributions with Lm / Lopt / Lmega / Lmean
indicator_plot <- NULL
if (nrow(length_indicators)) {
  ref <- length_indicators |>
    select(scientific_name, Lm, Lopt, Lmega, Lmean) |>
    pivot_longer(-scientific_name, names_to = "ref", values_to = "x")
  indicator_plot <- cleaned |>
    filter(scientific_name %in% ind_species, length_source == "measured") |>
    ggplot(aes(length_true_cm)) +
    geom_histogram(binwidth = BINWIDTH, fill = BREAM$gold, colour = BREAM$gold_dark, linewidth = 0.2) +
    geom_vline(data = ref, aes(xintercept = x, colour = ref), linewidth = 0.8) +
    facet_wrap(~scientific_name, scales = "free_y",
               labeller = labeller(scientific_name = sci2com)) +
    scale_colour_manual(values = REF_COLS, name = NULL,
                        breaks = c("Lm", "Lopt", "Lmega", "Lmean"),
                        labels = c("Maturity (Lm)", "Optimum (Lopt)",
                                   "Mega-spawner (Lmega)", "Catch mean")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Are caught fish mostly mature, and are big spawners present?",
         subtitle = "Size structure of the best-sampled species against biological reference points",
         x = "Fork length (cm)", y = "Number of fish",
         caption = "Maturity = size at first breeding; Optimum = size giving most yield per fish; Mega-spawner = large, highly\nfecund individuals. A healthy catch sits at or above maturity, with some fish near optimum and beyond.")
}

# =====================================================================
# 10. LB-SPR ASSESSMENT  + bootstrap & parameter-MC 95% CIs   (LBSPR only)
# =====================================================================
# LBSPR returns SPR, F/M and selectivity (SL50/SL95) from the length composition,
# given Linf, L50/L95 and M/K, under an equilibrium assumption.
# Reference points: SPR 0.40 healthy / 0.20 overfished; F/M = 1 means F = M.
# Section 10b then re-weights SPR for sequential hermaphrodites WITHOUT changing the fit.

make_lf <- function(x, Linf, binwidth = BINWIDTH) {
  x <- x[is.finite(x)]
  br <- seq(0, ceiling(max(c(x, Linf), na.rm = TRUE)) + binwidth, by = binwidth)
  mids <- br[-length(br)] + binwidth / 2
  counts <- as.numeric(table(cut(x, breaks = br, right = FALSE)))
  list(mids = mids, counts = counts)
}

# ── 10b. Sex-structured SPR (sequential hermaphroditism) ─────────────
# Standard LBSPR sets eggs-per-recruit = maturity(L) * L^FecB, i.e. it treats EVERY
# mature fish as an egg-producer. For sequential hermaphrodites that is wrong:
#   protogyny  -> the large fish are MALE (zero eggs);
#   protandry  -> the egg-producers ARE the large fish.
# This function recomputes SPR with a proportion-female-at-length weighting psi_f(L),
# faithfully re-implementing the LBSPR GTG forward model (LBSPRsim_) so that, for a
# gonochore (psi_f = 0.5, constant), SPR_gono reproduces the package SPR.
#
# DEFENSIBILITY: maturity/sex never enter the LBSPR likelihood, so the FIT (F/M,
# SL50, SL95) is unchanged. This takes the FITTED parameters and only reweights egg
# output - a transparent post-processing step, not a change to the estimator.
# Check SPR_gono_check ≈ fit@SPR to validate the replica.
#
#   psi_f(L): gonochore -> 0.5 (cancels) ; protogyny -> descending 1->0 ;
#             protandry -> ascending 0->1.  Ogive uses LD50/LD95 (same scale as Linf).
#   SPR_fem  : female (egg) SPR.   SPR_male: male-CAPACITY ratio (mature-male biomass-
#             per-recruit; a precautionary floor for protogyny, NOT a fertilisation
#             estimate - the mating function is unknown).   SPR_bind: status to report
#             = SPR_fem, except protogyny where it is min(SPR_fem, SPR_male).
spr_sex_structured <- function(FM, SL50, SL95,
                               Linf, MK, L50, L95, FecB = 3, CVLinf = 0.10,
                               sex_system = "gonochoristic", LD50 = NA, LD95 = NA,
                               anchor_SPR = NA_real_,
                               MaleExp = 3, binwidth = BINWIDTH, maxsd = 2, ngtg = 13) {
  raw <- tolower(trimws(as.character(sex_system)))
  if (length(raw) < 1 || is.na(raw) || raw == "") raw <- "gonochoristic"
  # Normalise to protandry / protogyny / gonochore. RUDIMENTARY hermaphrodites have NO
  # functional sex change -> treated as gonochore (no correction) but labelled as their
  # true mode in the output. grepl tolerates "protandrous"/"protandry" etc.
  mode <- if (grepl("protandr", raw)) "protandry" else
    if (grepl("protogyn", raw)) "protogyny" else "gonochore"  # incl. rudimentary/gonochoristic
  # a sex-changer without a documented sex-change length also falls back to gonochore
  ss_eff <- if (mode %in% c("protandry", "protogyny") && !is.finite(LD50)) "gonochore" else mode
  if (is.finite(LD50) && !is.finite(LD95)) LD95 <- LDELTA_FACTOR * LD50
  
  SDLinf <- CVLinf * Linf
  LBins  <- seq(0, 1.3 * Linf, by = binwidth)   # match LBSPRsim_ bin edges exactly
  LMids  <- LBins[-length(LBins)] + 0.5 * binwidth
  ngtg   <- max(ngtg, ceiling((2 * maxsd * SDLinf + 1) / binwidth))
  gtgLinfs <- seq(Linf - maxsd * SDLinf, Linf + maxsd * SDLinf, length = ngtg)
  recP <- dnorm(gtgLinfs, Linf, SDLinf); recP <- recP / sum(recP)
  
  # maturity-at-length per GTG (relative-size invariant, exactly as in LBSPRsim_)
  L50G <- L50 / Linf * gtgLinfs; L95G <- L95 / Linf * gtgLinfs; DeltaG <- L95G - L50G
  Mat <- sapply(seq_along(gtgLinfs), function(g)
    1 / (1 + exp(-log(19) * (LMids - L50G[g]) / DeltaG[g])))
  
  # proportion-female-at-length psi_f
  if (ss_eff == "gonochore") {
    PsiF <- matrix(0.5, nrow = length(LMids), ncol = ngtg)   # constant -> cancels
  } else {
    LD50G <- LD50 / Linf * gtgLinfs; LD95G <- LD95 / Linf * gtgLinfs; DeltaDG <- LD95G - LD50G
    Asc <- sapply(seq_along(gtgLinfs), function(g)
      1 / (1 + exp(-log(19) * (LMids - LD50G[g]) / DeltaDG[g])))
    PsiF <- if (ss_eff == "protandry") Asc else 1 - Asc       # protogyny = descending
  }
  
  Fec_gono <- Mat * LMids^FecB                  # original (gonochore-implicit)
  Fec_fem  <- Mat * PsiF * LMids^FecB           # female-only eggs
  Fec_male <- Mat * (1 - PsiF) * LMids^MaleExp  # male reproductive capacity
  
  # mortality ratios at bin edges, then the GTG survival recursion
  # selectivity at bin EDGES, with LBSPRsim_'s half-bin centering: SL50/SL95 shifted by 0.5*bin.
  # The +0.5*binwidth is exactly what makes SPR_gono_check reproduce the LBSPR package SPR (verified
  # against LBSPRsim_); without it the gonochore check ran ~5-12% low.
  selEdge <- 1 / (1 + exp(-log(19) * (LBins - (SL50 + 0.5 * binwidth)) / (SL95 - SL50)))
  MKMat  <- matrix(MK, nrow = length(LBins), ncol = ngtg)
  ZKLMat <- MKMat + matrix(FM * MK * selEdge, nrow = length(LBins), ncol = ngtg)
  NPRf <- NPRuf <- matrix(0, nrow = length(LBins), ncol = ngtg)
  NPRf[1, ] <- NPRuf[1, ] <- recP
  for (L in 2:length(LBins)) {
    ratio <- (gtgLinfs - LBins[L]) / (gtgLinfs - LBins[L - 1])
    ratio[!is.finite(ratio) | ratio < 0] <- 0
    NPRuf[L, ] <- NPRuf[L - 1, ] * ratio^MKMat[L - 1, ]
    NPRf[L, ]  <- NPRf[L - 1, ]  * ratio^ZKLMat[L - 1, ]
    ind <- gtgLinfs < LBins[L]; NPRf[L, ind] <- 0; NPRuf[L, ind] <- 0
  }
  NPRf[is.nan(NPRf) | NPRf < 0] <- 0; NPRuf[is.nan(NPRuf) | NPRuf < 0] <- 0
  
  # numbers AT length = (enter - leave) / rate
  NatLUF <- NatLF <- matrix(0, nrow = length(LMids), ncol = ngtg)
  for (L in seq_along(LMids)) {
    NatLUF[L, ] <- (NPRuf[L, ] - NPRuf[L + 1, ]) / MKMat[L, ]
    NatLF[L, ]  <- (NPRf[L, ]  - NPRf[L + 1, ])  / ZKLMat[L, ]
  }
  
  epr <- function(W) sum(NatLF * W) / sum(NatLUF * W)
  SPR_gono <- epr(Fec_gono)
  SPR_fem  <- epr(Fec_fem)
  SPR_male <- if (sum(NatLUF * Fec_male) > 0) epr(Fec_male) else NA_real_
  SPR_bind <- if (ss_eff == "protogyny") suppressWarnings(min(SPR_fem, SPR_male, na.rm = TRUE)) else SPR_fem

  # Package anchor (Option B): if the fitted package SPR is supplied, express the sex correction as
  # a RATIO applied to it (SPR_package x SPR_fem/SPR_gono), so any residual replica-vs-package gap
  # cancels. With the selectivity fix above SPR_gono ~ SPR_package, so anchor_factor ~ 1 (a no-op
  # safety net). anchor_SPR = NA (e.g. the OM in section Z, which has no package fit) -> factor 1.
  anchor_factor <- if (is.finite(anchor_SPR) && SPR_gono > 0) anchor_SPR / SPR_gono else 1
  SPR_fem_a  <- SPR_fem  * anchor_factor
  SPR_male_a <- SPR_male * anchor_factor
  SPR_bind_a <- if (ss_eff == "protogyny") suppressWarnings(min(SPR_fem_a, SPR_male_a, na.rm = TRUE)) else SPR_fem_a

  tibble(sex_system = raw, sex_applied = ss_eff, SPR_gono_check = SPR_gono,
         SPR_fem = SPR_fem_a, SPR_male = SPR_male_a, SPR_bind = SPR_bind_a,   # reported (anchored)
         anchor_factor = anchor_factor,
         SPR_fem_raw = SPR_fem, SPR_male_raw = SPR_male, SPR_bind_raw = SPR_bind)  # pre-anchor, for audit
}

fit_lbspr_one <- function(L, lh, binwidth = BINWIDTH, cvlinf = NULL) {
  if (!have_lbspr) return(NULL)
  # CVLinf: use the per-species value from the file unless the parameter Monte-Carlo
  # supplies a drawn one. coalesce-style fallback keeps it robust if the column is absent.
  if (is.null(cvlinf)) cvlinf <- if (!is.null(lh$CVLinf) && is.finite(lh$CVLinf)) lh$CVLinf else 0.10
  L <- L[is.finite(L)]; if (length(L) < 5) return(NULL)
  pars <- suppressMessages(methods::new("LB_pars"))
  pars@Species  <- lh$scientific_name
  pars@Linf     <- lh$Linf
  pars@L50      <- lh$L50
  pars@L95      <- lh$L95
  pars@MK       <- lh$MK
  pars@CVLinf   <- cvlinf
  pars@Walpha   <- lh$a
  pars@Wbeta    <- lh$b
  pars@BinWidth <- binwidth
  pars@L_units  <- "cm"
  lf <- make_lf(L, lh$Linf, binwidth)
  lenobj <- suppressMessages(methods::new("LB_lengths"))
  lenobj@LMids  <- lf$mids
  lenobj@LData  <- matrix(lf$counts, ncol = 1)
  lenobj@Years  <- 1L
  lenobj@NYears <- 1L
  fit <- tryCatch(suppressMessages(LBSPR::LBSPRfit(pars, lenobj, verbose = FALSE, useCPP = TRUE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  base <- tibble(method = "LBSPR", SPR = fit@SPR[1], FM = fit@FM[1],
                 SL50 = fit@SL50[1], SL95 = fit@SL95[1])
  if (isTRUE(SEX_STRUCTURED_SPR)) {
    sx <- spr_sex_structured(
      FM = fit@FM[1], SL50 = fit@SL50[1], SL95 = fit@SL95[1],
      Linf = lh$Linf, MK = lh$MK, L50 = lh$L50, L95 = lh$L95,
      FecB = lh$FecB, CVLinf = cvlinf,         # per-species FecB from file; cvlinf as used by LBSPRfit
      sex_system = lh$sex_system, LD50 = lh$LD50, LD95 = lh$LD95,
      anchor_SPR = fit@SPR[1],                 # Option B: anchor the sex correction to the package SPR
      MaleExp = lh$b, binwidth = binwidth)     # male capacity ~ mature-male biomass (L^b)
    base <- bind_cols(base, sx)
  }
  base
}

boot_lbspr <- function(L, lh, nboot = N_BOOT) {
  L <- L[is.finite(L)]; if (length(L) < 5) return(NULL)
  map_dfr(seq_len(nboot), function(i) {
    # resampled fits throw expected, non-actionable optimiser warnings; silence the
    # bootstrap only (the single point-estimate fit above keeps its warnings visible)
    r <- suppressWarnings(fit_lbspr_one(sample(L, replace = TRUE), lh))
    if (is.null(r)) return(NULL)
    cols <- intersect(c("SPR", "FM", "SPR_fem", "SPR_bind"), names(r))
    dplyr::mutate(r[, cols, drop = FALSE], rep = i, .before = 1)
  })
}

# percentile-CI summary shared by the bootstrap and the parameter Monte-Carlo (below).
# (Factored out so the two paths cannot drift apart; inline it back if you prefer.)
pct_ci <- function(df, metrics) {
  df |>
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") |>
    group_by(scientific_name, metric) |>
    summarise(median = median(value, na.rm = TRUE),
              lo95 = quantile(value, .025, na.rm = TRUE),
              hi95 = quantile(value, .975, na.rm = TRUE),
              lo80 = quantile(value, .10,  na.rm = TRUE),
              hi80 = quantile(value, .90,  na.rm = TRUE),
              .groups = "drop")
}

# ── Parameter Monte-Carlo: life-history uncertainty (Addition 1) ──────
# boot_lbspr() above resamples LENGTHS only, holding the life-history inputs fixed, so it
# sees sampling noise but NOT the uncertainty that dominates a data-poor LBSPR fit: M/K,
# Linf, CVLinf, L50 and the sex-change ogive LD50/LD95 (Hordyk 2014; Medeiros-Leal 2023;
# Magnusson 2013 - the length bootstrap is the most overconfident of the standard methods).
# This draws those inputs from distributions, refits, and pushes each draw through the same
# fit_lbspr_one() -> spr_sex_structured() path, giving a CI on the quantity that matters.
#   - independent multiplicative (lognormal) draws by default, spreads set by PARAM_MC_CV.
#     A correlated multivariate draw (Nadon 2016) is the refinement, but needs a covariance
#     we do not have, so independent draws are the honest default here.
#   - PARAM_MC_NESTED = TRUE also resamples the lengths inside each draw, so the interval
#     reflects parameter + sampling uncertainty together (total uncertainty).
#   - ordering kept consistent: L95 = L95_FACTOR*L50 (as in the main pipeline), LD95 > LD50.
# The bootstrap is retained alongside as the bias diagnostic, exactly as before.
rlnorm_mult <- function(x, cv) x * exp(stats::rnorm(1L, 0, cv))   # multiplicative jitter, keeps sign

draw_lh <- function(lh, cv = PARAM_MC_CV) {
  d <- lh
  # estimated-M species carry a wider M draw CV from §7b; reported M uses the global default
  mk_cv <- if (!is.null(lh$MK_cv) && is.finite(lh$MK_cv)) lh$MK_cv else cv[["MK"]]
  R <- fl_corr[[as.character(lh$scientific_name)]]     # FishLife correlation over (Linf,K,M,L50), or NULL
  if (!is.null(R)) {
    # correlated multiplicative draw: curated CENTRES, our marginal SIGMAS, FishLife DEPENDENCE.
    # M and K are drawn jointly, so M/K respects their correlation rather than being jittered as a ratio.
    sig  <- c(Linf = cv[["Linf"]], K = cv[["K"]], M = mk_cv, L50 = cv[["L50"]])
    mult <- exp(sig * rmvn_corr(R))                    # FishLife-correlated lognormal multipliers
    d$Linf <- lh$Linf * mult[["Linf"]]
    d$L50  <- lh$L50  * mult[["L50"]]
    Kd     <- lh$K    * mult[["K"]]
    Md     <- lh$M    * mult[["M"]]
  } else {
    # independent fallback (FishLife absent / lookup failed / toggle off)
    d$Linf <- rlnorm_mult(lh$Linf, cv[["Linf"]])
    d$L50  <- rlnorm_mult(lh$L50,  cv[["L50"]])
    Kd     <- rlnorm_mult(lh$K,    cv[["K"]])
    Md     <- rlnorm_mult(lh$M,    mk_cv)
  }
  d$K <- Kd; d$M <- Md; d$MK <- Md / Kd                # M/K formed from the joint draw
  d$L95 <- lh$L95 * (d$L50 / lh$L50)                   # preserve the adopted (workbook) L95/L50 ratio
  if (is.finite(lh$LD50)) {                            # only sex-changers carry an LD ogive
    d$LD50 <- rlnorm_mult(lh$LD50, cv[["LD50"]])
    gap    <- max(lh$LD95 - lh$LD50, 1e-6)             # keep the ogive ordered: LD95 > LD50
    d$LD95 <- d$LD50 + rlnorm_mult(gap, cv[["LD95"]])
  }
  cvlinf_centre <- if (!is.null(lh$CVLinf) && is.finite(lh$CVLinf)) lh$CVLinf else 0.10
  cvlinf <- rlnorm_mult(cvlinf_centre, cv[["CVLinf"]])  # centre from file; spread from PARAM_MC_CV
  list(lh = d, cvlinf = cvlinf)
}

param_mc_lbspr <- function(L, lh, ndraw = N_PARAM_MC, nested = PARAM_MC_NESTED) {
  L <- L[is.finite(L)]; if (length(L) < 5) return(NULL)
  map_dfr(seq_len(ndraw), function(i) {
    dr   <- draw_lh(lh)
    Luse <- if (isTRUE(nested)) sample(L, replace = TRUE) else L   # nested -> add sampling noise
    r <- suppressWarnings(fit_lbspr_one(Luse, dr$lh, cvlinf = dr$cvlinf))
    if (is.null(r)) return(NULL)
    cols <- intersect(c("SPR", "FM", "SPR_fem", "SPR_bind"), names(r))
    dplyr::mutate(r[, cols, drop = FALSE], draw = i, .before = 1)
  })
}

assess_species <- readiness |> filter(ready_lbspr) |> pull(scientific_name)

# Reproducible headline CIs: seed the length-bootstrap + parameter-Monte-Carlo draws with the
# same seed as the OM block (42) so the reported LBSPR 95% intervals are identical run-to-run.
set.seed(42)
assess_points <- list(); boot_store <- list(); param_store <- list()
for (sp in assess_species) {
  lh <- readiness |> filter(scientific_name == sp)
  L  <- measured_lengths(sp)
  
  r_lbspr <- fit_lbspr_one(L, lh)
  if (!is.null(r_lbspr)) {
    reliable_fit <- r_lbspr$FM <= FM_CEIL                      # implausibly high F/M -> unstable data-poor fit
    assess_points[[sp]] <- r_lbspr |>
      mutate(scientific_name = sp, n = length(L), reliable = reliable_fit, .before = 1)
    bl <- boot_lbspr(L, lh, N_BOOT)
    if (!is.null(bl) && nrow(bl)) boot_store[[sp]] <- bl |> mutate(scientific_name = sp)
    pm <- param_mc_lbspr(L, lh, N_PARAM_MC)
    if (!is.null(pm) && nrow(pm)) param_store[[sp]] <- pm |> mutate(scientific_name = sp)
  }
}

assess_results <- bind_rows(assess_points)
boot_results   <- bind_rows(boot_store)
param_results  <- bind_rows(param_store)

# Bootstrap CIs (LBSPR) for SPR and F/M  -- sampling uncertainty (bias diagnostic)
if (nrow(boot_results)) {
  boot_metrics <- intersect(c("SPR", "FM", "SPR_fem", "SPR_bind"), names(boot_results))
  boot_ci <- pct_ci(boot_results, boot_metrics) |> mutate(method = "LBSPR")
  write_csv(boot_ci, here("results", "lbspr_bootstrap_ci.csv"))
} else boot_ci <- tibble()

# Parameter Monte-Carlo CIs (LBSPR) -- life-history uncertainty, reported ALONGSIDE the bootstrap
if (nrow(param_results)) {
  pmc_metrics <- intersect(c("SPR", "FM", "SPR_fem", "SPR_bind"), names(param_results))
  param_ci <- pct_ci(param_results, pmc_metrics) |>
    mutate(method = "LBSPR",
           uncertainty = if (isTRUE(PARAM_MC_NESTED)) "param+sampling" else "parameter")
  write_csv(param_ci, here("results", "lbspr_param_mc_ci.csv"))
} else param_ci <- tibble()

if (nrow(assess_results)) {
  assess_display <- assess_results
  if (nrow(boot_ci)) {
    spr_ci <- boot_ci |> filter(metric == "SPR") |>
      transmute(scientific_name, SPR_lo95 = lo95, SPR_hi95 = hi95)
    fm_ci  <- boot_ci |> filter(metric == "FM") |>
      transmute(scientific_name, FM_lo95 = lo95, FM_hi95 = hi95)
    assess_display <- assess_display |>
      left_join(spr_ci, by = "scientific_name") |>
      left_join(fm_ci,  by = "scientific_name")
  }
  if (nrow(param_ci)) {        # parameter-MC interval reported next to the bootstrap one
    spr_pmc <- param_ci |> filter(metric == "SPR") |>
      transmute(scientific_name, SPR_pmc_lo95 = lo95, SPR_pmc_hi95 = hi95)
    fm_pmc  <- param_ci |> filter(metric == "FM") |>
      transmute(scientific_name, FM_pmc_lo95 = lo95, FM_pmc_hi95 = hi95)
    assess_display <- assess_display |>
      left_join(spr_pmc, by = "scientific_name") |>
      left_join(fm_pmc,  by = "scientific_name")
    message("CIs below: *_lo95/_hi95 = length bootstrap (sampling only); ",
            "*_pmc_lo95/_pmc_hi95 = parameter Monte-Carlo (life-history). The latter is ",
            "usually wider and is the more honest interval in data-poor LBSPR.")
  }
  print(assess_display)        # SPR/FM with bootstrap + parameter-MC 95% CIs where available
  write_csv(assess_results, here("results", "assessment_results.csv"))
}

# ---- 10b output: sex-structured SPR vs standard LBSPR ----------------
if (isTRUE(SEX_STRUCTURED_SPR) && nrow(assess_results) &&
    "SPR_bind" %in% names(assess_results)) {
  sex_compare <- assess_results |>
    transmute(
      scientific_name, n, sex_system, sex_applied,
      SPR_package = round(SPR, 3),            # standard LBSPR (gonochore-implicit)
      SPR_check   = round(SPR_gono_check, 3), # our replica of the same (validates the engine)
      SPR_female  = round(SPR_fem, 3),        # sex-structured egg-SPR (package-anchored)
      SPR_male    = round(SPR_male, 3),       # male-capacity ratio (protogyny floor, anchored)
      SPR_binding = round(SPR_bind, 3),       # status to report (anchored)
      anchor      = round(anchor_factor, 3),  # package/replica ratio applied; ~1.000 = replica matches package
      shift       = round(SPR_bind - SPR, 3)  # how far the correction moves the headline
    ) |>
    arrange(sex_applied, scientific_name)
  message("\n============ SEX-STRUCTURED SPR  (vs standard LBSPR) ============")
  message("Validation: SPR_package and SPR_check should now match to ~3 dp (selectivity fixed);")
  message("the 'anchor' factor should sit at ~1.000. SPR_female/_male/_binding are package-anchored.")
  message("Expectation: protandry -> SPR_female usually < package (safe, conservative);")
  message("protogyny -> SPR_female often > package, so SPR_binding takes the male floor.")
  print(sex_compare, n = Inf)
  message("================================================================\n")
  write_csv(sex_compare, here("results", "spr_sex_structured.csv"))
}

# ---- 10c output: LD50/LD95 sensitivity of the binding SPR (post-fit; no refit) -------
# The amendment is post-fit, so the sex-change ogive can be perturbed WITHOUT refitting LBSPR. Each
# assessed functional hermaphrodite is re-evaluated at LD50 x {0.8, 1.0, 1.2} (LD95 carried to hold
# the ogive width; the base row is the headline fit itself). Writes results/ld_sensitivity.csv.
if (isTRUE(SEX_STRUCTURED_SPR) && nrow(assess_results) && "SPR_bind" %in% names(assess_results)) {
  ld_zone <- function(s) as.character(cut(s, c(-Inf, 0.20, 0.40, Inf),
                                          labels = c("overfished", "cautionary", "healthy")))
  ld_rows <- list()
  for (sp in assess_results$scientific_name) {
    fit_sp <- assess_results |> filter(scientific_name == sp)
    lhr    <- readiness     |> filter(scientific_name == sp)
    if (!(lhr$sex_system %in% c("protandry", "protogyny")) || !is.finite(lhr$LD50)) next
    if (!isTRUE(fit_sp$reliable)) next                                                    # skip flagged (uninterpretable) fits
    width <- if (is.finite(lhr$LD95)) lhr$LD95 - lhr$LD50 else (LDELTA_FACTOR - 1) * lhr$LD50
    for (m in c(0.8, 1.0, 1.2)) {
      LD50m <- lhr$LD50 * m; LD95m <- LD50m + width
      if (m == 1.0) {
        SPR_fem <- fit_sp$SPR_fem; SPR_male <- fit_sp$SPR_male; SPR_bind <- fit_sp$SPR_bind
      } else {
        sx <- spr_sex_structured(
          FM = fit_sp$FM, SL50 = fit_sp$SL50, SL95 = fit_sp$SL95,
          Linf = lhr$Linf, MK = lhr$MK, L50 = lhr$L50, L95 = lhr$L95,
          FecB = lhr$FecB, CVLinf = lhr$CVLinf,
          sex_system = lhr$sex_system, LD50 = LD50m, LD95 = LD95m,
          anchor_SPR = fit_sp$SPR, MaleExp = lhr$b, binwidth = BINWIDTH)
        SPR_fem <- sx$SPR_fem; SPR_male <- sx$SPR_male; SPR_bind <- sx$SPR_bind
      }
      ld_rows[[length(ld_rows) + 1]] <- tibble(
        scientific_name = sp, sex_system = lhr$sex_system, LD_mult = m,
        LD50 = LD50m, LD95 = LD95m,
        SPR_fem = SPR_fem, SPR_male = SPR_male, SPR_bind = SPR_bind, zone = ld_zone(SPR_bind))
    }
  }
  if (length(ld_rows)) {
    ld_sensitivity <- bind_rows(ld_rows)
    write_csv(ld_sensitivity, here("results", "ld_sensitivity.csv"))
    message("\n============ LD50/LD95 SENSITIVITY (binding SPR at +/-20% sex-change ogive) ============")
    print(ld_sensitivity, n = Inf)
    message("========================================================================================\n")
  }
}

# ---- assessment figures (LBSPR-only): SPR stoplight + F/M ------------
spr_status <- NULL; fishing_pressure <- NULL
if (nrow(assess_results)) {
  asd <- assess_results |> left_join(sp_lookup, by = "scientific_name")
  if (nrow(boot_ci)) {
    asd <- asd |>
      left_join(boot_ci |> filter(metric == "SPR") |>
                  select(scientific_name, spr_lo = lo95, spr_hi = hi95), by = "scientific_name") |>
      left_join(boot_ci |> filter(metric == "FM") |>
                  select(scientific_name, fm_lo = lo95, fm_hi = hi95),  by = "scientific_name")
  }
  
  # (a) SPR stoplight: the headline status figure
  asd_spr <- asd |> mutate(common_name = fct_reorder(common_name, SPR))
  spr_err <- if (all(c("spr_lo", "spr_hi") %in% names(asd_spr))) {
    geom_segment(data = asd_spr,
                 aes(x = spr_lo, xend = spr_hi, y = common_name, yend = common_name),
                 colour = BREAM$ink, linewidth = 0.7)
  } else NULL
  # sex-structured SPR overlay: white diamond at the binding SPR + dotted shift line
  spr_sex_layer <- if (isTRUE(SEX_STRUCTURED_SPR) && "SPR_bind" %in% names(asd_spr)) {
    list(
      geom_segment(data = asd_spr,
                   aes(x = SPR, xend = SPR_bind, y = common_name, yend = common_name),
                   colour = BREAM$ink, linewidth = 0.4, linetype = "dotted"),
      geom_point(data = asd_spr, aes(x = SPR_bind, y = common_name),
                 size = 3.4, shape = 23, fill = "white", colour = BREAM$ink, stroke = 0.8)
    )
  } else NULL
  spr_status <- ggplot(asd_spr, aes(SPR, common_name)) +
    annotate("rect", xmin = 0,    xmax = 0.20, ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Overfished"]]) +
    annotate("rect", xmin = 0.20, xmax = 0.40, ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Cautionary"]]) +
    annotate("rect", xmin = 0.40, xmax = Inf,  ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Healthy"]]) +
    geom_vline(xintercept = c(0.20, 0.40), colour = "white", linewidth = 0.8) +
    spr_err +
    geom_point(size = 4.2, shape = 21, fill = BREAM$gold, colour = BREAM$ink, stroke = 0.8) +
    spr_sex_layer +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(0)) +
    coord_cartesian(xlim = c(0, 1.04), clip = "off") +   # room at SPR=1 so markers/CIs aren't sliced
    labs(title = "Spawning Potential Ratio (SPR): how much breeding capacity remains?",
         subtitle = "Red < 0.20 overfished  ·  amber 0.20–0.40 cautionary  ·  green > 0.40 healthy",
         x = "SPR   (0 = none left  →  1 = pristine)", y = NULL,
         caption = "SPR = eggs from the fished stock ÷ eggs an unfished stock would make. Bars = bootstrap 95% confidence\n(wide = few fish + borrowed growth values). Read as INDICATIVE: LBSPR assumes a textbook fishery, so where\nthe catch lacks large fish (e.g. black seabream) SPR can be pushed down artificially.\nGold = standard LBSPR; white diamond = sex-structured SPR (eggs from females only; binding for protogynous species).")
  
  # (b) fishing pressure F/M, on a log axis (clipped) with the F = M reference line
  asd_fm <- asd |> filter(is.finite(FM)) |> mutate(common_name = fct_reorder(common_name, FM))
  if (nrow(asd_fm)) {
    fm_err <- if (all(c("fm_lo", "fm_hi") %in% names(asd_fm))) {
      geom_segment(data = asd_fm,
                   aes(x = pmax(fm_lo, 0.05), xend = fm_hi, y = common_name, yend = common_name),
                   colour = BREAM$ink, linewidth = 0.7)
    } else NULL
    fishing_pressure <- ggplot(asd_fm, aes(FM, common_name)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
      fm_err +
      geom_point(size = 4.2, shape = 21, fill = BREAM$teal, colour = BREAM$ink, stroke = 0.8) +
      scale_x_log10(labels = scales::label_number(accuracy = 0.1)) +
      coord_cartesian(xlim = c(0.05, 20)) +
      labs(title = "Fishing pressure relative to natural mortality (F/M)",
           subtitle = "Below 1 = fishing removes fewer fish than nature does  ·  above 1 = fishing dominates",
           x = "F / M   (log scale)", y = NULL,
           caption = "Dashed line = F equals M. Bars = bootstrap 95% confidence and are very wide at these sample sizes — treat\nF/M as a rough signal only. The black seabream estimate is especially unstable (its catch lacks large fish).")
  }
}

# =====================================================================
# 10c. CATCH-CURVE CROSS-CHECK  (TropFishR length-converted catch curve)
# =====================================================================
# An INDEPENDENT read on total mortality Z, hence on F/M, from the descending limb of the
# length-converted catch curve (Pauly). LBSPR infers F/M from the whole length composition
# under equilibrium + logistic selectivity; the catch curve uses only the fully-selected
# right limb, so agreement is corroborating and divergence is diagnostic (commonly dome-shaped
# selectivity or non-equilibrium). It does NOT alter the LBSPR result; it sits beside it.
# F = Z - M with the workbook M; Linf/K are on the FL scale, matching the catch lengths.
# t0 is an additive constant on relative age and does not affect the slope Z, so t0 = 0.
cc_one <- function(L, lh, fm_lbspr) {
  if (!have_tropfishr) return(NULL)
  L <- L[is.finite(L)]; if (length(L) < MIN_N_ASSESS) return(NULL)
  br     <- seq(0, ceiling(max(c(L, lh$Linf))) + BINWIDTH, by = BINWIDTH)
  mids   <- br[-length(br)] + BINWIDTH / 2
  counts <- as.numeric(table(cut(L, breaks = br, right = FALSE)))
  peak <- which.max(counts)
  hi   <- suppressWarnings(max(which(counts >= 2)))
  if (!is.finite(hi) || hi <= peak + 1) return(NULL)        # need a usable descending limb
  lfq <- list(midLengths = mids, catch = counts, Linf = lh$Linf, K = lh$K, t0 = 0)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  cc <- tryCatch(TropFishR::catchCurve(lfq, reg_int = c(peak + 1, hi)),
                 error = function(e) NULL)
  if (is.null(cc) || is.null(cc$Z)) return(NULL)
  Z <- as.numeric(cc$Z)[1]; M <- lh$M; Fv <- Z - M
  tibble(scientific_name = lh$scientific_name, n = length(L),
         Z_cc = round(Z, 3), M = round(M, 3), F_cc = round(Fv, 3),
         FM_cc = round(Fv / M, 3), FM_lbspr = round(fm_lbspr, 3),
         reg_lo_cm = mids[peak + 1], reg_hi_cm = mids[hi])
}

catch_curve_cc <- NULL
if (have_tropfishr && exists("assess_results") && nrow(assess_results)) {
  fm_lookup <- assess_results |> select(scientific_name, FM) |> tibble::deframe()
  cc_list <- lapply(assess_species, function(sp) {
    lh <- readiness |> filter(scientific_name == sp)
    fm <- if (sp %in% names(fm_lookup)) fm_lookup[[sp]] else NA_real_
    cc_one(measured_lengths(sp), lh, fm)
  })
  catch_curve_cc <- bind_rows(cc_list)
  if (nrow(catch_curve_cc)) {
    message("\n========= CATCH-CURVE CROSS-CHECK (TropFishR Z vs LBSPR F/M) =========")
    print(catch_curve_cc, n = Inf)
    message("FM_cc (catch curve) close to FM_lbspr = mutually corroborating;")
    message("FM_cc >> FM_lbspr often signals dome-shaped selectivity (the §Z2 stressor).")
    message("=====================================================================\n")
    write_csv(catch_curve_cc, here("results", "catch_curve_crosscheck.csv"))
  }
} else if (!have_tropfishr) {
  message("[skip] TropFishR not installed - catch-curve cross-check skipped.")
}

# =====================================================================
# 11. DIURNAL vs NOCTURNAL   (auto-activates once night records exist)
# =====================================================================
# Day/night is judged in LOCAL time against an explicit Gibraltar SUMMER daylight
# window. Gibraltar runs on CEST (UTC+2) in summer, when the sun rises ~07:00-07:45
# and sets ~21:00-21:45 across Jun-Aug (e.g. 21 Jun: 07:04 / 21:41). We bracket that
# with a small twilight margin: a catch counts as "night" only if its LOCAL clock
# time is before 06:30 or after 21:45. Timestamps are stored in UTC, so they are
# converted to Europe/Gibraltar FIRST - skipping that conversion is exactly what made
# early-morning daytime catches look like "night" in earlier runs. The 11-12 July 2026 overnight
# tournament (~20:30-06:30 local) is an exception the clock cannot see: its catches were batch-logged
# in one block near dawn, so their per-catch timestamps do NOT reflect the true catch time. Every
# record from those two local dates is therefore forced to "night"; all other catches take the normal
# window-based classification, so night = the tournament plus any genuinely overnight-logged catch.
DAY_START <- 6.5     # 06:30 local  (before earliest summer sunrise ~07:00, with margin)
DAY_END   <- 21.75   # 21:45 local  (after latest summer sunset  ~21:41, with margin)
OVERNIGHT_DATES <- as.Date(c("2026-07-11", "2026-07-12"))   # batch-logged overnight tournament -> night

classify_daynight <- function(ts_utc, tz = LOCAL_TZ) {
  local    <- lubridate::with_tz(ts_utc, tz)               # UTC instant -> local clock
  hour_dec <- lubridate::hour(local) + lubridate::minute(local) / 60
  ifelse(hour_dec >= DAY_START & hour_dec < DAY_END, "day", "night")
}
cleaned <- cleaned |>
  mutate(diel = classify_daynight(timestamp_utc),
         # Gibraltar-local calendar date: as.Date() ignores a POSIXct's tzone attribute, so the
         # date must be taken with tz = LOCAL_TZ (with_tz alone is discarded by as.Date) to match
         # the overnight-tournament dates in local time rather than UTC.
         diel = if_else(as.Date(timestamp_utc, tz = LOCAL_TZ) %in% OVERNIGHT_DATES,
                        "night", diel))

# Diagnostic: print the LOCAL clock-time span of the catch so the day/night split is
# unambiguous and easy to sanity-check against the summer daylight window above.
loc_times <- lubridate::with_tz(cleaned$timestamp_utc, LOCAL_TZ)
loc_times <- loc_times[!is.na(loc_times)]
n_night   <- sum(cleaned$diel == "night", na.rm = TRUE)
if (length(loc_times)) {
  message(sprintf(
    "[diel] local catch times span %s-%s (%s); daylight window %02d:%02d-%02d:%02d. Night records: %d of %d.",
    format(min(loc_times), "%H:%M"), format(max(loc_times), "%H:%M"), LOCAL_TZ,
    floor(DAY_START), round((DAY_START %% 1) * 60),
    floor(DAY_END),   round((DAY_END   %% 1) * 60),
    n_night, nrow(cleaned)))
}

daynight_size_plot <- NULL; daynight_tests <- NULL; daynight_composition <- NULL
if (n_night >= NIGHT_MIN_N) {
  message(sprintf("[diel] %d night records present - running day/night comparison.", n_night))
  daynight_composition <- cleaned |> count(diel, scientific_name) |>
    pivot_wider(names_from = diel, values_from = n, values_fill = 0)
  write_csv(daynight_composition, here("results", "daynight_composition.csv"))
  
  daynight_size_plot <- cleaned |>
    filter(length_source == "measured") |>
    ggplot(aes(diel, length_true_cm, fill = diel)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey30", alpha = 0.85) +
    geom_jitter(width = 0.14, alpha = 0.45, size = 1.7, colour = "grey25") +
    facet_wrap(~scientific_name, scales = "free_y",
               labeller = labeller(scientific_name = sci2com)) +
    scale_fill_manual(values = DIEL_COLS, guide = "none") +
    labs(title = "Do day and night catches differ in size?",
         subtitle = "Night = the 11-12 July overnight tournament (batch-logged; see caption)",
         x = NULL, y = "Fork length (cm)",
         caption = "Exploratory only: a single overnight event is confounded with date, moon phase and tide.")
  
  # per-species Mann-Whitney on length, where both periods have >= 5 fish
  daynight_tests <- cleaned |>
    filter(length_source == "measured") |>
    group_by(scientific_name) |>
    filter(sum(diel == "day") >= 5, sum(diel == "night") >= 5) |>
    summarise(
      p_value = tryCatch(wilcox.test(length_true_cm ~ diel)$p.value, error = function(e) NA_real_),
      median_day   = median(length_true_cm[diel == "day"]),
      median_night = median(length_true_cm[diel == "night"]),
      .groups = "drop"
    )
  if (nrow(daynight_tests)) write_csv(daynight_tests, here("results", "daynight_size_tests.csv"))
} else {
  message(sprintf("[diel] only %d night records (< %d) - day/night section will activate after the overnight tournament.",
                  n_night, NIGHT_MIN_N))
}

# =====================================================================
# 12. LENGTH-WEIGHT RELATIONSHIP + CONDITION   (auto-activates with weights)
# =====================================================================
lw_results <- NULL; lw_plot <- NULL
n_weight <- sum(is.finite(cleaned$weight_kg))
if (n_weight >= WEIGHT_MIN_N) {
  message(sprintf("[LWR] %d weights present - fitting length-weight relationships.", n_weight))
  lw_dat <- cleaned |> filter(is.finite(weight_kg), is.finite(length_true_cm))
  # per-species fits where n >= 10, plus a pooled fit
  fit_lw <- function(d, label) {
    if (nrow(d) < 10) return(NULL)
    m <- lm(log10(weight_kg) ~ log10(length_true_cm), data = d)
    tibble(group = label, n = nrow(d),
           a = 10^coef(m)[1], b = coef(m)[2], r2 = summary(m)$r.squared)
  }
  lw_results <- bind_rows(
    lw_dat |> group_split(scientific_name) |>
      map_dfr(\(d) fit_lw(d, unique(as.character(d$scientific_name)))),
    fit_lw(lw_dat, "ALL Sparidae pooled")
  )
  if (!is.null(lw_results)) write_csv(lw_results, here("results", "length_weight.csv"))
  lw_plot <- lw_dat |>
    ggplot(aes(length_true_cm, weight_kg)) +
    geom_point(alpha = 0.5, size = 2.2, colour = BREAM$teal) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = BREAM$gold_dark, linewidth = 0.8) +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~scientific_name, scales = "free",
               labeller = labeller(scientific_name = sci2com)) +
    labs(title = "Length-weight relationship",
         subtitle = "Activates once enough weights are recorded",
         x = "Fork length (cm, log)", y = "Weight (kg, log)",
         caption = "A straight line on log-log axes; a slope near 3 indicates near-isometric growth.")
} else {
  message(sprintf("[LWR] only %d weights (< %d) - length-weight section will activate as weights accumulate.",
                  n_weight, WEIGHT_MIN_N))
}

# =====================================================================
# 13. TAGGING  — live effort summary + (COMMENTED) movement analysis
# =====================================================================
tag_summary <- cleaned |>
  summarise(
    n_total          = n(),
    n_newly_tagged   = sum(was_tagged, na.rm = TRUE),
    n_existing_tag   = sum(had_existing_tag, na.rm = TRUE),
    n_recaptures     = sum(!is.na(recapture_of_fish_id)),
    tag_return_rate  = ifelse(sum(was_tagged, na.rm = TRUE) > 0,
                              sum(!is.na(recapture_of_fish_id)) / sum(was_tagged, na.rm = TRUE), NA)
  )
print(tag_summary)
write_csv(tag_summary, here("results", "tagging_summary.csv"))

# -------- MOVEMENT / RECAPTURE ANALYSIS (kept commented until returns) -------
# Un-comment this block once recapture_of_fish_id / days_at_liberty start to fill.
# It needs the geosphere package for great-circle displacement.
#
# recaps <- cleaned |>
#   filter(had_existing_tag | !is.na(recapture_of_fish_id))
#
# if (nrow(recaps) > 0 && have_geo) {
#   # join each recapture back to its original tagging record by tag id
#   originals <- cleaned |>
#     filter(was_tagged) |>
#     transmute(new_tag_id,
#               orig_lat = latitude, orig_lon = longitude,
#               orig_date = date, orig_len = length_true_cm)
#   movement <- recaps |>
#     left_join(originals, by = c("existing_tag_id" = "new_tag_id")) |>
#     mutate(
#       days_at_liberty = as.numeric(date - orig_date),
#       displacement_km = geosphere::distHaversine(
#         cbind(orig_lon, orig_lat), cbind(longitude, latitude)) / 1000,
#       growth_increment_cm = length_true_cm - orig_len
#     )
#   write_csv(movement, here("results", "movement.csv"))
#
#   # movement vectors on the map
#   movement_map <- leaflet() |>
#     addProviderTiles("CartoDB.Positron") |>
#     addCircleMarkers(data = movement, lng = ~orig_lon, lat = ~orig_lat,
#                      radius = 4, color = "#2166ac", label = ~new_tag_id) |>
#     addCircleMarkers(data = movement, lng = ~longitude, lat = ~latitude,
#                      radius = 4, color = "#b2182b", label = ~existing_tag_id)
#   # growth increments give a direct check on the VBGF for recaptured fish
#   growth_check <- ggplot(movement, aes(days_at_liberty, growth_increment_cm)) +
#     geom_point() + labs(title = "Growth increment vs days at liberty") + theme_bream()
# }

# ===================================================================
# §Z  OPERATING-MODEL SIMULATION  (validation of the SPR amendment)
#     Self-contained. Depends on spr_sex_structured() from §10b.
#     Does NOT run on source() unless RUN_OM_SIM is TRUE (set at top).
# ===================================================================

# ---- Stage 1: operating model — TRUE state of a known stock ----
om_truth <- function(F_true, Linf = 35, K = 0.20, t0 = 0, M = 0.30,
                     L50 = 18, L95 = 22, a_lw = 0.01, b_lw = 3.0,
                     sex = "protogyny", LD50 = 25, LD95 = 28,
                     sel50 = 16, sel95 = 21, amax = NULL, da = 0.1,
                     # ---- stress-test hooks (defaults reproduce the baseline) ----
                     sel_type = "logistic", desc50 = NULL, desc95 = NULL,  # dome selectivity
                     kappa = 0) {                                          # compensatory sex change
  if (is.null(amax)) amax <- ceiling(-log(0.01) / M) + t0
  ages <- seq(0, amax, by = da)
  La <- Linf * (1 - exp(-K * (ages - t0))); La[La < 0] <- 0
  Wa <- a_lw * La^b_lw
  logistic <- function(L, p50, p95) 1 / (1 + exp(-log(19) * (L - p50) / (p95 - p50)))

  # selectivity: logistic (default) OR dome-shaped (stressor 1)
  sel_type <- match.arg(sel_type, c("logistic", "dome"))
  if (sel_type == "dome") {
    if (is.null(desc50)) desc50 <- sel95 + 0.30 * Linf   # location of the descending limb
    if (is.null(desc95)) desc95 <- desc50 + 0.15 * Linf
    asc  <- logistic(La, sel50, sel95)
    desc <- 1 / (1 + exp( log(19) * (La - desc50) / (desc95 - desc50)))  # +sign -> falls 1 to 0
    selA <- asc * desc                                   # the dome
  } else {
    selA <- logistic(La, sel50, sel95)
  }

  # sex-at-length: optionally compensatory (stressor 2). kappa = 0 -> fixed unfished
  # schedule (baseline). kappa > 0 -> the FISHED stock changes sex earlier; the UNFISHED
  # reference keeps the baseline schedule, so SPR compares like with like.
  u     <- F_true / (F_true + M)            # exploitation pressure, in [0,1)
  LD50e <- LD50 * (1 - kappa * u)           # earlier sex change under fishing
  LD95e <- LD50e + (LD95 - LD50)            # keep the ogive width
  mat   <- logistic(La, L50, L95)
  psi_of <- function(d50, d95) switch(sex,
                                      "gonochore"   = ,
                                      "rudimentary" = rep(0.5, length(ages)),
                                      "protandry"   = logistic(La, d50, d95),     # ASCENDING in length
                                      "protogyny"   = 1 - logistic(La, d50, d95), # DESCENDING in length
                                      stop("sex must be gonochore, rudimentary, protandry, or protogyny"))
  psi_f_F <- psi_of(LD50e, LD95e)           # fished schedule (shifted if kappa > 0)
  psi_f_0 <- psi_of(LD50,  LD95)            # unfished baseline schedule

  Z  <- M + F_true * selA
  Z0 <- M + 0 * selA                               # unfished Z as a per-age vector
  surv_F <- c(1, head(exp(-cumsum(Z  * da)), -1))  # fished survivorship (age 0 = 1)
  surv_0 <- c(1, head(exp(-cumsum(Z0 * da)), -1))  # unfished
  SPR_true <- sum(surv_F * mat * psi_f_F * Wa) / sum(surv_0 * mat * psi_f_0 * Wa)
  list(SPR_true = SPR_true, La = La, catch_at_age = surv_F * selA, Linf = Linf)
}

# ---- Stage 2: observation — simulate n sampled lengths ----
make_length_sample <- function(om, n, CVlen = 0.10, binwidth = 1) {
  Lbins <- seq(0, 1.3 * om$Linf, by = binwidth)
  Lmids <- Lbins[-length(Lbins)] + binwidth / 2
  probLA <- sapply(Lmids, function(L) dnorm(L, om$La, pmax(CVlen * om$La, 1e-3)))
  pLen <- as.numeric(om$catch_at_age %*% probLA)
  if (!any(is.finite(pLen)) || sum(pLen, na.rm = TRUE) <= 0) return(NULL)  # degenerate stock -> no sample
  pLen <- pLen / sum(pLen)
  counts <- as.numeric(rmultinom(1, n, pLen))      # THE sampling draw
  list(Lmids = Lmids, counts = counts, binwidth = binwidth, pLen = pLen)
}

# ---- Stage 3: estimation — fit from lengths ONLY, then apply the amendment ----
# CHANGED (validated == deployed): the estimator now routes through the SAME
# LBSPR::LBSPRfit path the real assessment uses (via fit_lbspr_one), instead of a
# bespoke single-Linf optim. The operating model therefore tests the estimator you
# actually report from — a growth-type-group fit carrying CVLinf — with the sex
# reweight applied on top exactly as in production. The assumed life-history is still
# NOT the OM truth; only the fitting machinery is now shared.
#   - bin counts are expanded to their midpoints and handed to fit_lbspr_one, which
#     re-bins them on the same grid, so the round-trip adds and loses no information;
#   - cvlinf is the assumed length CV given to BOTH the fit and the reweight, and
#     defaults to 0.10 to match the observation CV baked into make_length_sample();
#   - fit_lbspr_one() returns NULL on non-convergence or a package error, which the
#     replicate loops already treat as a dropped replicate.
# NOTE: est also carries SPR_fem (the female-only egg SPR). It is the quantity that
# should recover the OM's psi_f-weighted SPR_true; run_scenario()/run_stress() below
# currently summarise SPR_std and SPR_bind only. Carry SPR_fem through if you want the
# validation table to show the unbiased female estimator beside the naive SPR_std and
# the precautionary SPR_bind (recommended when reporting the validation).
estimate_from_sample <- function(samp, Linf, MK, L50, L95,
                                 sex, LD50, LD95, b_lw = 3, cvlinf = 0.10) {
  L <- rep(samp$Lmids, round(samp$counts))          # bin counts -> pseudo-lengths at midpoints
  if (length(L) < 5 || length(unique(L)) < 2) return(NULL)
  lh_om <- tibble::tibble(
    scientific_name = "OM_species",
    Linf = Linf, L50 = L50, L95 = L95, MK = MK,
    CVLinf = cvlinf, a = 0.01, b = b_lw, FecB = b_lw,
    sex_system = sex, LD50 = LD50, LD95 = LD95)
  fitres <- fit_lbspr_one(L, lh_om, binwidth = samp$binwidth, cvlinf = cvlinf)
  if (is.null(fitres)) return(NULL)                 # non-convergence / package error -> dropped
  sx <- spr_sex_structured(
    FM = fitres$FM, SL50 = fitres$SL50, SL95 = fitres$SL95,
    Linf = Linf, MK = MK, L50 = L50, L95 = L95, FecB = b_lw, CVLinf = cvlinf,
    sex_system = sex, LD50 = LD50, LD95 = LD95,
    anchor_SPR = fitres$SPR, MaleExp = b_lw, binwidth = samp$binwidth)
  list(FMhat = fitres$FM, SPR_std = sx$SPR_gono_check,
       SPR_fem = sx$SPR_fem, SPR_bind = sx$SPR_bind)
}

# ---- Stage 4: Monte-Carlo over replicates for ONE scenario ----
run_scenario <- function(F_true, n, R = 500, om_args = list(), est_args = list()) {
  om <- do.call(om_truth, c(list(F_true = F_true), om_args))
  zone <- function(s) cut(s, c(-Inf, 0.20, 0.40, Inf),
                          labels = c("overfished", "cautionary", "healthy"))
  true_zone <- zone(om$SPR_true)
  out <- replicate(R, {
    samp <- make_length_sample(om, n)
    if (is.null(samp) || sum(samp$counts > 0) < 3) return(c(NA, NA, NA))   # NULL -> degenerate stock
    est <- tryCatch(do.call(estimate_from_sample, c(list(samp = samp), est_args)),
                    error = function(e) NULL)
    if (is.null(est)) return(c(NA, NA, NA))        # non-convergence / error -> dropped replicate
    c(est$SPR_std, est$SPR_fem, est$SPR_bind)
  })
  out <- t(out); colnames(out) <- c("SPR_std", "SPR_fem", "SPR_bind")
  out <- out[complete.cases(out), , drop = FALSE]
  n_conv <- nrow(out)                              # replicates that converged AND gave finite SPR
  if (n_conv == 0) {                               # nothing usable -> all-NA row, still report n & R
    na6 <- c(median = NA, lo = NA, hi = NA, cv = NA, bias = NA, p_correct = NA)
    return(data.frame(n = n, SPR_true = om$SPR_true, R = R, n_converged = 0L,
                      metric = c("SPR_std", "SPR_fem", "SPR_bind"),
                      rbind(na6, na6, na6), row.names = NULL))
  }
  summ <- function(x) c(median = median(x), lo = unname(quantile(x, .1)),
                        hi = unname(quantile(x, .9)), cv = sd(x) / mean(x),
                        bias = median(x) - om$SPR_true,
                        p_correct = mean(zone(x) == true_zone))
  data.frame(n = n, SPR_true = om$SPR_true, R = R, n_converged = n_conv,
             metric = c("SPR_std", "SPR_fem", "SPR_bind"),
             rbind(summ(out[, "SPR_std"]), summ(out[, "SPR_fem"]), summ(out[, "SPR_bind"])), row.names = NULL)
}

# ---- driver (expensive; gated) ----
if (RUN_OM_SIM) {
  set.seed(42)
  est_args <- list(Linf = 35, MK = 1.5, L50 = 18, L95 = 22,
                   sex = "protogyny", LD50 = 25, LD95 = 28)
  om_args  <- list(Linf = 35, K = 0.20, M = 0.30, L50 = 18, L95 = 22,
                   sex = "protogyny", LD50 = 25, LD95 = 28, sel50 = 16, sel95 = 21)
  grid <- expand.grid(n = c(10, 20, 30, 50, 75, 100, 200, 500), F_true = 0.45)
  res  <- do.call(rbind, Map(function(n, f)
    run_scenario(f, n, R = 500, om_args = om_args, est_args = est_args),
    grid$n, grid$F_true))
  readr::write_csv(res, here::here("results", "om_simulation.csv"))
  print(res)
}

# =====================================================================
# Histogram of Tagged Fish by Species
# =====================================================================

tagged <- cleaned |>
  filter(was_tagged) |>
  count(scientific_name) |>
  mutate(common_name = sci2com[scientific_name]) |>
  arrange(desc(n)) |>
  mutate(common_name = fct_reorder(common_name, n))
if (interactive()) view(tagged)

tagged_histogram <- ggplot(tagged, aes(n, common_name)) +
  geom_col(fill = BREAM$teal) +
  scale_x_continuous(breaks = seq(0, max(tagged$n), by = 1)) +
  labs(title = "Number of tagged fish by species",
       x = "Count of tagged fish", y = NULL) +
  theme_bream()
# =====================================================================
# 14. SAVE FIGURES -> graphs and maps/   ;  tables already saved to results/
# =====================================================================
plots <- list(
  catchesspecies        = catchesspecies,
  lengthfrequency       = lengthfrequency,
  lengthspecies         = lengthspecies,
  lengthwater           = lengthwater,
  histogramlength       = histogramlength,
  catchdate             = catchdate,
  rankabundance         = rankabund,
  size_vs_reference     = indicator_plot,
  spr_status            = spr_status,
  fishing_pressure      = fishing_pressure,
  daynight_size         = daynight_size_plot,
  length_weight         = lw_plot,
  tagged_histogram      = tagged_histogram
)
plots <- compact(plots)   # drop any that are NULL because their section was skipped

# per-figure dimensions (inches); anything not listed falls back to the default
fig_size <- list(
  lengthfrequency       = c(11, 7.5), size_vs_reference  = c(10, 5.5),
  daynight_size         = c(11, 7.5), length_weight      = c(11, 7.5),
  catchesspecies        = c(8, 5.5),  lengthspecies      = c(8.5, 6),
  spr_status            = c(10, 5.5), fishing_pressure   = c(9, 4.5),
  rankabundance         = c(9.5, 6.5), histogramlength   = c(9, 5.5),
  catchdate             = c(9, 5),    lengthwater        = c(9, 6)
)
iwalk(plots, function(p, name) {
  d <- fig_size[[name]]; if (is.null(d)) d <- c(9, 6)
  tryCatch(
    ggsave(here("graphs and maps", paste0(name, ".png")), p,
           width = d[1], height = d[2], dpi = 300, bg = "white"),
    error = function(e) message(sprintf("[save] could not save %s: %s", name, conditionMessage(e))))
})

# Remove stale figures from earlier runs whose section is INACTIVE this time, so an
# old PNG (e.g. the pre-fix day/night plot, or the retired LBB/LIME comparison) does
# not linger in the folder and cause confusion.
conditional_figs <- c("daynight_size", "length_weight", "spr_status",
                      "fishing_pressure", "size_vs_reference", "assessment_comparison",
                      "diversity_rarefaction", "diversity_coverage")
for (f in setdiff(conditional_figs, names(plots))) {
  stale <- here("graphs and maps", paste0(f, ".png"))
  if (file.exists(stale)) {
    unlink(stale)
    message(sprintf("[clean] removed stale %s.png (not produced this run).", f))
  }
}

htmlwidgets::saveWidget(catch_map, here("graphs and maps", "leaflet.html"),
                        selfcontained = TRUE, title = "Bream catch map")

message("\nDone. Figures -> 'graphs and maps/' ; tables -> 'results/'.")