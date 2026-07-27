# 00_run_all.R: reproduce the full script pipeline in one R session.
# Run order matters: 03 uses functions defined in 02 and must follow it in the
# same session.
#
# The input folders below are supplied by you and are not distributed with this
# repository. See the README for the expected format.
library(here)

n_season <- length(list.files(here("all_records"),        pattern = "^bream_.*\\.csv$"))
n_hist   <- length(list.files(here("historical_tagging"), pattern = "\\.csv$"))

if (n_season == 0)
  stop("No files matching 'bream_*.csv' in all_records/. The assessment needs your own ",
       "catch records; see the README for the expected columns. To try the sex-change ",
       "amendment with no data at all, run:  cd sex-structured-lbspr && Rscript example.R",
       call. = FALSE)

# 01 merges the season catch with a historical tagging series. That series is
# optional: nothing else in this repository reads 01's output, so a missing
# historical_tagging/ folder skips the step instead of stopping the run.
if (n_hist > 0) {
  source(here("R", "01_combine_data.R"))   # -> combined_records/combined_dataset.csv
} else {
  message("Skipping 01_combine_data.R: no CSVs in historical_tagging/. ",
          "The assessment below does not depend on it.")
}

source(here("R", "02_analysis.R"))       # assessment, indicators, composition -> results/ + graphs and maps/
source(here("R", "03_stress_test.R"))    # misspecification battery + figures (uses 02's functions; slow)
source(here("R", "04_catalog.R"))        # photographic catalogue -> catalog/ (needs all_records/photos/)
message("\nPipeline complete. Outputs written to results/, graphs and maps/, and catalog/.")
