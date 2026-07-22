# 00_run_all.R: reproduce the full script pipeline in one R session.
# Run order matters: 03 uses functions defined in 02 and must follow it in the
# same session.
library(here)
source(here("R", "01_combine_data.R"))   # -> combined_records/combined_dataset.csv
source(here("R", "02_analysis.R"))       # assessment, indicators, diversity -> results/ + graphs and maps/
source(here("R", "03_stress_test.R"))    # misspecification battery + figures (uses 02's functions; slow)
source(here("R", "04_catalog.R"))        # photographic catalogue -> catalog/
message("\nPipeline complete. Outputs written to results/, graphs and maps/, and catalog/.")
