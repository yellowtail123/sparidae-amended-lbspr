# =====================================================================
# 00_run_all.R  —  reproduce the full script pipeline in one R session.
# =====================================================================
# RUN ORDER MATTERS.
#   01 writes combined_records/combined_dataset.csv, which 02 needs: the
#      cross-dataset de-duplication now gates the season analysis, so running 02
#      without it produces a season n one higher than the reported one and warns.
#   03 uses om_truth(), make_length_sample(), estimate_from_sample(),
#      spr_sex_structured(), BREAM, ZONE_FILL, theme_bream() and FINAL_RUN from
#      02, so it must follow 02 IN THE SAME SESSION. It is not standalone.
#   04 is standalone and slow (it re-encodes every photograph). Skip it with
#      RUN_CATALOG <- FALSE if you only want the analysis.
#   05 rebuilds the outputs that used to exist only inside the Rmd (the combined
#      assessment, the capture/recapture figure, the recovery figure and the
#      ogive-sensitivity figures). It needs 02's functions and two of 02's
#      results files, so it follows 02; it does not depend on 03 or 04.
#
# This run is long: 2000 bootstrap replicates and 2000 parameter Monte-Carlo
# draws per assessed species, a 1000-replicate recovery simulation, and a
# 1000-replicate stress battery at every severity level of four stressors.
# The preflight below exists so it cannot get forty minutes in and then stop on
# a missing file or an uninstalled package.
#
# Knitting the manuscript is a separate step; see README.
# =====================================================================

suppressPackageStartupMessages(library(here))

RUN_CATALOG <- TRUE    # FALSE to skip 04 (photo encoding; needs magick + jsonlite)

# ---- locate the scripts ---------------------------------------------
# Tolerates both layouts: scripts in R/ or scripts in the project root.
script_dir <- if (dir.exists(here("R"))) here("R") else here()
scripts <- c("01_combine_data.R", "02_analysis.R", "03_stress_test.R",
             if (isTRUE(RUN_CATALOG)) "04_catalog.R",
             "05_manuscript_outputs.R")
paths   <- file.path(script_dir, scripts)

# =====================================================================
# PREFLIGHT — fail now, loudly, rather than in forty minutes
# =====================================================================
fail <- character(0)
note <- function(...) fail <<- c(fail, sprintf(...))

missing_scripts <- scripts[!file.exists(paths)]
if (length(missing_scripts))
  note("Script(s) not found in %s: %s", script_dir, paste(missing_scripts, collapse = ", "))

# inputs the pipeline cannot run without
# NOTE: braces are required. At top level R closes an `if` at the end of its body, so a bare
# `else` on the next line is a parse error rather than a continuation.
if (!dir.exists(here("all_records"))) {
  note("Folder all_records/ not found (the bream_*.csv exports live there).")
} else if (!length(list.files(here("all_records"), pattern = "^bream_.*\\.csv$"))) {
  note("No bream_*.csv files found in all_records/.")
}

if (!dir.exists(here("historical_tagging"))) {
  note("Folder historical_tagging/ not found; 01_combine_data.R stops without it.")
} else if (!length(list.files(here("historical_tagging"), pattern = "\\.csv$"))) {
  note("No CSV files found in historical_tagging/.")
}

if (!file.exists(here("Sparid_LBSPR_LifeHistory_trinomial_Rready.csv")))
  note(paste("Sparid_LBSPR_LifeHistory_trinomial_Rready.csv not found in the project root.",
             "It supplies every Linf, L50, L95, M, K, M/K, CVLinf, FecB, sex_system, LD50 and",
             "LD95, so 02_analysis.R cannot produce a single assessment number without it."))

# packages: hard requirements halt, soft ones only skip their own section
hard <- c("tidyverse", "here", "janitor", "skimr", "sf", "leaflet",
          "LBSPR", "vegan", "iNEXT", "htmlwidgets", "scales", "patchwork")
if (isTRUE(RUN_CATALOG)) hard <- c(hard, "magick", "jsonlite")
soft <- c("FishLife", "fishmethods", "TropFishR", "rfishbase", "geosphere")

missing_hard <- hard[!vapply(hard, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_hard))
  note("Missing required package(s): %s\n    install.packages(c(%s))",
       paste(missing_hard, collapse = ", "),
       paste(sprintf('"%s"', missing_hard), collapse = ", "))

missing_soft <- soft[!vapply(soft, requireNamespace, logical(1), quietly = TRUE)]

if (length(fail)) {
  stop("Preflight failed. Fix these before running:\n  - ",
       paste(fail, collapse = "\n  - "), call. = FALSE)
}
if (length(missing_soft))
  message("[preflight] Optional package(s) absent, their sections will be skipped: ",
          paste(missing_soft, collapse = ", "),
          if ("FishLife" %in% missing_soft)
            "\n            FishLife is GitHub-only: pak::pkg_install(\"james-thorson/FishLife\")" else "")

# output folders (the scripts create their own, but make the intent explicit here)
for (d in c("results", "graphs and maps", "combined_records",
            if (isTRUE(RUN_CATALOG)) "catalog"))
  dir.create(here(d), showWarnings = FALSE, recursive = TRUE)

message("[preflight] passed. Scripts: ", script_dir)

# =====================================================================
# RUN
# =====================================================================
t0 <- Sys.time()
for (s in scripts) {
  message("\n", strrep("=", 70), "\n== ", s, "   (", format(Sys.time(), "%H:%M:%S"), ")\n",
          strrep("=", 70))
  ts <- Sys.time()
  source(file.path(script_dir, s), echo = FALSE)
  message(sprintf("== %s finished in %.1f min", s,
                  as.numeric(difftime(Sys.time(), ts, units = "mins"))))
}

message(sprintf("\nPipeline complete in %.1f min.",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message("Outputs: results/ (tables), 'graphs and maps/' (figures)",
        if (isTRUE(RUN_CATALOG)) ", catalog/ (photographic catalogue)" else "", ".")
message("Next: knit analysisfinal.Rmd (see README for switch settings).")

# Record what actually ran, so a figure or table can be traced to a session.
writeLines(c(paste("run completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             capture.output(sessionInfo())),
           here("results", "sessionInfo.txt"))

# rmarkdown::render(here("analysisfinal.Rmd"))   # uncomment to render the manuscript too
