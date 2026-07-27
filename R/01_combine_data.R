# =====================================================================
# BUILD combined_dataset.csv  —  own 2026 season + historical tagging records
# =====================================================================
# Merges the live season catch exports (all_records/bream_*.csv) with the historical
# tagging records (historical_tagging/*.csv) into one cleaned, length-reconciled file
# written to combined_records/combined_dataset.csv. That file is produced for onward
# analysis; nothing else in this repository reads it.
#
# WORKFLOW: whenever new season exports are added to all_records/, re-run this script.
# all_records/ is READ ONLY here: the script never writes to it or edits it.
#
# FOLDERS EXPECTED (under the project root; none is supplied with this repository):
#   all_records/         the bream_*.csv exports (and photos; only the CSVs are read)
#   historical_tagging/  the historical tagging CSV(s)
#   combined_records/    destination, created if absent; combined_dataset.csv written here
#
# LENGTH SCALE (important): the season lengths are fork length (FL); the historical
# records are total length (TL). Both are placed on the TL scale here (FL -> TL via
# / TL_FL_RATIO), so THE OUTPUT OF THIS SCRIPT IS ON THE TL SCALE. Any downstream
# analysis must convert the pooled series back to FL before assessing it, because the
# life-history parameters are applied on the FL scale. Do not change one side of that
# round-trip without the other.
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(here)
})

OWN_DIR     <- here("all_records")          # your 2026 bream_*.csv exports (READ ONLY)
HIST_DIR    <- here("historical_tagging")   # the historical tagging CSV(s)
OUT_DIR     <- here("combined_records")      # destination folder for the merged file
TL_FL_RATIO <- 0.92                          # FL = TL * ratio, for the length reconciliation

dir.create(OUT_DIR, showWarnings = FALSE)

# historical binomials -> the trinomials your app and life-history file use. Only the two
# subspecies cases need mapping; the other historical names already match the parameter file.
species_fix <- c(
  "Diplodus sargus"   = "Diplodus sargus sargus",
  "Diplodus cervinus" = "Diplodus cervinus cervinus"
)

# Read every CSV in a folder into one frame, tagged with its dataset label. Columns are read
# as text so files with slightly different schemas bind cleanly (the historical file carries a
# few columns the bream exports do not, e.g. location_name); proper types are restored once,
# on the pooled frame, below.
read_folder <- function(dir, pattern, label) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(files))
    stop("No files matching '", pattern, "' found in ", dir,
         " — check the folder name and contents.", call. = FALSE)
  files |>
    set_names(basename) |>
    map(\(f) read_csv(f, na = c("", "NA", "NR"), show_col_types = FALSE,
                      col_types = cols(.default = col_character()))) |>
    list_rbind(names_to = "source_file") |>
    mutate(dataset = label)
}

# Season: the SAME bream_*.csv pattern R/02_analysis.R uses, so the combined season rows
# are exactly the fish the base assessment sees.
own  <- read_folder(OWN_DIR,  "^bream_.*\\.csv$", "2026 season")
# Historical: every CSV in historical_tagging/ (drop additional historical files in later
# and they are picked up automatically).
hist <- read_folder(HIST_DIR, "\\.csv$",           "historical")

# Pool, restore column types, then reconcile species names and length scale.
pooled <- bind_rows(own, hist) |> clean_names()
pooled <- suppressMessages(type_convert(pooled))   # numeric / logical / date types on the pooled frame

combined <- pooled |>
  distinct(fish_id, .keep_all = TRUE) |>            # unique fish; a 2026 id wins any clash (bound first)
  mutate(
    scientific_name = coalesce(unname(species_fix[scientific_name]), scientific_name),
    length_true_cm  = if_else(length_type == "FL" & !is.na(length_true_cm),
                              length_true_cm / TL_FL_RATIO, length_true_cm),   # FL -> TL
    length_type     = if_else(!is.na(length_true_cm), "TL", length_type)
  )

# ---- cross-dataset de-duplication -----------------------------------
# The historical set and the 2026 exports were logged concurrently, so a fish caught by both carries a
# different fish_id in each (HUR-.../COM-... vs GOG-...) and slips through distinct(fish_id) above.
# Records sharing a TAG NUMBER on the SAME calendar day are one capture event logged twice; a shared tag
# on DIFFERENT days is a genuine tag -> recapture pair and is kept intact.
# KEEP RULE: default to the richer 2026 row (it carries a geolocation). EXCEPTION: a recapture of a tag
# first applied in the historical set is kept on the HISTORICAL row -- our lengths are fork length and the
# historical lengths are total (tail) length, so a growth increment is only valid with both ends measured
# the same way.
tag_origin <- combined |>
  dplyr::filter(was_tagged %in% TRUE, !is.na(new_tag_id)) |>
  group_by(.otag = as.character(new_tag_id)) |>
  summarise(orig_historical = any(dataset == "historical"), .groups = "drop")

combined <- combined |>
  mutate(.tag   = coalesce(as.character(new_tag_id), as.character(existing_tag_id)),
         .geo   = !is.na(latitude) & !is.na(longitude),
         .recap = had_existing_tag %in% TRUE) |>
  left_join(tag_origin, by = c(".tag" = ".otag")) |>
  mutate(.rank = if_else(.recap & orig_historical %in% TRUE,
                         2 * (dataset == "historical") + .geo,      # recapture of a historical tag -> keep historical
                         2 * .geo + (dataset == "2026 season")))    # else keep the richer own (geo) row

combined <- bind_rows(
  dplyr::filter(combined, is.na(.tag)),                            # untagged rows: all kept
  combined |> dplyr::filter(!is.na(.tag)) |>                       # tagged rows: one per tag-day capture event
    group_by(.tag, date) |> slice_max(.rank, n = 1, with_ties = FALSE) |> ungroup()
) |> dplyr::select(-.tag, -.geo, -.recap, -orig_historical, -.rank)

# ---- recapture linkage by tag number --------------------------------
# A record bearing an existing tag is a recapture of the record that first applied that tag number; link
# them and derive movement/growth from the pair. Displacement needs coordinates at BOTH ends, so a
# historical tagging point without a geolocation gives NA (intended, and flagged as such downstream).
hav_km <- function(lon1, lat1, lon2, lat2) {                       # great-circle km, NA-safe
  r <- 6371; d <- pi / 180
  a <- sin((lat2 - lat1) * d / 2)^2 + cos(lat1 * d) * cos(lat2 * d) * sin((lon2 - lon1) * d / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}
originals <- combined |>
  dplyr::filter(was_tagged %in% TRUE, !is.na(new_tag_id)) |>
  transmute(.mtag = as.character(new_tag_id), orig_fish_id = fish_id, orig_date = date,
            orig_lat = latitude, orig_lon = longitude, orig_len = length_true_cm)
combined <- combined |>
  mutate(.mtag = if_else(had_existing_tag %in% TRUE & !is.na(existing_tag_id),
                         as.character(existing_tag_id), NA_character_)) |>
  left_join(originals, by = ".mtag") |>
  mutate(
    recapture_of_fish_id = coalesce(orig_fish_id, recapture_of_fish_id),
    days_at_liberty  = if_else(!is.na(orig_date), as.numeric(as.Date(date) - as.Date(orig_date)), days_at_liberty),
    growth_length_cm = if_else(!is.na(orig_len) & !is.na(length_true_cm), length_true_cm - orig_len, growth_length_cm),
    displacement_km  = if_else(!is.na(orig_lat) & !is.na(latitude), hav_km(orig_lon, orig_lat, longitude, latitude), NA_real_)
  ) |> dplyr::select(-.mtag, -orig_fish_id, -orig_date, -orig_lat, -orig_lon, -orig_len)

out_path <- file.path(OUT_DIR, "combined_dataset.csv")
write_csv(combined, out_path)

message(sprintf("Wrote %d rows to %s  (%d from 2026 season, %d historical).",
                nrow(combined), out_path,
                sum(combined$dataset == "2026 season", na.rm = TRUE),
                sum(combined$dataset == "historical",  na.rm = TRUE)))
message(sprintf("  collapsed same-tag/same-day cross-dataset duplicates; %d recapture(s) linked by tag number.",
                sum(!is.na(combined$recapture_of_fish_id))))
