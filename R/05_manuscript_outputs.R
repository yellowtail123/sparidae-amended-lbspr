# =====================================================================
# 05_manuscript_outputs.R  —  the outputs that existed only inside the Rmd
# =====================================================================
# WHY THIS FILE EXISTS
#   Four manuscript items were generated only by analysisfinal.Rmd and had no
#   scripted equivalent: the combined-dataset assessment (Table 22 and its
#   figure), the capture/recapture figure, the closed-loop recovery figure, and
#   the sex-change-ogive sensitivity figure. Anything that exists in one place
#   only cannot be checked against anything, and the two paths had already
#   drifted once. Everything the manuscript reports is now reproducible from the
#   scripts alone, and the Rmd becomes a presentation layer rather than a second
#   analysis.
#
#   Each output below also carries a fix that the Rmd version did not have; the
#   FIX notes say which.
#
# NOT SELF-CONTAINED. Source 02_analysis.R first, in the same session: this file
# uses fit_lbspr_one(), spr_sex_structured(), life_history, sp_lookup, BREAM,
# ZONE_FILL, theme_bream() and the run-control constants. It also reads two files
# that 02 writes, results/om_simulation.csv (RUN_OM_SIM must have been TRUE) and
# results/ld_sensitivity.csv, rather than recomputing them. It does NOT depend on
# 03_stress_test.R, so it can run at any point after 02; it is numbered last only
# because it was added last.
#
# WRITES
#   results/combined_assessment_results.csv     (Table 22)
#   results/combined_tag_frequency.csv          (Table 23)
#   graphs and maps/combined_spr_status.png     (Figure 18)
#   graphs and maps/combined_recapture.png      (Figure 19)
#   graphs and maps/om_recovery.png             (Figure 20)
#   graphs and maps/om_zone_recovery.png        (companion to Figure 20)
#   graphs and maps/ld_sensitivity_female.png   (Figure 17)
#   graphs and maps/ld_sensitivity_male.png     (companion to Figure 17)
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
})

# ---- dependency guard: fail with a useful message, not a missing-object error ----
.needed <- c("fit_lbspr_one", "spr_sex_structured", "life_history", "BREAM",
             "ZONE_FILL", "theme_bream", "MIN_N_ASSESS", "FM_CEIL", "BINWIDTH",
             "TL_FL_RATIO", "COMBINED_CSV", "SEASON_LABEL")
.absent <- .needed[!vapply(.needed, exists, logical(1))]
if (length(.absent))
  stop("05_manuscript_outputs.R needs 02_analysis.R sourced first in this session.\n",
       "  Not found: ", paste(.absent, collapse = ", "), call. = FALSE)

if (!exists("bare_fig")) bare_fig <- function(p) p
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!have_patchwork)
  warning("[05] patchwork not installed - panel figures will not be assembled; the ",
          "individual panels are written instead. install.packages(\"patchwork\")",
          call. = FALSE)

FIG_DIR <- here("graphs and maps")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)

zone_of <- function(s) as.character(cut(s, c(-Inf, 0.20, 0.40, Inf),
                                        labels = c("overfished", "cautionary", "healthy")))
# Zone shading is for SPR and SPR_fem ONLY. It is never applied to a male-biomass
# quantity, which has no reference point; see the note in Section 10b of 02_analysis.R.
zone_bg_x <- list(
  annotate("rect", xmin = 0,    xmax = 0.20, ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Overfished"]]),
  annotate("rect", xmin = 0.20, xmax = 0.40, ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Cautionary"]]),
  annotate("rect", xmin = 0.40, xmax = Inf,  ymin = -Inf, ymax = Inf, fill = ZONE_FILL[["Healthy"]]),
  geom_vline(xintercept = c(0.20, 0.40), colour = "white", linewidth = 0.8)
)

# =====================================================================
# 1. COMBINED-DATASET ASSESSMENT   (Table 22, Figure 18)
# =====================================================================
# FIX (referee F2): the Rmd version never computed the reliability flag for this
# assessment, so Pagellus acarne appeared here with F/M = 12.11 and the
# `Reliable = no` it carries in Table 18 silently stripped off, and Diplodus
# sargus returned SPR = 1.00 with F/M = 0.00 -- both terms hard against their
# bounds -- presented as an estimate. Both flags are computed and reported now.
combined_assess_results <- NULL
combined_spr_status <- NULL
combined_recapture  <- NULL

if (!file.exists(COMBINED_CSV)) {
  warning("[05] ", basename(COMBINED_CSV), " not found - combined outputs skipped. ",
          "Run 01_combine_data.R first.", call. = FALSE)
} else {
  # Lengths in the combined file are on the TL scale (01_combine_data.R put them there so the
  # two series share one basis); the assessment runs on FL, so convert back exactly as the
  # life-history parameters were converted.
  combined <- read_csv(COMBINED_CSV, na = c("", "NA", "NR"), show_col_types = FALSE) |>
    clean_names() |>
    mutate(
      length_true_cm = if_else(length_type %in% "TL" & !is.na(length_true_cm),
                               length_true_cm * TL_FL_RATIO, length_true_cm),
      length_type    = if_else(!is.na(length_true_cm), "FL", length_type),
      dataset        = factor(dataset, levels = c(SEASON_LABEL, "historical")),
      was_tagged     = as.logical(was_tagged)
    )
  
  combined_sp_lookup <- combined |>
    distinct(scientific_name, common_name) |>
    filter(!is.na(scientific_name))
  
  combined_lengths <- function(sp) {
    v <- combined |>
      filter(scientific_name == sp, length_type == "FL", !is.na(length_true_cm)) |>
      pull(length_true_cm)
    v[is.finite(v)]
  }
  
  combined_n <- combined |>
    filter(length_type == "FL", !is.na(length_true_cm)) |>
    count(scientific_name, name = "n_measured")
  
  combined_ready <- life_history |>            # SAME FL-scaled parameters as the base assessment
    inner_join(combined_n, by = "scientific_name") |>
    filter(is.finite(Linf), is.finite(L50), n_measured >= MIN_N_ASSESS)
  
  pts <- list()
  for (sp in combined_ready$scientific_name) {
    lh <- combined_ready |> filter(scientific_name == sp)
    L  <- combined_lengths(sp)
    r  <- fit_lbspr_one(L, lh)                 # the SAME engine as the base assessment
    if (is.null(r)) next
    # how many of the pooled fish are historical, so the reader can see what pooling added
    n_hist <- combined |>
      filter(scientific_name == sp, dataset == "historical",
             length_type == "FL", !is.na(length_true_cm)) |>
      nrow()
    pts[[sp]] <- r |>
      mutate(scientific_name = sp, n = length(L), n_historical = n_hist,
             reliable = FM <= FM_CEIL,
             # A length-based fit that returns SPR at 1 and F/M at 0 has not estimated
             # anything; it has run to the edge of its parameter space. That is a different
             # failure from an implausibly high F/M and needs its own flag.
             at_bound = (SPR >= 0.999) | (FM <= 1e-6) | (FM >= 1e3),
             .before = 1)
  }
  combined_assess_results <- bind_rows(pts)
  
  if (nrow(combined_assess_results)) {
    combined_assess_results <- combined_assess_results |>
      mutate(status_note = case_when(
        at_bound  ~ "boundary fit: both terms at their limits, not an estimate",
        !reliable ~ sprintf("flagged: estimated F/M above the reliability ceiling (%s)", FM_CEIL),
        TRUE      ~ "reliable"
      ))
    write_csv(combined_assess_results, here("results", "combined_assessment_results.csv"))
    
    message("\n============ COMBINED ASSESSMENT (2026 season + historical) ============")
    print(combined_assess_results |>
            transmute(scientific_name, n, n_historical,
                      SPR = round(SPR, 3),
                      SPR_female = if ("SPR_fem" %in% names(combined_assess_results)) round(SPR_fem, 3) else NA,
                      FM = round(FM, 2), SL50 = round(SL50, 2), SL95 = round(SL95, 2),
                      status_note), n = Inf)
    n_bad <- sum(combined_assess_results$at_bound | !combined_assess_results$reliable)
    if (n_bad)
      message(sprintf("%d of %d pooled fits are flagged. The Methods (Section 2.4) describe the ",
                      n_bad, nrow(combined_assess_results)),
              "historical tagging as unsystematic with little usable baseline, so a flagged\n",
              "pooled fit is a result about poolability, not a stock status.")
    message("=======================================================================\n")
    
    # ---- Figure 18 -----------------------------------------------------
    has_sex <- "SPR_fem" %in% names(combined_assess_results)
    asd <- combined_assess_results |>
      left_join(combined_sp_lookup, by = "scientific_name") |>
      mutate(status = if (has_sex) SPR_fem else SPR,
             common_name = fct_reorder(common_name, status),
             flag = factor(if_else(at_bound | !reliable,
                                   "flagged: not interpreted", "reliable"),
                           levels = c("reliable", "flagged: not interpreted")))
    combined_spr_status <- ggplot(asd, aes(status, common_name)) + zone_bg_x +
      { if (has_sex) geom_segment(aes(x = SPR, xend = status, y = common_name, yend = common_name),
                                  colour = BREAM$ink, linewidth = 0.4, linetype = "dotted") } +
      { if (has_sex) geom_point(aes(x = SPR), size = 3.2, shape = 21, fill = BREAM$gold,
                                colour = BREAM$ink, stroke = 0.7) } +
      geom_point(aes(shape = flag, alpha = flag), size = 4.4, fill = BREAM$teal,
                 colour = BREAM$ink, stroke = 0.9) +
      scale_shape_manual(values = c("reliable" = 21, "flagged: not interpreted" = 24),
                         name = NULL, drop = FALSE) +
      scale_alpha_manual(values = c("reliable" = 1, "flagged: not interpreted" = 0.45),
                         guide = "none") +
      scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(0)) +
      coord_cartesian(xlim = c(0, 1.04), clip = "off") +
      labs(title = "Spawning potential ratio on the combined dataset",
           subtitle = "2026 season pooled with the historical tagging records, on the fork-length scale",
           x = "Spawning potential ratio (SPR)", y = NULL,
           caption = paste0(
             "Teal = the reported status", if (has_sex) " (sex-structured female, egg-based)" else "",
             "; gold = standard LBSPR. Open triangles are fits that are NOT interpreted:\neither the estimated F/M exceeds the reliability ceiling, or the fit has run to the edge of its parameter space (SPR at 1 with F/M at 0).\nPooling mixes two collection regimes and the historical records were logged without a standardised protocol, so this is a check on\nwhether the two series can be pooled at all, not a longer-run stock status."))
    ggsave(file.path(FIG_DIR, "combined_spr_status.png"), bare_fig(combined_spr_status),
           width = 10, height = 5.0, dpi = 300, bg = "white")
  }
  
  # =====================================================================
  # 2. CAPTURE / RECAPTURE BY SPECIES AND DATASET   (Table 23, Figure 19)
  # =====================================================================
  # FIX: the Rmd version faceted the two datasets with free x-scales under one
  # shared "Number of fish" axis title, so ten fish and twenty fish drew at the
  # same physical bar length and the single cross-panel comparison the layout
  # invites was false. The x scale is fixed here; only the species axis is free,
  # because the two series hold different species.
  combined_tag_freq <- combined |>
    filter(!is.na(scientific_name)) |>
    group_by(dataset, scientific_name) |>
    summarise(n_records    = n(),
              n_tagged     = sum(was_tagged, na.rm = TRUE),
              n_recaptured = sum(!is.na(recapture_of_fish_id)),
              .groups = "drop") |>
    left_join(combined_sp_lookup, by = "scientific_name") |>
    arrange(dataset, desc(n_records))
  write_csv(combined_tag_freq, here("results", "combined_tag_frequency.csv"))
  
  cr_long <- combined_tag_freq |>
    select(dataset, common_name, n_tagged, n_recaptured) |>
    pivot_longer(c(n_tagged, n_recaptured), names_to = "event", values_to = "n") |>
    filter(n > 0) |>
    mutate(event = recode(event, n_tagged = "Tagged", n_recaptured = "Recaptured"),
           event = factor(event, levels = c("Tagged", "Recaptured")))
  
  if (nrow(cr_long)) {
    cr_long <- cr_long |> mutate(common_name = fct_reorder(common_name, n, .fun = max))
    combined_recapture <- ggplot(cr_long, aes(n, common_name, fill = event)) +
      geom_col(position = position_dodge2(preserve = "single", padding = 0.1),
               colour = BREAM$ink, linewidth = 0.25) +
      facet_wrap(~dataset, scales = "free_y") +          # y free (different species), x FIXED
      scale_fill_manual(values = c(Tagged = BREAM$gold, Recaptured = BREAM$teal), name = NULL) +
      scale_x_continuous(breaks = scales::breaks_pretty()) +
      labs(title = "Fish tagged and recaptured, by species and series",
           subtitle = "Both panels share one horizontal scale, so bar lengths are directly comparable",
           x = "Number of fish", y = NULL,
           caption = "Recaptures are rare in both series, so the bars are dominated by first-capture tagging. Historical positions are recorded as\nanecdotal location names rather than satellite fixes, so no movement distance is derived from them. Counts are in Table 23.")
    ggsave(file.path(FIG_DIR, "combined_recapture.png"), bare_fig(combined_recapture),
           width = 11, height = 6.0, dpi = 300, bg = "white")
  }
}

# =====================================================================
# 3. CLOSED-LOOP RECOVERY   (Figure 20 + zone-recovery companion)
# =====================================================================
# FIX (referee, Figure 20): the Rmd version drew ONE dashed truth line for TWO
# estimators that target DIFFERENT quantities, and captioned it as the female-egg
# truth while it actually sat at the binding truth. A reader saw the female series
# apparently overshooting several-fold when its real bias is small. Each estimator
# is drawn against its own truth here, in its own colour.
om <- if (file.exists(here("results", "om_simulation.csv")))
  read_csv(here("results", "om_simulation.csv"), show_col_types = FALSE) else NULL

if (is.null(om)) {
  message("[05] results/om_simulation.csv not found - recovery figures skipped ",
          "(set RUN_OM_SIM <- TRUE in 02_analysis.R and re-run).")
} else {
  EST_LAB <- c(SPR_fem = "Sex-structured female (egg) SPR - the reported status",
               SPR_bind = "Precautionary binding SPR")
  EST_COL <- c("Sex-structured female (egg) SPR - the reported status" = BREAM$teal,
               "Precautionary binding SPR" = BREAM$gold_dark)
  omp <- om |>
    filter(metric %in% names(EST_LAB)) |>
    mutate(metric = factor(unname(EST_LAB[metric]), levels = unname(EST_LAB)))
  
  # SPR is on the vertical axis here, so the zone shading is horizontal bands.
  om_recovery <- ggplot(omp, aes(n, colour = metric)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 0.20, fill = ZONE_FILL[["Overfished"]]) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.20, ymax = 0.40, fill = ZONE_FILL[["Cautionary"]]) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.40, ymax = Inf,  fill = ZONE_FILL[["Healthy"]]) +
    geom_hline(yintercept = c(0.20, 0.40), colour = "white", linewidth = 0.7) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.7, alpha = 0.8) +
    geom_line(aes(y = SPR_true, linetype = "Its own true value"), linewidth = 0.7) +
    geom_line(aes(y = median, linetype = "Replicate median"), linewidth = 0.9) +
    geom_point(aes(y = median, fill = metric), shape = 21, colour = BREAM$ink,
               size = 2.9, stroke = 0.7, show.legend = FALSE) +
    scale_colour_manual(values = EST_COL, name = NULL) +
    scale_fill_manual(values = EST_COL, guide = "none") +
    scale_linetype_manual(values = c("Replicate median" = "solid",
                                     "Its own true value" = "dashed"), name = NULL) +
    scale_x_log10(breaks = sort(unique(omp$n))) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = "Closed-loop recovery of a stock whose true status is known",
         subtitle = "A protogynous operating stock sampled at a range of measured-sample sizes and re-estimated from the lengths alone",
         x = "Number of measured fish in the sample (log scale)",
         y = "Spawning potential ratio (SPR)",
         caption = "EACH estimator is drawn against ITS OWN truth, in its own colour: the two target different quantities and a single shared\nreference line would misrepresent both. Bars span the 10-90% replicate range. Status zones apply to the egg-based ratio.")
  ggsave(file.path(FIG_DIR, "om_recovery.png"), bare_fig(om_recovery),
         width = 10, height = 6.0, dpi = 300, bg = "white")
  
  # Companion: how often the categorical call is right. With the female (egg) ratio now the
  # reported status, its zone-recovery curve is the number the Results need, and the manuscript
  # has only ever reported the binding ratio's.
  om_zone <- ggplot(omp, aes(n, p_correct, colour = metric)) +
    geom_hline(yintercept = c(0.5, 0.95), colour = "grey80", linetype = "dotted") +
    geom_line(linewidth = 0.9) +
    geom_point(aes(fill = metric), shape = 21, colour = "white", size = 2.9, stroke = 0.7,
               show.legend = FALSE) +
    scale_colour_manual(values = EST_COL, name = NULL) +
    scale_fill_manual(values = EST_COL, guide = "none") +
    scale_x_log10(breaks = sort(unique(omp$n))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(title = "How often the status call is correct, as the sample grows",
         subtitle = "Proportion of replicates whose estimate falls in the stock's true status zone",
         x = "Number of measured fish in the sample (log scale)", y = "P(correct zone)",
         caption = sprintf("Computed over converged replicates only (convergence %.0f-%.0f%% across cells). The species assessed here carry between\n%d and %d measured fish, so the relevant part of the curve is its left half.",
                           100 * min(om$conv_rate, na.rm = TRUE), 100 * max(om$conv_rate, na.rm = TRUE),
                           27L, 50L))
  ggsave(file.path(FIG_DIR, "om_zone_recovery.png"), bare_fig(om_zone),
         width = 10, height = 5.5, dpi = 300, bg = "white")
  
  # the numbers Section 3.6 needs, printed rather than left to be read off a figure
  message("\n=========== RECOVERY AT THE SAMPLE SIZES ACTUALLY ASSESSED ===========")
  print(om |>
          filter(metric %in% c("SPR_fem", "SPR_bind"), n %in% c(20, 30, 50)) |>
          transmute(n, metric, SPR_true = round(SPR_true, 3), median = round(median, 3),
                    bias = round(bias, 3), p_correct = round(p_correct, 3),
                    conv_rate = round(conv_rate, 3)), n = Inf)
  message("=====================================================================\n")
}

# =====================================================================
# 4. SEX-CHANGE OGIVE SENSITIVITY   (Figure 17 + male companion)
# =====================================================================
# The Rmd version showed one bar per species for a +/-20% shift in LD50. The
# sensitivity now has TWO axes -- the ogive's LOCATION and its WIDTH -- because no
# species in the assemblage ships a measured LD95 and the width is applied as a
# constant to all twelve hermaphrodites. The status quantity has also changed.
# Two panels per figure, and the male arm is drawn WITHOUT zone shading.
ld <- if (file.exists(here("results", "ld_sensitivity.csv")))
  read_csv(here("results", "ld_sensitivity.csv"), show_col_types = FALSE) else NULL

if (is.null(ld) || !nrow(ld)) {
  message("[05] results/ld_sensitivity.csv not found or empty - ogive figures skipped.")
} else {
  AXIS_LAB <- c(location = "Ogive LOCATION\n(LD50 shifted, width held)",
                width    = "Ogive WIDTH\n(LD95/LD50 varied, LD50 held)")
  ld_lookup <- if (exists("sp_lookup")) {
    sp_lookup |> mutate(across(everything(), as.character))
  } else {
    tibble(scientific_name = unique(ld$scientific_name),
           common_name     = unique(ld$scientific_name))
  }
  ld <- ld |>
    left_join(ld_lookup, by = "scientific_name") |>
    # keep the raw axis key BEFORE it is replaced by the display label. Section 6 filters on it,
    # and the label cannot be used for that: "Ogive LOCATION (LD50 shifted, width held)" contains
    # the word "width", so a case-insensitive match on the label caught both axes at once.
    mutate(common_name = coalesce(common_name, scientific_name),
           axis_raw = as.character(axis),
           axis = factor(unname(AXIS_LAB[axis]), levels = unname(AXIS_LAB)))
  
  # The assessment estimate is a property of the SPECIES, not of the axis being swept, so it is
  # taken once from the location axis rather than looked for on each axis separately. That matters:
  # the workbook LD95 values are rounded, so a species whose actual LD95/LD50 is 1.1021 rather than
  # exactly 1.10 (Diplodus sargus is one) has NO is_base row on the width axis, and asking each axis
  # for its own base silently returned NA and dropped the point from the figure.
  base_vals <- ld |>
    dplyr::filter(is_base) |>
    dplyr::distinct(common_name, .keep_all = TRUE) |>
    transmute(common_name, fem_base = SPR_fem, male_base = SPR_male)
  missing_base <- setdiff(unique(ld$common_name), base_vals$common_name)
  if (length(missing_base))
    warning("[05] no base row found for: ", paste(missing_base, collapse = ", "),
            " - their assessment estimate will not be marked on the ogive figures.", call. = FALSE)
  
  span <- ld |>
    group_by(common_name, axis) |>
    summarise(fem_lo = min(SPR_fem), fem_hi = max(SPR_fem),
              male_lo = min(SPR_male), male_hi = max(SPR_male),
              zone_changes = n_distinct(zone) > 1, .groups = "drop") |>
    left_join(base_vals, by = "common_name") |>
    mutate(common_name = fct_reorder(common_name, fem_base),
           robust = factor(if_else(zone_changes, "crosses a status boundary",
                                   "status robust across the range"),
                           levels = c("status robust across the range",
                                      "crosses a status boundary")))
  write_csv(span, here("results", "ld_sensitivity_span.csv"))
  
  ld_female <- ggplot(span, aes(y = common_name)) + zone_bg_x +
    geom_segment(aes(x = fem_lo, xend = fem_hi, yend = common_name),
                 colour = BREAM$ink, linewidth = 0.9) +
    geom_point(aes(x = fem_base, shape = robust), size = 4.2,
               fill = BREAM$teal, colour = BREAM$ink, stroke = 0.9) +
    scale_shape_manual(values = c("status robust across the range" = 21,
                                  "crosses a status boundary" = 23),
                       name = NULL, drop = FALSE) +
    facet_wrap(~axis) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(0)) +
    coord_cartesian(xlim = c(0, 1.04), clip = "off") +
    labs(title = "Sensitivity of the reported status to the assumed sex-change ogive",
         subtitle = "Sex-structured female (egg) SPR; the point is the assessment estimate and the bar spans the range tested",
         x = "Spawning potential ratio (female, egg-based)", y = NULL,
         caption = "The amendment is a post-processing step, so the ogive can be varied without refitting the length data. LOCATION shifts LD50 and\ncarries the fitted width; WIDTH holds LD50 and varies LD95/LD50. The width axis matters because no species in this assemblage has a\nmeasured LD95: every hermaphrodite takes the same assumed width, so it is an assumption applied universally rather than a fallback.")
  ggsave(file.path(FIG_DIR, "ld_sensitivity_female.png"), bare_fig(ld_female),
         width = 11, height = 5.0, dpi = 300, bg = "white")
  
  # Male arm: SAME two axes, NO zone shading, log axis, read against nothing.
  ld_male <- ggplot(span, aes(y = common_name)) +
    geom_segment(aes(x = male_lo, xend = male_hi, yend = common_name),
                 colour = BREAM$ink, linewidth = 0.9) +
    geom_point(aes(x = male_base), size = 4.2, shape = 21,
               fill = BREAM$gold_dark, colour = BREAM$ink, stroke = 0.9) +
    facet_wrap(~axis) +
    labs(title = "Sensitivity of the male-capacity diagnostic to the same assumption",
         subtitle = "Mature-male biomass per recruit, fished relative to unfished  ·  no status zones apply to this quantity",
         x = "SPR (male-capacity diagnostic)", y = NULL,
         caption = "Shown on the same two axes as the female panel and deliberately without zone shading: this is a mature-male BIOMASS ratio, the\n0.20 and 0.40 reference points are defined for spawning output, and they have no established meaning here. It is a screening\ndiagnostic for whether one sex is being depleted faster than the other, not a fertilisation estimate.")
  ggsave(file.path(FIG_DIR, "ld_sensitivity_male.png"), bare_fig(ld_male),
         width = 11, height = 5.0, dpi = 300, bg = "white")
  
  message("\n============ OGIVE SENSITIVITY SPAN (both axes) ============")
  print(span, n = Inf)
  message("===========================================================\n")
}

# =====================================================================
# 5. MALE-CAPACITY DIAGNOSTIC, REDRAWN   (replaces spr_male_diagnostic.png)
# =====================================================================
# FIX: as first drawn, the flagged Pagellus acarne fit returned a deficit ratio of ~213, because
# its female ratio is 0.0007 and the ratio is near-undefined when the denominator approaches zero.
# One uninterpreted point forced the log axis across three orders of magnitude and compressed the
# entire real signal, the contrast between roughly 0.16 and 1.21, into a sliver. Flagged fits are
# excluded here rather than clipped, since the species is not interpreted anywhere else either.
ar <- if (file.exists(here("results", "assessment_results.csv")))
  read_csv(here("results", "assessment_results.csv"), show_col_types = FALSE) else NULL

if (is.null(ar) || !all(c("SPR_male", "SPR_fem", "reliable") %in% names(ar))) {
  message("[05] results/assessment_results.csv missing or lacks the sex-structured columns - ",
          "male diagnostic not redrawn.")
} else {
  ar_m <- ar |>
    left_join(if (exists("sp_lookup")) sp_lookup else
      transmute(ar, scientific_name, common_name = scientific_name),
      by = "scientific_name") |>
    mutate(common_name = coalesce(as.character(common_name), scientific_name)) |>
    filter(reliable %in% TRUE, is.finite(SPR_male), is.finite(SPR_fem), SPR_fem > 0.01) |>
    mutate(deficit = SPR_male / SPR_fem,
           common_name = fct_reorder(common_name, deficit))
  n_drop <- nrow(ar) - nrow(ar_m)
  
  if (nrow(ar_m)) {
    spr_male_diag2 <- ggplot(ar_m, aes(deficit, common_name)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
      geom_segment(aes(x = 1, xend = deficit, yend = common_name),
                   colour = BREAM$ink, linewidth = 0.6) +
      geom_point(size = 4.6, shape = 21, fill = BREAM$gold_dark,
                 colour = BREAM$ink, stroke = 0.9) +
      # Explicit breaks: the default log10 breaks put nothing between 0.1 and 1, so the two
      # protogynous species at ~0.16 and ~0.18 sat on bare axis and could not be read off.
      scale_x_log10(breaks = c(0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1, 1.5),
                    labels = c("0.1", "0.15", "0.2", "0.3", "0.5", "0.7", "1", "1.5")) +
      labs(title = "Male-capacity diagnostic: is one sex being depleted faster than the other?",
           subtitle = "Mature-male biomass per recruit relative to female egg production per recruit  \u00b7  no status zones apply to this quantity",
           x = "SPR (male) / SPR (female)   (log scale; 1 = both arms equally depleted)", y = NULL,
           caption = paste0(
             "Below 1, male capacity is the depleted arm, which is what protogyny predicts because the largest fish are male and are taken\nfirst. Above 1, the egg producers are, which is what protandry predicts. This is a screening diagnostic and not a fertilisation\nestimate: the mating function relating sex ratio to realised reproduction is unknown for these species. It carries NO reference point\nand is deliberately not read against the 0.20 / 0.40 zones, which are defined for spawning output.",
             if (n_drop > 0) sprintf("\n%d fit(s) flagged unreliable are excluded: the ratio is near-undefined where the female ratio approaches zero.", n_drop) else ""))
    ggsave(file.path(FIG_DIR, "spr_male_diagnostic.png"), bare_fig(spr_male_diag2),
           width = 9.5, height = 4.8, dpi = 300, bg = "white")
    message(sprintf("[05] male diagnostic redrawn over %d reliable fit(s); %d flagged fit(s) excluded.",
                    nrow(ar_m), n_drop))
    print(ar_m |> transmute(common_name, SPR_fem = round(SPR_fem, 3),
                            SPR_male = round(SPR_male, 3), deficit = round(deficit, 3)))
  }
}

message("05_manuscript_outputs.R complete: combined assessment, recapture, recovery and ogive outputs written.")

# =====================================================================
# 6. TABLES 21a AND 21b   (reshape only; nothing is recomputed here)
# =====================================================================
# The ogive sensitivity has two axes and both belong in the appendix, but a NEW numbered table
# would renumber everything after it in Appendix F, and the document's in-text cross-references
# are typed literals rather than REF fields, so they would NOT follow. 21a and 21b renumber
# nothing. Rows are species x arm: the female (egg) ratio is the reported status, the male ratio
# is the diagnostic, and keeping both in one column block keeps each table to six rows.
if (exists("ld") && !is.null(ld) && nrow(ld)) {
  ld_tab <- function(which_axis, file, label) {
    # exact match on the raw key, not a substring match on the display label
    d <- ld |> dplyr::filter(axis_raw == which_axis)
    if (!nrow(d)) return(invisible(NULL))
    long <- dplyr::bind_rows(
      d |> transmute(common_name, sex_system, arm = "female (egg), reported status",
                     mult, value = SPR_fem, zone),
      d |> transmute(common_name, sex_system, arm = "male capacity, diagnostic",
                     mult, value = SPR_male, zone)
    )
    wide <- long |>
      mutate(mult = sprintf("%.2f", mult)) |>
      select(-zone) |>
      pivot_wider(names_from = mult, values_from = value) |>
      arrange(sex_system, common_name, arm)
    # zone_changes applies to the STATUS quantity only; the male arm has no zones to cross
    zc <- d |> group_by(common_name) |>
      summarise(status_zone_changes = dplyr::n_distinct(zone) > 1, .groups = "drop")
    wide <- wide |> left_join(zc, by = "common_name") |>
      mutate(status_zone_changes = if_else(arm == "female (egg), reported status",
                                           status_zone_changes, NA))
    write_csv(wide, here("results", file))
    message("\n---- ", label, " (results/", file, ") ----")
    print(wide |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf, width = Inf)
    invisible(wide)
  }
  ld_tab("location", "table_21a_ogive_location.csv",
         "Table 21a: ogive LOCATION (LD50 shifted, fitted width carried)")
  ld_tab("width",    "table_21b_ogive_width.csv",
         "Table 21b: ogive WIDTH (LD95/LD50 varied, LD50 held)")
}

# =====================================================================
# 7. PANEL FIGURES   (assembly only; the panels themselves are unchanged)
# =====================================================================
# Three new visuals have no figure number, and inserting them as numbered figures would shift
# every number after them. The document's caption numbers are SEQ fields and would renumber
# themselves, but its in-text cross-references are typed literals and would NOT, leaving the
# prose pointing at the wrong figures silently. Pairing each new visual as panel (b) of an
# existing figure renumbers nothing, and in each case the pairing is the argument: the male
# diagnostic only means anything read against the egg-based status, the ogive arms share two
# axes and three species, and recovery and zone recovery are one experiment.
#
# Panels are stripped and combined here; they are also written individually above, so if
# patchwork is unavailable the document can still take two images under one caption.
if (have_patchwork) {
  library(patchwork)
  tag <- patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
  
  panel <- function(a, b, file, w, h, heights = c(1, 1), drop_top_x = FALSE) {
    if (!exists(a, inherits = TRUE) || !exists(b, inherits = TRUE)) {
      message("[panel] skipped ", file, ": ", a, " or ", b, " not available.")
      return(invisible(NULL))
    }
    pa <- bare_fig(get(a)); pb <- bare_fig(get(b))
    if (is.null(pa) || is.null(pb)) {
      message("[panel] skipped ", file, ": one panel is NULL."); return(invisible(NULL))
    }
    # where the two panels share an x variable, the upper axis title is redundant
    if (drop_top_x) pa <- pa + theme(axis.title.x = element_blank())
    out <- (pa / pb) + patchwork::plot_layout(heights = heights) + tag
    ggsave(file.path(FIG_DIR, file), out, width = w, height = h, dpi = 300, bg = "white")
    message("[panel] wrote ", file)
    invisible(out)
  }
  
  # Figure 15  (a) status  (b) male-capacity diagnostic
  panel("spr_status", "spr_male_diag2", "fig15_status_and_male_panel.png",
        w = 10, h = 9.0, heights = c(1, 0.85))
  # Figure 17  (a) female arm, zone-shaded  (b) male arm, unshaded
  panel("ld_female", "ld_male", "fig17_ogive_sensitivity_panel.png",
        w = 11, h = 9.5)
  # Figure 20  (a) recovery  (b) zone recovery   -- shared x, so drop the upper axis title
  panel("om_recovery", "om_zone", "fig20_recovery_panel.png",
        w = 10, h = 10.0, heights = c(1, 0.8), drop_top_x = TRUE)
} else {
  message("[panel] patchwork absent - individual panels written; combine them in the document.")
}

message("\nPanel figures and Tables 21a/21b written. Images are bare; legend text for the ",
        "02-generated figures is in results/figure_legends.txt.")