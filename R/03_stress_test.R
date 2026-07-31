# =====================================================================
# §Z2  STRESS TESTS  —  misspecification battery for the LB-SPR amendment
#      NOT self-contained: requires om_truth(), make_length_sample(),
#      estimate_from_sample() and spr_sex_structured() from 02_analysis.R, which
#      00_run_all.R sources first. The RUN_OM_* switches are defined BELOW and
#      therefore override any set in 02, so the whole battery runs under
#      00_run_all.R; set them FALSE here to skip it. The figure block at the foot
#      draws from whatever om_stress_*.csv are on disk, so delete stale CSVs before
#      re-running with changed settings rather than trusting the figures to refresh.
#
#      The principle: each stressor builds a TRUE stock that breaks ONE
#      LB-SPR assumption, then estimates it with the unchanged, naive
#      estimator (which always fits a logistic and uses the assumed
#      life-history). The bias that results is the cost of that wrong
#      assumption. We only ever modify the TRUTH side.
# =====================================================================


# ----------------------------------------------------------------------
# (A) CONFIG — put these next to RUN_OM_SIM in your run-control block.
# ----------------------------------------------------------------------
RUN_OM_DOME   <- TRUE   # stressor 1: dome-shaped (non-asymptotic) selectivity
RUN_OM_COMP   <- TRUE   # stressor 2: compensatory (plastic) sex change
RUN_OM_LHBIAS <- TRUE   # stressor 3: biased (not just uncertain) life-history
RUN_OM_RECVAR <- TRUE   # stressor 4: non-equilibrium recruitment
N_STRESS_REP  <- if (exists("FINAL_RUN") && isTRUE(FINAL_RUN)) 1000L else 500L  # replicates per severity level
# Each replicate now calls LBSPR::LBSPRfit (the deployed estimator), which is an order
# of magnitude heavier than the former hand-optim, so the whole battery at 9999 would
# run for hours. 1000 already pins the 10-90% bands and the zone probabilities; raise
# it for camera-ready if you want, but 9999 buys negligible extra stability at large cost.


# ----------------------------------------------------------------------
# (B) om_truth() — inherited from 02_analysis.R, not redefined here.
#     The canonical om_truth() (with the sel_type / kappa stress hooks)
#     lives in 02_analysis.R and is already in scope: 00_run_all.R sources
#     02 before 03, which is required anyway because run_stress() below
#     also uses make_length_sample(), estimate_from_sample() and
#     spr_sex_structured() from 02. The former drop-in copy here was
#     byte-for-byte identical to that canonical version and was removed so
#     the operating model has a single definition.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# (C) run_stress() — Monte-Carlo runner for ONE stressed scenario.
#     Mirrors run_scenario(), but (i) adds an optional per-replicate
#     recruitment distortion (stressor 4) and (ii) also reports
#     truth_in_band, a quick check of whether the 10-90% band covers
#     the truth. estimate_from_sample() is called UNCHANGED — that is
#     the firewall: the estimator never learns which assumption broke.
# ----------------------------------------------------------------------
run_stress <- function(F_true, n, R = 500, om_args = list(), est_args = list(),
                       rec_sigmaR = NULL) {
  om <- do.call(om_truth, c(list(F_true = F_true), om_args))
  zone <- function(s) cut(s, c(-Inf, 0.20, 0.40, Inf),
                          labels = c("overfished", "cautionary", "healthy"))
  true_zone <- zone(om$SPR_true)
  # Year-class index over the age grid. A recruitment deviation belongs to a COHORT, so it must
  # be constant across the age slices of the same year; one draw per slice (da = 0.1 yr) makes ten
  # independent draws inside each year class, which average out before they reach the length
  # composition and understate the stressor by roughly sqrt(da).
  cohort <- floor(om$ages)
  coh_id <- match(cohort, sort(unique(cohort)))
  n_coh  <- max(coh_id)
  
  out <- replicate(R, {
    om_use <- om
    if (!is.null(rec_sigmaR)) {                          # stressor 4: non-equilibrium recruitment
      # mean-1 log-normal deviation per YEAR CLASS (average recruitment unchanged, only its
      # year-to-year variability introduced). The distortion hits the SAMPLE only;
      # om$SPR_true stays the equilibrium reference.
      Rdev_y <- exp(rnorm(n_coh, mean = -0.5 * rec_sigmaR^2, sd = rec_sigmaR))
      om_use$catch_at_age <- om$catch_at_age * Rdev_y[coh_id]
    }
    samp <- make_length_sample(om_use, n)
    if (is.null(samp) || sum(samp$counts > 0) < 3) return(c(NA, NA, NA))
    est <- tryCatch(do.call(estimate_from_sample, c(list(samp = samp), est_args)),
                    error = function(e) NULL)
    if (is.null(est)) return(c(NA, NA, NA))
    c(est$SPR_std, est$SPR_fem, est$SPR_bind)
  })
  
  out <- t(out); colnames(out) <- c("SPR_std", "SPR_fem", "SPR_bind")
  out <- out[complete.cases(out), , drop = FALSE]
  n_conv <- nrow(out)
  
  truth_of <- c(SPR_std = om$SPR_std_true, SPR_fem = om$SPR_fem_true, SPR_bind = om$SPR_bind_true)
  summ <- function(x, truth) {
    lo <- unname(quantile(x, .1)); hi <- unname(quantile(x, .9))
    c(median        = median(x),
      lo            = lo,
      hi            = hi,
      cv            = sd(x) / mean(x),
      bias          = median(x) - truth,
      p_correct     = mean(zone(x) == zone(truth)),
      truth_in_band = as.numeric(truth >= lo && truth <= hi))
  }
  if (n_conv == 0) {
    na <- c(median = NA, lo = NA, hi = NA, cv = NA, bias = NA,
            p_correct = NA, truth_in_band = NA)
    return(data.frame(n = n,
                      SPR_true = c(om$SPR_std_true, om$SPR_fem_true, om$SPR_bind_true),
                      R = R, n_converged = 0L,
                      metric = c("SPR_std", "SPR_fem", "SPR_bind"),
                      rbind(na, na, na), row.names = NULL))
  }
  data.frame(n = n, SPR_true = unname(truth_of), R = R, n_converged = n_conv,
             metric = c("SPR_std", "SPR_fem", "SPR_bind"),
             rbind(summ(out[, "SPR_std"],  truth_of[["SPR_std"]]),
                   summ(out[, "SPR_fem"],  truth_of[["SPR_fem"]]),
                   summ(out[, "SPR_bind"], truth_of[["SPR_bind"]])),
             row.names = NULL)
}


# ----------------------------------------------------------------------
# (D) SHARED BASELINE — the same protogynous stock as RUN_OM_SIM.
#     base_om M/K = 0.30/0.20 = 1.5, matching base_est MK = 1.5, so the
#     "no-stress" end of every sweep should reproduce your §Z result.
# ----------------------------------------------------------------------
base_est <- if (exists("OM_BASE_EST")) OM_BASE_EST else
  list(Linf = 35, MK = 1.5, L50 = 18, L95 = 22, sex = "protogyny",
       LD50 = 25, LD95 = if (exists("LDELTA_FACTOR")) LDELTA_FACTOR * 25 else 27.5)
base_om  <- if (exists("OM_BASE_OM")) OM_BASE_OM else
  list(Linf = 35, K = 0.20, M = 0.30, L50 = 18, L95 = 22, sex = "protogyny",
       LD50 = 25, LD95 = if (exists("LDELTA_FACTOR")) LDELTA_FACTOR * 25 else 27.5,
       sel50 = 16, sel95 = 21)
if (!exists("OM_BASE_OM"))
  warning("[stress] OM_BASE_OM not found - 02_analysis.R has not been sourced in this session. ",
          "Falling back to a local baseline; source 02 first so the recovery test and the ",
          "stressor battery share one stock.", call. = FALSE)
message(sprintf("[stress] baseline sex-change ogive: LD50 = %.1f, LD95 = %.2f (width = %.2f x LD50).",
                base_om$LD50, base_om$LD95, base_om$LD95 / base_om$LD50))
F_STRESS <- 0.45     # the fishing pressure each stressor is evaluated at


# ----------------------------------------------------------------------
# (E1) STRESSOR 1 — dome-shaped selectivity   [highest priority]
#      Sweep the descending-limb location desc50 from far above Linf
#      (no dome -> should match baseline) down into the catch range
#      (strong dome). Estimator keeps fitting a logistic.
#      Expect SPR_bind biased DOWNWARD as the dome bites.
# ----------------------------------------------------------------------
if (RUN_OM_DOME) {
  set.seed(42)
  # The sweep needs a TRUE zero, and desc50 = 50 is not one. The descending logistic is never
  # truncated, so at desc50 = 50 (desc95 = 55) it still returns 0.934 at the top of the length
  # grid (1.3 * Linf = 45.5), i.e. it removes 6.6% of selectivity there. Claiming that the
  # zero-severity end reproduces the unstressed baseline was therefore not exactly true for this
  # stressor, and nothing tested it. The first row is now sel_type = "logistic", which is the
  # unstressed case by construction, and desc50 = 50 is retained as the first DOME level so the
  # residual is visible rather than assumed away.
  #
  # Note also that severity is far from linear in desc50: between 50 and 34 the deselection at
  # the top of the grid goes from 6.6% to 99.9%. Intermediate levels are added so the onset can
  # be located rather than straddled.
  dome_grid <- tibble::tibble(
    sel_type = c("logistic", rep("dome", 7L)),
    desc50   = c(NA_real_,  50, 42, 38, 34, 30, 28, 26))
  res_dome <- do.call(rbind, lapply(seq_len(nrow(dome_grid)), function(i) {
    g <- dome_grid[i, ]
    om_extra <- if (g$sel_type == "dome")
      list(sel_type = "dome", desc50 = g$desc50, desc95 = g$desc50 + 5) else
      list(sel_type = "logistic")
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args  = modifyList(base_om, om_extra),
                     est_args = base_est)                 # <- still logistic
    cbind(stressor = "dome", sel_type = g$sel_type, desc50 = g$desc50, sc)
  }))
  readr::write_csv(res_dome, here::here("results", "om_stress_dome.csv"))
  # The stated design requirement is that the zero-severity end reproduces the baseline. State it
  # as a check that runs, rather than as a claim in the appendix that nothing verifies.
  dome_zero <- subset(res_dome, sel_type == "logistic" & metric == "SPR_bind")
  dome_mild <- subset(res_dome, sel_type == "dome" & desc50 == 50 & metric == "SPR_bind")
  if (nrow(dome_zero) && nrow(dome_mild))
    message(sprintf(
      "[dome] zero-severity check: logistic median %.4f vs desc50=50 median %.4f (difference %.4f). %s",
      dome_zero$median[1], dome_mild$median[1], dome_mild$median[1] - dome_zero$median[1],
      "The untruncated descending limb still removes ~6.6% of selectivity at the top of the grid at desc50 = 50."))
  print(res_dome)
}


# ----------------------------------------------------------------------
# (E2) STRESSOR 2 — compensatory (plastic) sex change
#      Truth shifts the sex-change schedule earlier under fishing
#      (kappa > 0); estimator keeps the fixed unfished LD50/LD95.
#      kappa = 0 = baseline.
#
#      DIRECTION (verified by simulation): the fixed-schedule estimator does
#      NOT see the females that turn male early, so it OVER-counts female egg
#      output. Expect SPR_bind biased HIGH (optimistic) relative to the
#      compensated truth, and the true SPR itself to fall as kappa rises.
#      Note this is the OPPOSITE of the intuition that the male floor makes
#      the estimate over-conservative: the collapse of true female output
#      dominates. This is an important, slightly counter-intuitive result —
#      read it carefully before quoting a direction in the manuscript.
# ----------------------------------------------------------------------
if (RUN_OM_COMP) {
  set.seed(42)
  kappa_grid <- c(0, 0.1, 0.2, 0.3, 0.5)
  res_comp <- do.call(rbind, lapply(kappa_grid, function(k) {
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args  = modifyList(base_om, list(kappa = k)),
                     est_args = base_est)                 # <- fixed schedule
    cbind(stressor = "compensation", kappa = k, sc)
  }))
  readr::write_csv(res_comp, here::here("results", "om_stress_comp.csv"))
  print(res_comp)
}


# ----------------------------------------------------------------------
# (E2b) CROSSING OF THE FEMALE AND MALE ARMS UNDER COMPENSATION
#       Both true arms are analytic inside om_truth(), so the kappa at which SPR_fem_true
#       and SPR_male_true cross needs no fitting and no Monte-Carlo: it is a one-dimensional
#       root find on their difference, costing milliseconds. This is the maximum of the
#       binding truth and the point at which the sign of the bias in the reported ratio
#       turns, so it is reported in the Results and needs a reproducible source. Evaluated
#       at the stressed fishing mortality and at one value either side, because the crossing
#       moves with exploitation. Writes results/om_arm_crossing.csv; overwrites nothing else.
# ----------------------------------------------------------------------
RUN_OM_CROSS <- TRUE
if (RUN_OM_CROSS) {
  arm_gap <- function(k, F_true) {
    om <- do.call(om_truth, c(list(F_true = F_true),
                              modifyList(base_om, list(kappa = k))))
    om$SPR_fem_true - om$SPR_male_true          # + below the crossing, - above it
  }
  crossing_one <- function(F_true, interval = c(0, 0.9)) {
    r <- try(stats::uniroot(arm_gap, interval = interval, F_true = F_true, tol = 1e-8),
             silent = TRUE)
    if (inherits(r, "try-error"))
      return(data.frame(F_true = F_true, kappa_cross = NA_real_, SPR_at_cross = NA_real_))
    om <- do.call(om_truth, c(list(F_true = F_true),
                              modifyList(base_om, list(kappa = r$root))))
    data.frame(F_true = F_true, kappa_cross = round(r$root, 3),
               SPR_at_cross = round(mean(c(om$SPR_fem_true, om$SPR_male_true)), 3))
  }
  arm_crossing <- do.call(rbind, lapply(c(0.20, F_STRESS, 0.60), crossing_one))
  readr::write_csv(arm_crossing, here::here("results", "om_arm_crossing.csv"))
  message("\n===== FEMALE / MALE ARM CROSSING UNDER COMPENSATORY SEX CHANGE =====")
  print(arm_crossing)
  message("kappa_cross = compensation strength at which the binding truth peaks;")
  message("below it the male arm binds, above it the female arm binds.")
  message("====================================================================\n")
}


# ----------------------------------------------------------------------
# (E3) STRESSOR 3 — biased life-history (wrong, not just uncertain)
#      Truth built at M/K = MK_true (varying M, K fixed); estimator told
#      MK_assumed. The diagonal (MK_true == MK_assumed) is the zero-bias
#      baseline. This is distinct from your parameter Monte-Carlo, which
#      propagates uncertainty around a CORRECT centre.
# ----------------------------------------------------------------------
if (RUN_OM_LHBIAS) {
  set.seed(42)
  K_FIXED <- 0.20
  mk_grid <- c(1.2, 1.5, 1.8)
  grid <- expand.grid(MK_true = mk_grid, MK_assumed = mk_grid)
  res_lh <- do.call(rbind, Map(function(mkt, mka) {
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args  = modifyList(base_om,
                                           list(M = mkt * K_FIXED, K = K_FIXED)),
                     est_args = modifyList(base_est, list(MK = mka)))
    cbind(stressor = "lh_bias", MK_true = mkt, MK_assumed = mka, sc)
  }, grid$MK_true, grid$MK_assumed))
  readr::write_csv(res_lh, here::here("results", "om_stress_lhbias.csv"))
  print(res_lh)
}


# ----------------------------------------------------------------------
# (E5) STRESSOR 5 — MISSPECIFIED SEX-CHANGE OGIVE WIDTH   [gated OFF]
#      Truth built with LD95 = width_true * LD50; estimator always told
#      LDELTA_FACTOR * LD50, which is what it is told for every real
#      species. The diagonal case (width_true == LDELTA_FACTOR) is the
#      zero-bias baseline and reproduces the unstressed run exactly.
#
#      WHY THIS EXISTS. The assessment applies LD95 = LDELTA_FACTOR * LD50
#      to all twelve hermaphrodites and not one of them ships a measured
#      LD95, so the width is the single most exposed input in the whole
#      analysis. Section 10c of 02_analysis.R already sweeps it on the
#      ASSESSMENT side and answers "how far does the reported number move?".
#      This asks the strictly harder question the assessment side cannot:
#      "is the estimator BIASED when the true width is not the assumed one?"
#      It is the same shape as stressor 3, which does this for M/K.
#
#      COST OF TURNING IT ON: a fourth panel in the stressor figures, a
#      further block in the misspecification table, and a paragraph to
#      write. Left FALSE so it does not silently reshape figures that are
#      about to be finalised. Worth running ONCE regardless of whether it
#      goes in the paper: if the bias is small it is a free, bounded line in
#      the limitations, and if it is not, that is something to know before
#      submission rather than after.
# ----------------------------------------------------------------------
RUN_OM_OGIVE <- FALSE   # TRUE -> adds an ogive-width misspecification sweep (see note above)

if (isTRUE(RUN_OM_OGIVE)) {
  set.seed(42)
  ld_factor  <- if (exists("LDELTA_FACTOR")) LDELTA_FACTOR else 1.10   # what the estimator is told
  width_grid <- sort(unique(c(ld_factor, 1.20, 1.30, 1.40, 1.50)))     # what the stock really is
  res_ogive <- do.call(rbind, lapply(width_grid, function(w) {
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args  = modifyList(base_om,
                                           list(LD95 = w * base_om$LD50)),
                     est_args = base_est)     # <- estimator keeps LDELTA_FACTOR, as in production
    cbind(stressor = "ogive_width", width_true = w,
          LD95_true = w * base_om$LD50, LD95_assumed = base_est$LD95, sc)
  }))
  readr::write_csv(res_ogive, here::here("results", "om_stress_ogive.csv"))
  base_row <- subset(res_ogive, abs(width_true - ld_factor) < 1e-9 & metric == "SPR_fem")
  if (nrow(base_row))
    message(sprintf("[ogive] zero-misspecification check: bias in SPR_fem at width_true = %.2f is %+.4f (should be ~0).",
                    ld_factor, base_row$bias[1]))
  print(res_ogive)
}


# ----------------------------------------------------------------------
# (E4) STRESSOR 4 — non-equilibrium recruitment
#      Distort the OBSERVED age structure each replicate with mean-1
#      log-normal recruitment deviations; the equilibrium om$SPR_true is
#      the reference (per-recruit SPR is recruitment-invariant, so the
#      bias enters through the SAMPLE, not the truth). sigmaR = 0 =
#      baseline. Expect scatter/drift that grows with sigmaR.
# ----------------------------------------------------------------------
if (RUN_OM_RECVAR) {
  set.seed(42)
  sigma_grid <- c(0, 0.3, 0.6, 0.9)
  res_rec <- do.call(rbind, lapply(sigma_grid, function(s) {
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args = base_om, est_args = base_est,
                     rec_sigmaR = if (s == 0) NULL else s)
    cbind(stressor = "rec_var", sigmaR = s, sc)
  }))
  readr::write_csv(res_rec, here::here("results", "om_stress_recvar.csv"))
  print(res_rec)
}
# =====================================================================
# §Z3  STRESS-TEST FIGURES (tidied)  —  visualise the §Z2 CSVs
#      Reads results/om_stress_*.csv, writes PNGs to "graphs and maps/".
#      Run AFTER analysis.R (uses theme_bream / BREAM / ZONE_FILL).
# =====================================================================
library(ggplot2)

# Exports honour 02_analysis.R's FIG_TEXT switch: bare images, description in the Word
# legend. Fallback keeps this file usable if it is ever run without 02 in session.
if (!exists("bare_fig")) bare_fig <- function(p) p
STRESS_DIR <- here::here("results")
FIG_DIR    <- here::here("graphs and maps")

read_stress <- function(f) {
  p <- file.path(STRESS_DIR, f)
  if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE) else NULL
}

# short TWO-LINE panel headers (name + parameter & direction) so nothing truncates.
# Every panel now runs LEAST severe on the left to MOST severe on the right; see severity_plot
# below for how the dome panel is made to do that while keeping its true desc50 labels.
STRESS_LABS <- ggplot2::as_labeller(c(
  dome         = "Dome selectivity\n(desc50, cm; lower = stronger)",
  compensation = "Compensatory sex change\n(kappa; higher = stronger)",
  rec_var      = "Recruitment variability\n(sigma_R; higher = stronger)",
  ogive_width  = "Sex-change ogive width\n(true LD95/LD50; higher = stronger)"
))

zone_bg <- list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 0.20, fill = ZONE_FILL[["Overfished"]]),
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.20, ymax = 0.40, fill = ZONE_FILL[["Cautionary"]]),
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.40, ymax = Inf,  fill = ZONE_FILL[["Healthy"]]),
  geom_hline(yintercept = c(0.20, 0.40), colour = "white", linewidth = 0.7)
)

# Both estimators are carried, not just the binding one. SPR_fem is the REPORTED status and it
# targets the operating model's own female-egg truth; SPR_bind targets a different quantity. Each
# is therefore plotted against ITS OWN truth, and a single unlabelled reference line for two
# estimators that target different values is precisely the confusion to avoid here.
to_long <- function(df, sev) {
  if (is.null(df)) return(NULL)
  d <- df[df$metric %in% c("SPR_fem", "SPR_bind"), ]
  if (!nrow(d)) return(NULL)
  # n_converged / R arrived with the convergence-reporting fix. data.frame() errors on a NULL
  # column rather than dropping it, so a results file written before that fix would take this
  # function down; fill instead, and let the caption say nothing rather than say something wrong.
  nconv <- if (!is.null(d$n_converged)) d$n_converged else NA_real_
  reps  <- if (!is.null(d$R))           d$R           else NA_real_
  out <- data.frame(stressor = d$stressor, severity = d[[sev]], metric = d$metric,
                    SPR_true = d$SPR_true, median = d$median, lo = d$lo, hi = d$hi,
                    p_correct = d$p_correct,
                    n_converged = nconv, R = reps)
  # the dome sweep now carries a sel_type = "logistic" row with no desc50; it exists to verify
  # the zero-severity baseline claim numerically and has no position on a desc50 axis
  out[!is.na(out$severity), , drop = FALSE]
}
one_d <- do.call(rbind, list(
  to_long(read_stress("om_stress_dome.csv"),   "desc50"),
  to_long(read_stress("om_stress_comp.csv"),   "kappa"),
  to_long(read_stress("om_stress_recvar.csv"), "sigmaR"),
  # present only if RUN_OM_OGIVE was TRUE; read_stress returns NULL otherwise and rbind drops it,
  # so the figures stay at three panels unless the stressor has actually been run
  to_long(read_stress("om_stress_ogive.csv"),  "width_true")
))

if (!is.null(one_d)) {
  one_d$stressor <- factor(one_d$stressor,
                           levels = intersect(c("dome", "compensation", "rec_var", "ogive_width"),
                                              unique(one_d$stressor)))
  one_d$metric   <- factor(one_d$metric, levels = c("SPR_fem", "SPR_bind"),
                           labels = c("Estimated (female, egg SPR)", "Estimated (binding SPR)"))
  # SEVERITY DIRECTION. For the dome, severity rises as desc50 FALLS, so plotting desc50 directly
  # ran that panel backwards relative to the other two: its left edge was the most stressed case
  # while kappa = 0 and sigma_R = 0 sit at the left of theirs. A single caption cannot describe
  # both conventions, and the one that was there described the wrong one. Negating desc50 makes
  # every panel read least-severe to most-severe left to right; abs() in the axis labels restores
  # the true desc50 values, and leaves the other two panels untouched since they are positive.
  one_d$severity_plot <- ifelse(one_d$stressor == "dome", -one_d$severity, one_d$severity)
  sev_labels <- function(x) formatC(abs(x), format = "g")
  conv_note <- {
    cr <- one_d$n_converged / one_d$R
    if (any(cr < 0.99, na.rm = TRUE))
      sprintf("\nReplicate convergence ranges %.0f-%.0f%%; P(correct zone) is conditional on convergence, and non-convergence is not independent of severity.",
              100 * min(cr, na.rm = TRUE), 100 * max(cr, na.rm = TRUE)) else ""
  }
  
  # FIG 1 — estimated SPR (median + 10-90% bars) vs true SPR, per estimator, with status zones.
  # Each estimator gets its OWN truth line because they target different quantities: the female
  # (egg) ratio targets the operating model's female-egg truth, the binding ratio targets the
  # lesser of the female and male-capacity truths. One line serving both is what made the
  # recovery figure misreadable.
  EST_COLS <- c("Estimated (female, egg SPR)" = BREAM$teal,
                "Estimated (binding SPR)"     = BREAM$gold_dark)
  TRU_COLS <- c("True SPR (female, egg)" = BREAM$teal,
                "True SPR (binding)"     = BREAM$gold_dark)
  one_d$truth_lab <- ifelse(grepl("female", as.character(one_d$metric)),
                            "True SPR (female, egg)", "True SPR (binding)")

  stress_spr <- ggplot(one_d, aes(severity_plot)) + zone_bg +
    geom_linerange(aes(ymin = lo, ymax = hi, colour = metric),
                   linewidth = 0.7, alpha = 0.75) +
    geom_line(aes(y = SPR_true, colour = truth_lab, linetype = "truth"), linewidth = 0.7) +
    geom_point(aes(y = SPR_true, colour = truth_lab), shape = 4, size = 2, show.legend = FALSE) +
    geom_line(aes(y = median, colour = metric, linetype = "estimate"), linewidth = 0.85) +
    geom_point(aes(y = median, fill = metric), colour = BREAM$ink,
               shape = 21, size = 2.8, stroke = 0.7, show.legend = FALSE) +
    facet_wrap(~stressor, scales = "free_x", labeller = STRESS_LABS) +
    scale_colour_manual(values = c(EST_COLS, TRU_COLS), name = NULL) +
    scale_fill_manual(values = EST_COLS, guide = "none") +
    scale_linetype_manual(values = c(estimate = "solid", truth = "dashed"), name = NULL) +
    scale_x_continuous(labels = sev_labels) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(title = "Stress test: estimated SPR against the truth as each assumption breaks",
         subtitle = "Severity increases left to right in every panel  \u00b7  bars = 10\u201390% of replicates  \u00b7  zones: red <0.2 \u00b7 amber 0.2\u20130.4 \u00b7 green >0.4",
         x = "Stressor control parameter (see panel header)", y = "SPR",
         caption = paste0(
           "Each estimator is scored against its own truth: the female (egg) ratio against the operating model's female-egg SPR, the binding\nratio against the lesser of the female and male-capacity truths. Dashed lines are those truths, solid lines the replicate medians.\nThe dome panel's axis is drawn so that severity increases rightwards as it does in the other two; its tick labels are the true desc50 in cm.",
           conv_note))
  ggsave(file.path(FIG_DIR, "stress_spr_vs_severity.png"), bare_fig(stress_spr),
         width = 11.5, height = 5.4, dpi = 300, bg = "white")

  # FIG 2 — status-zone accuracy, per estimator
  stress_zone <- ggplot(one_d, aes(severity_plot, p_correct, colour = metric)) +
    geom_hline(yintercept = c(0.5, 0.95), colour = "grey80", linetype = "dotted") +
    geom_line(linewidth = 0.9) +
    geom_point(aes(fill = metric), colour = "white", shape = 21, size = 2.8, stroke = 0.7,
               show.legend = FALSE) +
    facet_wrap(~stressor, scales = "free_x", labeller = STRESS_LABS) +
    scale_colour_manual(values = EST_COLS, name = NULL) +
    scale_fill_manual(values = EST_COLS, guide = "none") +
    scale_x_continuous(labels = sev_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
    labs(title = "Stress test: how often the status call stays correct",
         subtitle = "Probability the estimate lands in the true status zone  \u00b7  severity increases left to right in every panel",
         x = "Stressor control parameter (see panel header)", y = "P(correct zone)",
         caption = paste0(
           "Computed over converged replicates only.",
           if (nzchar(conv_note)) conv_note else
             "\nAll cells converged in at least 99% of replicates."))
  ggsave(file.path(FIG_DIR, "stress_zone_accuracy.png"), bare_fig(stress_zone),
         width = 11.5, height = 5.4, dpi = 300, bg = "white")
}

# FIG 3 — life-history bias heatmap (legend on the right, wrapped subtitle).
# Reports bias in SPR_fem, the quantity the assessment now reports as status. SPR_bind is still
# in the results file for anyone who wants it; it is not what the paper reads against the zones.
lh <- read_stress("om_stress_lhbias.csv")
if (!is.null(lh)) {
  lhb <- lh[lh$metric == "SPR_fem", ]
  if (!nrow(lhb)) lhb <- lh[lh$metric == "SPR_bind", ]   # older results file: fall back, do not fail
  stress_lh <- ggplot(lhb, aes(factor(MK_assumed), factor(MK_true), fill = bias)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = sprintf("%+.2f", bias)), colour = BREAM$ink, size = 4) +
    geom_tile(data = subset(lhb, MK_true == MK_assumed), fill = NA, colour = BREAM$ink, linewidth = 1.4) +
    scale_fill_gradient2(low = BREAM$teal, mid = "white", high = "#B5532E", midpoint = 0, name = "SPR\nbias") +
    coord_equal() +
    labs(title = "Stress test: life-history bias in M/K",
         subtitle = paste0("Bias in ", unique(as.character(lhb$metric))[1],
                           " (estimate \u2212 truth).\nBoxed diagonal = correctly specified M/K (\u2248 0); off-diagonal = the cost of a wrong M/K."),
         x = "Assumed M/K (what the estimator is told)",
         y = "True M/K (what the stock really is)") +
    theme(legend.position = "right")
  ggsave(file.path(FIG_DIR, "stress_lhbias_heatmap.png"), bare_fig(stress_lh),
         width = 7.2, height = 5.6, dpi = 300, bg = "white")
}

message("Tidied stress-test figures written to 'graphs and maps/'.")

