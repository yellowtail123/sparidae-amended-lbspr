# =====================================================================
# §Z2  STRESS TESTS  —  misspecification battery for the LB-SPR amendment
#      Self-contained. Depends on om_truth(), make_length_sample(),
#      estimate_from_sample() and spr_sex_structured() from §Z / §10b.
#      Nothing here runs on source() unless its RUN_OM_* switch is TRUE.
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
N_STRESS_REP  <- 1000     # Monte-Carlo replicates per severity level
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
  nage <- length(om$catch_at_age)
  
  out <- replicate(R, {
    om_use <- om
    if (!is.null(rec_sigmaR)) {                          # stressor 4: non-equilibrium recruitment
      # log-normal deviations with mean 1 (so average recruitment is unchanged,
      # only its year-to-year variability is introduced). The distortion hits the
      # SAMPLE only; om$SPR_true stays the equilibrium reference.
      Rdev <- exp(rnorm(nage, mean = -0.5 * rec_sigmaR^2, sd = rec_sigmaR))
      om_use$catch_at_age <- om$catch_at_age * Rdev
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
  
  summ <- function(x) {
    lo <- unname(quantile(x, .1)); hi <- unname(quantile(x, .9))
    c(median        = median(x),
      lo            = lo,
      hi            = hi,
      cv            = sd(x) / mean(x),
      bias          = median(x) - om$SPR_true,
      p_correct     = mean(zone(x) == true_zone),
      truth_in_band = as.numeric(om$SPR_true >= lo && om$SPR_true <= hi))
  }
  if (n_conv == 0) {
    na <- c(median = NA, lo = NA, hi = NA, cv = NA, bias = NA,
            p_correct = NA, truth_in_band = NA)
    return(data.frame(n = n, SPR_true = om$SPR_true, R = R, n_converged = 0L,
                      metric = c("SPR_std", "SPR_fem", "SPR_bind"),
                      rbind(na, na, na), row.names = NULL))
  }
  data.frame(n = n, SPR_true = om$SPR_true, R = R, n_converged = n_conv,
             metric = c("SPR_std", "SPR_fem", "SPR_bind"),
             rbind(summ(out[, "SPR_std"]), summ(out[, "SPR_fem"]), summ(out[, "SPR_bind"])),
             row.names = NULL)
}


# ----------------------------------------------------------------------
# (D) SHARED BASELINE — the same protogynous stock as RUN_OM_SIM.
#     base_om M/K = 0.30/0.20 = 1.5, matching base_est MK = 1.5, so the
#     "no-stress" end of every sweep should reproduce your §Z result.
# ----------------------------------------------------------------------
base_est <- list(Linf = 35, MK = 1.5, L50 = 18, L95 = 22,
                 sex = "protogyny", LD50 = 25, LD95 = 28)
base_om  <- list(Linf = 35, K = 0.20, M = 0.30, L50 = 18, L95 = 22,
                 sex = "protogyny", LD50 = 25, LD95 = 28, sel50 = 16, sel95 = 21)
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
  desc_grid <- c(50, 34, 30, 28, 26)     # >1.3*Linf(=45.5) is effectively no dome
  res_dome <- do.call(rbind, lapply(desc_grid, function(d) {
    sc <- run_stress(F_true = F_STRESS, n = 200, R = N_STRESS_REP,
                     om_args  = modifyList(base_om,
                                           list(sel_type = "dome",
                                                desc50 = d, desc95 = d + 5)),
                     est_args = base_est)                 # <- still logistic
    cbind(stressor = "dome", desc50 = d, sc)
  }))
  readr::write_csv(res_dome, here::here("results", "om_stress_dome.csv"))
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
#      read it carefully before quoting a direction in any write-up.
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

STRESS_DIR <- here::here("results")
FIG_DIR    <- here::here("graphs and maps")

read_stress <- function(f) {
  p <- file.path(STRESS_DIR, f)
  if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE) else NULL
}

# short TWO-LINE panel headers (name + parameter & direction) so nothing truncates
STRESS_LABS <- ggplot2::as_labeller(c(
  dome         = "Dome selectivity\n(desc50 \u2193 = stronger)",
  compensation = "Compensatory sex change\n(kappa \u2191 = stronger)",
  rec_var      = "Recruitment variability\n(sigma_R \u2191 = stronger)"
))

zone_bg <- list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 0.20, fill = ZONE_FILL[["Overfished"]]),
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.20, ymax = 0.40, fill = ZONE_FILL[["Cautionary"]]),
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.40, ymax = Inf,  fill = ZONE_FILL[["Healthy"]]),
  geom_hline(yintercept = c(0.20, 0.40), colour = "white", linewidth = 0.7)
)

to_long <- function(df, sev) {
  if (is.null(df)) return(NULL)
  d <- df[df$metric == "SPR_bind", ]
  data.frame(stressor = d$stressor, severity = d[[sev]], SPR_true = d$SPR_true,
             median = d$median, lo = d$lo, hi = d$hi, p_correct = d$p_correct)
}
one_d <- do.call(rbind, list(
  to_long(read_stress("om_stress_dome.csv"),   "desc50"),
  to_long(read_stress("om_stress_comp.csv"),   "kappa"),
  to_long(read_stress("om_stress_recvar.csv"), "sigmaR")
))

if (!is.null(one_d)) {
  one_d$stressor <- factor(one_d$stressor, levels = c("dome", "compensation", "rec_var"))
  
  # FIG 1 — estimated SPR (median + clean 10-90% BARS) vs true SPR, with status zones
  stress_spr <- ggplot(one_d, aes(severity)) + zone_bg +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = BREAM$gold_dark, linewidth = 0.7) +
    geom_line(aes(y = SPR_true, colour = "True SPR", linetype = "True SPR"), linewidth = 0.7) +
    geom_point(aes(y = SPR_true), colour = BREAM$ink, shape = 4, size = 2) +
    geom_line(aes(y = median, colour = "Estimated (SPR_bind)", linetype = "Estimated (SPR_bind)"), linewidth = 0.85) +
    geom_point(aes(y = median), fill = BREAM$gold, colour = BREAM$ink, shape = 21, size = 2.8, stroke = 0.7) +
    facet_wrap(~stressor, scales = "free_x", labeller = STRESS_LABS) +
    scale_colour_manual(values = c("Estimated (SPR_bind)" = BREAM$gold_dark, "True SPR" = BREAM$ink), name = NULL) +
    scale_linetype_manual(values = c("Estimated (SPR_bind)" = "solid", "True SPR" = "dashed"), name = NULL) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(title = "Stress test: estimated SPR vs the truth as each assumption breaks",
         subtitle = "Bars = 10\u201390% of replicates  \u00b7  zones: red <0.2 \u00b7 amber 0.2\u20130.4 \u00b7 green >0.4",
         x = "Stressor control parameter (see panel header)", y = "SPR")
  ggsave(file.path(FIG_DIR, "stress_spr_vs_severity.png"), stress_spr,
         width = 11, height = 5.0, dpi = 300, bg = "white")
  
  # FIG 2 — status-zone accuracy
  stress_zone <- ggplot(one_d, aes(severity, p_correct)) +
    geom_hline(yintercept = c(0.5, 0.95), colour = "grey80", linetype = "dotted") +
    geom_line(colour = BREAM$teal, linewidth = 0.9) +
    geom_point(fill = BREAM$teal, colour = "white", shape = 21, size = 2.8, stroke = 0.7) +
    facet_wrap(~stressor, scales = "free_x", labeller = STRESS_LABS) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent) +
    labs(title = "Stress test: how often the status call stays correct",
         subtitle = "Probability the estimate lands in the true status zone, across replicates",
         x = "Stressor control parameter (see panel header)", y = "P(correct zone)")
  ggsave(file.path(FIG_DIR, "stress_zone_accuracy.png"), stress_zone,
         width = 11, height = 5.0, dpi = 300, bg = "white")
}

# FIG 3 — life-history bias heatmap (legend on the right, wrapped subtitle)
lh <- read_stress("om_stress_lhbias.csv")
if (!is.null(lh)) {
  lhb <- lh[lh$metric == "SPR_bind", ]
  stress_lh <- ggplot(lhb, aes(factor(MK_assumed), factor(MK_true), fill = bias)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = sprintf("%+.2f", bias)), colour = BREAM$ink, size = 4) +
    geom_tile(data = subset(lhb, MK_true == MK_assumed), fill = NA, colour = BREAM$ink, linewidth = 1.4) +
    scale_fill_gradient2(low = BREAM$teal, mid = "white", high = "#B5532E", midpoint = 0, name = "SPR\nbias") +
    coord_equal() +
    labs(title = "Stress test: life-history bias in M/K",
         subtitle = "Bias in SPR_bind (estimate \u2212 truth).\nBoxed diagonal = correctly specified M/K (\u2248 0); off-diagonal = the cost of a wrong M/K.",
         x = "Assumed M/K (what the estimator is told)",
         y = "True M/K (what the stock really is)") +
    theme(legend.position = "right")
  ggsave(file.path(FIG_DIR, "stress_lhbias_heatmap.png"), stress_lh,
         width = 7.2, height = 5.6, dpi = 300, bg = "white")
}

message("Tidied stress-test figures written to 'graphs and maps/'.")

