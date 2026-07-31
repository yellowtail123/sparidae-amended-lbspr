# =====================================================================
# Sex-structured LB-SPR  —  the sequential-hermaphroditism amendment
# =====================================================================
# A self-contained, sourceable copy of the sex-structured spawning-potential-ratio
# amendment used in the BGTW Sparidae stock assessment. Standard LB-SPR treats EVERY
# mature fish as an egg-producer; for sequential hermaphrodites that is wrong (protogyny
# -> the large fish are MALE and produce no eggs; protandry -> the egg-producers ARE the
# large fish). This reweights egg output by a proportion-female-at-length ogive and adds
# a mature-male capacity diagnostic, applied as a POST-PROCESSING step that leaves the
# LBSPR fit (F/M, selectivity, package SPR) untouched. See README.md for the mathematics.
#
# PROVENANCE: these functions are copied VERBATIM from the assessment pipeline,
# `R/02_analysis.R` section 10b, which remains the canonical source. This folder is a
# standalone extract so the method can be read and reused on its own. If the canonical
# version changes, re-copy from there.
#
# USAGE:
#   source("sex_structured_lbspr.R")
#   sx <- spr_sex_structured(FM, SL50, SL95, Linf, MK, L50, L95,
#                            sex_system = "protogyny", LD50 = ..., LD95 = ...)
#   # FM, SL50, SL95 come from a standard LBSPR fit; fit_lbspr_one() below wraps that.
# See example.R for a runnable demonstration.
# =====================================================================

suppressPackageStartupMessages({
  library(tibble)   # tibble()
  library(dplyr)    # bind_cols()
})

# ---- constants the amendment needs (mirrors the pipeline's run-control block) --------
BINWIDTH           <- 1     # length-bin width (cm)
LDELTA_FACTOR      <- 1.10  # LD95 = 1.10 * LD50 when only LD50 is given
SEX_STRUCTURED_SPR <- TRUE  # master switch; FALSE -> behaves like standard LBSPR

# ---- optional-package guard (LBSPR is only needed for the standard fit wrapper) ------
need <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) message(sprintf("[optional] package '%s' not installed; fit_lbspr_one() will return NULL.", pkg))
  ok
}
have_lbspr <- need("LBSPR")

# ---- length-frequency helper (for the LBSPR fit; base R only) ------------------------
make_lf <- function(x, Linf, binwidth = BINWIDTH) {
  x <- x[is.finite(x)]
  br <- seq(0, ceiling(max(c(x, Linf), na.rm = TRUE)) + binwidth, by = binwidth)
  mids <- br[-length(br)] + binwidth / 2
  counts <- as.numeric(table(cut(x, breaks = br, right = FALSE)))
  list(mids = mids, counts = counts)
}

# ── Sex-structured SPR (sequential hermaphroditism) ──────────────────
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
#             per-recruit; a DIAGNOSTIC read as SPR_male / SPR_fem against parity, NOT a
#             fertilisation estimate and NOT read against the 0.20 / 0.40 zones - those are
#             defined for spawning output and none exists for a male biomass ratio).
#             SPR_fem: the status to report.   SPR_bind: precautionary minimum, retained
#             for continuity and audit; not the status
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

# ── Standard LB-SPR fit (LBSPR package) + the amendment applied to it ─
# Fits the standard length-based SPR to a vector of measured lengths `L` given a
# life-history row `lh`, then (if SEX_STRUCTURED_SPR) feeds the FITTED F/M and
# selectivity straight into spr_sex_structured() - showing exactly how the amendment
# is layered on top of an unchanged fit. Returns NULL if LBSPR is unavailable or the
# fit fails. `lh` needs: scientific_name, Linf, L50, L95, MK, a, b, and (for the
# amendment) FecB, CVLinf, sex_system, LD50, LD95.
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
