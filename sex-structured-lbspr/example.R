# =====================================================================
# example.R  —  runnable demonstration of the sex-structured LB-SPR amendment
# =====================================================================
# Applies the amendment to three real BGTW sparids, one per sex system, using
# ILLUSTRATIVE fitted selectivity + fishing pressure (F/M). In a real assessment those
# three numbers come from a standard LBSPR fit to measured lengths (see fit_lbspr_one()
# and the note at the bottom); here they are fixed so the demo is deterministic and needs
# no catch data. Run from this folder:  Rscript example.R
#
# What to look for in the printed table:
#   gonochore  -> SPR_fem == SPR_std                 (the amendment reduces to standard LBSPR)
#   protandry  -> the large egg-producers are fished (SPR_fem shifts vs standard)
#   protogyny  -> SPR_bind = min(SPR_fem, SPR_male)  (the precautionary male-capacity floor binds)
# =====================================================================

# resolve this script's own folder, so it runs from anywhere
.args <- commandArgs(trailingOnly = FALSE)
.self <- sub("^--file=", "", .args[grep("^--file=", .args)])
here_dir <- if (length(.self)) dirname(normalizePath(.self)) else getwd()

source(file.path(here_dir, "sex_structured_lbspr.R"))

# Life-history parameters are the thesis values on the FORK-LENGTH scale (published
# total lengths x 0.92); LD50/LD95 are NA for the gonochore, where psi_f = 1/2.
lh <- read.csv(file.path(here_dir, "example_life_history.csv"),
               stringsAsFactors = FALSE)

# --- illustrative fitted parameters (would come from LBSPRfit) -----------------------
ILLUSTRATIVE_FM <- 1.5    # fishing mortality = 1.5 x natural mortality (F/M)
# selectivity taken at the maturity ogive for the illustration (SL50 = L50, SL95 = L95)

rows <- lapply(seq_len(nrow(lh)), function(i) {
  r <- lh[i, ]
  sx <- spr_sex_structured(
    FM = ILLUSTRATIVE_FM, SL50 = r$L50, SL95 = r$L95,
    Linf = r$Linf, MK = r$MK, L50 = r$L50, L95 = r$L95,
    FecB = r$FecB, CVLinf = r$CVLinf,
    sex_system = r$sex_system, LD50 = r$LD50, LD95 = r$LD95,
    MaleExp = r$b)
  data.frame(
    species  = r$scientific_name,
    sex      = r$sex_system,
    SPR_std  = round(sx$SPR_gono_check, 3),   # standard LBSPR (the replica of the package SPR)
    SPR_fem  = round(sx$SPR_fem, 3),          # female (egg-based) SPR
    SPR_male = round(sx$SPR_male, 3),         # male reproductive-capacity ratio (floor)
    SPR_bind = round(sx$SPR_bind, 3),         # status reported = min for protogyny, else SPR_fem
    shift    = round(sx$SPR_bind - sx$SPR_gono_check, 3),
    stringsAsFactors = FALSE)
})
out <- do.call(rbind, rows)

cat("\nSex-structured LB-SPR amendment  —  illustrative F/M =", ILLUSTRATIVE_FM,
    ", selectivity at the maturity ogive\n")
cat("(SPR_std = standard LBSPR; SPR_bind = the status the amendment reports)\n\n")
print(out, row.names = FALSE)
cat("\nNote: gonochore SPR_fem equals SPR_std (the amendment is a generalisation, not a",
    "\ncompeting method); for protogyny the male floor SPR_male can bind below SPR_fem.\n")

# --- With real data instead of the illustration --------------------------------------
# Given a vector `L` of measured fork lengths and a life-history row `lh_row`, the whole
# fit -> reweight path is one call (needs the LBSPR package installed):
#
#   res <- fit_lbspr_one(L, lh_row)   # res$SPR is standard LBSPR; res$SPR_bind is amended
#
# fit_lbspr_one() runs the standard LBSPR fit and passes the fitted F/M, SL50, SL95
# straight into spr_sex_structured(), so the estimator itself is never changed.
