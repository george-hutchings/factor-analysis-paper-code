# build_clinical_analysis_inputs.R
# Build and save the two shared inputs used downstream:
#   Y           - MRI matrix (subjects x voxels, scan order, complete-baseline only)
#   subjects_df - USUBJID for each row of Y
#
# Run once on the cluster to produce inputs/ before running any analysis.
# Saves: inputs/Y.rds, inputs/subjects.rds

TRAINING_DATA <- "/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds"
CORE_LONGITUDINAL <- "/data/ms/processed/clinical/projects/longitudinal/longitudinal-table-CORE.csv"
SCAN_DIR <- "/data/ms/processed/mri/MS_Share/4George/T2lesion_ind"
OUTPUT_DIR <- "/data/users/uu85g9/factor-analysis/noms/mri/inputs"

clinical_vars <- c("EDSS", "VISFNC", "BRNFNC", "BOWFNC", "CLRFNC", "CRBFNC",
                   "PYRFNC", "SENFNC", "T25FWM", "HPT9M", "PASAT")

rule <- function() cat(strrep("-", 72), "\n")
section <- function(title) { cat("\n"); rule(); cat(" ", title, "\n"); rule(); cat("\n") }
info <- function(x) cat("  [INFO]", x, "\n")

get_scan_ids <- function(scan_dir) {
  files <- list.files(scan_dir, pattern = "\\.nii\\.gz$")
  files <- grep("ses-V1x", files, value = TRUE)
  d2301 <- sort(grep("CFTY720D2301", files, value = TRUE))
  d2309 <- sort(grep("CFTY720D2309", files, value = TRUE))

  scan_id <- function(x) {
    m <- regmatches(x, regexpr("sub-.+?_ses", x))
    gsub("x", "_", sub("_ses$", "", sub("^sub-", "", m)))
  }

  c(sapply(d2301, scan_id, USE.NAMES = FALSE),
    sapply(d2309, scan_id, USE.NAMES = FALSE))
}

section("Build clinical analysis inputs")
for (f in c(TRAINING_DATA, CORE_LONGITUDINAL))
  if (!file.exists(f)) stop("File not found: ", f)
if (!dir.exists(SCAN_DIR)) stop("Scan dir not found: ", SCAN_DIR)

Y_full <- readRDS(TRAINING_DATA)
scan_ids <- get_scan_ids(SCAN_DIR)

if (nrow(Y_full) != length(scan_ids))
  stop("Full MRI rows do not match scan-order IDs: ", nrow(Y_full), " vs ", length(scan_ids))

core <- read.csv(CORE_LONGITUDINAL)
subjects_df <- core[
  core$STUDY %in% c("CFTY720D2301", "CFTY720D2309") & core$DAY == 1,
  c("USUBJID", clinical_vars),
  drop = FALSE
]

dup_ids <- unique(subjects_df$USUBJID[duplicated(subjects_df$USUBJID)])
if (length(dup_ids) > 0)
  stop("Duplicate baseline rows for USUBJID: ", paste(dup_ids, collapse = ", "))

subjects_df <- subjects_df[match(scan_ids, subjects_df$USUBJID), , drop = FALSE]
keep <- complete.cases(subjects_df[, clinical_vars, drop = FALSE])

subjects_df <- subjects_df[keep, "USUBJID", drop = FALSE]
Y <- Y_full[keep, , drop = FALSE]

info(sprintf("Built Y: %d x %d", nrow(Y), ncol(Y)))
info(sprintf("Built subjects_df: %d rows", nrow(subjects_df)))
info(sprintf("Dropped subjects: %d", sum(!keep)))
info(sprintf("Variables used to filter subjects: %s", paste(clinical_vars, collapse = ", ")))
info(sprintf("First 5 IDs: %s", paste(head(subjects_df$USUBJID, 5), collapse = ", ")))
info(sprintf("Last 5 IDs:  %s", paste(tail(subjects_df$USUBJID, 5), collapse = ", ")))

dir.create(OUTPUT_DIR, showWarnings = FALSE)
saveRDS(Y,           file.path(OUTPUT_DIR, "Y.rds"))
saveRDS(subjects_df, file.path(OUTPUT_DIR, "subjects.rds"))
info(sprintf("Saved Y.rds and subjects.rds to: %s", OUTPUT_DIR))
