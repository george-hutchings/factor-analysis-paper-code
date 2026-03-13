# check_current_clinical_inputs.R
# Rebuild the shared inputs from raw sources and compare them to the files
# currently used by compression / coxph / EDSS correlation.
# Read-only: prints diagnostics, saves nothing.

TRAINING_DATA <- "/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds"
CORE_LONGITUDINAL <- "/data/ms/processed/clinical/projects/longitudinal/longitudinal-table-CORE.csv"
SCAN_DIR <- "/data/ms/processed/mri/MS_Share/4George/T2lesion_ind"
MRI_CLINICAL <- "/data/users/uu85g9/factor-analysis-old/mri-clinical/mri-clinical_processed.rds"
SUBJECTS_FILE <- "/data/users/uu85g9/factor-analysis-old/mri-clinical/mri-clinical_clinical-only_processed.rds"

clinical_vars <- c("EDSS", "VISFNC", "BRNFNC", "BOWFNC", "CLRFNC", "CRBFNC",
                   "PYRFNC", "SENFNC", "T25FWM", "HPT9M", "PASAT")

rule <- function() cat(strrep("-", 72), "\n")
section <- function(title) { cat("\n"); rule(); cat(" ", title, "\n"); rule(); cat("\n") }
info <- function(x) cat("  [INFO]", x, "\n")
pass <- function(x) cat("  [PASS]", x, "\n")
fail <- function(x) cat("  [FAIL]", x, "\n")

n_checks <- 0
n_pass <- 0
check <- function(ok, yes, no) {
  n_checks <<- n_checks + 1
  if (ok) {
    n_pass <<- n_pass + 1
    pass(yes)
  } else {
    fail(no)
  }
}

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

section("1. Load data")
for (f in c(TRAINING_DATA, CORE_LONGITUDINAL, MRI_CLINICAL, SUBJECTS_FILE))
  if (!file.exists(f)) stop("File not found: ", f)
if (!dir.exists(SCAN_DIR)) stop("Scan dir not found: ", SCAN_DIR)

Y_full <- readRDS(TRAINING_DATA)
core <- read.csv(CORE_LONGITUDINAL)
current_dat <- readRDS(MRI_CLINICAL)
current_subjects <- readRDS(SUBJECTS_FILE)
scan_ids <- get_scan_ids(SCAN_DIR)

if (!"USUBJID" %in% colnames(current_subjects)) stop("USUBJID column missing from subjects file")

current_Y <- as.matrix(current_dat[, grepl("^X\\d+", colnames(current_dat))])
current_ids <- as.character(current_subjects$USUBJID)

info(sprintf("Full MRI Y: %d x %d", nrow(Y_full), ncol(Y_full)))
info(sprintf("Current MRI-clinical Y: %d x %d", nrow(current_Y), ncol(current_Y)))
info(sprintf("Current subjects_df: %d rows", nrow(current_subjects)))
info(sprintf("Scan-order IDs: %d", length(scan_ids)))

section("2. Rebuild from raw sources")
built_subjects <- core[
  core$STUDY %in% c("CFTY720D2301", "CFTY720D2309") & core$DAY == 1,
  c("USUBJID", clinical_vars),
  drop = FALSE
]

dup_ids <- unique(built_subjects$USUBJID[duplicated(built_subjects$USUBJID)])
if (length(dup_ids) > 0)
  stop("Duplicate baseline rows for USUBJID: ", paste(dup_ids, collapse = ", "))

built_subjects <- built_subjects[match(scan_ids, built_subjects$USUBJID), , drop = FALSE]
keep <- complete.cases(built_subjects[, clinical_vars, drop = FALSE])

built_ids <- as.character(built_subjects$USUBJID[keep])
built_subjects <- built_subjects[keep, , drop = FALSE]
built_Y <- Y_full[keep, , drop = FALSE]
dropped_ids <- scan_ids[!keep]

info(sprintf("Built subjects_df: %d rows", nrow(built_subjects)))
info(sprintf("Built Y: %d x %d", nrow(built_Y), ncol(built_Y)))
info(sprintf("Dropped subjects from raw baseline data: %d", length(dropped_ids)))
if (length(dropped_ids) <= 30)
  info(sprintf("Dropped IDs: %s", paste(dropped_ids, collapse = ", ")))

section("3. Genuine ID alignment")
check(length(scan_ids) == nrow(Y_full),
      sprintf("Scan-order IDs match full MRI rows: %d", length(scan_ids)),
      sprintf("Scan-order IDs=%d vs full MRI rows=%d", length(scan_ids), nrow(Y_full)))
check(length(unique(scan_ids)) == length(scan_ids),
      "Scan-order IDs are unique",
      "Duplicate scan-order IDs found")
check(identical(built_ids, current_ids),
      "Raw-built subject IDs exactly match current subjects file",
      "Raw-built subject IDs do not match current subjects file")
check(identical(current_ids, scan_ids[keep]),
      "Current subject IDs are exactly the kept scan-order IDs",
      "Current subject IDs are not exactly the kept scan-order IDs")

section("4. MRI matrix")
check(ncol(built_Y) == ncol(current_Y),
      sprintf("MRI column counts match: %d", ncol(built_Y)),
      sprintf("Built=%d vs current=%d", ncol(built_Y), ncol(current_Y)))
check(identical(unname(built_Y), unname(current_Y)),
      "Raw-built Y exactly matches current MRI-clinical MRI matrix",
      "Raw-built Y does not match current MRI-clinical MRI matrix")

section("Summary")
cat(sprintf("  Checks passed: %d / %d\n", n_pass, n_checks))
if (n_pass == n_checks) cat("\n  ALL CHECKS PASSED\n") else
  cat(sprintf("\n  WARNING: %d check(s) FAILED\n", n_checks - n_pass))
cat("\n  No files were written.\n")
