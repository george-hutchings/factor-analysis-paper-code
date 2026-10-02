# make_pca_loading_maps.R
#
# Computes the first NPC principal component loading vectors of the NO.MS T2
# lesion matrix and writes them into the brain mask as a 4D NIfTI (one volume
# per PC), for the PCA figures (Figures 5, S8 and S9).
#
# Deterministic by construction: PCA via eigendecomposition of the small
# (n x n) gram matrix (p = 51595 >> n = 2093), PCs ordered by eigenvalue, and
# each PC signed so its dominant loading is positive. No RNG anywhere, so
# repeated runs on the same data give byte-identical output.
#
# Parameters come from environment variables (set by make_pca_supp_figure.sh)
# with sensible defaults, so the script also runs standalone: Rscript make_pca_loading_maps.R

suppressMessages(library(RNifti))

getenv <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) v else default
}

data_file <- getenv("DATA", "/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds")
mask_file <- getenv("MASK", "/data/ms/processed/mri/MS_Share/4George/final_mask.nii.gz")
outdir    <- getenv("OUTDIR", ".")
npc       <- as.integer(getenv("NPC", "3"))
doscale   <- as.logical(getenv("DOSCALE", "TRUE"))  # repo PCA convention; FALSE for raw

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

Y <- readRDS(data_file)
mask_nifti <- readNifti(mask_file)
mask <- as.logical(mask_nifti > 0)
stopifnot(ncol(Y) == sum(mask))

# PCA via the small n x n gram matrix (p >> n). The loading vectors returned
# here are identical to prcomp(Y, scale. = doscale)$rotation up to column sign.
Yc <- scale(Y, center = TRUE, scale = doscale)
G  <- tcrossprod(Yc)                                 # n x n
e  <- eigen(G, symmetric = TRUE)
U  <- e$vectors[, seq_len(npc), drop = FALSE]
d  <- sqrt(pmax(e$values[seq_len(npc)], 0))
L  <- crossprod(Yc, U) %*% diag(1 / d, npc)          # p x npc loading vectors

# deterministic sign convention: dominant loading positive (as for the MoCos in Fig 4)
for (k in seq_len(npc)) if (L[which.max(abs(L[, k])), k] < 0) L[, k] <- -L[, k]

# map each PC loading vector into the brain, write a 4D NIfTI (one vol per PC)
arr <- array(0, c(dim(mask_nifti), npc))
for (k in seq_len(npc)) {
  v <- rep(0, length(mask))
  v[mask] <- L[, k]
  arr[, , , k] <- v
}
writeNifti(arr, file.path(outdir, "pca_loadings.nii.gz"), template = mask_nifti)

# variance explained, recorded for the figure caption
ve <- 100 * e$values / sum(e$values)
writeLines(
  c("PC\tvar_explained_percent",
    sprintf("%d\t%.2f", seq_len(npc), ve[seq_len(npc)])),
  file.path(outdir, "pca_variance_explained.txt"))

cat(sprintf("wrote %s (%d volumes)\n", file.path(outdir, "pca_loadings.nii.gz"), npc))
cat("variance explained:",
    paste0(sprintf("%.1f%%", ve[seq_len(npc)]), collapse = ", "), "\n")
