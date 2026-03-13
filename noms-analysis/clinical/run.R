## Habib V2 PX clinical analysis (cluster version)
library(tidyverse)
library(haven)

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
r = slurm_id
verbose = 0

# Load base clinical data
base = read.csv('/data/ms/processed/mri/MS_Share/4George/03_baseline_results.csv')
source_saspath <- '/data/ms/unprocessed/clinical/NOMS_Version2_20211022-OXF-ANALYTICS/'

studies <- c("CFTY720D2201", "CFTY720D2301","CFTY720D2302", "CFTY720D2309", "CFTY720D2312",
             "CFTY720D2306", "CBAF312A2304", "COMB157G2301", "COMB157G2302")

# Load and process SDMT data
sdmt <- read_sas(paste0(source_saspath, "sdmt.sas7bdat")) %>% as.data.frame()
sdmt_mod <- sdmt %>%
  dplyr::select(c(USUBJID, STUDY, DATE, DAY, SDMT)) %>%
  filter(!is.na(DAY), DAY <= 1, !is.na(SDMT), STUDY %in% studies) %>%
  arrange(USUBJID, desc(DAY)) %>%
  filter(!duplicated(.$USUBJID)) %>%
  bind_rows(., sdmt[which((sdmt$DAY > 1) & (!is.na(sdmt$SDMT)) & (sdmt$STUDY %in% studies) & (!is.na(sdmt$DAY))),
                    c("USUBJID", "STUDY", "DATE", "DAY", "SDMT")]) %>%
  arrange(USUBJID, DATE, desc(SDMT)) %>%
  filter(!duplicated(.[c("USUBJID", "STUDY", "DATE")])) %>%
  rename(T25DAY = DAY) %>%
  dplyr::select(c(USUBJID, T25DAY, STUDY, SDMT)) %>%
  mutate(T25DAY = ifelse(T25DAY <= 1, 1, T25DAY)) %>%
  filter(T25DAY == 1) %>%
  dplyr::select(c("USUBJID", "STUDY", "SDMT"))

# Join SDMT to base
base = left_join(base, sdmt_mod, by = c('USUBJID', 'STUDY'))

Y_sdmt = base %>%
  dplyr::select(c("EDSS", "T25FWM", "HPT9M", "PASAT", "SDMT", "VOLT2", "NBV2", "NUMGDT1", "RELAPSE")) %>%
  as.matrix()

# Setup
D = ncol(Y_sdmt)
K = 10
set.seed(1234 + r)
lambda_init = matrix(rnorm(D * K), D, K)

v0s = c(1, .8, .4, .2, .1, .05)

# Source algorithm
source('/home/uu85g9/factor-analysis/mcem_fa_algorithm.R')

set.seed(1234 + r)
res = dpe_func(
  v0s = v0s,
  Y = Y_sdmt,
  K = K,
  v1 = 10,
  do_stick_breaking = TRUE,
  tol = 0.0000001,
  conv_param = TRUE,
  maxiter = 101,
  verbose = verbose,
  lambda = lambda_init,
  nburn = 100,
  n_mc = 100,
  px_sigma = TRUE,
  px_rotate = TRUE,
  varimax_every = -1,
  ibp_alpha = 2
)

out_dir = '/data/users/uu85g9/factor-analysis/noms/clinical'
dir.create(out_dir, showWarnings = FALSE)
saveRDS(res, file = paste0(out_dir, '/results_', r, '.rds'))
