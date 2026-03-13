# dpe_random_pxsigma_restart_dpe3

Restart experiment branching from the original `dpe_random_pxsigma` run.

## Rationale

DPE2 iteration 70 was selected by examining the convergence plots - the algorithm had stabilized but not yet completed DPE2. By branching here and jumping directly to DPE3, we test whether the additional DPE2 iterations affect the final factor structure.

## Setup

- **Source checkpoint:** `dpe2_iter70.rds` (v0=0.01)
- **Original schedule:** [0.1, 0.01, 0.001, 0.0001]
- **Restart schedule:** [0.001, 0.0001] (DPE3 and DPE4 only)
