# Projected Effects of Naloxone and Buprenorphine Scale-Up on Opioid Overdose Mortality in Six Pennsylvania Counties

R code for a Bayesian-calibrated Markov cohort model that estimates
**county-specific** effects of harm-reduction (naloxone) and treatment
(buprenorphine) interventions on opioid overdose mortality in six Pennsylvania
counties: Allegheny and Philadelphia (large urban), Dauphin and Erie
(intermediate-sized), and Clearfield and Columbia (rural).

The model is calibrated separately to each county's annual overdose death rates
and used to project overdose mortality over a five-year horizon (2025–2029)
under proportional increases (+10 % / +20 % / +30 %) in naloxone and
buprenorphine dispensing.

---

## Overview

The analysis runs in two stages:

1. **Calibration** (`calibration_PA.R`). A seven-state Markov cohort model
   (monthly cycles) representing progression through opioid use, treatment, and
   overdose death is calibrated, one county at a time, to observed annual
   overdose death rates (2018–2022) using Incremental Mixture Importance
   Sampling (IMIS). Seven parameters are estimated: three logistic intercepts
   and four covariate slopes. The script writes a posterior resample per county.

2. **Projection / trajectory** (`intervention_trajectory_analysis.R`).
   The posterior resamples are loaded, the model is re-run over the warm-up
   years, and death rates are projected to 2029 under naloxone, buprenorphine,
   and combined intervention scenarios. The script saves the observed-and-
   projected trajectory figures (Figures 2–4 of the manuscript).

### Model states

Non-use (NU) → prescription-opioid use (PU) → misuse (MU) → opioid use disorder
(OUD) → treatment (MOUD), with overdose death and non-overdose death as
absorbing states.

### Calibrated transition probabilities

Three monthly transition probabilities are logistic functions of standardized
county-level covariates (eAppendix 2 of the manuscript):

```
logit(p1) = b0 + b1 * opioid_dispensing_rate
logit(p2) = b2 + b3 * buprenorphine_dispensing_rate
logit(p3) = b4 + b5 * pct_seizures_fentanyl  -  b6 * naloxone
```

- `p1` — non-user → prescription-opioid user
- `p2` — OUD → MOUD (treatment)
- `p3` — OUD → overdose death

All covariates are z-standardized using county-specific means and standard
deviations computed over the 2018–2022 calibration window.

### Naloxone as a depreciating stock

Because naloxone kits can remain available across multiple years rather than
only in the year dispensed, naloxone enters the model as a depreciating
**stock** rather than a single-year flow:

```
stock[t] = delta * stock[t-1] + flow[t]
```

with `delta = 0.45`. Each year's stock carries forward 45 % of the previous
year's stock plus new dispensing; the remainder is lost from circulation.
This annual carryover corresponds to the ~6.5 % monthly kit depreciation
estimated by Zang et al. (2022). The stock is warmed up from 2007 using the
full dispensing history before the calibration window, and `delta` must be set
identically in both scripts.

---

## Repository structure

```
.
├── calibration_PA.R                            # Stage 1: per-county IMIS calibration
├── intervention_trajectory_analysis.R # Stage 2: projection + trajectory figures
└── README.md
```

Not included (see [Data files](#data-files)):

- `PA_OUD_Covariates_2007_2024.xlsx` — the covariate workbook (licensed data;
  user-supplied).
- The **IMIS** package source — obtained separately; not redistributed.
- `routput/` and `trajectory_results/` — created at runtime for calibration
  posteriors/diagnostics and the figures, respectively.

---

## Software and package requirements

- **R** (>= 4.1 recommended).
- CRAN packages:
  - `readxl` — read the covariate workbook
  - `lhs` — Latin hypercube prior sampling
  - `ggplot2`, `reshape2` — plotting / reshaping (calibration)
  - `dplyr`, `gridExtra`, `grid`, `cowplot` — plotting / layout (trajectory)
  - `devtools`, `mvtnorm`, `matrixStats` — calibration support
- **IMIS**: the Incremental Mixture Importance Sampling routine (Raftery & Bao,
  2010) is an external package that is **not** redistributed here. Obtain it
  separately, then point the `load_all("path\\IMIS")` line in the IMIS section
  of `calibration_PA.R` at your local copy of the package source and call
  `library(IMIS)`.

Install the CRAN packages with:

```r
install.packages(c("readxl","lhs","ggplot2","reshape2","dplyr",
                   "gridExtra","cowplot","devtools","mvtnorm","matrixStats"))
```

---

## Data files

Two kinds of input are needed. One is embedded in the code; the other is
user-supplied.

### Observed overdose death targets (in-code)

The calibration targets — annual opioid-involved overdose death rates per
100,000 — are hard-coded in the `TARGETS` list in `calibration_PA.R` (and
mirrored as `obs=` in the trajectory script). For 2018 onward these are drawn
from **CDC WONDER Multiple Cause of Death data**, which is publicly available.
Deaths are identified with ICD-10 codes for opioid-related poisoning:

- **Underlying cause of death** (drug poisoning): `X40`–`X44`, `X60`–`X64`,
  `X85`, `Y10`–`Y14`
- **Multiple cause of death** (opioid involvement): `T40.0`, `T40.1`, `T40.2`,
  `T40.3`, `T40.4`, `T40.6`

`T40.4` alone identifies synthetic-opioid (fentanyl-type) involvement. Crude
rates are used. Small county/year cells are suppressed by CDC WONDER when fewer
than 10 deaths are observed (e.g. Columbia and Clearfield have no 2024 value)
and are therefore absent from `TARGETS`. Calibration uses 2018–2022; 2023–2024
are reserved for out-of-sample validation.

### Covariate workbook (user-supplied)

`PA_OUD_Covariates_2007_2024.xlsx` is **not** included, because it is built
from licensed data (opioid, buprenorphine, and naloxone dispensing from IQVIA;
fentanyl prevalence from NFLIS). You must supply this workbook yourself. Its
structure, which both scripts expect, is:

- **One worksheet per county**, named exactly: `Philly`, `Allegheny`,
  `Dauphin`, `Erie`, `Columbia`, `Clearfield`.
- The first **3 rows** of each sheet are headers and are skipped
  (`skip = 3`). Data begin on row 4, one row per year (2007–2024).
- **10 columns, in this order** (no header row is read; order matters):

  | # | Column | Meaning |
  |---|--------|---------|
  | 1 | year | calendar year (integer, 2007–2024) |
  | 2 | op_rx_rate | opioid dispensing rate (prescriptions per 100,000/yr) |
  | 3 | op_rx_stand | opioid dispensing, standardized (not used; re-standardized internally) |
  | 4 | nal_rate | naloxone dispensing rate (raw annual **flow**, per 100,000/yr) |
  | 5 | nal_stand | naloxone, standardized (not used; stock is built and standardized internally) |
  | 6 | bup_rate | buprenorphine dispensing rate (per 100,000/yr) |
  | 7 | bup_stand | buprenorphine, standardized (not used; re-standardized internally) |
  | 8 | fent_rate | fentanyl prevalence — a **proportion** (share of drug-seizure items; NFLIS, state-level), **not** a rate |
  | 9 | fent_stand | fentanyl, standardized (not used; re-standardized internally) |
  | 10 | source | free-text data-source note |

Notes:

- The scripts read the **raw** columns (2, 4, 6, 8) and standardize internally
  over the calibration window, so the pre-standardized columns (3, 5, 7, 9) are
  placeholders and need not be populated for the code to run.
- Column 8 (`fent_rate`) is a **proportion**, not a dispensing rate, despite the
  column name. Because NFLIS reports at the state level, this measure is common
  to all counties within a state.
- Missing fentanyl values in the most recent years are forward-filled from the
  last observed value; naloxone flow is converted to a depreciating stock
  before standardization.

Approximate covariate ranges across the six counties (manuscript Table 1):
naloxone dispensing 39–59,893; buprenorphine 1,548–252,863; opioid
38,357–689,553 (per 100,000/yr); fentanyl share of seizures 0.001–0.263.

---

## How to run

1. Place `PA_OUD_Covariates_2007_2024.xlsx` in the working directory.
2. Point the IMIS `load_all()` path in `calibration_PA.R` at your local IMIS
   package source.
3. **Calibrate each county.** In `calibration_PA.R`, set `COUNTY` to one of
   `Philly | Allegheny | Dauphin | Erie | Columbia | Clearfield` and run the
   whole script. Repeat for all six counties. Each run writes
   `routput/imis_resample_<county>2.csv` (plus diagnostic CSVs and plots).
4. **Project trajectories.** Run `intervention_trajectory_analysis.R`.
   It reads the six posterior files from `routput/` and writes the figures to
   `trajectory_results/`.

Configuration lives in the `CONFIG` block near the top of each script
(calibration window, naloxone stock parameters, projection method, intervention
percentages). Defaults reproduce the main analysis: calibration window
2018–2022, 2023–2024 held out for validation, naloxone stock `delta = 0.45`
(identical in both scripts), covariate projection by linear trend fit to
2018–2024, and intervention scenarios of +10 % / +20 % / +30 %. Overdose deaths
averted relative to the no-increase baseline are computed sample-wise across 250
posterior draws and summarized as means with 95 % uncertainty intervals.

A sensitivity analysis for Allegheny County calibrates to 2018–2024 (rather than
2018–2022); to reproduce it, set `CAL_END`/calibration window accordingly for
that county (manuscript eAppendix 3).

---

## Outputs

- **Calibration** (`routput/`): `imis_resample_<county>2.csv` (posterior
  resample, read by the trajectory script), plus `psa_data.csv`,
  `model_pred.csv`, `posterior_summary.csv`, and calibration / prior-posterior
  plots.
- **Trajectory** (`trajectory_results/`): PNG figures of observed and projected
  overdose death rates per 100,000, by county, for the naloxone, buprenorphine,
  and combined intervention scenarios.

---

## Data availability

- **Overdose mortality** targets are derived from the CDC Wide-ranging Online
  Data for Epidemiologic Research (CDC WONDER), Multiple Cause of Death data,
  which are publicly accessible (https://wonder.cdc.gov). Deaths are identified
  with ICD-10 codes for opioid-related poisoning (`X40–X44`, `X60–X64`, `X85`,
  `Y10–Y14`, `T40.0–T40.4`, `T40.6`); the resulting rates are embedded in the
  code.
- **Dispensing** rates for prescription opioids, naloxone, and buprenorphine
  were obtained from the IQVIA Xponent database (retail, mail-order, and
  long-term-care pharmacy channels). These data are **not** publicly available;
  access can be requested from IQVIA at
  https://www.iqvia.com/insights/the-iqvia-institute/available-iqvia-data.
- **Fentanyl** seizure data were obtained from the National Forensic Laboratory
  Information System (NFLIS, https://www.nflis.deadiversion.usdoj.gov).
- The IMIS routine is obtained separately and is not included here.

The expected covariate-workbook structure is documented above so that users
with access to equivalent sources can reconstruct the input file.
<!--
---

## Funding and acknowledgments

This work was funded by contract 75D30121C12574 from the Centers for Disease
Control and Prevention (CDC) and by the National Science Foundation under grant
agreement CMMI-2240408. The findings and conclusions are those of the authors
and do not necessarily represent the official position of the CDC or the views
of the U.S. National Science Foundation.

The University of Pittsburgh Institutional Review Board determined that this
study did not constitute human-subjects research.
-->

---

## Citation

<!--
This manuscript is currently under peer review. The citation below is a
placeholder and will be updated upon publication.

Ahmed, A. A., Rahimian, M. A., Chen, Q., & Kumar, P. County-Level Heterogeneity
in Opioid Harm Reduction and Treatment Effects: A Simulation Modeling Analysis.
(Under review.)

BibTeX:

@article{ahmed_county_heterogeneity,
  title   = {County-Level Heterogeneity in Opioid Harm Reduction and Treatment
             Effects: A Simulation Modeling Analysis},
  author  = {Ahmed, Abdulrahman A. and Rahimian, M. Amin and Chen, Qiushi and
             Kumar, Praveen},
  note    = {Manuscript under review},
  year    = {2026}
}
-->

<!--
*This manuscript is currently under review; full citation details will be
provided once the process is finalized.*
-->

If you use this code, please cite:

Ahmed, A. A., Rahimian, M. A., Chen, Q., & Kumar, P. (2026). Projected Effects of Naloxone and Buprenorphine Scale-Up on Opioid Overdose Mortality in Six Pennsylvania Counties. *medRxiv* [Preprint]. https://doi.org/10.64898/2026.05.11.26352900

```bibtex
@article{ahmed2026projected,
  title   = {Projected Effects of Naloxone and Buprenorphine Scale-Up on Opioid
             Overdose Mortality in Six {Pennsylvania} Counties},
  author  = {Ahmed, Abdulrahman A. and Rahimian, M. Amin and Chen, Qiushi and
             Kumar, Praveen},
  journal = {medRxiv},
  year    = {2026},
  doi     = {10.64898/2026.05.11.26352900},
  url     = {https://doi.org/10.64898/2026.05.11.26352900},
  note    = {Preprint}
}


