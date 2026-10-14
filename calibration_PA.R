# =============================================================
# PA COUNTY OUD CALIBRATION
#
# Bayesian calibration (IMIS) of a 7-state Markov cohort model of
# opioid use and overdose mortality, for one Pennsylvania county.
#
# WHAT THIS SCRIPT DOES
#   1. Reads annual county covariates from an Excel workbook
#      (one sheet per county) and the observed overdose death
#      targets hard-coded in the TARGETS list below.
#   2. Standardizes covariates over the calibration window.
#   3. Calibrates seven parameters (three logistic intercepts and
#      four covariate slopes) with Incremental Mixture Importance
#      Sampling (IMIS).
#   4. Writes the posterior resample and diagnostic/plot outputs
#      to OUTPUT_DIR.
#
# HOW TO RUN
#   This script calibrates ONE county per run. Set COUNTY below,
#   then run the whole file. Repeat for all six counties; each run
#   writes imis_resample_<county>2.csv, which the trajectory
#   script reads. See the README for data-file structure and the
#   full run order.
#
# DATA / DEPENDENCIES
#   - PA_OUD_Covariates_2007_2024.xlsx is user-supplied; its
#     covariates derive from licensed sources (see README) and are
#     not redistributed here.
#   - The IMIS routine is an external package obtained separately
#     (see the load_all() line in the IMIS section) and is not
#     redistributed here.
# =============================================================
set.seed(1234)
library(readxl); library(lhs); library(ggplot2); library(reshape2)

# =============================================================
# CONFIGURATION (edit here only)
# =============================================================
# County lookup - all 6 counties with population and sheet name
COUNTY_INFO <- list(
  Philly     = list(sheet = "Philly",     pop = 1583219),
  Allegheny  = list(sheet = "Allegheny",  pop = 1248340),
  Dauphin    = list(sheet = "Dauphin",    pop = 276645),
  Erie       = list(sheet = "Erie",       pop = 279609),
  Columbia   = list(sheet = "Columbia",   pop = 66647),
  Clearfield = list(sheet = "Clearfield", pop = 81708)
)

COUNTY     <- "Clearfield"   # change only this: Philly|Allegheny|Dauphin|Erie|Columbia|Clearfield
SHEET_NAME <- COUNTY_INFO[[COUNTY]]$sheet
COUNTY_POP <- COUNTY_INFO[[COUNTY]]$pop

COV_FILE   <- "PA_OUD_Covariates_2007_2024.xlsx"

# Calibration window: choose any contiguous years in 2007-2024
CAL_START  <- 2018               # first calibration year
CAL_END    <- 2022               # last  calibration year

# -- NALOXONE ACCUMULATION CONFIG -----------------------------
# Treat naloxone as a STOCK that persists across years instead of a yearly FLOW.
#   "flow"  : original behaviour (this year's dispensing only)
#   "cumA"  : Option A - pure cumulative sum (delta = 1, never depreciates)
#   "cumB"  : Option B - depreciating stock: stock[t] = delta*stock[t-1] + flow[t]
# delta = 0.45 corresponds to ~6.5%/month loss of circulating kits
# (0.935^12 ~ 0.45); see README for the source.
NAL_MODE         <- "cumB"       # "flow" | "cumA" | "cumB"   (default: cumB)
NAL_DELTA        <- 0.45         # carryover fraction for cumB (0-1); shelf-life of kits
NAL_WARMUP_START <- 2007         # year to begin accumulating from (uses full history)
# -------------------------------------------------------------

# All observed overdose death rates (per 100k) by county and year
# 2015-2017: original calibration values (different/older definition - see seam caveat)
# 2018-2024: CDC WONDER D157, T40.0-T40.4 + T40.6 opioid-involved, crude rate
TARGETS <- list(
  Philly     = c("2015"=48.33, "2016"=61.96, "2017"=88.70,
                 "2018"=47.0,  "2019"=51.3,  "2020"=57.8, "2021"=59.6,
                 "2022"=65.7,  "2023"=60.8,  "2024"=44.9),
  Allegheny  = c("2015"=41.76, "2016"=63.09, "2017"=77.14,
                 "2018"=32.1,  "2019"=36.9,  "2020"=47.0, "2021"=45.8,
                 "2022"=45.7,  "2023"=43.3,  "2024"=25.7),
  Dauphin    = c("2015"=24.93, "2016"=34.07, "2017"=38.27,
                 "2018"=20.6,  "2019"=20.1,  "2020"=12.5, "2021"=12.2,
                 "2022"=14.5,  "2023"=12.4,  "2024"=6.1),
  Erie       = c("2015"=23.30, "2016"=34.20, "2017"=48.68,
                 "2018"=20.6,  "2019"=21.5,  "2020"=24.2, "2021"=29.7,
                 "2022"=34.7,  "2023"=35.9,  "2024"=19.4),
  Columbia   = c("2015"=20.37, "2016"=24.29, "2017"=26.37,
                 "2018"=16.8,  "2019"=15.4,  "2020"=23.1, "2021"=32.4,
                 "2022"=20.0,  "2023"=27.5),   # no 2024 in CDC pull
  Clearfield = c("2015"=18.23, "2016"=18.29,   # 2017 suppressed (small N)
                 "2018"=18.9,  "2019"=13.9,  "2020"=17.8, "2021"=31.2,
                 "2022"=19.3,  "2023"=29.8)    # no 2024 in CDC pull
)

# Auto-slice targets for the chosen county and calibration window
target <- unname(TARGETS[[COUNTY]][as.character(CAL_START:CAL_END)])
if (anyNA(target)) stop(paste("No target data for", COUNTY, "in years", CAL_START, "-", CAL_END,
                              "- check TARGETS list"))

# Prior bounds for 7 parameters:
# [1] NU->PU intercept  [2] OUD->MOUD intercept  [3] OUD->DeathOD intercept
# [4] opioid Rx slope   [5] buprenorphine slope   [6] fentanyl slope   [7] naloxone slope
PRIORS <- list(
  #          [1]      [2]     [3]    [4]   [5]   [6]   [7]
  Philly     = list(lower = c(-15.40, -15.5,  -9.0,  0.2,  0.1,  0.0,  0.0),
                    upper = c(  5.16, -4.19, -4.00,  1.0,  1.4,  1.9,  0.6)),
  Allegheny  = list(lower = c(-9.40, -19.5,  -12.5,  0.15,  0.2,  0.0,  0.0),
                    upper = c(  2.16,   -2.19,   -3.0,  2.5,  2.2,  0.4,  0.6)),
  Dauphin    = list(lower = c(-10.40,  -7.5,  -8.0,  0.3,  0.0,  0.0,  0.0),
                    upper = c(  0.16,   1.19,   0.53, 0.8,  0.5,  0.9,  1.5)),
  Erie       = list(lower = c(-10.40,  -8.5,  -9.0,  0.3,  0.0,  0.1,  0.05),
                    upper = c(  0.16,   0.19,  -6.53, 0.8,  0.7,  0.9,  0.3)),
  Columbia   = list(lower = c(-10.40, -10.5, -8.0,  0.3,  0.0,  0.0,  0.0),
                    upper = c( 0.16,   2.19,   -3.53, 2.8,  0.5,  0.6,  0.4)),
  Clearfield = list(lower = c(-8.40,  -9.5, -7.0,  0.3,  0.0,  0.3,  0.0),
                    upper = c( -1.16,   -1.19,   -3.53, 1.0,  0.9,  1.3,  0.9))
)

v1_lower <- matrix(PRIORS[[COUNTY]]$lower, nrow=1, ncol=7)
v1_upper <- matrix(PRIORS[[COUNTY]]$upper, nrow=1, ncol=7)

PSA_SIZE   <- 250
PSA_SEED   <- 23
IMIS_B     <- 1500               # IMIS initial sample size
IMIS_B_RE  <- 3000               # IMIS resample size
IMIS_K     <- 150                # IMIS iterations
OUTPUT_DIR <- "routput"
# -------------------------------------------------------------

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
stopifnot(length(target) == (CAL_END - CAL_START + 1))

# =============================================================
# LOAD & STANDARDIZE COVARIATES
# =============================================================
raw <- read_excel(COV_FILE, sheet = SHEET_NAME, skip = 3, col_names = FALSE)
# Columns: year, op_rx_rate, op_rx_stand, nal_rate, nal_stand,
#          bup_rate, bup_stand, fent_rate, fent_stand, source
names(raw) <- c("year","op_rx_rate","op_rx_stand","nal_rate","nal_stand",
                "bup_rate","bup_stand","fent_rate","fent_stand","source")
raw$year <- as.integer(raw$year)
raw <- raw[!is.na(raw$year), ]
raw <- raw[order(raw$year), ]   # ensure year-sorted for the stock recursion

# Forward-fill fentanyl NAs (2023-2024 missing in source data)
# Use last known value (2022) as flat projection
last_fent <- tail(raw$fent_rate[!is.na(raw$fent_rate)], 1)
raw$fent_rate[is.na(raw$fent_rate)] <- last_fent

# -- Convert naloxone flow -> stock on the FULL raw history ----
# Done BEFORE windowing so the stock warms up from NAL_WARMUP_START.
build_nal_stock <- function(years, flow, mode, delta, warm_start) {
  keep <- years >= warm_start
  yy <- years[keep]; ff <- flow[keep]
  ff[is.na(ff)] <- 0
  stock <- numeric(length(ff))
  if (mode == "flow") {
    stock <- ff
  } else if (mode == "cumA") {
    stock <- cumsum(ff)
  } else if (mode == "cumB") {
    s <- 0
    for (i in seq_along(ff)) { s <- delta * s + ff[i]; stock[i] <- s }
  } else stop("Unknown NAL_MODE: ", mode)
  out <- rep(NA_real_, length(years))
  out[keep] <- stock
  out
}
raw$nal_input <- build_nal_stock(raw$year, raw$nal_rate,
                                 NAL_MODE, NAL_DELTA, NAL_WARMUP_START)

cat(sprintf("Naloxone mode: %s%s | warm-up from %d\n",
            NAL_MODE,
            if (NAL_MODE=="cumB") sprintf(" (delta=%.2f)", NAL_DELTA) else "",
            NAL_WARMUP_START))
cat("Naloxone flow vs stock (full history):\n")
print(data.frame(year=raw$year,
                 flow=round(raw$nal_rate,1),
                 stock=round(raw$nal_input,1)))

# Standardize using the calibration window as reference
REF_START <- CAL_START; REF_END <- CAL_END
ref_data  <- raw[raw$year >= REF_START & raw$year <= REF_END, ]
stopifnot(nrow(ref_data) >= 2)

standardize_with <- function(x, ref) {
  m <- mean(ref, na.rm=TRUE); s <- sd(ref, na.rm=TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

# Filter to calibration window for model input
cal_data <- raw[raw$year >= CAL_START & raw$year <= CAL_END, ]
stopifnot(nrow(cal_data) == (CAL_END - CAL_START + 1))

cal_data$op_rx_std <- standardize_with(cal_data$op_rx_rate, ref_data$op_rx_rate)
cal_data$nal_std   <- standardize_with(cal_data$nal_input,  ref_data$nal_input)  # STOCK, not flow
cal_data$bup_std   <- standardize_with(cal_data$bup_rate,   ref_data$bup_rate)
cal_data$fent_std  <- standardize_with(cal_data$fent_rate,  ref_data$fent_rate)

# Extract covariate vectors (one value per calibration year, in order)
n_years       <- nrow(cal_data)
opioid_rx     <- cal_data$op_rx_std
naloxone_rate <- cal_data$nal_std
bupre_rate    <- cal_data$bup_std
fent_rate_vec <- cal_data$fent_std

cat(sprintf("\nCounty: %s | Window: %d-%d (%d years)\n",
            COUNTY, CAL_START, CAL_END, n_years))
cat("Standardized covariates (naloxone = standardized STOCK):\n")
print(data.frame(year=cal_data$year,
                 op_rx=round(opioid_rx,4), nal=round(naloxone_rate,4),
                 bup=round(bupre_rate,4), fent=round(fent_rate_vec,4)))

# Standard error for likelihood
std_target <- sqrt(target * COUNTY_POP / 1e5) / COUNTY_POP * 1e5

# =============================================================
# MARKOV MODEL
# =============================================================
logi_fun <- function(x) 1 / (1 + exp(-x))

oud_model <- function(params) {
  n_cycles <- 12 * n_years
  pr <- numeric(3)

  v_m_init <- c(NU=0.909, PU=0.006, MU=0.04, OUD=0.04, OUD_treat=0.005, DeathOd=0, Death=0)
  m_M <- matrix(NA, nrow=n_cycles+1, ncol=7,
                dimnames=list(0:n_cycles, names(v_m_init)))
  m_M[1,] <- v_m_init

  # Fixed transition matrix
  m_P <- matrix(0, 7, 7, dimnames=list(names(v_m_init), names(v_m_init)))
  m_P["PU",]       <- c(0.373, 0, 0.05, 0.0012778, 0, 0, 0.000238)
  m_P["PU","PU"]   <- max(0, 1 - sum(m_P["PU",]))
  m_P["MU",]       <- c(0.039, 0, 0, 0.023, 0, 1.3e-7, 0.000238)
  m_P["MU","MU"]   <- max(0, 1 - sum(m_P["MU",]))
  m_P["OUD",]      <- c(0.0162, 0, 0.013, 0, 0, 0, 0.00258)
  m_P["OUD_treat",]<- c(0.009, 0, 0, 0.0089, 0, 0, 0.000238)
  m_P["OUD_treat","OUD_treat"] <- max(0, 1 - sum(m_P["OUD_treat",]))
  m_P["DeathOd","DeathOd"] <- 1
  m_P["Death","Death"]     <- 1

  # Year-change cycles: update matrix at start of each year
  year_change <- (0:(n_years-1)) * 12 + 1   # cycles 1, 13, 25, ...

  update_matrix <- function(yr_idx) {
    pr1 <- logi_fun(params[1] + params[4] * opioid_rx[yr_idx])
    pr2 <- logi_fun(params[2] + params[5] * bupre_rate[yr_idx])
    pr3 <- logi_fun(params[3] + params[6] * fent_rate_vec[yr_idx]
                    - params[7] * naloxone_rate[yr_idx])
    m_P["NU","PU"]        <<- pr1
    m_P["NU","NU"]        <<- max(0, 1 - (m_P["NU","MU"] + m_P["NU","OUD"] + pr1))
    m_P["OUD","OUD_treat"]<<- pr2
    m_P["OUD","DeathOd"]  <<- pr3
    m_P["OUD","OUD"]      <<- max(0, 1 - sum(m_P["OUD",]) + m_P["OUD","OUD"])
  }

  update_matrix(1)   # initialise before first cycle

  for (t in 1:n_cycles) {
    m_M[t+1,] <- m_M[t,] %*% m_P
    if (t %in% year_change) {
      yr_idx <- which(year_change == t)
      update_matrix(yr_idx)
    }
  }

  m_M <- cbind(m_M, alive = rowSums(m_M[,1:5]))

  rate_years <- numeric(n_years)
  od_prev    <- 0
  for (i in 1:n_years) {
    k <- (i-1) * 12
    rate_years[i] <- (m_M[i*12, 6] - od_prev) /
      mean(m_M[k:(i*12), 8]) * 1e5
    od_prev <- m_M[i*12, 6]
  }
  return(rate_years)
}

# =============================================================
# PRIOR SAMPLING
# =============================================================
n_calibrated_params <- 7

sample.prior.lhs <- function(n) {
  A <- randomLHS(n=n, k=n_calibrated_params)
  draws <- data.frame(
    NU_PU_intercept       = qunif(A[,1], v1_lower[1], v1_upper[1]),
    OUD_Rx_intercept      = qunif(A[,2], v1_lower[2], v1_upper[2]),
    OUD_DeathOD_intercept = qunif(A[,3], v1_lower[3], v1_upper[3]),
    NU_PU_slope           = qunif(A[,4], v1_lower[4], v1_upper[4]),
    OUD_Rx_bup            = qunif(A[,5], v1_lower[5], v1_upper[5]),
    OUD_DeathOD_fent      = qunif(A[,6], v1_lower[6], v1_upper[6]),
    OUD_DeathOD_nal       = qunif(A[,7], v1_lower[7], v1_upper[7])
  )
  return(as.matrix(draws))
}

sample.prior <- sample.prior.lhs

# =============================================================
# PSA
# =============================================================
psa_run <- function(par_vector) {
  if (is.null(dim(par_vector))) par_vector <- t(par_vector)
  psa <- matrix(0, nrow=nrow(par_vector), ncol=n_years)
  for (j in seq_len(nrow(par_vector))) {
    psa[j,] <- oud_model(as.numeric(par_vector[j,]))
  }
  return(psa)
}

set.seed(1234)
input_params <- sample.prior.lhs(PSA_SIZE)
psa_data     <- psa_run(input_params)
colnames(psa_data) <- as.character(CAL_START:CAL_END)
write.csv(psa_data, file.path(OUTPUT_DIR,"psa_data.csv"), row.names=FALSE)

# Quick PSA plot
psa_long <- melt(as.data.frame(psa_data), variable.name="Year", value.name="Rate")
stats <- aggregate(Rate ~ Year, psa_long, function(x)
  c(mean=mean(x), lower=quantile(x,.025), upper=quantile(x,.975)))
stats <- do.call(data.frame, stats)
names(stats) <- c("Year","Mean","Lower","Upper")
stats$Target <- target

p_psa <- ggplot(stats, aes(x=Year, y=Mean, group=1)) +
  geom_point(color="blue", size=3) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), color="blue", width=0.2) +
  geom_point(aes(y=Target), color="red", size=3, shape=18) +
  labs(title=paste(COUNTY,"- PSA Coverage"), x="Year", y="Deaths per 100k")
print(p_psa)

# =============================================================
# LIKELIHOOD
# =============================================================
l_likelihood <- function(par_vector) {
  if (is.null(dim(par_vector))) par_vector <- t(par_vector)
  llik <- numeric(nrow(par_vector))
  for (j in seq_len(nrow(par_vector))) {
    res <- tryCatch(oud_model(as.numeric(par_vector[j,])), error=function(e) NA)
    if (anyNA(res)) { llik[j] <- -Inf; next }
    llik[j] <- sum(dnorm(target, mean=res, sd=std_target, log=TRUE))
  }
  return(llik)
}

# =============================================================
# PRIOR DENSITY (uniform over bounds)
# =============================================================
l_prior <- function(par_vector) {
  if(is.null(dim(par_vector))) par_vector <- t(par_vector)
  lprior <- rep(0,nrow(par_vector))
  lprior <- lprior+dunif(par_vector[,1],v1_lower[1],v1_upper[1],log=TRUE)
  lprior <- lprior+dunif(par_vector[,2],v1_lower[2],v1_upper[2],log=TRUE)
  lprior <- lprior+dunif(par_vector[,3],v1_lower[3],v1_upper[3],log=TRUE)
  lprior <- lprior+dunif(par_vector[,4],v1_lower[4],v1_upper[4],log=TRUE)
  lprior <- lprior+dunif(par_vector[,5],v1_lower[5],v1_upper[5],log=TRUE)
  lprior <- lprior+dunif(par_vector[,6],v1_lower[6],v1_upper[6],log=TRUE)
  lprior <- lprior+dunif(par_vector[,7],v1_lower[7],v1_upper[7],log=TRUE)
  return(lprior)
}

prior      <- function(pv) exp(l_prior(pv))
likelihood <- function(pv) exp(l_likelihood(pv))

# =============================================================
# IMIS
# =============================================================
# The IMIS routine is an external package obtained separately and
# is NOT redistributed in this repository. Point load_all() at your
# local copy of the IMIS package source. See the README for where
# to obtain it.
library(devtools)
library(mvtnorm)
load_all("path\\IMIS")   # <-- set to your local IMIS package directory
library(IMIS)

set.seed(123)
imis_res <- IMIS(B=IMIS_B, B.re=IMIS_B_RE, number_k=IMIS_K, D=0)

cat(sprintf("Unique samples: %d\n", length(unique(imis_res$resample[,1]))))
cat(sprintf("Effective N:    %.1f\n",
            sum(table(imis_res$resample[,1]))^2 / sum(table(imis_res$resample[,1])^2)))
cat(sprintf("Max weight:     %.4f\n",
            max(table(imis_res$resample[,1])) / sum(table(imis_res$resample[,1]))))

write.csv(imis_res$stat,     file.path(OUTPUT_DIR,"imis_stat.csv"))
write.csv(imis_res$resample, file.path(OUTPUT_DIR,
                                       sprintf("imis_resample_%s2.csv", tolower(COUNTY))))
write.csv(imis_res$center,   file.path(OUTPUT_DIR,"imis_center.csv"))

# =============================================================
# POSTERIOR PREDICTIONS & PLOTS
# =============================================================
resample_mat <- imis_res$resample
if (ncol(resample_mat) == n_calibrated_params + 1) {
  resample_mat <- resample_mat[, -1]   # drop leading index column if present
}
posterior_data <- as.data.frame(resample_mat)
colnames(posterior_data) <- colnames(sample.prior(1))

set.seed(PSA_SEED)
idx        <- sample(nrow(posterior_data), size=min(PSA_SIZE, nrow(posterior_data)))
samp_post  <- as.matrix(posterior_data[idx,])
model_pred <- psa_run(samp_post)
colnames(model_pred) <- as.character(CAL_START:CAL_END)
write.csv(model_pred, file.path(OUTPUT_DIR,"model_pred.csv"), row.names=FALSE)

# Calibration plot
mp_long <- melt(as.data.frame(model_pred), variable.name="Year", value.name="Rate")
mp_stats <- aggregate(Rate ~ Year, mp_long, function(x)
  c(Mean=mean(x), Lower=quantile(x,.05), Upper=quantile(x,.95)))
mp_stats <- do.call(data.frame, mp_stats)
names(mp_stats) <- c("Year","Mean","Lower","Upper")
mp_stats$Target    <- target
mp_stats$std_lower <- target - 1.96*std_target
mp_stats$std_upper <- target + 1.96*std_target

p_cal <- ggplot(mp_stats, aes(x=Year, group=1)) +
  geom_errorbar(aes(ymin=Lower,     ymax=Upper),     color="blue",  width=0.2) +
  geom_errorbar(aes(ymin=std_lower, ymax=std_upper), color="cyan3", width=0.2) +
  geom_point(aes(y=Mean),   color="black", size=2) +
  geom_point(aes(y=Target), color="red",   size=3, shape=18) +
  labs(title=paste(COUNTY,"- Model vs Observed"), x="Year", y="Deaths per 100k")
print(p_cal)
ggsave(file.path(OUTPUT_DIR,"model_rate_with_targets.pdf"),
       p_cal, width=8, height=4, dpi=400)

# Prior-posterior overlay
prior_samp          <- as.data.frame(sample.prior.lhs(PSA_SIZE))
prior_samp$type     <- "prior"
posterior_data$type <- "post"
pp_df  <- melt(rbind(prior_samp, posterior_data[,1:8]), id.vars="type")
p_pp   <- ggplot(pp_df, aes(value, fill=type, colour=type)) +
  geom_density(alpha=0.1) +
  facet_wrap(~variable, scales="free", ncol=4) +
  theme(axis.text=element_text(size=6), strip.text=element_text(size=8))
ggsave(file.path(OUTPUT_DIR,"prior_post.pdf"), p_pp, width=10, height=8, dpi=400)

# Posterior summary
library(matrixStats)
pm <- as.matrix(posterior_data[, sapply(posterior_data, is.numeric)])
post_summary <- data.frame(
  variable = colnames(pm),
  mean     = colMeans(pm),
  median   = colMedians(pm),
  sd       = apply(pm, 2, sd),
  q2.5     = apply(pm, 2, quantile, .025),
  q97.5    = apply(pm, 2, quantile, .975)
)
print(post_summary)
write.csv(post_summary, file.path(OUTPUT_DIR,"posterior_summary.csv"), row.names=FALSE)
cat("\nDone. Outputs in:", OUTPUT_DIR, "\n")
