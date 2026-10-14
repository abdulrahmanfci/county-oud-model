# =============================================================================
# INTERVENTION TRAJECTORY ANALYSIS (STOCK)
#
# Generates the observed-and-projected overdose-death trajectory figures.
# For each county this script loads the calibration posteriors, re-runs the
# Markov model over the warm-up years, projects 2025-2029 under naloxone,
# buprenorphine and combined intervention scenarios, and saves the figures
# to OUTPUT_DIR.
#
# RUN calibration_PA.R for all six counties FIRST: this script reads the
# imis_resample_<county>2.csv posterior files it produces (in POSTERIOR_DIR).
# Naloxone is standardized as a depreciating STOCK (delta must match
# calibration_PA.R). See the README for the full run order.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2)
  library(gridExtra); library(grid); library(cowplot)
})
options(warn = -1)

# =============================================================================
# CONFIG
# =============================================================================
SELECTED_COUNTIES <- NULL          # NULL = all 6; or c("Allegheny","Erie","Clearfield")
COV_PROJ_METHOD   <- "trend"       # "constant" or "trend"
FIT_YEARS         <- 2018:2022

COV_FILE      <- "PA_OUD_Covariates_2007_2024.xlsx"
POSTERIOR_DIR <- "routput"
OUTPUT_DIR    <- "trajectory_results"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ── NALOXONE STOCK CONFIG (must match calibration_PA.R) ──────────────────────
NAL_MODE         <- "cumB"         # "flow" | "cumA" | "cumB"
NAL_DELTA        <- 0.45           # carryover fraction for cumB
NAL_WARMUP_START <- 2007           # begin accumulating here
# ─────────────────────────────────────────────────────────────────────────────

NAL_PCTS <- c(0, 10, 20, 30)
BUP_PCTS <- c(0, 10, 20, 30)
COMPOUND_COMBOS <- list(
  c(nal=10, bup=10), c(nal=10, bup=30),
  c(nal=30, bup=10), c(nal=30, bup=30)
)
COMPOUND_SINGLE_NAL <- c(10, 30)
COMPOUND_SINGLE_BUP <- c(10, 30)

WU_YEARS   <- 2018:2024
CAL_YEARS  <- 2018:2022
CAL_END    <- 2022
PROJ_YEARS <- 2025:2029
N_SAMPLES  <- 250
PSA_SEED   <- 23

# =============================================================================
# COUNTY CONFIG - obs = CDC WONDER D157 T40-opioid crude rates (per 100k)
#   2015-2017 kept from older source; 2018-2024 from CDC WONDER.
#   (Only 2018+ is plotted; small counties have no 2024.)
# =============================================================================
COUNTIES <- list(
  Philly     = list(sheet="Philly",    type="Large Urban", pop=1583219,
                    obs=c("2015"=48.33,"2016"=61.96,"2017"=88.70,
                          "2018"=47.0,"2019"=51.3,"2020"=57.8,"2021"=59.6,
                          "2022"=65.7,"2023"=60.8,"2024"=44.9),
                    st_fent=1.133),
  Allegheny  = list(sheet="Allegheny", type="Large Urban", pop=1248340,
                    obs=c("2015"=41.76,"2016"=63.09,"2017"=77.14,
                          "2018"=32.1,"2019"=36.9,"2020"=47.0,"2021"=45.8,
                          "2022"=45.7,"2023"=43.3,"2024"=25.7),
                    st_fent=1.133),
  Dauphin    = list(sheet="Dauphin",   type="Medium", pop=276645,
                    obs=c("2015"=24.93,"2016"=34.07,"2017"=38.27,
                          "2018"=20.6,"2019"=20.1,"2020"=12.5,"2021"=12.2,
                          "2022"=14.5,"2023"=12.4,"2024"=6.1),
                    st_fent=1.133),
  Erie       = list(sheet="Erie",      type="Medium", pop=279609,
                    obs=c("2015"=23.30,"2016"=34.20,"2017"=48.68,
                          "2018"=20.6,"2019"=21.5,"2020"=24.2,"2021"=29.7,
                          "2022"=34.7,"2023"=35.9,"2024"=19.4),
                    st_fent=1.133),
  Columbia   = list(sheet="Columbia",  type="Rural", pop=66647,
                    obs=c("2015"=20.37,"2016"=24.29,"2017"=26.37,
                          "2018"=16.8,"2019"=15.4,"2020"=23.1,"2021"=32.4,
                          "2022"=20.0,"2023"=27.5),
                    st_fent=1.133),
  Clearfield = list(sheet="Clearfield",type="Rural", pop=81708,
                    obs=c("2015"=18.23,"2016"=18.29,
                          "2018"=18.9,"2019"=13.9,"2020"=17.8,"2021"=31.2,
                          "2022"=19.3,"2023"=29.8),
                    st_fent=1.133)
)

# Validation triangle colour (distinct from intervention oranges/greens)
VAL_COLOR <- "#b8860b"   # dark goldenrod - separate from +10% light orange
if (!is.null(SELECTED_COUNTIES)) COUNTIES <- COUNTIES[SELECTED_COUNTIES]

# =============================================================================
# NALOXONE STOCK BUILDER (identical logic to calibration_PA.R)
# =============================================================================
build_nal_stock <- function(years, flow, mode = NAL_MODE,
                            delta = NAL_DELTA, warm_start = NAL_WARMUP_START) {
  keep <- years >= warm_start
  ff <- flow[keep]; ff[is.na(ff)] <- 0
  stock <- numeric(length(ff))
  if (mode == "flow") {
    stock <- ff
  } else if (mode == "cumA") {
    stock <- cumsum(ff)
  } else if (mode == "cumB") {
    s <- 0
    for (i in seq_along(ff)) { s <- delta * s + ff[i]; stock[i] <- s }
  } else stop("Unknown NAL_MODE: ", mode)
  out <- rep(NA_real_, length(years)); out[keep] <- stock
  out
}

# Continue a stock recursion forward given a starting stock and future flows.
# start_stock = stock value at the last observed year (e.g. 2023);
# future_flow = vector of projected annual flows (2025..2029, or with gap year).
continue_nal_stock <- function(start_stock, future_flow, delta = NAL_DELTA) {
  s <- start_stock; out <- numeric(length(future_flow))
  for (i in seq_along(future_flow)) { s <- delta * s + future_flow[i]; out[i] <- s }
  out
}

# =============================================================================
# MARKOV MODEL
# =============================================================================
run_model <- function(params, st_op, st_fent, st_nal, st_bup, n_years,
                      init_state = c(0.909,0.006,0.04,0.04,0.005,0.0,0.0)) {
  n_cycles <- 12L * n_years
  m_M <- matrix(0, n_cycles+1L, 7L); m_M[1L,] <- init_state
  P <- matrix(0,7,7)
  P[2,1]<-0.373; P[2,3]<-0.05; P[2,4]<-0.0012778; P[2,7]<-0.000238
  P[2,2]<-max(0,1-P[2,1]-P[2,3]-P[2,4]-P[2,7])
  P[3,1]<-0.039; P[3,4]<-0.023; P[3,6]<-1.3e-7; P[3,7]<-0.000238
  P[3,3]<-max(0,1-P[3,1]-P[3,4]-P[3,6]-P[3,7])
  P[4,1]<-0.0162; P[4,3]<-0.013; P[4,7]<-0.00258
  P[4,4]<-max(0,1-P[4,1]-P[4,3]-P[4,7])
  P[5,1]<-0.009; P[5,4]<-0.0089; P[5,7]<-0.000238
  P[5,5]<-max(0,1-P[5,1]-P[5,4]-P[5,7])
  P[6,6]<-1; P[7,7]<-1
  logi <- function(x) 1/(1+exp(-pmin(pmax(x,-500),500)))
  yc   <- (seq_len(n_years)-1L)*12L
  upd <- function(i) {
    pr1<-logi(params[1]+params[4]*st_op[i])
    pr2<-logi(params[2]+params[5]*st_bup[i])
    pr3<-logi(params[3]+params[6]*st_fent[i]-params[7]*st_nal[i])
    fs<-0.0162+0.013+0.00258; mp23<-1-fs
    if(!is.na(pr2+pr3)&&pr2+pr3>mp23){sf<-mp23/(pr2+pr3);pr2<-pr2*sf;pr3<-pr3*sf}
    P[1,2]<<-pr1; P[1,1]<<-max(0,1-P[1,3]-P[1,4]-pr1)
    P[4,5]<<-pr2; P[4,6]<<-pr3
    P[4,4]<<-max(0,1-P[4,1]-P[4,3]-pr2-pr3-P[4,7])
  }
  upd(1L)
  for (t in seq_len(n_cycles)) {
    m_M[t+1L,] <- m_M[t,] %*% P
    if ((t-1L) %in% yc && t>1L) upd(which(yc==(t-1L)))
  }
  alive <- rowSums(m_M[,1:5]); od <- init_state[6]
  rates <- numeric(n_years)
  for (i in seq_len(n_years)) {
    k<-(i-1L)*12L
    rates[i]<-(m_M[i*12L,6]-od)/mean(alive[(k+1L):(i*12L)])*1e5
    od<-m_M[i*12L,6]
  }
  list(rates=rates, end_state=m_M[n_cycles+1L,])
}

# =============================================================================
# COVARIATES
# =============================================================================
load_covariates <- function(sheet) {
  raw <- read_excel(COV_FILE, sheet=sheet, skip=3, col_names=FALSE)
  names(raw) <- c("year","op_raw","op_std","nal_raw","nal_std",
                  "bup_raw","bup_std","fent_raw","fent_std","source")
  raw$year <- suppressWarnings(as.integer(raw$year))
  raw <- raw[!is.na(raw$year),]
  raw <- raw[order(raw$year),]
  lf <- NA_real_
  for (i in seq_len(nrow(raw))) {
    if (!is.na(raw$fent_raw[i])) lf <- raw$fent_raw[i] else raw$fent_raw[i] <- lf
  }
  # Build naloxone STOCK on full history (2007+), used everywhere below
  raw$nal_stock <- build_nal_stock(raw$year, raw$nal_raw)
  raw
}

sw <- function(x, ref) {
  ref <- ref[!is.na(ref)]; m <- mean(ref); s <- sd(ref)
  if (is.na(s)||s==0) rep(0,length(x)) else (ifelse(is.na(x),0,x)-m)/s
}

build_wu_covs <- function(sheet) {
  raw <- load_covariates(sheet)
  cal <- raw[raw$year %in% CAL_YEARS,]
  wu  <- raw[raw$year %in% WU_YEARS,]; wu <- wu[order(wu$year),]
  last_yr_row <- wu[wu$year == max(WU_YEARS), ]

  list(
    # naloxone standardized on STOCK (matches calibration); others on flow
    st_nal      = sw(wu$nal_stock, cal$nal_stock),
    st_bup      = sw(wu$bup_raw,   cal$bup_raw),
    st_fent     = sw(wu$fent_raw,  cal$fent_raw),
    st_op       = sw(ifelse(is.na(wu$op_raw),0,wu$op_raw),
                     ifelse(is.na(cal$op_raw),0,cal$op_raw)),

    # standardization params for projection
    cal_bup_m   = mean(cal$bup_raw,   na.rm=TRUE),
    cal_bup_s   = sd(cal$bup_raw,     na.rm=TRUE),
    cal_nal_m   = mean(cal$nal_stock, na.rm=TRUE),   # STOCK mean/sd
    cal_nal_s   = sd(cal$nal_stock,   na.rm=TRUE),

    # last warm-up year values (2024)
    bup_last       = last_yr_row$bup_raw,
    nal_flow_last  = last_yr_row$nal_raw,    # last-year flow (for extrapolating future flow)
    nal_stock_last = last_yr_row$nal_stock,  # last-year stock (recursion starting point)

    st_bup_last = (last_yr_row$bup_raw   - mean(cal$bup_raw,na.rm=T))   / sd(cal$bup_raw,na.rm=T),
    st_nal_last = (last_yr_row$nal_stock - mean(cal$nal_stock,na.rm=T)) / sd(cal$nal_stock,na.rm=T),
    raw         = raw
  )
}

# =============================================================================
# PROJECTION COVARIATE BUILDERS
# =============================================================================
# Buprenorphine (and any flow covariate): extrapolate flow, apply pct, standardize.
build_proj_flow <- function(raw_base, cal_m, cal_s, pct,
                            raw_df=NULL, method=COV_PROJ_METHOD) {
  n_proj <- length(PROJ_YEARS)
  if (method == "constant") {
    raw_proj <- rep(raw_base * (1 + pct/100), n_proj)
  } else {
    fit_data <- raw_df[raw_df$year %in% FIT_YEARS,]
    if (nrow(fit_data) < 2) {
      raw_proj <- rep(raw_base * (1 + pct/100), n_proj)
    } else {
      extrap <- predict(lm(value ~ year, data=fit_data),
                        newdata=data.frame(year=PROJ_YEARS))
      raw_proj <- pmax(extrap, 0) * (1 + pct/100)
    }
  }
  (raw_proj - cal_m) / cal_s
}

# Naloxone: extrapolate FLOW, apply pct to flow, run stock recursion forward
# from the last warm-up year's stock, then standardize the projected STOCK.
# Advances through any years between the last warm-up year and the projection
# years (a gap only exists if PROJ_YEARS does not start immediately after WU).
build_proj_nal_stock <- function(wu, pct, method=COV_PROJ_METHOD) {
  raw_nal_df <- data.frame(year=wu$raw$year, value=wu$raw$nal_raw)
  # Advance stock from (last warm-up year + 1) through the end of projection.
  future_years <- seq(max(WU_YEARS) + 1L, max(PROJ_YEARS))
  if (method == "constant") {
    flow_future <- rep(wu$nal_flow_last * (1 + pct/100), length(future_years))
  } else {
    fit_data <- raw_nal_df[raw_nal_df$year %in% FIT_YEARS,]
    if (nrow(fit_data) < 2) {
      flow_future <- rep(wu$nal_flow_last * (1 + pct/100), length(future_years))
    } else {
      extrap <- predict(lm(value ~ year, data=fit_data),
                        newdata=data.frame(year=future_years))
      flow_future <- pmax(extrap, 0) * (1 + pct/100)
    }
  }
  stock_future <- continue_nal_stock(wu$nal_stock_last, flow_future, NAL_DELTA)
  # Keep only PROJ_YEARS
  proj_idx   <- future_years %in% PROJ_YEARS
  stock_proj <- stock_future[proj_idx]
  (stock_proj - wu$cal_nal_m) / wu$cal_nal_s
}

# =============================================================================
# POSTERIOR LOADING
# =============================================================================
load_posterior <- function(county_lower) {
  fp   <- file.path(POSTERIOR_DIR, paste0("imis_resample_",county_lower,"2.csv"))
  raw  <- read.csv(fp)
  post <- if (ncol(raw) > 7) as.matrix(raw[,-1]) else as.matrix(raw)
  ok   <- !is.na(post[,6]) & !is.na(post[,7]) & post[,6]>0 & post[,7]>0
  post[ok,,drop=FALSE]
}

# =============================================================================
# PSA CORE
# =============================================================================
run_psa_core <- function(samp, wu, scenarios) {
  n_wu   <- length(WU_YEARS)
  n_proj <- length(PROJ_YEARS)

  wu_mat     <- matrix(NA_real_, nrow(samp), n_wu)
  end_states <- vector("list", nrow(samp))
  for (i in seq_len(nrow(samp))) {
    r <- tryCatch(run_model(samp[i,], wu$st_op, wu$st_fent,
                            wu$st_nal, wu$st_bup, n_wu), error=function(e) NULL)
    if (!is.null(r) && all(is.finite(r$rates)) && all(r$rates>=0) && all(r$rates<1000)) {
      wu_mat[i,]      <- r$rates
      end_states[[i]] <- r$end_state
    }
  }
  ok_wu <- !is.na(wu_mat[,1])
  cat(sprintf("    Warm-up: %d/%d valid\n", sum(ok_wu), nrow(samp)))

  results <- list()
  for (key in names(scenarios)) {
    sc <- scenarios[[key]]
    proj_mat <- matrix(NA_real_, nrow(samp), n_proj)
    for (i in which(ok_wu)) {
      r <- tryCatch(
        run_model(samp[i,], rep(0,n_proj), rep(sc$st_fent,n_proj),
                  sc$st_nal_v, sc$st_bup_v, n_proj, end_states[[i]]),
        error=function(e) NULL)
      if (!is.null(r) && all(is.finite(r$rates)) && all(r$rates>=0) && all(r$rates<1000))
        proj_mat[i,] <- r$rates
    }
    ok_proj  <- ok_wu & !is.na(proj_mat[,1])
    full_mat <- cbind(wu_mat[ok_proj,], proj_mat[ok_proj,])
    results[[key]] <- data.frame(
      year = c(WU_YEARS, PROJ_YEARS),
      mid  = apply(full_mat, 2, median),
      lo   = apply(full_mat, 2, quantile, 0.025),
      hi   = apply(full_mat, 2, quantile, 0.975)
    )
    cat(sprintf("    %s: %d/%d valid (%.0f%% kept)\n",
                key, sum(ok_proj), sum(ok_wu), 100*sum(ok_proj)/sum(ok_wu)))
  }
  results
}

# =============================================================================
# BUILD SCENARIOS - naloxone uses stock projection; bup uses flow projection
# =============================================================================
build_nal_scenarios <- function(wu, cfg) {
  raw_bup_df <- data.frame(year=wu$raw$year, value=wu$raw$bup_raw)
  setNames(lapply(NAL_PCTS, function(pct) {
    list(
      st_nal_v = build_proj_nal_stock(wu, pct),
      st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, 0, raw_bup_df),
      st_fent  = cfg$st_fent
    )
  }), as.character(NAL_PCTS))
}

build_bup_scenarios <- function(wu, cfg) {
  raw_bup_df <- data.frame(year=wu$raw$year, value=wu$raw$bup_raw)
  setNames(lapply(BUP_PCTS, function(pct) {
    list(
      st_nal_v = build_proj_nal_stock(wu, 0),
      st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, pct, raw_bup_df),
      st_fent  = cfg$st_fent
    )
  }), as.character(BUP_PCTS))
}

build_cmp_scenarios <- function(wu, cfg) {
  raw_bup_df <- data.frame(year=wu$raw$year, value=wu$raw$bup_raw)
  scenarios  <- list()

  scenarios[["0"]] <- list(
    st_nal_v = build_proj_nal_stock(wu, 0),
    st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, 0, raw_bup_df),
    st_fent  = cfg$st_fent
  )
  for (combo in COMPOUND_COMBOS) {
    np <- combo["nal"]; bp <- combo["bup"]
    key <- paste0("N",np,"B",bp)
    scenarios[[key]] <- list(
      st_nal_v = build_proj_nal_stock(wu, np),
      st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, bp, raw_bup_df),
      st_fent  = cfg$st_fent
    )
  }
  for (np in COMPOUND_SINGLE_NAL) {
    key <- paste0("NAL",np)
    scenarios[[key]] <- list(
      st_nal_v = build_proj_nal_stock(wu, np),
      st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, 0, raw_bup_df),
      st_fent  = cfg$st_fent
    )
  }
  for (bp in COMPOUND_SINGLE_BUP) {
    key <- paste0("BUP",bp)
    scenarios[[key]] <- list(
      st_nal_v = build_proj_nal_stock(wu, 0),
      st_bup_v = build_proj_flow(wu$bup_last, wu$cal_bup_m, wu$cal_bup_s, bp, raw_bup_df),
      st_fent  = cfg$st_fent
    )
  }
  scenarios
}

# =============================================================================
# COLORS & LABELS
# =============================================================================
gen_palette <- function(pcts, base_dark) {
  pcts_chr <- as.character(pcts)
  cols <- c("steelblue", colorRampPalette(c("#a8d5b8", base_dark))(length(pcts) - 1))
  setNames(cols, pcts_chr)
}

# Baseline = steelblue. Intervention ramps chosen to stay clearly distinct
# from steelblue at every level (avoid dark-navy that reads as blue).
NAL_PAL  <- setNames(
  c("steelblue", colorRampPalette(c("#7fc97f", "#1a7a1a"))(length(NAL_PCTS) - 1)),
  as.character(NAL_PCTS)
)  # greens for naloxone
NAL_LABS <- setNames(
  c("Baseline (no intervention)", paste0("+", NAL_PCTS[-1], "% Naloxone")),
  as.character(NAL_PCTS)
)

BUP_PAL <- setNames(
  c("steelblue", colorRampPalette(c("#f4a460", "#8b0000"))(length(BUP_PCTS) - 1)),
  as.character(BUP_PCTS)
)  # orange -> pure dark red for buprenorphine (distinct from blue)
BUP_LABS <- setNames(
  c("Baseline (no intervention)", paste0("+", BUP_PCTS[-1], "% Buprenorphine")),
  as.character(BUP_PCTS)
)

CMP_KEYS <- c(
  "0",
  sapply(COMPOUND_COMBOS, function(c) paste0("N",c["nal"],"B",c["bup"])),
  paste0("NAL", COMPOUND_SINGLE_NAL),
  paste0("BUP", COMPOUND_SINGLE_BUP)
)
n_combo <- length(COMPOUND_COMBOS)
n_solo  <- length(COMPOUND_SINGLE_NAL) + length(COMPOUND_SINGLE_BUP)
CMP_PAL <- setNames(c(
  "steelblue",
  colorRampPalette(c("#c7a8d9", "#3d1f6e"))(n_combo),
  colorRampPalette(c("#52b69a", "#1a7a4a"))(length(COMPOUND_SINGLE_NAL)),
  colorRampPalette(c("#e07b39", "#a03010"))(length(COMPOUND_SINGLE_BUP))
), CMP_KEYS)
CMP_LABS <- setNames(c(
  "Baseline (no intervention)",
  sapply(COMPOUND_COMBOS, function(c) paste0("NAL+",c["nal"],"% + BUP+",c["bup"],"%")),
  paste0("NAL+", COMPOUND_SINGLE_NAL, "% only"),
  paste0("BUP+", COMPOUND_SINGLE_BUP, "% only")
), CMP_KEYS)
CMP_LTY <- setNames(
  c("solid", rep("solid", n_combo), rep("dashed", n_solo)),
  CMP_KEYS
)

# =============================================================================
# PLOT ONE COUNTY
# =============================================================================
plot_county <- function(cname, cfg, scenario_data, pal, labs,
                        pcts, lty_map=NULL, int_title="") {
  obs_all <- data.frame(year=as.integer(names(cfg$obs)), deaths=as.numeric(cfg$obs))
  obs_all <- obs_all[obs_all$year >= 2018,]

  # Poisson-style 95% CI on the OBSERVED targets (belongs to the data, drawn black)
  # se(rate) = sqrt(rate * pop / 1e5) / pop * 1e5 ; CI = rate +/- 1.96*se
  obs_all$se     <- sqrt(obs_all$deaths * cfg$pop / 1e5) / cfg$pop * 1e5
  obs_all$obs_lo <- pmax(obs_all$deaths - 1.96 * obs_all$se, 0)
  obs_all$obs_hi <- obs_all$deaths + 1.96 * obs_all$se

  # Split: calibration period (<=2022) vs validation period (2023-2024)
  obs_cal <- obs_all[obs_all$year <= 2022, ]                    # black line
  obs_val <- obs_all[obs_all$year %in% c(2023, 2024), ]         # black triangles

  ylim_top <- max(obs_all$obs_hi, na.rm=TRUE) * 1.7

  p <- ggplot() +
    # --- OBSERVED CALIBRATION DATA (2018-2022): black line + black CI ---
    geom_errorbar(data=obs_cal, aes(x=year, ymin=obs_lo, ymax=obs_hi),
                  color="black", width=0.25, linewidth=0.5, inherit.aes=FALSE) +
    geom_line(data=obs_cal, aes(x=year, y=deaths),
              color="black", lwd=1.3, inherit.aes=FALSE) +
    geom_point(data=obs_cal, aes(x=year, y=deaths),
               color="black", size=2.0, inherit.aes=FALSE) +
    # --- OBSERVED VALIDATION DATA (2023-2024): black triangles + black CI ---
    geom_errorbar(data=obs_val, aes(x=year, ymin=obs_lo, ymax=obs_hi),
                  color="black", width=0.25, linewidth=0.5, inherit.aes=FALSE) +
    geom_point(data=obs_val, aes(x=year, y=deaths),
               color="black", shape=17, size=3.0, inherit.aes=FALSE) +

    {layers <- list()
    for (k in pcts) {
      d <- scenario_data[[k]]; if (is.null(d)) next
      col     <- pal[k]
      is_base <- (k == "0")
      lwd_k   <- if (is_base) 1.6 else 1.0
      lty_k   <- if (!is.null(lty_map) && k %in% names(lty_map)) lty_map[k] else "solid"

      d_hist   <- d[d$year %in% WU_YEARS,]
      d_proj   <- d[d$year %in% PROJ_YEARS,]
      # Branch point: last warm-up year (2024). At this year all scenarios
      # equal the baseline (intervention has not started yet), so intervention
      # lines connect smoothly from the 2024 baseline point into 2025-2029.
      branch_year <- max(WU_YEARS)
      base_branch <- scenario_data[["0"]][scenario_data[["0"]]$year == branch_year, ]

      if (is_base) {
        layers[[paste0(k,"_hr")]] <- geom_ribbon(
          data=d_hist, aes(x=year,ymin=lo,ymax=hi),
          fill=col, alpha=0.18, inherit.aes=FALSE)
        layers[[paste0(k,"_hl")]] <- geom_line(
          data=d_hist, aes(x=year,y=mid),
          color=col, lwd=lwd_k, inherit.aes=FALSE)
      }

      # Baseline: warm-up already drawn above; its projection continues from branch year.
      # Interventions: connect from the 2024 baseline point, then diverge 2025-2029.
      d_pa <- if (is_base) rbind(base_branch, d_proj) else rbind(base_branch, d_proj)
      layers[[paste0(k,"_pr")]] <- geom_ribbon(
        data=d_pa, aes(x=year,ymin=lo,ymax=hi),
        fill=col, alpha=if(lty_k=="dashed") 0.07 else 0.13, inherit.aes=FALSE)
      layers[[paste0(k,"_pl")]] <- geom_line(
        data=d_pa, aes(x=year,y=mid),
        color=col, lwd=lwd_k, linetype=lty_k, inherit.aes=FALSE)
      layers[[paste0(k,"_pp")]] <- geom_point(
        data=d_proj, aes(x=year,y=mid),
        color=col, shape=if(lty_k=="dashed") 17 else 16,
        size=1.8, inherit.aes=FALSE)
    }
    layers} +

    geom_vline(xintercept=2022.5, lty="dotted", color="darkorange", lwd=0.6, alpha=0.7) +
    geom_vline(xintercept=2024.5, lty="dashed",  color="#c0392b",   lwd=0.8, alpha=0.8) +

    annotate("text",x=2018.1,y=ylim_top*0.99,hjust=0,vjust=1,size=2.3,
             color="steelblue", label="Cal.\n2018-22") +
    annotate("text",x=2023.1,y=ylim_top*0.99,hjust=0,vjust=1,size=2.3,
             color="darkorange",label="Val.\n2023") +
    annotate("text",x=2025.1,y=ylim_top*0.99,hjust=0,vjust=1,size=2.3,
             color="#c0392b",  label="Proj.\n2025-29") +

    scale_x_continuous(breaks=c(2018,2020,2022,2024,2026,2028), limits=c(2018,2029)) +
    scale_y_continuous(limits=c(0,ylim_top), expand=expansion(mult=c(0,.02))) +
    labs(title=paste0(cname," (",cfg$type,")"), x=NULL, y=NULL) +
    theme_minimal(base_size=12) +
    theme(plot.title      = element_text(face="bold",hjust=0.5,size=14),
          panel.grid.minor= element_blank(),
          axis.text.x     = element_text(angle=30,hjust=1,size=12),
          axis.text.y     = element_text(size=12),
          axis.title      = element_text(size=13),
          legend.position = "none")
  p
}

# =============================================================================
# LEGEND
# =============================================================================
make_combined_legend <- function(pal, labs, pcts, lty_map=NULL, int_title="") {
  lab_vals <- unname(labs[pcts])
  col_vals <- unname(pal[pcts])
  lty_vals <- if (!is.null(lty_map)) unname(lty_map[pcts]) else rep("solid", length(pcts))

  sc_df <- data.frame(
    x   = rep(1:2, length(pcts)),
    y   = rep(seq_along(pcts), each=2),
    grp = rep(pcts, each=2),
    lbl = factor(rep(lab_vals, each=2), levels=lab_vals)
  )

  sc_leg <- ggplot(sc_df, aes(x=x, y=y, color=lbl, linetype=lbl, group=grp)) +
    geom_line(lwd=1.4) +
    scale_color_manual(values=setNames(col_vals, lab_vals), name=int_title) +
    scale_linetype_manual(values=setNames(lty_vals, lab_vals), name=int_title) +
    guides(color=guide_legend(ncol=2, byrow=TRUE,
                              override.aes=list(linewidth=1.2)),
           linetype=guide_legend(ncol=2, byrow=TRUE)) +
    theme_void() +
    theme(legend.position ="right",
          legend.key.width = unit(1.2,"cm"),
          legend.title     = element_text(size=12,face="bold"),
          legend.text      = element_text(size=11))

  ref_df <- data.frame(
    x   = c(1,2, 1,2, 1,2, 1,2),
    y   = c(1,1, 2,2, 3,3, 4,4),
    grp = c("a","a","b","b","c","c","d","d"),
    lbl = factor(rep(c(
      "Observed deaths, calibration 2018-2022 (CDC WONDER)",
      "Observed deaths, validation 2023-2024 (\u25b2)",
      "Observed 95% CI (black bars)",
      "Model fit / baseline (blue band = 95% posterior interval)"), each=2),
      levels=c("Observed deaths, calibration 2018-2022 (CDC WONDER)",
               "Observed deaths, validation 2023-2024 (\u25b2)",
               "Observed 95% CI (black bars)",
               "Model fit / baseline (blue band = 95% posterior interval)"))
  )
  val_pt <- ref_df[ref_df$grp=="b" & ref_df$x==1,]

  ref_leg <- ggplot(ref_df, aes(x=x,y=y,color=lbl,linetype=lbl,group=grp)) +
    geom_line(lwd=1.3) +
    geom_point(data=val_pt, aes(x=x,y=y), shape=17, size=2.5,
               color="black", inherit.aes=FALSE) +
    scale_color_manual(
      values=c("Observed deaths, calibration 2018-2022 (CDC WONDER)"="black",
               "Observed deaths, validation 2023-2024 (\u25b2)"="black",
               "Observed 95% CI (black bars)"="black",
               "Model fit / baseline (blue band = 95% posterior interval)"="steelblue"),
      name="Reference") +
    scale_linetype_manual(
      values=c("Observed deaths, calibration 2018-2022 (CDC WONDER)"="solid",
               "Observed deaths, validation 2023-2024 (\u25b2)"="blank",
               "Observed 95% CI (black bars)"="solid",
               "Model fit / baseline (blue band = 95% posterior interval)"="solid"),
      name="Reference") +
    guides(color=guide_legend(ncol=1, byrow=TRUE,
                              override.aes=list(linewidth=1.2)),
           linetype="none") +
    theme_void() +
    theme(legend.position ="right",
          legend.key.width = unit(1.2,"cm"),
          legend.title     = element_text(size=12,face="bold"),
          legend.text      = element_text(size=11))

  list(sc=cowplot::get_legend(sc_leg),
       ref=cowplot::get_legend(ref_leg))
}

# =============================================================================
# SAVE
# =============================================================================
save_plot <- function(plot_list, pal, labs, pcts, lty_map=NULL,
                      int_title="", main_title="", suffix="") {
  legs  <- make_combined_legend(pal, labs, pcts, lty_map, int_title)
  ncols <- if (length(plot_list)==1) 1 else 3
  panels <- gridExtra::arrangeGrob(grobs=plot_list, ncol=ncols)

  leg_row <- gridExtra::arrangeGrob(legs$sc, legs$ref, ncol=2, widths=c(0.7, 0.3))

  title_grob <- grid::textGrob(main_title,
                               gp=grid::gpar(fontface="bold",fontsize=14))

  x_label <- grid::textGrob("Year",
                            gp=grid::gpar(fontsize=15, fontface="plain"))
  y_label <- grid::textGrob("Annual OD Deaths per 100,000 people",
                            rot=90,
                            gp=grid::gpar(fontsize=15, fontface="plain"))

  panels_with_y <- gridExtra::arrangeGrob(y_label, panels,
                                          ncol=2, widths=c(0.02, 0.98))

  final <- gridExtra::arrangeGrob(
    title_grob,
    panels_with_y,
    x_label,
    leg_row,
    nrow=4, heights=c(0.05, 0.73, 0.03, 0.19)
  )

  out_w <- if (length(plot_list)==1) 7 else 16
  out_h <- if (length(plot_list)==1) 7 else 11

  out_file <- file.path(OUTPUT_DIR, paste0("trajectory", suffix, ".png"))
  ggsave(out_file, final, width=out_w, height=out_h, dpi=200, bg="white")
  cat(sprintf("Saved: %s\n", out_file))
}

# =============================================================================
# RUN
# =============================================================================
cat(sprintf("\n=== Trajectory v7 STOCK | Method: %s | NAL_MODE: %s (delta=%.2f) ===\n",
            COV_PROJ_METHOD, NAL_MODE, NAL_DELTA))

county_data <- list()
for (cname in names(COUNTIES)) {
  cfg <- COUNTIES[[cname]]
  cat(sprintf("\n[%s]\n", cname))
  post <- tryCatch(load_posterior(tolower(cname)),
                   error=function(e){message("  ERROR: ",e$message); NULL})
  if (is.null(post)||nrow(post)<10) next
  set.seed(PSA_SEED)
  samp <- post[sample(nrow(post),min(N_SAMPLES,nrow(post))),,drop=FALSE]
  wu   <- build_wu_covs(cfg$sheet)

  cat(sprintf("  2023 actual: nal_flow=%.1f, nal_stock=%.1f, bup=%.1f, st_nal=%.3f, st_bup=%.3f\n",
              wu$nal_flow_last, wu$nal_stock_last, wu$bup_last, wu$st_nal_last, wu$st_bup_last))

  cat("  Building NAL scenarios...\n");      nal_sc <- build_nal_scenarios(wu, cfg)
  cat("  Building BUP scenarios...\n");      bup_sc <- build_bup_scenarios(wu, cfg)
  cat("  Building Compound scenarios...\n"); cmp_sc <- build_cmp_scenarios(wu, cfg)

  cat("  Running NAL PSA...\n")
  nal_res <- tryCatch(run_psa_core(samp, wu, nal_sc), error=function(e){message(e$message);NULL})
  cat("  Running BUP PSA...\n")
  bup_res <- tryCatch(run_psa_core(samp, wu, bup_sc), error=function(e){message(e$message);NULL})
  cat("  Running Compound PSA...\n")
  cmp_res <- tryCatch(run_psa_core(samp, wu, cmp_sc), error=function(e){message(e$message);NULL})

  county_data[[cname]] <- list(cfg=cfg, nal=nal_res, bup=bup_res, cmp=cmp_res)
}

# =============================================================================
# OUTPUT
# =============================================================================
method_str <- if (COV_PROJ_METHOD=="trend") {
  paste0("Trend (",min(FIT_YEARS),"-",max(FIT_YEARS),")")
} else "Constant baseline"

county_suffix <- if (!is.null(SELECTED_COUNTIES)) paste0("_",paste(tolower(SELECTED_COUNTIES),collapse="_")) else "_all"
cov_suffix    <- paste0("_", COV_PROJ_METHOD,
                        if(COV_PROJ_METHOD=="trend") paste0("_",min(FIT_YEARS),"_",max(FIT_YEARS)) else "")

cat("\n=== Saving NAL plot ===\n")
nal_plots <- lapply(names(county_data), function(cn) {
  d <- county_data[[cn]]
  if (is.null(d$nal)) return(NULL)
  plot_county(cn, d$cfg, d$nal, NAL_PAL, NAL_LABS, as.character(NAL_PCTS))
})
names(nal_plots) <- names(county_data)
nal_plots <- Filter(Negate(is.null), nal_plots)
save_plot(nal_plots, NAL_PAL, NAL_LABS, as.character(NAL_PCTS),
          int_title="Naloxone scenarios",
          main_title=paste0("Overdose Deaths per 100,000 people: Observed & Projected\n",
                            "Harm reduction (naloxone) Dispensing Intervention"),
          suffix=paste0("_NAL", county_suffix, cov_suffix))

cat("\n=== Saving BUP plot ===\n")
bup_plots <- lapply(names(county_data), function(cn) {
  d <- county_data[[cn]]
  if (is.null(d$bup)) return(NULL)
  plot_county(cn, d$cfg, d$bup, BUP_PAL, BUP_LABS, as.character(BUP_PCTS))
})
names(bup_plots) <- names(county_data)
bup_plots <- Filter(Negate(is.null), bup_plots)
save_plot(bup_plots, BUP_PAL, BUP_LABS, as.character(BUP_PCTS),
          int_title="Buprenorphine scenarios",
          main_title=paste0("Overdose Deaths per 100,000 people: Observed & Projected\n",
                            "Treatment (Buprenorphine) Intervention"),
          suffix=paste0("_BUP", county_suffix, cov_suffix))

cat("\n=== Saving Compound plot ===\n")
cmp_plots <- lapply(names(county_data), function(cn) {
  d <- county_data[[cn]]
  if (is.null(d$cmp)) return(NULL)
  plot_county(cn, d$cfg, d$cmp, CMP_PAL, CMP_LABS, CMP_KEYS, lty_map=CMP_LTY)
})
names(cmp_plots) <- names(county_data)
cmp_plots <- Filter(Negate(is.null), cmp_plots)
save_plot(cmp_plots, CMP_PAL, CMP_LABS, CMP_KEYS, lty_map=CMP_LTY,
          int_title="Combined scenarios",
          main_title=paste0("Overdose Deaths per 100,000 people: Observed & Projected\n",
                            "Combined harm reduction (naloxone) + treatment (buprenorphine) Intervention"),
          suffix=paste0("_CMP", county_suffix, cov_suffix))

cat("\n=== DONE ===\n")
