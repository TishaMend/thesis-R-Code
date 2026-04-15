###############################################################################
##                                                                           ##
##   MODELLING AND FORECASTING THE UNDERWRITING CYCLE IN GENERAL INSURANCE   ##
##   Using Hidden Markov Switching Time Series Models                        ##
##                                                                           ##
##   Author       : Tisha Alpeshbhai Mendpara                                ##
##   Student ID   : 230014950                                                ##
##   Supervisor   : Dr Iqbal Owadally                                        ##
##   BSc Actuarial Science - Bayes Business School                           ##
##   City St George's, University of London                                  ##
##                                                                           ##
##   Data    : Canadian P&C loss ratios, 1985-2022 (n = 38)                  ##
##   Lines   : Personal Property (PP) and Auto                               ##
##                                                                           ##
###############################################################################
##
##  SCRIPT STRUCTURE
##  ────────────────────────────────────────────────────────────────────────
##  0.  Setup and data preparation
##  1.  Descriptive statistics
##  2.  Stationarity tests  (ADF + KPSS, levels and first differences)
##  3.  ARIMA mean models   (PP levels, PP I(1), Auto I(1))
##  4.  GARCH volatility    (PP levels, PP diff., Auto diff.)
##  5.  Spectral analysis   (smoothed periodogram)
##  6.  Wavelet analysis    (Morlet CWT)
##  7.  Hidden Markov Models
##        7a. PP  - 2-state levels HMM    [primary]
##        7b. PP  - 2-state residual HMM  [robustness - degenerate]
##        7c. PP  - 3-state levels HMM    [extension]
##        7d. Auto - 2-state GARCH-residual HMM [primary]
##        7e. Auto - 3-state HMM          [extension - fails]
##  8.  HMM Model election - BIC Comparison
##  9.  Expected regime durations and stationary distribution
##  10. Out-of-sample rolling forecast evaluation (2011–2022)
##  11. Regime-conditioned forecasts (2023–2025)
## 
##  ────────────────────────────────────────────────────────────────────────

###############################################################################
## 0.  SETUP AND DATA PREPARATION
###############################################################################

# ── 0.1 Load required packages ────────────────────────────────────────────────
# Install any missing packages with:  install.packages("package_name")

library(urca)        # ADF unit-root test (ur.df)
library(tseries)     # KPSS stationarity test
library(forecast)    # ARIMA modelling (Arima, accuracy, checkresiduals)
library(rugarch)     # GARCH modelling (ugarchspec, ugarchfit, ugarchforecast)
library(depmixS4)    # Hidden Markov Models (depmix, fit, posterior, getpars)
library(WaveletComp) # Continuous wavelet transform (analyze.wavelet, wt.image)
library(dplyr)       # Data manipulation (filter, mutate, etc.)
library(tidyr)       # Data reshaping (pivot_longer)
library(ggplot2)     # Grammar of graphics for plotting and visualisation
library(patchwork)   # Combining multiple ggplot objects into composite figures
library(readxl)      # Importing Excel (.xlsx) data files

# Import Excel 
# The raw loss ratio data are stored in an Excel file.
# To replicate the analysis, place the file in the working directory

pc_raw <- read_excel("P&C.Canada.xlsx")

# ── 0.2 Extract Personal Property (PP) loss ratio series ─────────────────────
# Column 10 = PP net incurred claims / net written premiums (loss ratio).
# Column  1 = Year.  We keep only rows where Year is a valid number.

rows_keep <- !is.na(as.numeric(P_C_Canada[[1]]))

df_pp <- data.frame(
  Year  = as.numeric(pc_raw[[1]])[rows_keep],   # 1985–2022
  PP_LR = as.numeric(pc_raw[[10]])[rows_keep]   # PP loss ratio
)

# Convert to annual time-series object (start = 1985, frequency = 1)
pp_ts <- ts(df_pp$PP_LR, start = df_pp$Year[1], frequency = 1)


# ── 0.3 Extract Auto loss ratio series ───────────────────────────────────────
# Column 9 = Auto loss ratio

df_auto <- data.frame(
  Year    = as.numeric(pc_raw[[1]])[rows_keep],
  Auto_LR = as.numeric(pc_raw[[9]])[rows_keep]
)

auto_ts <- ts(df_auto$Auto_LR, start = df_auto$Year[1], frequency = 1)


# ── 0.4 First differences ─────────────────────────────────────────────────────
# ΔYₜ = Yₜ − Yₜ₋₁   (removes one unit root if present)

pp_ts_d1   <- diff(pp_ts)
auto_ts_d1 <- diff(auto_ts)


# ── 0.5 Quick sanity checks ───────────────────────────────────────────────────
# Check lengths
cat("PP observations  :", length(pp_ts),   "\n")   # expect 38
cat("Auto observations:", length(auto_ts), "\n")   # expect 38

# Common plotting settings

par(mfrow = c(1, 1),
    mar = c(4.5, 4.5, 3, 1))       

# 1) PP loss ratio
plot(pp_ts,  main = "PP loss ratio, 1985-2022",
     type = "l",
     lwd  = 2,
     col  = "#1B4F72",
     xlab = "Year",
     ylab = "Loss Ratio",
     ylim = range(pp_ts, auto_ts))
grid(col = "grey85")

# 2) Auto loss ratio
plot(auto_ts, main = "Auto loss ratio, 1985-2022",
     type = "l",
     lwd  = 2,
     col  = "#B03A2E",
     xlab = "Year",
     ylab = "Loss Ratio",
     ylim = range(pp_ts, auto_ts))
grid(col = "grey85")


###############################################################################
## 1.  DESCRIPTIVE STATISTICS
###############################################################################

cat("\n── PP descriptive statistics ──────────────────────────────────────────\n")
print(summary(pp_ts))
cat("SD:", round(sd(pp_ts), 4), "\n")

cat("\n── Auto descriptive statistics ─────────────────────────────────────────\n")
print(summary(auto_ts))
cat("SD:", round(sd(auto_ts), 4), "\n")


###############################################################################
## 2.  STATIONARITY TESTS
###############################################################################
##
##  ADF  (Augmented Dickey-Fuller): H₀ = unit root (non-stationary)
##  KPSS (Kwiatkowski-Phillips-Schmidt-Shin): H₀ = level-stationary
##

# ── 2.1 PP levels (trend specification, 3 lags) ────────────────────────────
cat("\n── ADF: PP levels ──────────────────────────────────────────────────────\n")
adf_pp_lvl <- ur.df(pp_ts, type = "trend", lags = 3)
summary(adf_pp_lvl)
# Test statistic: -2.621  |  5% c.v.: -3.50  -  fail to reject H₀  (non-stationary)

cat("\n── KPSS: PP levels ─────────────────────────────────────────────────────\n")
kpss_pp_lvl <- kpss.test(pp_ts, null = "Level")
print(kpss_pp_lvl)
# Statistic: 0.392  |  p = 0.081 - fail to reject H₀  (stationary)
# CONCLUSION: Conflicting - PP is NEAR-STATIONARY. 

# ── 2.2 PP first differences (drift specification, 3 lags) ─────────────────
cat("\n── ADF: PP first differences ───────────────────────────────────────────\n")
adf_pp_d1 <- ur.df(pp_ts_d1, type = "drift", lags = 3)
summary(adf_pp_d1)
# Test statistic: -4.070  |  1% c.v.: -3.58  -  REJECT H₀  (stationary)

cat("\n── KPSS: PP first differences ──────────────────────────────────────────\n")
kpss_pp_d1 <- kpss.test(pp_ts_d1, null = "Level")
print(kpss_pp_d1)
# Statistic: 0.070  |  p > 0.10  -  fail to reject H₀  (stationary)
# CONCLUSION: ΔPP is unambiguously stationary - conservative I(1) adopted

# ── 2.3 Auto levels ─────────────────────────────────────────────────────────
cat("\n── ADF: Auto levels ────────────────────────────────────────────────────\n")
adf_auto_lvl <- ur.df(auto_ts, type = "trend", lags = 3)
summary(adf_auto_lvl)
# Expected: fail to reject H₀

cat("\n── KPSS: Auto levels ───────────────────────────────────────────────────\n")
kpss_auto_lvl <- kpss.test(auto_ts, null = "Level")
print(kpss_auto_lvl)
# Expected: reject H₀  -  Auto is clearly non-stationary in levels

# ── 2.4 Auto first differences ─────────────────────────────────────────────
cat("\n── ADF: Auto first differences ─────────────────────────────────────────\n")
adf_auto_d1 <- ur.df(auto_ts_d1, type = "drift", lags = 3)
summary(adf_auto_d1)

cat("\n── KPSS: Auto first differences ────────────────────────────────────────\n")
kpss_auto_d1 <- kpss.test(auto_ts_d1, null = "Level")
print(kpss_auto_d1)
# CONCLUSION: ΔAuto is stationary - Auto modelled as I(1)


###############################################################################
## 3.  ARIMA MEAN MODELS
###############################################################################
##
##  Model selection criteria:
##    • AIC  (Akaike Information Criterion) — lower is better
##    • Ljung-Box test p-value — want p > 0.05 (no residual autocorrelation)
##    • In-sample RMSE — lower is better
##
##  Three candidate orders are compared for each series.
##

# ── 3.1 PP levels ──────────────────────────────────────────────────────────
m1_pp_lvl <- Arima(pp_ts, order = c(1, 0, 0))   # AR(1)
m2_pp_lvl <- Arima(pp_ts, order = c(0, 0, 1))   # MA(1)
m3_pp_lvl <- Arima(pp_ts, order = c(1, 0, 1))   # ARMA(1,1)

cat("\n── PP levels — AIC comparison ──────────────────────────────────────────\n")
cat("AR(1)   AIC:", round(AIC(m1_pp_lvl), 3), "\n")
cat("MA(1)   AIC:", round(AIC(m2_pp_lvl), 3), "\n")
cat("ARMA    AIC:", round(AIC(m3_pp_lvl), 3), "\n")
# AR(1) has lowest AIC - selected as primary levels benchmark

cat("\n── PP levels — ARIMA coefficients ─────────────────────────────────────\n")
print(coef(m1_pp_lvl))           # AR(1) coefficient and intercept

cat("\n── PP levels — RMSE ────────────────────────────────────────────────────\n")
rmse_pp_lvl <- c(
  AR1   = accuracy(m1_pp_lvl)["Training set", "RMSE"],
  MA1   = accuracy(m2_pp_lvl)["Training set", "RMSE"],
  ARMA  = accuracy(m3_pp_lvl)["Training set", "RMSE"]
)
print(round(rmse_pp_lvl, 4))

# Selected model: AR(1) on levels
best_pp_lvl <- m1_pp_lvl

# Residuals From selected model
res_pp_lvl <- residuals(best_pp_lvl)

Box.test(res_pp_lvl,
         lag  = 8,          
         type = "Ljung-Box")

res_pp <- residuals(best_pp_lvl)
df_res <- data.frame(res = as.numeric(res_pp))

# common theme
my_theme <- theme_minimal(base_size = 11)

# Top: residual time plot
p_top <- autoplot(res_pp) +
  labs(title = "Residuals from AR(1) for PP",
       x = "Year", y = "Residual") +
  my_theme

# Bottom-left: ACF
p_acf <- ggAcf(res_pp) +
  labs(title = "ACF of residuals", x = "Lag", y = "ACF") +
  my_theme

# Bottom-right: histogram with normal curve
p_hist <- ggplot(df_res, aes(x = res)) +
  geom_histogram(aes(y = ..density..),
                 bins = 10, colour = "grey40", fill = "skyblue3") +
  stat_function(fun = dnorm,
                args = list(mean = mean(res_pp, na.rm = TRUE),
                            sd   = sd(res_pp,   na.rm = TRUE)),
                colour = "red3", linewidth = 0.9) +
  labs(title = "Residual distribution", x = "Residual", y = "Density") +
  my_theme

(p_top) / (p_acf + p_hist + plot_layout(widths = c(1, 1)))


# ACF / PACF correlograms
par(mfrow = c(2, 1))
acf(pp_ts,  main = "ACF - PP loss ratio (levels)")
pacf(pp_ts, main = "PACF - PP loss ratio (levels)")
par(mfrow = c(1, 1))

# ── 3.2 PP first differences ────────────────────────────────────────────────
m1_pp_d1 <- Arima(pp_ts_d1, order = c(1, 0, 0))   # AR(1) on ΔPP
m2_pp_d1 <- Arima(pp_ts_d1, order = c(0, 0, 1))   # MA(1) on ΔPP
m3_pp_d1 <- Arima(pp_ts_d1, order = c(1, 0, 1))   # ARMA(1,1) on ΔPP

cat("\n── PP first diff — AIC comparison ──────────────────────────────────────\n")
cat("AR(1)   AIC:", round(AIC(m1_pp_d1), 3), "\n")
cat("MA(1)   AIC:", round(AIC(m2_pp_d1), 3), "\n")
cat("ARMA    AIC:", round(AIC(m3_pp_d1), 3), "\n")
# MA(1) on differences has lowest AIC - selected

cat("\n── PP first diff — MA(1) coefficients ──────────────────────────────────\n")
print(coef(m2_pp_d1))
# MA(1) coefficient 

cat("\n── PP first diff — RMSE ─────────────────────────────────────────────────\n")
rmse_pp_d1 <- c(
  AR1   = accuracy(m1_pp_d1)["Training set", "RMSE"],
  MA1   = accuracy(m2_pp_d1)["Training set", "RMSE"],
  ARMA  = accuracy(m3_pp_d1)["Training set", "RMSE"]
)
print(round(rmse_pp_d1, 4))

best_pp_d1 <- m2_pp_d1

res_pp_d1 <- residuals(best_pp_d1)

Box.test(res_pp_d1,
         lag  = 7,          
         type = "Ljung-Box")

# residuals from PP I(1) model
res_pp_d1  <- residuals(best_pp_d1)  
df_res_d1  <- data.frame(res = as.numeric(res_pp_d1))

my_theme <- theme_minimal(base_size = 11)

# Top: residual time plot
p_top_d1 <- autoplot(res_pp_d1) +
  labs(title = "Residuals from MA(1) on ΔPP",
       x = "Year", y = "Residual") +
  my_theme

# Bottom-left: ACF
p_acf_d1 <- ggAcf(res_pp_d1) +
  labs(title = "ACF of residuals", x = "Lag", y = "ACF") +
  my_theme

# Bottom-right: histogram
p_hist_d1 <- ggplot(df_res_d1, aes(x = res)) +
  geom_histogram(aes(y = ..density..),
                 bins = 10, colour = "grey40", fill = "skyblue3") +
  stat_function(fun = dnorm,
                args = list(mean = mean(res_pp_d1, na.rm = TRUE),
                            sd   = sd(res_pp_d1,   na.rm = TRUE)),
                colour = "red3", linewidth = 0.9) +
  labs(title = "Residual distribution", x = "Residual", y = "Density") +
  my_theme

# Same layout
p_top_d1 / (p_acf_d1 + p_hist_d1 + plot_layout(widths = c(1, 1)))

par(mfrow = c(2, 1))
acf(pp_ts_d1,  main = "ACF - ΔPP loss ratio")
pacf(pp_ts_d1, main = "PACF - ΔPP loss ratio")
par(mfrow = c(1, 1))

# ── 3.3 Auto first differences (I(1) confirmed) ────────────────────────────
m1_auto <- Arima(auto_ts, order = c(1, 1, 0))   # ARIMA(1,1,0)
m2_auto <- Arima(auto_ts, order = c(0, 1, 1))   # ARIMA(0,1,1)
m3_auto <- Arima(auto_ts, order = c(1, 1, 1))   # ARIMA(1,1,1)

cat("\n── Auto — AIC comparison ───────────────────────────────────────────────\n")
cat("ARIMA(1,1,0) AIC:", round(AIC(m1_auto), 3), "\n")
cat("ARIMA(0,1,1) AIC:", round(AIC(m2_auto), 3), "\n")
cat("ARIMA(1,1,1) AIC:", round(AIC(m3_auto), 3), "\n")
# ARIMA(0,1,1) has lowest AIC - selected

cat("\n── Auto — ARIMA(0,1,1) coefficients ───────────────────────────────────\n")
print(coef(m2_auto))

cat("\n── Auto — RMSE ─────────────────────────────────────────────────────────\n")
rmse_auto <- c(
  ARIMA110 = accuracy(m1_auto)["Training set", "RMSE"],
  ARIMA011 = accuracy(m2_auto)["Training set", "RMSE"],
  ARIMA111 = accuracy(m3_auto)["Training set", "RMSE"]
)
print(round(rmse_auto, 4))

best_auto <- m2_auto

res_auto <- residuals(best_auto)

Box.test(res_auto,
         lag   = 7,         
         type  = "Ljung-Box")

# residuals from your chosen Auto I(1) model, e.g. ARIMA(0,1,1)
res_auto_d1 <- residuals(best_auto)   
df_auto_d1  <- data.frame(res = as.numeric(res_auto_d1))

my_theme <- theme_minimal(base_size = 11)

# Top: residual time plot
p_top_auto <- autoplot(res_auto_d1) +
  labs(title = "Residuals from ARIMA(0,1,1) for ΔAuto",
       x = "Year", y = "Residual") +
  my_theme

# Bottom-left: ACF
p_acf_auto <- ggAcf(res_auto_d1) +
  labs(title = "ACF of residuals", x = "Lag", y = "ACF") +
  my_theme

# Bottom-right: histogram
p_hist_auto <- ggplot(df_auto_d1, aes(x = res)) +
  geom_histogram(aes(y = ..density..),
                 bins = 10, colour = "grey40", fill = "skyblue3") +
  stat_function(fun = dnorm,
                args = list(mean = mean(res_auto_d1, na.rm = TRUE),
                            sd   = sd(res_auto_d1,   na.rm = TRUE)),
                colour = "red3", linewidth = 0.9) +
  labs(title = "Residual distribution", x = "Residual", y = "Density") +
  my_theme

# Same layout 
p_top_auto / (p_acf_auto + p_hist_auto + plot_layout(widths = c(1, 1)))

par(mfrow = c(2, 1))
acf(diff(auto_ts),  main = "ACF - ΔAuto loss ratio")
pacf(diff(auto_ts), main = "PACF - ΔAuto loss ratio")
par(mfrow = c(1, 1))

# ── 3.4 ARIMA summary table ───────────────────────────────────────────────────
arima_tbl <- data.frame(
  Series   = c("PP levels", "", "", "PP diff.", "", "", "Auto I(1)", "", ""),
  Model    = c("AR(1)", "MA(1)", "ARMA(1,1)",
               "AR(1)", "MA(1)", "ARMA(1,1)",
               "ARIMA(1,1,0)", "ARIMA(0,1,1)", "ARIMA(1,1,1)"),
  AIC      = round(c(AIC(m1_pp_lvl), AIC(m2_pp_lvl), AIC(m3_pp_lvl),
                     AIC(m1_pp_d1),  AIC(m2_pp_d1),  AIC(m3_pp_d1),
                     AIC(m1_auto),   AIC(m2_auto),   AIC(m3_auto)), 2),
  RMSE     = round(c(rmse_pp_lvl, rmse_pp_d1, rmse_auto), 4),
  Selected = c("✓", "", "", "", "✓", "", "", "✓", "")
)


###############################################################################
## 4.  GARCH VOLATILITY MODELS
###############################################################################
##
##  GARCH(1,1): σₜ² = ω + α·εₜ₋₁² + β·σₜ₋₁²
##    α = shock sensitivity (ARCH effect)
##    β = persistence (GARCH effect)
##    α + β → 1 implies near-integrated volatility (IGARCH-like)
##
##  For PP and Auto, α is near zero and β close to 1, indicating very persistent but weakly shock‑driven volatility
##  GARCH offers only limited improvement over a constant‑variance assumption.
##

# ── 4.1 PP levels — AR(1)–GARCH(1,1) ────────────────────────────────────────

spec_garch_pp_lvl <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_garch_pp_lvl <- ugarchfit(spec = spec_garch_pp_lvl, data = pp_ts)

cat("\n── PP levels — AR(1)–GARCH(1,1) results ───────────────────────────────\n")
show(fit_garch_pp_lvl)

# Extract key parameters
alpha_pp_lvl <- coef(fit_garch_pp_lvl)["alpha1"]
beta_pp_lvl  <- coef(fit_garch_pp_lvl)["beta1"]
cat("α:", round(alpha_pp_lvl, 4),
    "  β:", round(beta_pp_lvl, 4),
    "  α+β:", round(alpha_pp_lvl + beta_pp_lvl, 4), "\n")
# Expected: α ≈ 0, β ≈ 1 - no ARCH effect, near-constant variance

# RMSE 
pp_garch_fitted_lvl <- as.numeric(fitted(fit_garch_pp_lvl))
rmse_garch_pp_lvl   <- sqrt(mean((as.numeric(pp_ts) - pp_garch_fitted_lvl)^2))
cat("PP levels GARCH RMSE:", round(rmse_garch_pp_lvl, 4), "\n")

# Standardised residual diagnostics
resid_pp_lvl <- residuals(fit_garch_pp_lvl, standardize = TRUE)
par(mfrow = c(2, 1))
acf(resid_pp_lvl,    main = "ACF - PP levels GARCH standardised residuals")
acf(resid_pp_lvl^2, main = "ACF - PP levels GARCH squared standardised residuals")
par(mfrow = c(1, 1))

# ── 4.2 PP first differences — MA(1)–GARCH(1,1) ──────────────────────────────
spec_garch_pp_d1 <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 1), include.mean = TRUE),
  distribution.model = "norm"
)

fit_garch_pp_d1 <- ugarchfit(spec = spec_garch_pp_d1, data = pp_ts_d1)

cat("\n── PP first diff — MA(1)–GARCH(1,1) results ───────────────────────────\n")
show(fit_garch_pp_d1)

# Extract key parameters
alpha_pp_d1 <- coef(fit_garch_pp_d1)["alpha1"]
beta_pp_d1  <- coef(fit_garch_pp_d1)["beta1"]
cat("α:", round(alpha_pp_d1, 4),
    "  β:", round(beta_pp_d1, 4),
    "  α+β:", round(alpha_pp_d1 + beta_pp_d1, 4), "\n")
# Similar pattern to levels - GARCH adds no value for PP
# Note: levels spec is retained as primary; differenced results are robustness only.

# RMSE
rmse <- function(actual, fitted) {
  sqrt(mean((actual - fitted)^2, na.rm = TRUE))
}


# Actual and fitted values (difference scale)
pp_d1_actual <- as.numeric(pp_ts_d1)
pp_d1_fitted <- as.numeric(fitted(fit_garch_pp_d1))

rmse_garch_pp_d1 <- rmse(pp_d1_actual, pp_d1_fitted)
cat("PP Δ loss ratio GARCH RMSE:", round(rmse_garch_pp_d1, 4), "\n")

rmse_garch_pp_d1 <- rmse(pp_d1_actual, pp_d1_fitted)
rmse_garch_pp_d1

# Standardised residual diagnostics
resid_pp_d1 <- residuals(fit_garch_pp_d1, standardize = TRUE)
par(mfrow = c(2, 1))
acf(resid_pp_d1,    main = "ACF - ΔPP GARCH standardised residuals")
acf(resid_pp_d1^2, main = "ACF — ΔPP GARCH squared standardised residuals")
par(mfrow = c(1, 1))


# ── 4.3 Auto — ARIMA(0,1,1)–GARCH(1,1) ───────────────────────────────────────
spec_garch_auto <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 1), include.mean = TRUE),
  distribution.model = "norm"
)

fit_garch_auto <- ugarchfit(spec = spec_garch_auto, data = auto_ts_d1)

cat("\n── Auto — ARIMA(0,1,1)–GARCH(1,1) results ─────────────────────────────\n")
show(fit_garch_auto)

# Extract key parameters
alpha_auto <- coef(fit_garch_auto)["alpha1"]
beta_auto  <- coef(fit_garch_auto)["beta1"]
cat("α:", round(alpha_auto, 4),
    "  β:", round(beta_auto, 4),
    "  α+β:", round(alpha_auto + beta_auto, 4), "\n")
#Expected: α ≈ 0, β ≈ 1 - no ARCH effect, near-constant variance

# RMSE on differenced scale
auto_garch_fitted <- as.numeric(fitted(fit_garch_auto))
rmse_garch_auto   <- sqrt(mean((as.numeric(auto_ts_d1) - auto_garch_fitted)^2))
cat("Auto GARCH RMSE (diff scale):", round(rmse_garch_auto, 4), "\n")


# Standardised residual diagnostics
resid_auto <- residuals(fit_garch_auto, standardize = TRUE)
par(mfrow = c(2, 1))
acf(resid_auto,    main = "ACF - Auto GARCH standardised residuals")
acf(resid_auto^2, main = "ACF - Auto GARCH squared standardised residuals")
par(mfrow = c(1, 1))

# ── 4.4 GARCH summary table ───────────────────────────────────────────────────
garch_tbl <- data.frame(
  Specification   = c("AR(1)–GARCH(1,1) - PP levels",
                      "MA(1)–GARCH(1,1) - PP diff.",
                      "ARIMA(0,1,1)–GARCH(1,1) — Auto"),
  alpha           = round(c(alpha_pp_lvl, alpha_pp_d1, alpha_auto), 4),
  beta            = round(c(beta_pp_lvl,  beta_pp_d1,  beta_auto),  4),
  alpha_plus_beta = round(c(alpha_pp_lvl + beta_pp_lvl,
                            alpha_pp_d1  + beta_pp_d1,
                            alpha_auto   + beta_auto),  4)
)
cat("\n── GARCH parameter summary ─────────────────────────────────────────────\n")
print(garch_tbl, row.names = FALSE)


###############################################################################
## 5.  SPECTRAL ANALYSIS
###############################################################################
##
##  Smoothed power spectrum: Daniell kernel with spans c(3,3)
##  Dominant frequency f* → implied cycle period = 1/f* years
##
##  PP:   dominant f* ≈ 0.025  -  period ≈ 40 years (no stable 6–8 yr cycle)
##  Auto: dominant f* ≈ 0.100  -  period ≈ 10 years (decadal cycle)
##

# ── 5.1 PP levels ──────────────────────────────────────────────────────────
spec_pp_lvl <- spectrum(pp_ts, spans = c(3, 3), log = "no" ,
                        main = "Smoothed spectrum - PP loss ratio")
dom_f_pp_lvl <- spec_pp_lvl$freq[which.max(spec_pp_lvl$spec)]
cat("\n── PP levels spectrum ──────────────────────────────────────────────────\n")
cat("Dominant frequency:", dom_f_pp_lvl, "\n")
cat("Implied period (yrs):", round(1 / dom_f_pp_lvl, 1), "\n")

# ── 5.2 Auto levels ─────────────────────────────────────────────────────────
spec_auto <- spectrum(auto_ts, spans = c(3, 3), log = "no" ,
                      main = "Smoothed spectrum - Auto loss ratio")
dom_f_auto <- spec_auto$freq[which.max(spec_auto$spec)]
cat("\n── Auto levels spectrum ────────────────────────────────────────────────\n")
cat("Dominant frequency:", dom_f_auto, "\n")
cat("Implied period (yrs):", round(1 / dom_f_auto, 1), "\n")


###############################################################################
## 6.  WAVELET ANALYSIS
###############################################################################
##
##  Continuous Morlet wavelet transform — time-frequency decomposition.
##  Unlike spectral analysis, wavelets reveal HOW cyclical content
##  changes over time (non-stationary cycles).
##
##  Parameters:
##    dt = 1 (annual data)
##    dj = 1/20 (scale resolution)
##    lowerPeriod = 2, upperPeriod = 16 (plausible cycle range in years)
##    n.sim = 200 (Monte Carlo simulations for significance testing)
##

# ── 6.1 PP levels ──────────────────────────────────────────────────────────
pp_df_lvl <- data.frame(
  time = as.numeric(time(pp_ts)),
  pp   = as.numeric(pp_ts)
)


wt_pp_lvl <- analyze.wavelet(
  pp_df_lvl, "pp",
  loess.span  = 0,       # no pre-detrending
  dt          = 1,
  dj          = 1 / 20,
  lowerPeriod = 2,
  upperPeriod = 16,
  make.pval   = TRUE,
  n.sim       = 200
)

wt.image(wt_pp_lvl,
         main      = "Wavelet power spectrum - PP loss ratio",
         legend.params = list(lab = "Power"),
         periodlab = "Period (years)")


# ── 6.2 Auto levels ─────────────────────────────────────────────────────────
auto_df <- data.frame(
  time = as.numeric(time(auto_ts)),
  auto = as.numeric(auto_ts)
)

wt_auto <- analyze.wavelet(
  auto_df, "auto",
  loess.span  = 0,
  dt          = 1,
  dj          = 1 / 20,
  lowerPeriod = 2,
  upperPeriod = 16,
  make.pval   = TRUE,
  n.sim       = 200
)

wt.image(wt_auto,
         main      = "Wavelet power spectrum - Auto loss ratio",
         legend.params = list(lab = "Power"),
         periodlab = "Period (years)")


###############################################################################
## 7.  HIDDEN MARKOV MODELS
###############################################################################
##
##  Model: Yₜ | Sₜ = k ~ N(μₖ, σₖ²)
##  Sₜ evolves as a first-order Markov chain with transition matrix P.
##  Estimated by EM (Baum-Welch) algorithm via depmixS4.
##
##  set.seed(123) ensures reproducibility of EM starting values.
##

# ── 7a. PP - 2-state HMM on LEVELS [PRIMARY SPECIFICATION] ───────────────────
#
#  PP is near-stationary, and the level series carries the
#  economically meaningful regime structure (hard/soft market phases).
#  The decoded state path maps directly onto underwriting conditions.


y_pp <- as.numeric(pp_ts)
set.seed(123)
mod_hmm_pp_lvl <- depmix(
  y_pp ~ 1,
  data    = data.frame(y_pp = y_pp),
  nstates = 2,
  family  = gaussian()
)

fit_hmm_pp_lvl <- fit(mod_hmm_pp_lvl)

cat("\n── PP levels 2-state HMM - parameter summary ───────────────────────────\n")
summary(fit_hmm_pp_lvl)
# State 1: higher mean LR - Hard market
# State 2: lower  mean LR - Soft market

# Extract transition matrix and state parameters
pars_pp_lvl    <- getpars(fit_hmm_pp_lvl)
tmat_pp_lvl    <- matrix(pars_pp_lvl[3:6], nrow = 2, byrow = TRUE)
means_pp_lvl   <- pars_pp_lvl[c(7, 9)]
sds_pp_lvl     <- pars_pp_lvl[c(8, 10)]

cat("\nTransition matrix:\n"); print(round(tmat_pp_lvl, 3))
cat("\nState means:", round(means_pp_lvl, 3), "\n")
cat("State SDs:  ", round(sds_pp_lvl,   3), "\n")

# Decoded state sequence (Viterbi)
post_pp_lvl <- posterior(fit_hmm_pp_lvl)

# Plot: loss ratio coloured by decoded HMM state
plot(time(pp_ts), y_pp, type = "l", col = "grey50", lwd = 1.5,
     main = "PP loss ratio - HMM regime classification (2-state)",
     xlab = "Year", ylab = "Loss ratio")
points(time(pp_ts), y_pp,
       col = c("tomato", "steelblue")[post_pp_lvl$state], pch = 19, cex = 1.1)
legend("topright",
       legend = c("State 1 - Hard market", "State 2 - Soft market"),
       col    = c("tomato", "steelblue"), pch = 19, bty = "n")

set.seed(123)
mod_hmm_pp_3st <- depmix(
  y_pp ~ 1,
  data    = data.frame(y_pp = y_pp),
  nstates = 3,
  family  = gaussian()
)

fit_hmm_pp_3st <- fit(mod_hmm_pp_3st, verbose = TRUE)
loglik_pp_3state <- as.numeric(logLik(fit_hmm_pp_3st))


# ── 7b. PP — 2-state HMM on AR(1) RESIDUALS [ROBUSTNESS — DEGENERATE] ────────
#
#  This asks whether DEVIATIONS from the AR(1) trend show regime structure.
#  Result: State 2 has P(stay) ≈ 0.002 - degenerate (duration ≈ 1 year).
#  Conclusion: insufficient regime structure in residuals - LEVELS HMM preferred.

# AR(1) residuals from the linear mean model
eps_pp_lvl <- as.numeric(residuals(best_pp_lvl))  

set.seed(123)
mod_hmm_pp_res <- depmix(
  eps_pp_lvl ~ 1,
  data    = data.frame(eps_pp_lvl = eps_pp_lvl),
  nstates = 2,
  family  = gaussian()
)

fit_hmm_pp_res <- fit(mod_hmm_pp_res)

cat("\n── PP residuals 2-state HMM — summary (ROBUSTNESS CHECK) ─────────────\n")
summary(fit_hmm_pp_res)
# Expected: degenerate - P(stay) ≈ 0 for one state

# Count regime entries to confirm instability
states_pp_res <- posterior(fit_hmm_pp_res)$state
cat("Number of entries into State 1:", sum(diff(c(0, states_pp_res)) == 1), "\n")
# High count across only 38 obs - not a genuine persistent regime

# ── 7c. PP — 3-state HMM on LEVELS [EXTENSION] ────────────────────────────────
#
#  Tests whether a transitional state exists between hard and soft markets.
#  BIC penalty = 6 × ln(38) ≈ 21.8 log-likelihood units (12 vs 6 parameters).
#  Included for qualitative insight; 2-state is primary specification.

set.seed(123)
mod_hmm_pp_3st <- depmix(
  y_pp ~ 1,
  data    = data.frame(y_pp = y_pp),
  nstates = 3,
  family  = gaussian()
)

fit_hmm_pp_3st <- fit(mod_hmm_pp_3st, verbose = TRUE)

cat("\n── PP levels 3-state HMM — summary ────────────────────────────────────\n")
summary(fit_hmm_pp_3st)

pars_pp_3st  <- getpars(fit_hmm_pp_3st)
tmat_pp_3st  <- matrix(pars_pp_3st[3:11], nrow = 3, byrow = TRUE)
means_pp_3st <- pars_pp_3st[c(7, 9, 11)]
sds_pp_3st   <- pars_pp_3st[c(8, 10, 12)]

cat("\nTransition matrix (3-state PP):\n"); print(round(tmat_pp_3st, 3))
cat("State means:", round(means_pp_3st, 3), "\n")

post_pp_3st <- posterior(fit_hmm_pp_3st)
loglik_pp_3state <- as.numeric(logLik(fit_hmm_pp_3st))

# Plot 3-state decoded path
plot(time(pp_ts), y_pp, type = "l", col = "grey50", lwd = 1.5,
     main = "PP loss ratio - HMM regime classification (3-state, levels)",
     xlab = "Year", ylab = "Loss ratio")
points(time(pp_ts), y_pp,
       col = c("tomato", "goldenrod", "steelblue")[post_pp_3st$state],
       pch = 19, cex = 1.1)
legend("topright",
       legend = c("State 1", "State 2 - Transitional", "State 3"),
       col    = c("tomato", "goldenrod", "steelblue"), pch = 19, bty = "n")

# 3-state HMM suggests some additional nuance in PP loss-ratio dynamics (three nearby mean levels), 
# but patterns are not strong enough to replace the simpler 2-state model.

# ── 7d. Auto — 2-state HMM on GARCH RESIDUALS [PRIMARY] ──────────────────────
#
#  Why residuals? Auto is I(1) — differencing is required for stationarity.
#  GARCH removes volatility clustering too. The HMM then identifies regimes
#  in what ARIMA-GARCH cannot explain: persistent above/below average years.
#
#  State 1 (benign):  near-zero residuals, high SD, high persistence
#  State 2 (stressed): strongly negative residuals, low SD, brief episodes

y_auto <- as.numeric(resid_auto)   # standardised ARIMA-GARCH residuals

set.seed(123)
mod_hmm_auto <- depmix(
  y_auto ~ 1,
  data    = data.frame(y_auto = y_auto),
  nstates = 2,
  family  = gaussian()
)

fit_hmm_auto <- fit(mod_hmm_auto)

cat("\n── Auto 2-state HMM (GARCH residuals) — summary ───────────────────────\n")
summary(fit_hmm_auto)

pars_auto_hmm  <- getpars(fit_hmm_auto)
tmat_auto_hmm  <- matrix(pars_auto_hmm[3:6], nrow = 2, byrow = TRUE)
means_auto_hmm <- pars_auto_hmm[c(7, 9)]
sds_auto_hmm   <- pars_auto_hmm[c(8, 10)]

cat("\nTransition matrix (Auto):\n"); print(round(tmat_auto_hmm, 3))
cat("State means:", round(means_auto_hmm, 3), "\n")

post_auto <- posterior(fit_hmm_auto)
loglik_auto_2state <- as.numeric(logLik(fit_hmm_auto))

# Plot decoded Auto residual regimes
plot(y_auto, type = "l", col = "grey50", lwd = 1.5,
     main = "Auto GARCH residuals - HMM regime classification (2-state)",
     xlab = "Year index", ylab = "Standardised residual")
points(seq_along(y_auto), y_auto,
       col = c("tomato", "steelblue")[post_auto$state], pch = 19, cex = 1.1)
abline(h = 0, lty = 2, col = "grey70")
legend("topright",
       legend = c("State 1 - Benign", "State 2 - Stressed"),
       col    = c("tomato", "steelblue"), pch = 19, bty = "n")

# ── 7e. Auto — 3-state HMM [EXTENSION — EXPECTED TO FAIL] ────────────────────
#
#  With n = 38 and p = 12 parameters, BIC penalty ≈ 21.8 units.
#  Diagnostic signs of failure:
#    • P(stay) = 0 for one state (junk-drawer state)
#    • P(stay) = 1 for another (absorbing state — never exits)
#    • Transition row sums ≠ 1 (numerical instability in EM)
#
#  This failure is reported as a methodological finding, not hidden.

y_auto_lvl <- as.numeric(auto_ts)   # NOTE: 3-state on levels for comparability

set.seed(123)
mod_hmm_auto_3st <- depmix(
  y_auto_lvl ~ 1,
  data    = data.frame(y_auto_lvl = y_auto_lvl),
  nstates = 3,
  family  = gaussian()
)

fit_hmm_auto_3st <- fit(mod_hmm_auto_3st, verbose = TRUE)

cat("\n── Auto 3-state HMM — summary (EXPECTED DEGENERACY) ───────────────────\n")
summary(fit_hmm_auto_3st)

pars_auto_3st <- getpars(fit_hmm_auto_3st)
tmat_auto_3st <- matrix(pars_auto_3st[3:11], nrow = 3, byrow = TRUE)
cat("\nTransition matrix (3-state Auto):\n")
print(round(tmat_auto_3st, 3))
cat("\nRow sums (should all equal 1.0):\n")
print(round(rowSums(tmat_auto_3st), 3))
# Row sums ≠ 1 confirms EM convergence failure
loglik_auto_3state <- as.numeric(logLik(fit_hmm_auto_3st))


# 3-state Auto HMM is numerically degenerate (row sums ≠ 1, one absorbing and one “junk” state),
# so it’s useful only as evidence that adding a third regime overfits this 38-year sample.



###############################################################################
## 8. HMM MODEL SELECTION — BIC COMPARISON
###############################################################################
##
## Bayesian Information Criterion (BIC):
##   BIC = -2 × logLik + p × ln(n)
##
## where:
##   p = K^2 + K
##       = K(K−1) transition probabilities
##       + K state-dependent means
##       + K state-dependent variances
##   n = number of time-series observations
##
## Lower BIC indicates preferred model (penalises over-parameterisation).
###############################################################################

library(depmixS4)

# Number of observations (annual data)
n <- length(pp_ts)  # expected to be 38

# ---------------------------------------------------------------------------
# Helper function to estimate a Gaussian HMM and extract log-likelihood
# ---------------------------------------------------------------------------
fit_hmm_loglik <- function(ts_data, k) {
  
  hmm_spec <- depmix(
    response = ts_data ~ 1,
    data     = data.frame(ts_data = ts_data),
    nstates  = k,
    family   = gaussian()
  )
  
  hmm_fit <- fit(hmm_spec, verbose = FALSE)
  
  return(as.numeric(logLik(hmm_fit)))
}

# ---------------------------------------------------------------------------
# Fit HMMs and extract log-likelihoods
# ---------------------------------------------------------------------------

# Property (PP) — levels
loglik_pp_2state   <- fit_hmm_loglik(pp_ts, k = 2)
loglik_pp_3state   <- fit_hmm_loglik(pp_ts, k = 3)

# Auto — residuals or levels (as used elsewhere in analysis)
loglik_auto_2state <- fit_hmm_loglik(auto_ts, k = 2)
loglik_auto_3state <- fit_hmm_loglik(auto_ts, k = 3)

# ---------------------------------------------------------------------------
# BIC calculation function
# ---------------------------------------------------------------------------
bic_fn <- function(loglik, k, n) {
  p   <- k^2 + k
  bic <- -2 * loglik + p * log(n)
  c(logLik = round(loglik, 2),
    params = p,
    BIC    = round(bic, 2))
}

# ---------------------------------------------------------------------------
# BIC comparison table
# ---------------------------------------------------------------------------
cat("\n── HMM BIC model comparison ────────────────────────────────────────────\n")

bic_tbl <- rbind(
  "PP 2-state (levels)"   = bic_fn(loglik_pp_2state,   k = 2, n = n),
  "PP 3-state (levels)"   = bic_fn(loglik_pp_3state,   k = 3, n = n),
  "Auto 2-state"          = bic_fn(loglik_auto_2state, k = 2, n = n),
  "Auto 3-state"          = bic_fn(loglik_auto_3state, k = 3, n = n)
)

print(bic_tbl)

# Interpretation:
# Given the small sample size (n = 38), BIC consistently favours the two-state
# specification for both PP and Auto, reflecting the strong penalty imposed on
# additional regime complexity.


###############################################################################
## 9.  EXPECTED REGIME DURATIONS AND STATIONARY DISTRIBUTION
###############################################################################
##
##  For a K-state Markov chain:
##    E[Duration of state k] = 1 / (1 - P_kk)
##
##  Stationary distribution (2-state):
##    π₁ = (1 - P₂₂) / (2 - P₁₁ - P₂₂)
##    π₂ = 1 - π₁
##

# ── PP transition probabilities ───────────────────────────────────────────────
p11_pp <- tmat_pp_lvl[1, 1]   # P(stay in State 1 | currently in State 1)
p22_pp <- tmat_pp_lvl[2, 2]   # P(stay in State 2 | currently in State 2)

dur1_pp <- 1 / (1 - p11_pp)   # expected duration of hard market (years)
dur2_pp <- 1 / (1 - p22_pp)   # expected duration of soft market (years)

pi1_pp  <- (1 - p22_pp) / (2 - p11_pp - p22_pp)   # long-run share in State 1
pi2_pp  <- 1 - pi1_pp

cat("\n── PP regime durations ─────────────────────────────────────────────────\n")
cat("P(stay in S1):", round(p11_pp, 3), "\n")
cat("P(stay in S2):", round(p22_pp, 3), "\n")
cat("E[Duration S1 — Hard market]:", round(dur1_pp, 1), "years\n")
cat("E[Duration S2 — Soft market]:", round(dur2_pp, 1), "years\n")
cat("Long-run share S1 (π₁):", round(pi1_pp * 100, 1), "%\n")
cat("Long-run share S2 (π₂):", round(pi2_pp * 100, 1), "%\n")

# ── Auto transition probabilities ─────────────────────────────────────────────
p11_auto <- tmat_auto_hmm[1, 1]
p22_auto <- tmat_auto_hmm[2, 2]

dur1_auto <- 1 / (1 - p11_auto)
dur2_auto <- 1 / (1 - p22_auto)

pi1_auto  <- (1 - p22_auto) / (2 - p11_auto - p22_auto)
pi2_auto  <- 1 - pi1_auto

cat("\n── Auto regime durations ───────────────────────────────────────────────\n")
cat("P(stay in S1):", round(p11_auto, 3), "\n")
cat("P(stay in S2):", round(p22_auto, 3), "\n")
cat("E[Duration S1 — Benign]:",   round(dur1_auto, 1), "years\n")
cat("E[Duration S2 — Stressed]:", round(dur2_auto, 1), "years\n")
cat("Long-run share S1 (π₁):", round(pi1_auto * 100, 1), "%\n")
cat("Long-run share S2 (π₂):", round(pi2_auto * 100, 1), "%\n")

# ── Regime duration summary table ─────────────────────────────────────────────
regime_tbl <- data.frame(
  Line          = c("Personal Property", "Personal Property", "Auto", "Auto"),
  State         = c("S1 — Hard", "S2 — Soft", "S1 - Benign", "S2 - Stressed"),
  P_stay        = round(c(p11_pp, p22_pp, p11_auto, p22_auto), 3),
  E_Duration_yr = round(c(dur1_pp, dur2_pp, dur1_auto, dur2_auto), 1),
  LongRun_pct   = round(c(pi1_pp, pi2_pp, pi1_auto, pi2_auto) * 100, 1)
)
cat("\n── Regime duration summary ─────────────────────────────────────────────\n")
print(regime_tbl, row.names = FALSE)
# NOTE: PP hard-market duration ≈ 10 years is consistent with Auto's
#       spectral peak at f = 0.1 (period = 10 years) — convergent evidence
#       from two independent analytical frameworks.


###############################################################################
## 10.  OUT-OF-SAMPLE ROLLING FORECAST EVALUATION
###############################################################################
##
##  Protocol: expanding window, train on 1985–2010 (26 obs), forecast 2011–2022
##  (12 one-step-ahead forecasts). Re-estimate all models each year.
##
##  Models compared:
##    PP:   ARIMA(0,1,0) vs 2-state HMM (levels)
##    Auto: ARIMA(0,1,1) vs ARIMA-GARCH(1,1) vs HMM on GARCH residuals
##

train_end  <- 2010
test_start <- 2011
test_end   <- 2022
n_test     <- 12

years_oos   <- test_start:test_end
actual_pp   <- as.numeric(window(pp_ts,   start = test_start, end = test_end))
actual_auto <- as.numeric(window(auto_ts, start = test_start, end = test_end))

# Storage
arima_pp_fc   <- numeric(n_test)
hmm_pp_fc     <- numeric(n_test)
arima_auto_fc <- numeric(n_test)
garch_auto_fc <- numeric(n_test)

# ── PP rolling forecasts ──────────────────────────────────────────
cat("Running PP rolling forecasts...\n")
for (i in seq_len(n_test)) {
  train_pp <- window(pp_ts, end = train_end + i - 1)
  y_tmp    <- as.numeric(train_pp)
  
  # Benchmark: AR(1) on levels (selected model)
  fit_arima_tmp    <- Arima(train_pp, order = c(1, 0, 0))
  arima_pp_fc[i]   <- as.numeric(forecast(fit_arima_tmp, h = 1)$mean)
  
  # HMM: 2-state on levels, regime-conditioned forecast
  # Forecast = weighted average of state means using transition probs
  set.seed(123)
  fit_hmm_tmp <- tryCatch({
    m <- depmix(y_tmp ~ 1,
                data    = data.frame(y_tmp = y_tmp),
                nstates = 2,
                family  = gaussian())
    fit(m, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(fit_hmm_tmp)) {
    pars_tmp     <- getpars(fit_hmm_tmp)
    tmat_tmp     <- matrix(pars_tmp[3:6], nrow = 2, byrow = TRUE)
    means_tmp    <- c(pars_tmp[7], pars_tmp[9])
    post_tmp     <- posterior(fit_hmm_tmp)
    last_state   <- tail(post_tmp$state, 1)
    
    # Regime-conditioned: E[Y_{t+1}] = sum_j P(s_{t+1}=j|s_t) * mu_j
    hmm_pp_fc[i] <- sum(tmat_tmp[last_state, ] * means_tmp)
  } else {
    # Fallback: use ARIMA forecast if HMM fails
    hmm_pp_fc[i] <- arima_pp_fc[i]
    cat("  PP HMM fallback at year", test_start + i - 1, "\n")
  }
  cat("PP", test_start + i - 1, "done |",
      "ARIMA:", round(arima_pp_fc[i], 4),
      "HMM:", round(hmm_pp_fc[i], 4),
      "Actual:", round(actual_pp[i], 4), "\n")
}

# ── Auto rolling forecasts ────────────────────────────────────────
cat("\nRunning Auto rolling forecasts...\n")
for (i in seq_len(n_test)) {
  train_auto    <- window(auto_ts, end = train_end + i - 1)
  train_auto_d1 <- diff(train_auto)
  
  # Benchmark: ARIMA(0,1,1)
  fit_arima_a      <- Arima(train_auto, order = c(0, 1, 1))
  arima_auto_fc[i] <- as.numeric(forecast(fit_arima_a, h = 1)$mean)
  
  # ARIMA(0,1,1)-GARCH(1,1) — point forecast only, no HMM adjustment
  spec_tmp <- ugarchspec(
    variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(0, 1), include.mean = TRUE),
    distribution.model = "norm"
  )
  
  fit_garch_tmp <- tryCatch(
    ugarchfit(spec_tmp, data = train_auto_d1,
              solver = "hybrid", solver.control = list(trace = 0)),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  
  if (!is.null(fit_garch_tmp)) {
    fc_g             <- ugarchforecast(fit_garch_tmp, n.ahead = 1)
    # Convert back to level: last level + forecast of difference
    garch_auto_fc[i] <- tail(as.numeric(train_auto), 1) +
      as.numeric(fitted(fc_g)[1])
    # Sanity check — loss ratio must be positive
    if (garch_auto_fc[i] < 0 | garch_auto_fc[i] > 2) {
      garch_auto_fc[i] <- arima_auto_fc[i]
      cat("  Auto GARCH sanity fallback at year", test_start + i - 1, "\n")
    }
  } else {
    garch_auto_fc[i] <- arima_auto_fc[i]
    cat("  Auto GARCH fit fallback at year", test_start + i - 1, "\n")
  }
  
  cat("Auto", test_start + i - 1, "done |",
      "ARIMA:", round(arima_auto_fc[i], 4),
      "GARCH:", round(garch_auto_fc[i], 4),
      "Actual:", round(actual_auto[i], 4), "\n")
}

# ── Accuracy metrics ──────────────────────────────────────────────
rmse <- function(a, f) sqrt(mean((a - f)^2))
mae  <- function(a, f) mean(abs(a - f))

oos_tbl <- data.frame(
  Line  = c("PP", "PP", "Auto", "Auto"),
  Model = c("AR(1) levels", "HMM 2-state levels",
            "ARIMA(0,1,1)", "ARIMA-GARCH(1,1)"),
  RMSE  = round(c(rmse(actual_pp,   arima_pp_fc),
                  rmse(actual_pp,   hmm_pp_fc),
                  rmse(actual_auto, arima_auto_fc),
                  rmse(actual_auto, garch_auto_fc)), 4),
  MAE   = round(c(mae(actual_pp,   arima_pp_fc),
                  mae(actual_pp,   hmm_pp_fc),
                  mae(actual_auto, arima_auto_fc),
                  mae(actual_auto, garch_auto_fc)), 4)
)

cat("\n── Out-of-sample accuracy (2011-2022) ──────────────────────\n")
print(oos_tbl, row.names = FALSE)

# ── Plot ──────────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# PP
ylim_pp <- range(c(actual_pp, arima_pp_fc, hmm_pp_fc)) + c(-0.03, 0.03)
plot(years_oos, actual_pp, type = "b", pch = 19, lwd = 1.5,
     col = "black", ylim = ylim_pp,
     xlab = "Year", ylab = "Loss ratio",
     main = "PP: Out-of-sample forecasts (2011-2022)")
lines(years_oos, arima_pp_fc, col = "steelblue", lty = 2, lwd = 1.5)
lines(years_oos, hmm_pp_fc,   col = "tomato",    lty = 3, lwd = 1.5)
legend("topright",
       legend = c("Actual", "AR(1)", "HMM"),
       col    = c("black", "steelblue", "tomato"),
       lty    = c(1, 2, 3), pch = c(19, NA, NA),
       bty = "n", cex = 0.85)

# Auto
ylim_auto <- range(c(actual_auto, arima_auto_fc, garch_auto_fc)) + c(-0.02, 0.02)
plot(years_oos, actual_auto, type = "b", pch = 19, lwd = 1.5,
     col = "black", ylim = ylim_auto,
     xlab = "Year", ylab = "Loss ratio",
     main = "Auto: Out-of-sample forecasts (2011-2022)")
lines(years_oos, arima_auto_fc, col = "steelblue", lty = 2, lwd = 1.5)
lines(years_oos, garch_auto_fc, col = "darkgreen", lty = 4, lwd = 1.5)
legend("topright",
       legend = c("Actual", "ARIMA", "ARIMA-GARCH"),
       col    = c("black", "steelblue", "darkgreen"),
       lty    = c(1, 2, 4), pch = c(19, NA, NA),
       bty = "n", cex = 0.85)

par(mfrow = c(1, 1))


###############################################################################
## 11.  REGIME-CONDITIONED FORECASTS (2023–2025)
###############################################################################
##
##  Given the current regime at end of 2022, forecast 3 years ahead.
##  95% prediction interval: forecast ± 1.96 × state SD
##

horizon      <- 3
forecast_yrs <- 2023:2025

# ── 11.1 PP regime-conditioned forecast ───────────────────────────────────────
post_pp_full     <- posterior(fit_hmm_pp_lvl)
current_state_pp <- tail(post_pp_full$state, 1)
fc_mean_pp       <- means_pp_lvl[current_state_pp]
fc_sd_pp         <- sds_pp_lvl[current_state_pp]

pp_fc_tbl <- data.frame(
  Year     = forecast_yrs,
  Forecast = round(rep(fc_mean_pp, horizon), 4),
  Lower95  = round(rep(fc_mean_pp - 1.96 * fc_sd_pp, horizon), 4),
  Upper95  = round(rep(fc_mean_pp + 1.96 * fc_sd_pp, horizon), 4),
  Regime   = ifelse(current_state_pp == 1, "Hard market (S1)", "Soft market (S2)")
)

cat("\n── PP regime-conditioned forecasts 2023–2025 ───────────────────────────\n")
print(pp_fc_tbl, row.names = FALSE)

# ── 11.2 Auto regime-conditioned forecast ─────────────────────────────────────
post_auto_full     <- posterior(fit_hmm_auto)
current_state_auto <- tail(post_auto_full$state, 1)

# 3-step ARIMA-GARCH forecast for Auto level
fc_garch_auto  <- ugarchforecast(fit_garch_auto, n.ahead = horizon)
arima_garch_fc <- tail(as.numeric(auto_ts), 1) +
  cumsum(as.numeric(fitted(fc_garch_auto)))

# Add HMM regime residual adjustment
fc_mean_auto <- means_auto_hmm[current_state_auto]
fc_sd_auto   <- sds_auto_hmm[current_state_auto]

auto_fc_tbl <- data.frame(
  Year     = forecast_yrs,
  Forecast = round(arima_garch_fc + fc_mean_auto, 4),
  Lower95  = round(arima_garch_fc + fc_mean_auto - 1.96 * fc_sd_auto, 4),
  Upper95  = round(arima_garch_fc + fc_mean_auto + 1.96 * fc_sd_auto, 4),
  Regime   = ifelse(current_state_auto == 1, "Benign (S1)", "Stressed (S2)")
)

cat("\n── Auto regime-conditioned forecasts 2023–2025 ─────────────────────────\n")
print(auto_fc_tbl, row.names = FALSE)

# ── 11.3 Forecast plots ───────────────────────────────────────────────────────
par(mfrow = c(1, 2))

# PP
plot(time(pp_ts), as.numeric(pp_ts), type = "l", col = "grey50", lwd = 1.5,
     xlim = c(1985, 2026), ylim = c(0.45, 0.85),
     xlab = "Year", ylab = "Loss ratio",
     main = "PP: Regime-conditioned forecast (2023–2025)")
polygon(c(forecast_yrs, rev(forecast_yrs)),
        c(pp_fc_tbl$Upper95, rev(pp_fc_tbl$Lower95)),
        col = rgb(0.8, 0.2, 0.2, 0.15), border = NA)
lines(forecast_yrs, pp_fc_tbl$Forecast, col = "tomato", lwd = 2)
points(forecast_yrs, pp_fc_tbl$Forecast, col = "tomato", pch = 17, cex = 1.2)
abline(v = 2022.5, lty = 3, col = "grey60")
legend("topleft",
       legend = c("Historical", "Regime forecast", "95% PI"),
       col    = c("grey50", "tomato", rgb(0.8, 0.2, 0.2, 0.4)),
       lty    = c(1, 1, NA), lwd = c(1.5, 2, NA),
       fill   = c(NA, NA, rgb(0.8, 0.2, 0.2, 0.15)),
       border = NA, bty = "n", cex = 0.85)

# Auto
plot(time(auto_ts), as.numeric(auto_ts), type = "l", col = "grey50", lwd = 1.5,
     xlim = c(1985, 2026), ylim = c(0.55, 1.00),
     xlab = "Year", ylab = "Loss ratio",
     main = "Auto: Regime-conditioned forecast (2023–2025)")
polygon(c(forecast_yrs, rev(forecast_yrs)),
        c(auto_fc_tbl$Upper95, rev(auto_fc_tbl$Lower95)),
        col = rgb(0.2, 0.4, 0.8, 0.15), border = NA)
lines(forecast_yrs, auto_fc_tbl$Forecast, col = "steelblue", lwd = 2)
points(forecast_yrs, auto_fc_tbl$Forecast, col = "steelblue", pch = 17, cex = 1.2)
abline(v = 2022.5, lty = 3, col = "grey60")
legend("topleft",
       legend = c("Historical", "Regime forecast", "95% PI"),
       col    = c("grey50", "steelblue", rgb(0.2, 0.4, 0.8, 0.4)),
       lty    = c(1, 1, NA), lwd = c(1.5, 2, NA),
       fill   = c(NA, NA, rgb(0.2, 0.4, 0.8, 0.15)),
       border = NA, bty = "n", cex = 0.85)

par(mfrow = c(1, 1))



