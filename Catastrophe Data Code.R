# ============================================================================
# COMPLETE MMPP ANALYSIS FOR CATASTROPHE INSURANCE

# Author: Subhaan Tahir
# Supervisor: Hong Duong
# ============================================================================

# Load required libraries
library(tidyverse)
library(lubridate)
library(MASS)
library(Matrix)

# Install and load additional packages if needed
if(!require(expm)) {
  install.packages("expm")
  library(expm)
}

if(!require(evd)) {
  install.packages("evd")
  library(evd)
}

# Set seed for reproducibility
set.seed(12345)

# ============================================================================
# DATA PREPROCESSING
# ============================================================================

cat("\n================== CATASTROPHE DATA PREPROCESSING ==================\n\n")

# Read catastrophe data
cat("Reading main catastrophe dataset...\n")
cat_data <- read.csv("C:/Users/subha/Downloads/Catastrophe_Data.csv", stringsAsFactors = FALSE)

# Clean column names
names(cat_data) <- gsub(" ", ".", names(cat_data))
names(cat_data) <- gsub("\\(|\\)|'", "", names(cat_data))

# Check for supplementary US data
us_data_path <- "C:/Users/subha/Downloads/events-US-1980-2024-Q4.csv"
us_data <- NULL

if(file.exists(us_data_path)) {
  cat("Reading supplementary US events data for validation...\n")
  us_data <- read.csv(us_data_path, stringsAsFactors = FALSE)
  names(us_data) <- gsub(" ", ".", names(us_data))
  names(us_data) <- gsub("\\(|\\)|'", "", names(us_data))
  cat("US supplementary data loaded -", nrow(us_data), "events\n")
  
  # IMPORTANT: Check if costs are in millions in supplementary data
  if("Total.Damages..Adjusted....Millions.of.US.." %in% names(us_data)) {
    cat("Note: Supplementary dataset costs are in MILLIONS of USD\n")
  }
}

# Detailed data preparation function
prepare_catastrophe_data <- function(data, threshold = 100000, start_year = 1980) {
  
  cat("Initial dataset size:", nrow(data), "events\n")
  
  # Create proper date column
  data$Start_Date <- as.Date(paste(
    data$Start.Year, 
    ifelse(is.na(data$Start.Month), 1, data$Start.Month),
    ifelse(is.na(data$Start.Day), 1, data$Start.Day), 
    sep = "-"
  ))
  
  # Apply filters step by step
  data_filtered <- data %>%
    filter(!is.na(Start_Date))
  cat("After date validation:", nrow(data_filtered), "events\n")
  
  data_filtered <- data_filtered %>%
    filter(Start.Year >= start_year)
  cat("After year filter (≥", start_year, "):", nrow(data_filtered), "events\n")
  
  # Find damage column
  damage_cols <- names(data_filtered)[grepl("Total.*Damage", names(data_filtered))]
  damage_col <- damage_cols[!grepl("Adjusted", damage_cols)][1]
  if(is.na(damage_col)) damage_col <- damage_cols[1]
  
  if(!is.null(damage_col) && !is.na(damage_col)) {
    cat("Using damage column:", damage_col, "\n")
    
    # Apply damage threshold
    data_filtered <- data_filtered %>%
      filter(!is.na(.data[[damage_col]]))
    cat("With damage data:", nrow(data_filtered), "events\n")
    
    data_filtered <- data_filtered %>%
      filter(.data[[damage_col]] >= threshold)
    cat("Above damage threshold ($", threshold/1000, "M):", nrow(data_filtered), "events\n")
  }
  
  # Sort by date
  data_filtered <- data_filtered %>%
    arrange(Start_Date)
  
  return(data_filtered)
}

# Prepare the data
events <- prepare_catastrophe_data(cat_data)

# Calculate arrival times in years since 1980
start_date <- as.Date("1980-01-01")
events$arrival_time <- as.numeric(difftime(events$Start_Date, start_date, units = "days")) / 365.25

# Remove any negative times
events <- events %>% filter(arrival_time >= 0)

# Extract key variables
arrival_times <- events$arrival_time
n_events <- length(arrival_times)
T_obs <- max(arrival_times)
inter_arrivals <- diff(c(0, arrival_times))

cat("\n--- Final Dataset Summary ---\n")
cat("Number of catastrophic events:", n_events, "\n")
cat("Observation period:", round(T_obs, 2), "years\n")
cat("Average rate:", round(n_events/T_obs, 3), "events/year\n")
cat("Mean inter-arrival time:", round(mean(inter_arrivals), 3), "years\n")
cat("SD of inter-arrival time:", round(sd(inter_arrivals), 3), "years\n\n")

# DATA UNITS WARNING
damage_col <- "Total.Damage...000.US.."
if(damage_col %in% names(events)) {
  cat("\n!!! DATA UNITS WARNING !!!\n")
  cat("The column name suggests damages are in '000s (thousands) of US$\n")
  cat("Mean damage value:", mean(events[[damage_col]], na.rm=TRUE), "\n")
  cat("If this is $", mean(events[[damage_col]], na.rm=TRUE), 
      "K, then actual mean is $", mean(events[[damage_col]], na.rm=TRUE)/1000, "M\n")
  cat("Please verify the correct units for damage values\n\n")
}

# ============================================================================
# HOMOGENEOUS POISSON PROCESS (HPP)
# ============================================================================

cat("================== HOMOGENEOUS POISSON PROCESS (HPP) ==================\n\n")

# Maximum Likelihood Estimation for HPP
hpp_rate <- n_events / T_obs

# Log-likelihood for HPP
hpp_loglik <- n_events * log(hpp_rate) - hpp_rate * T_obs

# Standard error using Fisher Information
hpp_se <- sqrt(hpp_rate / T_obs)

# 95% Confidence interval for rate
hpp_ci_lower <- hpp_rate - 1.96 * hpp_se
hpp_ci_upper <- hpp_rate + 1.96 * hpp_se

cat("HPP Parameter Estimates:\n")
cat("─────────────────────────\n")
cat("Rate (λ):", sprintf("%.6f", hpp_rate), "events/year\n")
cat("Standard Error:", sprintf("%.6f", hpp_se), "\n")
cat("95% CI: [", sprintf("%.6f", hpp_ci_lower), ", ", sprintf("%.6f", hpp_ci_upper), "]\n")
cat("Log-likelihood:", sprintf("%.2f", hpp_loglik), "\n\n")

# ============================================================================
# MARKOV-MODULATED POISSON PROCESS (MMPP)
# ============================================================================

cat("================== MARKOV-MODULATED POISSON PROCESS (MMPP) ==================\n\n")

# Function for 2-state MMPP using EM algorithm
fit_mmpp_em <- function(times, T_max, max_iter = 100, tol = 1e-6) {
  
  n <- length(times)
  if(n == 0) stop("No events to fit")
  
  inter_times <- diff(c(0, times))
  
  cat("Fitting 2-state MMPP using EM algorithm\n")
  cat("Number of events:", n, "\n")
  cat("Observation period:", round(T_max, 2), "years\n\n")
  
  # Initialize parameters using quantile-based approach
  q25 <- quantile(inter_times, 0.25)
  q75 <- quantile(inter_times, 0.75)
  
  # Initial rates (ensure λ1 < λ2)
  lambda1_init <- 1 / max(q75, 0.001)
  lambda2_init <- 1 / max(q25, 0.0001)
  
  if(lambda1_init >= lambda2_init) {
    lambda1_init <- lambda2_init * 0.5
  }
  
  lambda <- c(lambda1_init, lambda2_init)
  
  # Initial transition rates (continuous time)
  q12 <- 1.0
  q21 <- 0.5
  
  Q <- matrix(c(-q12, q12, q21, -q21), nrow = 2, byrow = TRUE)
  
  # Function to compute stationary distribution
  compute_stationary <- function(Q) {
    q12 <- Q[1,2]
    q21 <- Q[2,1]
    if(q12 + q21 == 0) return(c(0.5, 0.5))
    pi1 <- q21 / (q12 + q21)
    pi2 <- q12 / (q12 + q21)
    return(c(pi1, pi2))
  }
  
  loglik_old <- -Inf
  
  for(iter in 1:max_iter) {
    
    # E-step: Compute posterior probabilities
    # Forward pass
    alpha <- matrix(0, n, 2)
    pi_init <- compute_stationary(Q)
    
    # First observation
    alpha[1,] <- pi_init * lambda * exp(-lambda * inter_times[1])
    alpha[1,] <- alpha[1,] / sum(alpha[1,])
    
    # Subsequent observations
    if(n > 1) {
      for(i in 2:n) {
        dt <- inter_times[i]
        # Transition probabilities
        P_trans <- expm(Q * dt)
        # Update
        for(j in 1:2) {
          alpha[i,j] <- sum(alpha[i-1,] * P_trans[,j]) * lambda[j] * exp(-lambda[j] * dt)
        }
        # Normalize
        if(sum(alpha[i,]) > 0) {
          alpha[i,] <- alpha[i,] / sum(alpha[i,])
        }
      }
    }
    
    # Backward pass
    beta <- matrix(0, n, 2)
    beta[n,] <- 1
    
    if(n > 1) {
      for(i in (n-1):1) {
        dt <- inter_times[i+1]
        P_trans <- expm(Q * dt)
        for(j in 1:2) {
          beta[i,j] <- sum(P_trans[j,] * lambda * exp(-lambda * dt) * beta[i+1,])
        }
        if(sum(beta[i,]) > 0) {
          beta[i,] <- beta[i,] / sum(beta[i,])
        }
      }
    }
    
    # Compute gamma (posterior probabilities)
    gamma <- alpha * beta
    gamma <- gamma / rowSums(gamma + 1e-10)
    
    # M-step: Update parameters
    
    # Update lambdas
    for(j in 1:2) {
      expected_events <- sum(gamma[,j])
      expected_time <- sum(gamma[,j] * inter_times)
      if(expected_time > 0) {
        lambda[j] <- expected_events / expected_time
      }
    }
    
    # Ensure ordering
    if(lambda[1] > lambda[2]) {
      lambda <- rev(lambda)
      gamma <- gamma[,2:1]
    }
    
    # Update Q matrix (simplified approach)
    if(n > 1) {
      runs <- rle(apply(gamma, 1, which.max))
      
      if(length(runs$lengths[runs$values == 1]) > 0) {
        mean_sojourn_1 <- mean(runs$lengths[runs$values == 1]) * mean(inter_times)
        q12 <- 1 / max(mean_sojourn_1, 0.1)
      } else {
        q12 <- 1.0
      }
      
      if(length(runs$lengths[runs$values == 2]) > 0) {
        mean_sojourn_2 <- mean(runs$lengths[runs$values == 2]) * mean(inter_times)
        q21 <- 1 / max(mean_sojourn_2, 0.1)
      } else {
        q21 <- 0.5
      }
      
      Q <- matrix(c(-q12, q12, q21, -q21), nrow = 2, byrow = TRUE)
    }
    
    # Compute log-likelihood
    loglik <- 0
    pi_stat <- compute_stationary(Q)
    
    # Initial state
    loglik <- loglik + log(sum(pi_stat * lambda * exp(-lambda * inter_times[1])))
    
    # Subsequent events
    if(n > 1) {
      for(i in 2:n) {
        dt <- inter_times[i]
        P_trans <- expm(Q * dt)
        trans_prob <- gamma[i-1,] %*% P_trans
        event_prob <- sum(trans_prob * lambda * exp(-lambda * dt))
        loglik <- loglik + log(max(event_prob, 1e-10))
      }
    }
    
    # Survival probability for last interval
    dt_last <- T_max - times[n]
    if(dt_last > 0) {
      surv_prob <- sum(gamma[n,] * exp(-lambda * dt_last))
      loglik <- loglik + log(max(surv_prob, 1e-10))
    }
    
    # Check convergence
    if(abs(loglik - loglik_old) < tol) {
      cat("Converged after", iter, "iterations\n")
      break
    }
    
    loglik_old <- loglik
  }
  
  if(iter == max_iter) {
    cat("Maximum iterations reached\n")
  }
  
  # Final stationary distribution
  pi_stat <- compute_stationary(Q)
  
  # Build transition probability matrix for display
  P_display <- expm(Q * 1)  # One-year transition matrix
  
  result <- list(
    lambda = lambda,
    Q = Q,
    P = P_display,
    pi = pi_stat,
    gamma = gamma,
    loglik = loglik,
    converged = (iter < max_iter),
    n_iterations = iter
  )
  
  return(result)
}

# Fit MMPP
mmpp_result <- fit_mmpp_em(arrival_times, T_obs)

# Display results
cat("\nMMPP Parameter Estimates:\n")
cat("─────────────────────────\n")
cat("State 1 (Dormant) rate (λ₁):", sprintf("%.6f", mmpp_result$lambda[1]), "events/year\n")
cat("State 2 (Active) rate (λ₂):", sprintf("%.6f", mmpp_result$lambda[2]), "events/year\n")
cat("Rate ratio (λ₂/λ₁):", sprintf("%.2f", mmpp_result$lambda[2]/mmpp_result$lambda[1]), "\n\n")

cat("Generator Matrix Q:\n")
print(round(mmpp_result$Q, 4))

cat("\nOne-year Transition Matrix:\n")
print(round(mmpp_result$P, 4))

cat("\nStationary Distribution:\n")
cat("P(Dormant):", sprintf("%.3f", mmpp_result$pi[1]), "\n")
cat("P(Active):", sprintf("%.3f", mmpp_result$pi[2]), "\n")
cat("Log-likelihood:", sprintf("%.2f", mmpp_result$loglik), "\n")
cat("Iterations to convergence:", mmpp_result$n_iterations, "\n\n")

# ============================================================================
# LIKELIHOOD RATIO TEST
# ============================================================================

cat("================== LIKELIHOOD RATIO TEST ==================\n\n")

LR_statistic <- 2 * (mmpp_result$loglik - hpp_loglik)
df <- 3
p_value <- pchisq(LR_statistic, df = df, lower.tail = FALSE)

cat("H₀: HPP is adequate (no regime switching needed)\n")
cat("H₁: MMPP is preferred (regime switching present)\n\n")
cat("Test Statistic: ", sprintf("%.3f", LR_statistic), "\n")
cat("Degrees of Freedom: ", df, "\n")
cat("p-value: ", sprintf("%.6f", p_value), "\n")
cat("Decision at α=0.05: ", ifelse(p_value < 0.05, 
                                   "Reject H₀ - MMPP significantly better", 
                                   "Fail to reject H₀"), "\n\n")

# ============================================================================
# MODEL COMPARISON (INFORMATION CRITERIA)
# ============================================================================

cat("================== MODEL COMPARISON ==================\n\n")

k_hpp <- 1
k_mmpp <- 4

AIC_hpp <- -2 * hpp_loglik + 2 * k_hpp
AIC_mmpp <- -2 * mmpp_result$loglik + 2 * k_mmpp
BIC_hpp <- -2 * hpp_loglik + k_hpp * log(n_events)
BIC_mmpp <- -2 * mmpp_result$loglik + k_mmpp * log(n_events)

cat("Information Criteria:\n")
cat("─────────────────────\n")
cat(sprintf("%-10s AIC: %8.2f  BIC: %8.2f\n", "HPP", AIC_hpp, BIC_hpp))
cat(sprintf("%-10s AIC: %8.2f  BIC: %8.2f\n", "MMPP", AIC_mmpp, BIC_mmpp))
cat("\n")
cat("ΔAIC (HPP - MMPP): ", sprintf("%6.2f", AIC_hpp - AIC_mmpp), 
    " (", ifelse(AIC_mmpp < AIC_hpp, "MMPP preferred", "HPP preferred"), ")\n")
cat("ΔBIC (HPP - MMPP): ", sprintf("%6.2f", BIC_hpp - BIC_mmpp),
    " (", ifelse(BIC_mmpp < BIC_hpp, "MMPP preferred", "HPP preferred"), ")\n\n")

# ============================================================================
# DIAGNOSTIC TESTS
# ============================================================================

cat("================== DIAGNOSTIC TESTS ==================\n\n")

# 1. OVERDISPERSION TEST
cat("1. Overdispersion Analysis:\n")
cat("───────────────────────────\n")

monthly_bins <- seq(0, ceiling(T_obs * 12) + 1) / 12
monthly_counts <- hist(arrival_times, breaks = monthly_bins, plot = FALSE)$counts

if(length(monthly_counts) > floor(T_obs * 12)) {
  monthly_counts <- monthly_counts[1:floor(T_obs * 12)]
}

mean_count <- mean(monthly_counts)
var_count <- var(monthly_counts)
dispersion_index <- var_count / mean_count

cat("Mean monthly count: ", sprintf("%.3f", mean_count), "\n")
cat("Variance of monthly counts: ", sprintf("%.3f", var_count), "\n")
cat("Index of Dispersion: ", sprintf("%.3f", dispersion_index), "\n")
cat("Interpretation: ")
if(dispersion_index > 1.5) {
  cat("Strong overdispersion detected (supports MMPP)\n")
} else if(dispersion_index > 1.2) {
  cat("Moderate overdispersion detected\n")
} else {
  cat("No significant overdispersion (HPP may be adequate)\n")
}

# 2. KOLMOGOROV-SMIRNOV TEST
cat("\n2. Kolmogorov-Smirnov Test for Exponential Inter-arrivals:\n")
cat("──────────────────────────────────────────────────────────\n")

if(any(duplicated(inter_arrivals))) {
  inter_arrivals_test <- inter_arrivals + rnorm(length(inter_arrivals), 0, 1e-10)
  ks_test <- ks.test(inter_arrivals_test, "pexp", rate = hpp_rate)
  cat("Note: Small jitter added to handle tied values\n")
} else {
  ks_test <- ks.test(inter_arrivals, "pexp", rate = hpp_rate)
}

cat("Test Statistic: ", sprintf("%.4f", ks_test$statistic), "\n")
cat("p-value: ", sprintf("%.6f", ks_test$p.value), "\n")
cat("Decision at α=0.05: ", ifelse(ks_test$p.value < 0.05,
                                   "Reject exponentiality (supports MMPP)",
                                   "Cannot reject exponentiality"), "\n")

# 3. TEMPORAL CLUSTERING
cat("\n3. Temporal Clustering Detection:\n")
cat("──────────────────────────────────\n")

acf_result <- acf(monthly_counts, lag.max = 12, plot = FALSE)
significant_lags <- sum(abs(acf_result$acf[-1]) > 2/sqrt(length(monthly_counts)))

cat("Number of significant lags (out of 12): ", significant_lags, "\n")
cat("Clustering present: ", ifelse(significant_lags > 0, "Yes", "No"), "\n")

lb_test <- Box.test(monthly_counts, lag = 6, type = "Ljung-Box")
cat("Ljung-Box test p-value: ", sprintf("%.4f", lb_test$p.value), "\n")
cat("Decision: ", ifelse(lb_test$p.value < 0.05,
                         "Significant autocorrelation (clustering)",
                         "No significant autocorrelation"), "\n\n")

# ============================================================================
# INSURANCE RISK METRICS
# ============================================================================

cat("================== INSURANCE RISK METRICS ==================\n\n")

T_horizon <- 1

cat("Risk Metrics for ", T_horizon, "-Year Horizon:\n")
cat("────────────────────────────────────────\n\n")

# HPP metrics
hpp_expected <- hpp_rate * T_horizon
hpp_var90 <- qpois(0.90, hpp_expected)
hpp_var95 <- qpois(0.95, hpp_expected)
hpp_var99 <- qpois(0.99, hpp_expected)
hpp_var995 <- qpois(0.995, hpp_expected)

cat("HPP Model:\n")
cat("  Expected events: ", sprintf("%.2f", hpp_expected), "\n")
cat("  VaR(90%): ", hpp_var90, " events\n")
cat("  VaR(95%): ", hpp_var95, " events\n")
cat("  VaR(99%): ", hpp_var99, " events\n")
cat("  VaR(99.5%): ", hpp_var995, " events (Solvency II standard)\n\n")

# MMPP metrics - proper simulation
n_sim <- 10000
set.seed(12345)
mmpp_simulated <- numeric(n_sim)

for(i in 1:n_sim) {
  state <- sample(1:2, 1, prob = mmpp_result$pi)
  
  t_current <- 0
  n_events_sim <- 0
  
  while(t_current < T_horizon) {
    t_event <- rexp(1, mmpp_result$lambda[state])
    
    if(state == 1) {
      t_trans <- rexp(1, abs(mmpp_result$Q[1,2]))
    } else {
      t_trans <- rexp(1, abs(mmpp_result$Q[2,1]))
    }
    
    if(t_event < t_trans && t_current + t_event < T_horizon) {
      n_events_sim <- n_events_sim + 1
      t_current <- t_current + t_event
    } else if(t_trans < t_event && t_current + t_trans < T_horizon) {
      t_current <- t_current + t_trans
      state <- 3 - state
    } else {
      break
    }
  }
  
  mmpp_simulated[i] <- n_events_sim
}

mmpp_expected <- mean(mmpp_simulated)
mmpp_var90 <- quantile(mmpp_simulated, 0.90)
mmpp_var95 <- quantile(mmpp_simulated, 0.95)
mmpp_var99 <- quantile(mmpp_simulated, 0.99)
mmpp_var995 <- quantile(mmpp_simulated, 0.995)

cat("MMPP Model:\n")
cat("  Expected events: ", sprintf("%.2f", mmpp_expected), "\n")
cat("  VaR(90%): ", mmpp_var90, " events\n")
cat("  VaR(95%): ", mmpp_var95, " events\n")
cat("  VaR(99%): ", mmpp_var99, " events\n")
cat("  VaR(99.5%): ", mmpp_var995, " events (Solvency II standard)\n\n")

# Capital Requirements
cat("Capital Requirements Analysis:\n")
cat("─────────────────────────────────\n")

capital_increase_pct <- (mmpp_var995 - hpp_var995) / hpp_var995 * 100

cat("Solvency II Capital (99.5% VaR):\n")
cat("  HPP-based: ", hpp_var995, " events\n")
cat("  MMPP-based: ", mmpp_var995, " events\n")
cat("  Increase: ", sprintf("%.1f%%", capital_increase_pct), "\n\n")

# ============================================================================
# 5. SEVERITY INTEGRATION - REGIME-DEPENDENT SEVERITY ANALYSIS
# ============================================================================

cat("\n================== SEVERITY INTEGRATION BY REGIME ==================\n\n")

damage_col <- "Total.Damage...000.US.."
if(damage_col %in% names(events) && !is.null(mmpp_result$gamma)) {
  
  # Convert damages from thousands to millions
  damages_millions <- events[[damage_col]] / 1000
  regime_assignments <- apply(mmpp_result$gamma, 1, which.max)
  
  cat("Severity Analysis by Regime (in millions USD):\n")
  cat("────────────────────────────────────────────\n\n")
  
  regime_severity <- list()
  
  for(regime in 1:2) {
    regime_damages <- damages_millions[regime_assignments == regime]
    regime_damages <- regime_damages[!is.na(regime_damages)]
    
    if(length(regime_damages) > 10) {
      cat("Regime", regime, "(", ifelse(regime == 1, "Dormant", "Active"), "):\n")
      cat("  Number of events:", length(regime_damages), "\n")
      cat("  Mean severity: $", sprintf("%.2f", mean(regime_damages)), "M\n")
      cat("  Median severity: $", sprintf("%.2f", median(regime_damages)), "M\n")
      cat("  95th percentile: $", sprintf("%.2f", quantile(regime_damages, 0.95)), "M\n")
      cat("  Max severity: $", sprintf("%.2f", max(regime_damages)), "M\n\n")
      
      regime_severity[[regime]] <- regime_damages
    }
  }
  
  if(length(regime_severity) == 2) {
    wilcox_test <- wilcox.test(regime_severity[[1]], regime_severity[[2]])
    cat("Wilcoxon Test for Severity Difference:\n")
    cat("  Test statistic:", wilcox_test$statistic, "\n")
    cat("  p-value:", sprintf("%.6f", wilcox_test$p.value), "\n")
    cat("  Conclusion:", ifelse(wilcox_test$p.value < 0.05,
                                "Significant difference in severity between regimes",
                                "No significant difference in severity"), "\n\n")
  }
}

# ============================================================================
# 7. SENSITIVITY ANALYSIS
# ============================================================================

cat("\n================== SENSITIVITY ANALYSIS ==================\n\n")

cat("Sensitivity to Damage Threshold:\n")
cat("─────────────────────────────────\n")

thresholds <- c(50000, 100000, 200000, 500000)
sensitivity_results <- data.frame()

for(thresh in thresholds) {
  events_temp <- prepare_catastrophe_data(cat_data, threshold = thresh)
  
  if(nrow(events_temp) > 100) {
    events_temp$arrival_time <- as.numeric(
      difftime(events_temp$Start_Date, as.Date("1980-01-01"), units = "days")
    ) / 365.25
    
    arrival_times_temp <- events_temp$arrival_time
    n_events_temp <- length(arrival_times_temp)
    T_obs_temp <- max(arrival_times_temp)
    
    hpp_rate_temp <- n_events_temp / T_obs_temp
    
    inter_arrivals_temp <- diff(c(0, arrival_times_temp))
    q25 <- quantile(inter_arrivals_temp, 0.25)
    q75 <- quantile(inter_arrivals_temp, 0.75)
    
    lambda1_est <- 1 / q75
    lambda2_est <- 1 / q25
    rate_ratio_est <- lambda2_est / lambda1_est
    
    monthly_bins <- seq(0, ceiling(T_obs_temp * 12)) / 12
    monthly_counts_temp <- hist(arrival_times_temp, breaks = monthly_bins, plot = FALSE)$counts
    dispersion_index_temp <- var(monthly_counts_temp) / mean(monthly_counts_temp)
    
    sensitivity_results <- rbind(sensitivity_results, data.frame(
      Threshold = thresh / 1000,
      N_Events = n_events_temp,
      HPP_Rate = hpp_rate_temp,
      Rate_Ratio = rate_ratio_est,
      Dispersion = dispersion_index_temp
    ))
  }
}

print(sensitivity_results)

# ============================================================================
# 8. ECONOMIC IMPACT QUANTIFICATION (WITH CORRECTED UNITS)
# ============================================================================

cat("\n\n================== ECONOMIC IMPACT QUANTIFICATION ==================\n\n")

portfolio_size <- 10000  # $10 billion
annual_premium_rate <- 0.02
reinsurance_attachment <- 1000

cat("Portfolio Assumptions:\n")
cat("─────────────────────\n")
cat("Total Exposure: $", portfolio_size/1000, " billion\n")
cat("Annual Premium Income: $", portfolio_size * annual_premium_rate, " million\n")
cat("Reinsurance Attachment: $", reinsurance_attachment, " million\n\n")

cat("Solvency Capital Requirement (SCR) Impact:\n")
cat("──────────────────────────────────────────\n")

# Convert from thousands to millions
avg_severity_millions <- mean(events[[damage_col]], na.rm = TRUE) / 1000
median_severity_millions <- median(events[[damage_col]], na.rm = TRUE) / 1000

cat("Severity statistics (millions USD):\n")
cat("  Mean severity: $", sprintf("%.2f", avg_severity_millions), "M\n")
cat("  Median severity: $", sprintf("%.2f", median_severity_millions), "M\n\n")

# Use median for more realistic capital calculations
hpp_capital <- hpp_var995 * median_severity_millions
mmpp_capital <- mmpp_var995 * median_severity_millions
capital_difference <- mmpp_capital - hpp_capital

cat("Using median severity for capital calculations:\n")
cat("HPP-based SCR: $", sprintf("%.0f", hpp_capital), " million\n")
cat("MMPP-based SCR: $", sprintf("%.0f", mmpp_capital), " million\n")
cat("Additional Capital Required: $", sprintf("%.0f", capital_difference), " million\n")
cat("As % of Premium: ", sprintf("%.1f%%", 
                                 capital_difference / (portfolio_size * annual_premium_rate) * 100), "\n\n")

# ============================================================================
# 16. EXTREME VALUE THEORY INTEGRATION
# ============================================================================

cat("\n\n================== EXTREME VALUE THEORY ANALYSIS ==================\n\n")

damages_evt <- events[[damage_col]][!is.na(events[[damage_col]])]
# Convert to millions
damages_evt_millions <- damages_evt / 1000

cat("Threshold Selection for GPD:\n")
cat("────────────────────────────\n")

threshold_evt <- quantile(damages_evt_millions, 0.90)
cat("Selected threshold: $", sprintf("%.2f", threshold_evt), 
    " million (90th percentile)\n\n")

excesses <- damages_evt_millions[damages_evt_millions > threshold_evt] - threshold_evt
n_excesses <- length(excesses)

cat("Fitting Generalized Pareto Distribution:\n")
cat("────────────────────────────────────────\n")
cat("Number of exceedances: ", n_excesses, "\n")
cat("Exceedance rate: ", sprintf("%.3f", n_excesses/n_events), "\n\n")

if(n_excesses > 30) {
  gpd_fit <- fpot(damages_evt_millions, threshold = threshold_evt, std.err = TRUE)
  
  cat("GPD Parameters:\n")
  cat("  Scale (σ): $", sprintf("%.2f", gpd_fit$estimate[1]), " million\n")
  cat("  Shape (ξ): ", sprintf("%.4f", gpd_fit$estimate[2]), "\n\n")
  
  xi <- gpd_fit$estimate[2]
  
  if(xi > 0) {
    cat("Tail Interpretation: Heavy-tailed (ξ > 0)\n")
    cat("  Tail index α = 1/ξ = ", sprintf("%.2f", 1/xi), "\n\n")
  }
  
  cat("Extreme Quantiles (Return Levels in millions USD):\n")
  cat("────────────────────────────────────────────────\n")
  
  return_periods <- c(100, 250, 500, 1000)
  for(rp in return_periods) {
    p_exceed_threshold <- n_excesses / n_events
    p_rp <- 1 - 1/(rp * hpp_rate)
    
    if(xi != 0) {
      gpd_quantile <- threshold_evt + (gpd_fit$estimate[1]/xi) * 
        ((((1-p_rp)/p_exceed_threshold)^(-xi)) - 1)
    } else {
      gpd_quantile <- threshold_evt + gpd_fit$estimate[1] * 
        log((1-p_rp)/p_exceed_threshold)
    }
    
    cat(sprintf("%4d", rp), "-year return level: $", 
        sprintf("%.2f", gpd_quantile), " million\n")
  }
}

cat("\n\n================== GENERATING VISUALIZATIONS ==================\n\n")

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Plot 1: Cumulative event arrivals
plot(arrival_times, 1:n_events, type = "s", lwd = 2,
     xlab = "Years since 1980", ylab = "Cumulative events",
     main = "Cumulative Event Arrivals")
abline(a = 0, b = hpp_rate, col = "red", lty = 2, lwd = 2)
legend("topleft", c("Observed", "HPP expected"), 
       col = c("black", "red"), lty = c(1, 2), lwd = 2)

# Plot 2: Inter-arrival time distribution
hist(inter_arrivals, breaks = 30, freq = FALSE,
     xlab = "Inter-arrival time (years)", 
     main = "Inter-arrival Distribution",
     col = "lightgray")
curve(dexp(x, rate = hpp_rate), add = TRUE, col = "red", lwd = 2)

# Plot 3: Monthly event counts
hist(monthly_counts, breaks = seq(-0.5, max(monthly_counts) + 0.5, 1),
     xlab = "Events per month", main = "Monthly Event Distribution",
     col = "lightblue")
abline(v = mean_count, col = "red", lwd = 2, lty = 2)
abline(v = median(monthly_counts), col = "blue", lwd = 2, lty = 3)
legend("topright", c("Mean", "Median"), col = c("red", "blue"), lty = c(2, 3))

# Plot 4: Autocorrelation function
acf(monthly_counts, main = "Autocorrelation of Monthly Counts", lag.max = 12)

# Plot 5: Q-Q plot for exponential inter-arrivals
qqplot(qexp(ppoints(length(inter_arrivals)), rate = hpp_rate),
       inter_arrivals,
       main = "Q-Q Plot: Exponential Fit",
       xlab = "Theoretical quantiles",
       ylab = "Sample quantiles")
abline(0, 1, col = "red", lty = 2, lwd = 2)

# Plot 6: Model comparison
comparison_data <- matrix(c(AIC_hpp, AIC_mmpp, BIC_hpp, BIC_mmpp), 
                          nrow = 2, byrow = TRUE)
barplot(comparison_data, beside = TRUE,
        names.arg = c("HPP", "MMPP"),
        main = "Model Comparison",
        ylab = "Information Criterion",
        col = c("lightcoral", "lightblue"),
        legend.text = c("AIC", "BIC"),
        args.legend = list(x = "topright"))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n================== DISSERTATION SUMMARY ==================\n\n")

evidence_score <- 0
evidence_summary <- list()

if(p_value < 0.05) {
  evidence_score <- evidence_score + 1
  evidence_summary$lr_test <- "Significant (supports MMPP)"
} else {
  evidence_summary$lr_test <- "Not significant"
}

if(AIC_mmpp < AIC_hpp) {
  evidence_score <- evidence_score + 1
  evidence_summary$aic <- "MMPP preferred"
} else {
  evidence_summary$aic <- "HPP preferred"
}

if(BIC_mmpp < BIC_hpp) {
  evidence_score <- evidence_score + 1
  evidence_summary$bic <- "MMPP preferred"
} else {
  evidence_summary$bic <- "HPP preferred"
}

if(dispersion_index > 1.5) {
  evidence_score <- evidence_score + 1
  evidence_summary$overdispersion <- "Strong (supports MMPP)"
}

if(ks_test$p.value < 0.05) {
  evidence_score <- evidence_score + 1
  evidence_summary$exponential <- "Rejected (supports MMPP)"
}

if(significant_lags > 0) {
  evidence_score <- evidence_score + 1
  evidence_summary$clustering <- "Present (supports MMPP)"
}

cat("EVIDENCE SUMMARY (", evidence_score, "/6 criteria support MMPP):\n")
cat("────────────────────────────────────────────────────\n")
cat("1. Likelihood Ratio Test: ", evidence_summary$lr_test, "\n")
cat("2. AIC Selection: ", evidence_summary$aic, "\n")
cat("3. BIC Selection: ", evidence_summary$bic, "\n")
cat("4. Overdispersion: ", evidence_summary$overdispersion, "\n")
cat("5. Exponentiality: ", evidence_summary$exponential, "\n")
cat("6. Temporal Clustering: ", evidence_summary$clustering, "\n\n")

cat("KEY FINDINGS:\n")
cat("─────────────\n")
cat("• Rate intensification during active periods: ", 
    sprintf("%.1fx", mmpp_result$lambda[2]/mmpp_result$lambda[1]), "\n")
cat("• Time in active state: ", sprintf("%.1f%%", mmpp_result$pi[2] * 100), "\n")
cat("• Capital requirement increase: ", sprintf("%.1f%%", capital_increase_pct), "\n")

if(evidence_score >= 4) {
  cat("\nCONCLUSION: STRONG SUPPORT for Hypothesis 2\n")
  cat("The MMPP significantly outperforms the HPP in modeling catastrophe arrivals.\n")
}

cat("\n────────────────────────────────────────────────────────────────────────\n")
cat("Analysis complete with enhanced robustness checks and validation.\n")
cat("────────────────────────────────────────────────────────────────────────\n\n")

# Save results
results_for_dissertation <- list(
  hpp = list(
    rate = hpp_rate,
    se = hpp_se,
    ci = c(hpp_ci_lower, hpp_ci_upper),
    loglik = hpp_loglik
  ),
  mmpp = list(
    lambda = mmpp_result$lambda,
    Q = mmpp_result$Q,
    pi = mmpp_result$pi,
    loglik = mmpp_result$loglik,
    converged = mmpp_result$converged
  ),
  insurance_metrics = list(
    capital = list(
      hpp = hpp_var995,
      mmpp = mmpp_var995,
      increase_pct = capital_increase_pct,
      capital_difference = capital_difference
    )
  ),
  evidence_score = evidence_score
)

saveRDS(results_for_dissertation, "dissertation_mmpp_complete_analysis.rds")
