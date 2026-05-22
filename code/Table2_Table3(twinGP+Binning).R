data_path <- file.path("data")
terrain_path <- file.path("data", "weightedTerrainData.csv")
results_dir <- file.path("results","intermediate")

library(data.table)
library(twingp)
library(dplyr)
library(readr)
library(tidyr)



turbine_ids <- 1:66
testset_2018 <- c(1:46, 48:50, 52, 54:60, 62:66)
exclude_ids <- c(47, 51, 53, 61)


terrain_data <- read.csv(terrain_path)

scale_01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
terrain_data[, 2:4] <- lapply(terrain_data[, 2:4, drop = FALSE], scale_01)

terrain_mat <- as.matrix(terrain_data[1:66, 2:4])
storage.mode(terrain_mat) <- "double"

# -------------------------
# Binning (speed only)
# -------------------------
binning <- function(train.x, train.y, test.x, bin_width = 0.5){
  train.y <- as.numeric(train.y)
  train.x <- as.numeric(train.x)
  test.x  <- as.numeric(test.x)
  
  start <- 0
  end <- round(max(train.x, na.rm = TRUE))
  n_bins <- round((end - start)/bin_width, 0) + 1
  x_bin <- y_bin <- numeric(n_bins)
  
  for (n in 2:n_bins){
    bin_element <- which(train.x > (start + (n-1)*bin_width) & train.x < (start + n*bin_width))
    x_bin[n] <- mean(train.x[bin_element], na.rm = TRUE)
    y_bin[n] <- mean(train.y[bin_element], na.rm = TRUE)
  }
  
  binned_data <- data.frame(x_bin, y_bin)
  binned_data <- binned_data[!is.na(binned_data$y_bin), ]
  splinefit <- smooth.spline(x = binned_data$x_bin, y = binned_data$y_bin, all.knots = TRUE)
  y_pred <- predict(splinefit, test.x)$y
  y_pred[y_pred < 0] <- 0
  y_pred
}

# -------------------------
# Load all turbine-year CSVs ONCE (cache)
# -------------------------
cache_key <- function(id, year) sprintf("T%02d_%d", id, year)

data_cache <- vector("list", length = length(turbine_ids) * 2)
names(data_cache) <- as.vector(outer(sprintf("T%02d", turbine_ids), c("2017", "2018"), paste, sep = "_"))

for (id in turbine_ids) {
  for (yr in c(2017, 2018)) {
    f <- sprintf("%s/Turbine%d_%d.csv", data_path, id, yr)
    
    d <- tryCatch(
      fread(file = f, showProgress = FALSE),   # <-- IMPORTANT: file=
      error = function(e) NULL
    )
    
    if (!is.null(d)) {
      d[, wind_speed := as.numeric(wind_speed)]
      d[, temperature := as.numeric(temperature)]
      d[, power := as.numeric(power)]
    }
    
    data_cache[[cache_key(id, yr)]] <- d
  }
}

get_data <- function(id, year) data_cache[[cache_key(id, year)]]

# -------------------------
# Build full 2017 training pool ONCE (X, y, tid, and terrain S)
# -------------------------
X2017_list <- vector("list", 66)
y2017_list <- vector("list", 66)
tid2017_list <- vector("list", 66)

for (j in turbine_ids) {
  d <- get_data(j, 2017)
  if (!is.null(d)) {
    X2017_list[[j]] <- cbind(d$wind_speed, d$temperature)
    y2017_list[[j]] <- d$power
    tid2017_list[[j]] <- rep(j, nrow(d))
  } else {
    X2017_list[[j]] <- matrix(numeric(0), ncol = 2)
    y2017_list[[j]] <- numeric(0)
    tid2017_list[[j]] <- integer(0)
  }
}

X2017  <- do.call(rbind, X2017_list)
y2017  <- as.numeric(unlist(y2017_list, use.names = FALSE))
tid2017 <- as.integer(unlist(tid2017_list, use.names = FALSE))

stopifnot(nrow(X2017) == length(y2017), length(y2017) == length(tid2017))

# Terrain per training row (precomputed once)
S2017 <- terrain_mat[tid2017, , drop = FALSE]

# Helper: terrain for test turbine i repeated n times
S_test_for <- function(i, n) matrix(rep(terrain_mat[i, ], each = n), nrow = n, ncol = 3)

# Helper: fast LOO slice from full 2017 pool
get_train_loo_2017 <- function(leave_out_id) {
  idx <- tid2017 != leave_out_id
  list(
    X = X2017[idx, , drop = FALSE],
    y = y2017[idx],
    tid = tid2017[idx],
    S = S2017[idx, , drop = FALSE]
  )
}

# -------------------------
# Main results collector
# -------------------------
results_long <- data.frame()

# =========================================================
# LOOP 1: Binning + XGBoost (kept separate from TwinGP)
# =========================================================
for (year in c(2017, 2018)) {
  test_ids <- if (year == 2017) turbine_ids else testset_2018
  
  for (i in test_ids) {
    cat("[Loop 1] Evaluating Turbine", i, "Year", year, "\n")
    
    test_data <- get_data(i, year)
    if (is.null(test_data)) next
    
    test_speed <- test_data$wind_speed
    test_temp  <- test_data$temperature
    test_power <- test_data$power
    
    X_test  <- cbind(test_speed, test_temp)
    Xs_test <- cbind(X_test, S_test_for(i, nrow(X_test)))
    
    tr <- get_train_loo_2017(i)
    X_train <- tr$X
    y_train <- tr$y
    S_train <- tr$S
    
    # --- Binning(speed only) ---
    set.seed(i)
    t1 <- Sys.time()
    bin_pred <- binning(train.x = X_train[, 1], train.y = y_train, test.x = test_speed)
    t2 <- Sys.time()
    
    rmse_bin <- sqrt(mean((bin_pred - test_power)^2, na.rm = TRUE))
    
    results_long <- rbind(results_long, data.frame(
      Method = "Binning", Turbine = i, Year = year,
      RMSE = rmse_bin, NLPD = NA_real_,
      Runtime = round(as.numeric(difftime(t2, t1, units = "secs")), 4)
    ))
    
  }
}

# =========================================================
# LOOP 2: TwinGP only (separate loop, no jitter)
# =========================================================
for (year in c(2017, 2018)) {
  test_ids <- if (year == 2017) turbine_ids else testset_2018
  
  for (i in test_ids) {
    cat("[Loop 2] Evaluating Turbine", i, "Year", year, "\n")
    
    test_data <- get_data(i, year)
    if (is.null(test_data)) next
    
    test_speed <- test_data$wind_speed
    test_temp  <- test_data$temperature
    test_power <- test_data$power
    
    X_test  <- cbind(test_speed, test_temp)
    Xs_test <- cbind(X_test, S_test_for(i, nrow(X_test)))
    
    tr <- get_train_loo_2017(i)
    X_train <- tr$X
    y_train <- tr$y
    S_train <- tr$S
    
    # --- TwinGP(x) ---
    set.seed(i)
    t1 <- Sys.time()
    twin_out <- twingp(x = X_train, y = y_train, x_test = X_test)
    t2 <- Sys.time()
    
    pred <- as.numeric(twin_out$mu)
    pred_sd <- as.numeric(twin_out$sigma)
    
    rmse <- sqrt(mean((pred - test_power)^2, na.rm = TRUE))
    nlpd <- mean(0.5 * log(2 * pi * pred_sd^2) +
                   0.5 * ((test_power - pred)^2) / (pred_sd^2),
                 na.rm = TRUE)
    
    results_long <- rbind(results_long, data.frame(
      Method = "TwinGP(x)", Turbine = i, Year = year,
      RMSE = rmse, NLPD = nlpd,
      Runtime = round(as.numeric(difftime(t2, t1, units = "secs")), 4)
    ))
    
    # --- TwinGP(x+s) ---
  #  set.seed(i)
  #  Xs_train <- cbind(X_train, S_train)
    
  #  t1 <- Sys.time()
  #  twin_out_s <- twingp(x = Xs_train, y = y_train, x_test = Xs_test)
  #  t2 <- Sys.time()
    
  #  pred_s <- as.numeric(twin_out_s$mu)
  #  pred_sd_s <- as.numeric(twin_out_s$sigma)
    
  #  rmse_s <- sqrt(mean((pred_s - test_power)^2, na.rm = TRUE))
  #  nlpd_s <- mean(0.5 * log(2 * pi * pred_sd_s^2) +
  #                   0.5 * ((test_power - pred_s)^2) / (pred_sd_s^2),
  #                 na.rm = TRUE)
    
  #  results_long <- rbind(results_long, data.frame(
  #    Method = "TwinGP(x+s)", Turbine = i, Year = year,
  #    RMSE = rmse_s, NLPD = nlpd_s,
  #    Runtime = round(as.numeric(difftime(t2, t1, units = "secs")), 4)
  #  ))
  }
}

# -------------------------
# Final summary table (mean over turbines)
# -------------------------
summary_table <- results_long %>%
  group_by(Method, Year) %>%
  summarise(
    RMSE = mean(RMSE, na.rm = TRUE),
    NLPD = if (all(is.na(NLPD))) NA_real_ else mean(NLPD, na.rm = TRUE),
    Runtime = mean(Runtime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Year,
    values_from = c(RMSE, NLPD, Runtime),
    names_sep = ""
  ) %>%
  transmute(
    Method,
    rmse2017 = RMSE2017, rmse2018 = RMSE2018,
    nlpd2017 = NLPD2017, nlpd2018 = NLPD2018,
    runtime2017 = Runtime2017, runtime2018 = Runtime2018
  ) %>%
  arrange(factor(Method, levels = c("Binning", "TwinGP(x)", "TwinGP(x+s)")))

print(summary_table)

write_csv(
  results_long,
  file = file.path(results_dir, "loo_results_all_methods_long.csv")
)

write_csv(
  summary_table,
  file = file.path(results_dir, "loo_results_all_methods_summary.csv")
)


# ============================================================
# UPDATE final results.csv : Table 2 and Table 3 rows
# ============================================================
source(file.path("code", "update_final_results.R"))

# choose which summary_table row should represent the paper row "TwinGP"
#twingp_label <- "TwinGP(x+s)"
twingp_label <- "TwinGP(x)"

get_val <- function(method_name, col_name) {
  v <- summary_table[summary_table$Method == method_name, col_name]
  if (length(v) == 0) return(NA_real_)
  as.numeric(v[1])
}

# -------------------------
# Table 2 : 2017
# -------------------------
update_final_results(
  method   = "Binning",
  table_id = "Table 2",
  rmse     = get_val("Binning", "rmse2017"),
  nlpd     = get_val("Binning", "nlpd2017"),
  runtime  = get_val("Binning", "runtime2017")
)

update_final_results(
  method   = "TwinGP",
  table_id = "Table 2",
  rmse     = get_val(twingp_label, "rmse2017"),
  nlpd     = get_val(twingp_label, "nlpd2017"),
  runtime  = get_val(twingp_label, "runtime2017")
)

# -------------------------
# Table 3 : 2018
# -------------------------
update_final_results(
  method   = "Binning",
  table_id = "Table 3",
  rmse     = get_val("Binning", "rmse2018"),
  nlpd     = get_val("Binning", "nlpd2018"),
  runtime  = get_val("Binning", "runtime2018")
)

update_final_results(
  method   = "TwinGP",
  table_id = "Table 3",
  rmse     = get_val(twingp_label, "rmse2018"),
  nlpd     = get_val(twingp_label, "nlpd2018"),
  runtime  = get_val(twingp_label, "runtime2018")
)