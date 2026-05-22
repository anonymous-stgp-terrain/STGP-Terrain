# Table4(Binning-hetGP).R
# Runs both x-only and x+s versions of the hetGP model for Table 4.

data_path    <- file.path("data")
terrain_path <- file.path("data", "weightedTerrainData.csv")
results_dir  <- file.path("results", "intermediate")

library(dplyr)
library(hetGP)
library(readr)

testset  <- 38:44
trainset <- setdiff(1:66, testset)

# ── Terrain ────────────────────────────────────────────────────────────────────
scale_01 <- function(x) (x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
terrain_data <- read.csv(terrain_path)
scaled_terrain <- as.data.frame(lapply(terrain_data, scale_01))
terrain_features <- scaled_terrain[, c("weighted_slope","weighted_rix","weighted_ridge")]
terrain_features$turbine_id <- seq_len(nrow(terrain_features))

# ── Load turbine files ─────────────────────────────────────────────────────────
turbine_files <- list.files(data_path,
  pattern = "Turbine[0-9]+_2017\\.csv", full.names=TRUE)
turbine_files <- turbine_files[order(as.integer(
  sub(".*Turbine([0-9]+)_2017\\.csv","\\1", basename(turbine_files))))]

# ── Build training data ────────────────────────────────────────────────────────
train_df <- dplyr::bind_rows(lapply(trainset, function(i) {
  df <- read.csv(turbine_files[i])[, c("wind_speed","temperature","power")]
  df$turbine_id <- i
  df
})) %>%
  left_join(terrain_features, by="turbine_id") %>%
  select(-turbine_id)

cat("Train rows:", nrow(train_df), "\n")

best_n        <- 3.0
train_rounded <- train_df
train_rounded[, 1:2] <- round(train_rounded[, 1:2] / best_n) * best_n

run_hetgp <- function(X_train, Z_train, X_test, y_test, label) {
  cat("\nFitting hetGP:", label, "\n")
  t0 <- proc.time()
  model <- mleHetGP(X=X_train, Z=Z_train,
                    lower=rep(0.001, ncol(X_train)),
                    upper=rep(1000,  ncol(X_train)))
  train_time <- (proc.time()-t0)["elapsed"]
  cat("Training time:", round(train_time,2), "sec\n")

  results <- data.frame(Turbine=integer(), RMSE=numeric(), Test_Time_s=numeric())
  for (tid in testset) {
    test_df <- read.csv(turbine_files[tid])[, c("wind_speed","temperature","power")]
    test_df$turbine_id <- tid
    test_df <- left_join(test_df, terrain_features, by="turbine_id") %>% select(-turbine_id)

    if (ncol(X_train) == 2) {
      X_te <- as.matrix(test_df[, c("wind_speed","temperature")])
    } else {
      X_te <- as.matrix(test_df[, c("wind_speed","temperature",
                                     "weighted_slope","weighted_rix","weighted_ridge")])
    }
    y_te <- test_df$power

    t1 <- proc.time()
    preds <- predict(model, X_te)$mean
    test_time <- (proc.time()-t1)["elapsed"]

    rmse <- sqrt(mean((preds-y_te)^2))
    results <- rbind(results, data.frame(
      Turbine=tid, RMSE=round(rmse,4), Test_Time_s=round(test_time,2)
    ))
    cat("Turbine", tid, "| RMSE:", round(rmse,4), "\n")
  }
  list(results=results, train_time=train_time)
}

# ── x-only ─────────────────────────────────────────────────────────────────────
X_train_x <- as.matrix(train_rounded[, c("wind_speed","temperature")])
Z_train   <- train_rounded$power

out_x <- run_hetgp(X_train_x, Z_train, NULL, NULL, "x-only")
avg_rmse_x <- mean(out_x$results$RMSE)
cat("\nBinning-hetGP(x)  avg RMSE:", round(avg_rmse_x,4), "\n")

# ── x+s ────────────────────────────────────────────────────────────────────────
X_train_xs <- as.matrix(train_rounded[, c("wind_speed","temperature",
                                           "weighted_slope","weighted_rix","weighted_ridge")])

out_xs <- run_hetgp(X_train_xs, Z_train, NULL, NULL, "x+s")
avg_rmse_xs <- mean(out_xs$results$RMSE)
cat("\nBinning-hetGP(x+s) avg RMSE:", round(avg_rmse_xs,4), "\n")

# ── Save ───────────────────────────────────────────────────────────────────────
write_csv(out_x$results,  file.path(results_dir, "table4_hetgp_x_results.csv"))
write_csv(out_xs$results, file.path(results_dir, "table4_hetgp_xs_results.csv"))

# ── Update final results.csv ───────────────────────────────────────────────────
source(file.path("code", "update_final_results.R"))
update_final_results(method="Binning-hetGP", table_id="Table 4", rmse=avg_rmse_x,  version="x")
update_final_results(method="Binning-hetGP", table_id="Table 4", rmse=avg_rmse_xs, version="xs")
