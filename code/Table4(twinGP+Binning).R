data_path    <- file.path("data")
terrain_path <- file.path("data", "weightedTerrainData.csv")
results_dir  <- file.path("results", "intermediate")

library(data.table)
library(twingp)
library(dplyr)
library(readr)
library(tidyr)

testset     <- c(38:44)
trainset    <- setdiff(1:66, testset)
turbine_ids <- 1:66

terrain_data <- read.csv(terrain_path)
scale_01 <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
terrain_data[, 2:4] <- lapply(terrain_data[, 2:4, drop=FALSE], scale_01)
terrain_mat <- as.matrix(terrain_data[1:66, 2:4])
storage.mode(terrain_mat) <- "double"

# ── Binning ────────────────────────────────────────────────────────────────────
binning <- function(train.x, train.y, test.x, bin_width = 0.5){
  train.y <- as.numeric(train.y); train.x <- as.numeric(train.x); test.x <- as.numeric(test.x)
  start <- 0; end <- round(max(train.x, na.rm=TRUE))
  n_bins <- round((end-start)/bin_width, 0) + 1
  x_bin <- y_bin <- numeric(n_bins)
  for (n in 2:n_bins){
    bin_element <- which(train.x > (start+(n-1)*bin_width) & train.x < (start+n*bin_width))
    x_bin[n] <- mean(train.x[bin_element], na.rm=TRUE)
    y_bin[n] <- mean(train.y[bin_element], na.rm=TRUE)
  }
  binned_data <- data.frame(x_bin, y_bin)
  binned_data <- binned_data[!is.na(binned_data$y_bin), ]
  splinefit <- smooth.spline(x=binned_data$x_bin, y=binned_data$y_bin, all.knots=TRUE)
  y_pred <- predict(splinefit, test.x)$y
  y_pred[y_pred < 0] <- 0
  y_pred
}

# ── Cache data ─────────────────────────────────────────────────────────────────
cache_key <- function(id, year) sprintf("T%02d_%d", id, year)
data_cache <- vector("list", length=length(turbine_ids)*2)
names(data_cache) <- as.vector(outer(sprintf("T%02d", turbine_ids), c("2017","2018"), paste, sep="_"))

for (id in turbine_ids) {
  for (yr in c(2017, 2018)) {
    f <- sprintf("%s/Turbine%d_%d.csv", data_path, id, yr)
    d <- tryCatch(fread(file=f, showProgress=FALSE), error=function(e) NULL)
    if (!is.null(d)) {
      d[, wind_speed  := as.numeric(wind_speed)]
      d[, temperature := as.numeric(temperature)]
      d[, power       := as.numeric(power)]
    }
    data_cache[[cache_key(id, yr)]] <- d
  }
}
get_data <- function(id, year) data_cache[[cache_key(id, year)]]

# ── Training pool ──────────────────────────────────────────────────────────────
X2017_list <- y2017_list <- tid2017_list <- vector("list", 66)
for (j in turbine_ids) {
  d <- get_data(j, 2017)
  if (!is.null(d)) {
    X2017_list[[j]]   <- cbind(d$wind_speed, d$temperature)
    y2017_list[[j]]   <- d$power
    tid2017_list[[j]] <- rep(j, nrow(d))
  } else {
    X2017_list[[j]]   <- matrix(numeric(0), ncol=2)
    y2017_list[[j]]   <- numeric(0)
    tid2017_list[[j]] <- integer(0)
  }
}
X2017   <- do.call(rbind, X2017_list)
y2017   <- as.numeric(unlist(y2017_list, use.names=FALSE))
tid2017 <- as.integer(unlist(tid2017_list, use.names=FALSE))
S2017   <- terrain_mat[tid2017, , drop=FALSE]

S_test_for <- function(i, n) matrix(rep(terrain_mat[i,], each=n), nrow=n, ncol=3)

get_train <- function() {
  idx <- !(tid2017 %in% testset)
  list(X=X2017[idx,,drop=FALSE], y=y2017[idx], tid=tid2017[idx], S=S2017[idx,,drop=FALSE])
}

results_long <- data.frame()

# ── LOOP 1: IEC Binning (x only — terrain not applicable) ─────────────────────
for (i in testset) {
  cat("[Binning] Turbine", i, "\n")
  test_data <- get_data(i, 2017); if (is.null(test_data)) next
  tr <- get_train()
  set.seed(i); t1 <- Sys.time()
  bin_pred <- binning(train.x=tr$X[,1], train.y=tr$y, test.x=test_data$wind_speed)
  t2 <- Sys.time()
  results_long <- rbind(results_long, data.frame(
    Method="IEC Binning", Turbine=i, Year=2017, version="x",
    RMSE=sqrt(mean((bin_pred-test_data$power)^2, na.rm=TRUE)), NLPD=NA_real_,
    Runtime=round(as.numeric(difftime(t2,t1,units="secs")),4)
  ))
}

# ── LOOP 2: TwinGP(x) ─────────────────────────────────────────────────────────
for (i in testset) {
  cat("[TwinGP x] Turbine", i, "\n")
  test_data <- get_data(i, 2017); if (is.null(test_data)) next
  X_test  <- cbind(test_data$wind_speed, test_data$temperature)
  tr <- get_train()
  set.seed(i); t1 <- Sys.time()
  twin_out <- twingp::twingp(x=tr$X, y=tr$y, x_test=X_test)
  t2 <- Sys.time()
  pred <- as.numeric(twin_out$mu); pred_sd <- as.numeric(twin_out$sigma)
  rmse <- sqrt(mean((pred-test_data$power)^2, na.rm=TRUE))
  nlpd <- mean(0.5*log(2*pi*pred_sd^2)+0.5*((test_data$power-pred)^2)/(pred_sd^2), na.rm=TRUE)
  results_long <- rbind(results_long, data.frame(
    Method="TwinGP", Turbine=i, Year=2017, version="x",
    RMSE=rmse, NLPD=nlpd,
    Runtime=round(as.numeric(difftime(t2,t1,units="secs")),4)
  ))
}

# ── LOOP 3: TwinGP(x+s) ───────────────────────────────────────────────────────
for (i in testset) {
  cat("[TwinGP x+s] Turbine", i, "\n")
  test_data <- get_data(i, 2017); if (is.null(test_data)) next
  X_test  <- cbind(test_data$wind_speed, test_data$temperature)
  Xs_test <- cbind(X_test, S_test_for(i, nrow(X_test)))
  tr <- get_train()
  Xs_train <- cbind(tr$X, tr$S)
  set.seed(i); t1 <- Sys.time()
  twin_out <- twingp::twingp(x=Xs_train, y=tr$y, x_test=Xs_test)
  t2 <- Sys.time()
  pred <- as.numeric(twin_out$mu); pred_sd <- as.numeric(twin_out$sigma)
  rmse <- sqrt(mean((pred-test_data$power)^2, na.rm=TRUE))
  nlpd <- mean(0.5*log(2*pi*pred_sd^2)+0.5*((test_data$power-pred)^2)/(pred_sd^2), na.rm=TRUE)
  results_long <- rbind(results_long, data.frame(
    Method="TwinGP", Turbine=i, Year=2017, version="xs",
    RMSE=rmse, NLPD=nlpd,
    Runtime=round(as.numeric(difftime(t2,t1,units="secs")),4)
  ))
}

# ── Summary & save ─────────────────────────────────────────────────────────────
summary_table <- results_long %>%
  group_by(Method, Year, version) %>%
  summarise(RMSE=mean(RMSE,na.rm=TRUE), NLPD=mean(NLPD,na.rm=TRUE),
            Runtime=mean(Runtime,na.rm=TRUE), .groups="drop")

print(summary_table)
write_csv(results_long,   file.path(results_dir, "table4_twinGP_binning_long.csv"))
write_csv(summary_table,  file.path(results_dir, "table4_twinGP_binning_summary.csv"))

# ── Update final results.csv ───────────────────────────────────────────────────
source(file.path("code", "update_final_results.R"))

get_rmse <- function(m, v) {
  r <- summary_table$RMSE[summary_table$Method==m & summary_table$version==v]
  if (length(r)==0) NA_real_ else r
}

update_final_results(method="IEC Binning", table_id="Table 4", rmse=get_rmse("IEC Binning","x"),  version="x")
update_final_results(method="TwinGP",      table_id="Table 4", rmse=get_rmse("TwinGP","x"),        version="x")
update_final_results(method="TwinGP",      table_id="Table 4", rmse=get_rmse("TwinGP","xs"),       version="xs")
