testset=c(1:66)
testset1=c(1:46,48:50,52,54:60,62:66)


data_dir <- file.path("data")
data_path <- data_dir
input_folder <- file.path("data", "processed data")
terrain_data_path <- file.path("data", "weightedTerrainData.csv")
output_folder <- file.path("results","intermediate")



stopifnot(dir.exists(input_folder))

read_vec_csv <- function(path, colname = NULL) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(colname)) {
    if (!colname %in% names(df)) stop("Column '", colname, "' not found in: ", path)
    v <- df[[colname]]
  } else {
    if (ncol(df) != 1) stop("Expected 1 column in: ", path, " but found ", ncol(df))
    v <- df[[1]]
  }
  if (is.character(v)) suppressWarnings(vn <- as.numeric(v)) else vn <- v
  if (is.numeric(vn) && !all(is.na(vn))) return(vn)
  return(v)
}

read_mat_csv <- function(path) {
  df_try <- try(read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  
  if (!inherits(df_try, "try-error")) {
    df <- df_try
    rn <- rownames(df)
    if (any(is.na(rn)) || any(rn == "") || any(duplicated(rn))) {
      df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  } else {
    df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
  df_num <- df
  for (j in seq_along(df_num)) {
    suppressWarnings(df_num[[j]] <- as.numeric(df_num[[j]]))
  }
  
  na_ratio <- function(x) mean(is.na(x))
  before <- mean(vapply(df, na_ratio, numeric(1)))
  after  <- mean(vapply(df_num, na_ratio, numeric(1)))
  if (after <= before + 1e-12) df <- df_num
  m <- as.matrix(df)
  return(m)
}

temp_vector <- read_vec_csv(file.path(input_folder, "temp_vector.csv"), colname = "temp")

speed_vector <- read_vec_csv(file.path(input_folder, "speed_vector.csv"), colname = "speed")

power_matrix <- read_mat_csv(file.path(input_folder, "power_matrix.csv"))

sd_matrix <- read_mat_csv(file.path(input_folder, "sd_matrix.csv"))

sigma2_matrix <- read_mat_csv(file.path(input_folder, "sigma2_matrix.csv"))


library(dplyr)

turbine_files <- list.files(data_path, pattern = "Turbine[1-9]{1}_2017.csv|Turbine[1-6][0-9]_2017.csv", full.names = TRUE)
set.seed(15)
turbine_data_list <- list()

for (i in seq_along(turbine_files)) {
  dataset <- read.csv(turbine_files[i])
  turbine_data_list[[i]] <- dataset
}

all_speeds <- c()
for (i in seq_along(turbine_data_list)) {
  speed <- turbine_data_list[[i]]$wind_speed
  all_speeds <- c(all_speeds, speed)  }

all_temp <- c()
for (i in seq_along(turbine_data_list)) {
  temp <- turbine_data_list[[i]]$temperature
  all_temp <- c(all_temp, temp)  }


terrain_data <- read.csv(terrain_data_path)
scale_01 <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

scaled_terrain_data <- as.data.frame(lapply(terrain_data, scale_01))

##############

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)


rmse_list <- c(); rmse_list_18 <- c()
sigma2_list <- c(); sigma2_list_18 <- c()
tau2_list <- c(); tau2_list_18 <- c()
nlpd_list <- c(); nlpd_list_18 <- c()


runtimes_total <- c()   # total per i
runtimes_fit <- c()     # Stage A+B+C per i
runtimes_pred_17 <- c() # prediction time 2017 per i
runtimes_pred_18 <- c() # prediction time 2018 per i (only for i in testset1)

pred_list <- list(); pred_list_18 <- list()
pred_sd_list <- list(); pred_sd_list_18 <- list()



scale_01_speed <- function(x) (x - min(all_speeds)) / (max(all_speeds) - min(all_speeds))
scale_01_temp  <- function(x) (x - min(all_temp))   / (max(all_temp)   - min(all_temp))

closest_terrain_id <- function(i, scaled_terrain_data) {
  Xall <- as.matrix(scaled_terrain_data[, 2:4])
  u0 <- as.numeric(Xall[i, ])
  d2 <- rowSums((Xall - matrix(u0, nrow=nrow(Xall), ncol=3, byrow=TRUE))^2)
  d2[i] <- Inf
  which.min(d2)
}

dbg <- TRUE
eps_jit <- 1e-8
eps_inv <- 1e-8

# ---- OUTSIDE LOOP (once) ----
z_raw  <- scale_01_speed(speed_vector)
z1_raw <- scale_01_temp(temp_vector)

ord_global <- order(z_raw)
z_sorted  <- z_raw[ord_global]
z1_sorted <- z1_raw[ord_global]

Ez_global  <- as.matrix(dist(z_sorted,  diag=TRUE, upper=TRUE))
Ez1_global <- as.matrix(dist(z1_sorted, diag=TRUE, upper=TRUE))
m_global <- length(z_sorted)




runtimes_total   <- c()
runtimes_fit     <- c()
runtimes_pred_17 <- c()
runtimes_pred_18 <- c()

rmse_list_ok    <- c()   # NEW: OK RMSE 2017
rmse_list_ok_18 <- c()   # NEW: OK RMSE 2018 (only i in testset1)

pred_ok_list    <- list()  # NEW: OK preds 2017
pred_ok_list_18 <- list()  # NEW: OK preds 2018

runtimes_total   <- c()   # total per i
runtimes_fit     <- c()   # ML fitting per i
runtimes_pred_17 <- c()   # prediction time 2017 per i
runtimes_pred_18 <- c()   # prediction time 2018 per i (only for i in testset1)



for (i in testset) {
  
  t_total0 <- Sys.time()
  set.seed(i)

  scale_01_test <- function(x) (x - min(all_speeds)) / (max(all_speeds) - min(all_speeds))
  scale_01_test_temperature <- function(x) (x - min(all_temp)) / (max(all_temp) - min(all_temp))

  test_data_path <- file.path(
    data_path,
    paste0("Turbine", i, "_2017.csv")
  )
  test_data <- read.csv(test_data_path)
  
  test_data$scaled_wind_speed  <- scale_01_test(test_data$wind_speed)
  test_data$scaled_temperature <- scale_01_test_temperature(test_data$temperature)
  
  z.ev  <- test_data$scaled_wind_speed
  z.ev1 <- test_data$scaled_temperature
  f.ev  <- test_data$power

  do_2018 <- exists("testset1") && (i %in% testset1)
  if (do_2018) {
    test_data_path1 <- file.path(
      data_path,
      paste0("Turbine", i, "_2018.csv")
    )
    test_data1 <- read.csv(test_data_path1)
    
    test_data1$scaled_wind_speed  <- scale_01_test(test_data1$wind_speed)
    test_data1$scaled_temperature <- scale_01_test_temperature(test_data1$temperature)
    
    z.ev_18  <- test_data1$scaled_wind_speed
    z.ev1_18 <- test_data1$scaled_temperature
    f.ev_18  <- test_data1$power
  }
  
  # ----------------------------
  # training objects (-i)
  # ----------------------------
  x <- scaled_terrain_data[, 2:4]
  x <- x[-i, , drop = FALSE]
  
  z  <- scale_01_test(speed_vector)
  z1 <- scale_01_test_temperature(temp_vector)
  
  ord <- order(z)
  z  <- z[ord]
  z1 <- z1[ord]
  
  y <- power_matrix[, ord, drop = FALSE]
  y <- y[-i, , drop = FALSE]
  
  n <- nrow(x)
  m <- length(z)
  N <- n * m
  vec.y <- c(t(y))
  
  Ex1 <- as.matrix(dist(x[, 1], diag = TRUE, upper = TRUE))
  Ex2 <- as.matrix(dist(x[, 2], diag = TRUE, upper = TRUE))
  Ex3 <- as.matrix(dist(x[, 3], diag = TRUE, upper = TRUE))
  Ez  <- as.matrix(dist(z,  diag = TRUE, upper = TRUE))
  Ez1 <- as.matrix(dist(z1, diag = TRUE, upper = TRUE))
  
  # ----------------------------
  # nugget matrix + scalar nugget per test point (delta-style)
  # ----------------------------
  nug_all <- sigma2_matrix[, ord, drop = FALSE]
  nug <- nug_all[-i, , drop = FALSE]  # match training (-i)
  
  # nearest column in z for each z.ev (same as before)
  nearest_col <- vapply(z.ev, function(zz) which.min(abs(z - zz)), integer(1))
  
  # u is your target terrain vector (length 3 here)
  u <- c(scaled_terrain_data[i, 2], scaled_terrain_data[i, 3], scaled_terrain_data[i, 4])
  
  # find the single row index z0 whose x[row,] is closest to u
  z0 <- which.min(rowSums((x - matrix(u, nrow=nrow(x), ncol=ncol(x), byrow=TRUE))^2))
  # (that’s squared Euclidean distance; sqrt not needed for argmin)
  
  # now pick nugget from that row only, for each chosen column k
  nug_scalar <- vapply(nearest_col, function(k) nug[z0, k], numeric(1))
  
  if (do_2018) {
    nearest_col_18 <- vapply(z.ev_18, function(zz) which.min(abs(z - zz)), integer(1))
    z0 <- which.min(rowSums((x - matrix(u, nrow=nrow(x), ncol=ncol(x), byrow=TRUE))^2))
    nug_scalar_18 <- vapply(nearest_col_18, function(k) nug[z0, k], numeric(1))
  }
  
  # ----------------------------
  # ML objective (same as yours)
  # ----------------------------
  ML <- function(para) {
    theta <- para
    eps <- 1e-10
    
    Rx <- 1 / (1 + (Ex1 / theta[1])^2 + (Ex2 / theta[2])^2 + (Ex3 / theta[3])^2)
    Rx <- Rx + eps * diag(n)
    Rxinv <- solve(Rx)
    
    Rz <- exp(-Ez / theta[4] - Ez1 / theta[5])
    Rz <- Rz + eps * diag(m)
    Rzinv <- solve(Rz)
    
    a_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
    b_mat <- Rxinv %*% as.matrix(y) %*% Rzinv
    a <- c(t(a_mat))
    b <- c(t(b_mat))
    
    mu <- sum(b) / sum(a)
    sigma2 <- (1/(m*n)) * sum((vec.y - mu) * (b - mu * a))
    
    logdetRx <- determinant(Rx, logarithm = TRUE)$modulus[1]
    logdetRz <- determinant(Rz, logarithm = TRUE)$modulus[1]
    val <- m * n * log(sigma2) + m * logdetRx + n * logdetRz
    
    # gradients
    grad <- numeric(5)
    denom_a <- sum(a)
    denom_b <- sum(b)
    dmu <- function(da, db) (sum(db) * denom_a - denom_b * sum(da)) / denom_a^2
    
    for (k in 1:3) {
      Exk <- switch(k, Ex1, Ex2, Ex3)
      dRx  <- 2 * Rx^2 * (Exk^2) / theta[k]^3
      dKx  <- -Rxinv %*% dRx %*% Rxinv
      
      da_mat <- dKx %*% rep(1, n) %*% t(rep(1, m)) %*% Rzinv
      db_mat <- dKx %*% y %*% Rzinv
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      
      d_sigma2 <- mean(
        -d_mu * a * (b - mu * a)
        - mu  * da * (b - mu * a)
        + (vec.y - mu * a) * (db - d_mu * a - mu * da)
      )
      grad[k] <- (m * n / sigma2) * d_sigma2 + m * sum(Rxinv * dRx)
    }
    
    for (k in 4:5) {
      Ez_k <- if (k == 4) Ez else Ez1
      dRz  <- Rz * Ez_k / theta[k]^2
      dKz  <- -Rzinv %*% dRz %*% Rzinv
      
      da_mat <- Rxinv %*% rep(1, n) %*% t(rep(1, m)) %*% dKz
      db_mat <- Rxinv %*% y %*% dKz
      da <- c(t(da_mat)); db <- c(t(db_mat))
      d_mu <- dmu(da, db)
      
      d_sigma2 <- mean(
        -d_mu * a * (b - mu * a)
        - mu  * da * (b - mu * a)
        + (vec.y - mu * a) * (db - d_mu * a - mu * da)
      )
      grad[k] <- (m * n / sigma2) * d_sigma2 + n * sum(Rzinv * dRz)
    }
    
    attr(val, "gradient") <- grad
    val
  }
  
  ML_log <- function(phi) {
    theta <- exp(phi)
    out   <- ML(theta)
    g     <- attr(out, "gradient") * theta
    attr(out, "gradient") <- g
    out
  }
  
  phi0 <- log(c(0.1, 0.1, 0.1, 1, 1))
  
  # ----------------------------
  # FIT timing
  # ----------------------------
  t_fit0 <- Sys.time()
  
  fit <- nlminb(
    start     = phi0,
    objective = function(p) ML_log(p),
    gradient  = function(p) attr(ML_log(p), "gradient"),
    control   = list(iter.max = 200, rel.tol = 1e-8)
  )
  theta <- exp(fit$par)
  
  t_fit1 <- Sys.time()
  runtimes_fit <- c(runtimes_fit, as.numeric(difftime(t_fit1, t_fit0, units = "secs")))
  cat("Turbine", i, " nlminb convergence code:", fit$convergence, "\n")
  
  # ----------------------------
  # build Rx,Rz,a,b,mu,sigma2
  # ----------------------------
  Rx    <- 1/(1 + (Ex1/theta[1])^2 + (Ex2/theta[2])^2 + (Ex3/theta[3])^2)
  Rxinv <- solve(Rx + 1e-8 * diag(n))
  
  Rz    <- exp(-Ez/theta[4] - Ez1/theta[5])
  Rzinv <- solve(Rz + 1e-8 * diag(m))
  
  a <- c(t((Rxinv %*% rep(1, n)) %*% (t(rep(1, m)) %*% Rzinv)))
  b <- c(t(Rxinv %*% as.matrix(y) %*% Rzinv))
  
  mu     <- sum(b) / sum(a)
  sigma2 <- (1/(m*n)) * sum((vec.y - mu) * (b - mu * a))
  
  # LK matrices
  A <- matrix(a, nrow = n, ncol = m, byrow = TRUE)
  B <- matrix(b, nrow = n, ncol = m, byrow = TRUE)
  
  # ----------------------------
  # OK COEF for your requested OK formula:
  # yhat_ok = mu + t(rx_vec %*% COEF %*% rz_vec)
  # ----------------------------
  # IMPORTANT: define COEF consistently with byrow ordering:
  COEF <- Rxinv %*% (y - mu) %*% Rzinv
  
  # ----------------------------
  # correlation vectors
  # ----------------------------
  basis.x <- function(h) 1/(1 + sum((h/theta[1:3])^2))
  basis.z <- function(h, h1) exp(-abs(h)/theta[4] - abs(h1)/theta[5])
  
  r.x <- function(u) {
    Auu <- t(t(x) - u)
    apply(Auu, 1, basis.x)
  }
  r.z <- function(v, v1) basis.z(v - z, v1 - z1)
  

  # ==========================================================
  # Prediction 2017 (timed): LK + OK (your requested form)
  # ==========================================================
  t_p170 <- Sys.time()
  
  yhat_lk <- numeric(length(z.ev))
  yhat_ok <- numeric(length(z.ev))
  yhat_sd <- numeric(length(z.ev))
  
  rx_vec <- r.x(u)  # compute ONCE (u fixed)
  
  for (j in seq_along(z.ev)) {
    rz_vec <- r.z(z.ev[j], z.ev1[j])
    
    # LK mean
    denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec)
    num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec)
    yhat_lk[j] <- num / denom
    
    # OK mean (as you requested)
    yhat_ok[j] <- as.numeric(mu + t(rx_vec %*% COEF %*% rz_vec))
    
    # Joseph (2005) MSPE_LK
    alpha <- as.numeric(t(rx_vec) %*% Rxinv %*% rx_vec)
    beta  <- as.numeric(t(rz_vec) %*% Rzinv %*% rz_vec)
    alphaBeta <- alpha * beta
    
    mspe_pred <- sigma2 * (1 - alphaBeta + alphaBeta * (1 - denom)^2 / denom^2)
    
    # delta-style nugget add-on (keep this; remove "+ nug_scalar[j]" if you want latent only)
    yhat_sd[j] <- sqrt(pmax(mspe_pred + nug_scalar[j], 1e-12))
  }
  
  t_p171 <- Sys.time()
  runtimes_pred_17 <- c(runtimes_pred_17, as.numeric(difftime(t_p171, t_p170, units = "secs")))
  
  # ==========================================================
  # Prediction 2018 (if requested, timed): LK + OK
  # ==========================================================
  if (do_2018) {
    
    t_p180 <- Sys.time()
    
    yhat_lk_18 <- numeric(length(z.ev_18))
    yhat_ok_18 <- numeric(length(z.ev_18))
    yhat_sd_18 <- numeric(length(z.ev_18))
    
    # rx_vec same (u fixed)
    for (j in seq_along(z.ev_18)) {
      rz_vec_18 <- r.z(z.ev_18[j], z.ev1_18[j])
      
      denom <- as.numeric(t(rx_vec) %*% A %*% rz_vec_18)
      num   <- as.numeric(t(rx_vec) %*% B %*% rz_vec_18)
      yhat_lk_18[j] <- num / denom
      
      yhat_ok_18[j] <- as.numeric(mu + t(rx_vec %*% COEF %*% rz_vec_18))
      
      alpha <- as.numeric(t(rx_vec) %*% Rxinv %*% rx_vec)
      beta  <- as.numeric(t(rz_vec_18) %*% Rzinv %*% rz_vec_18)
      alphaBeta <- alpha * beta
      
      mspe_pred_18 <- sigma2 * (1 - alphaBeta + alphaBeta * (1 - denom)^2 / denom^2)
      yhat_sd_18[j] <- sqrt(pmax(mspe_pred_18 + nug_scalar_18[j], 1e-12))
    }
    
    t_p181 <- Sys.time()
    runtimes_pred_18 <- c(runtimes_pred_18, as.numeric(difftime(t_p181, t_p180, units = "secs")))
    
    # store 2018 preds (if you have these lists)
    if (exists("pred_list_18"))    pred_list_18[[as.character(i)]]    <- yhat_lk_18
    if (exists("pred_ok_list_18")) pred_ok_list_18[[as.character(i)]] <- yhat_ok_18
    if (exists("pred_sd_list_18")) pred_sd_list_18[[as.character(i)]] <- yhat_sd_18
    
    # RMSE 2018 (if you keep these vectors)
    if (exists("rmse_list_18")) {
      rmse_list_18    <- c(rmse_list_18,    sqrt(mean((f.ev_18 - yhat_lk_18)^2, na.rm = TRUE)))
      rmse_list_ok_18 <- c(rmse_list_ok_18, sqrt(mean((f.ev_18 - yhat_ok_18)^2, na.rm = TRUE)))
      sigma2_list_18  <- c(sigma2_list_18, sigma2)
    }
  }
  
  # ----------------------------
  # Save per turbine predictions (2017)
  # ----------------------------
  if (exists("pred_list"))    pred_list[[as.character(i)]]    <- yhat_lk
  if (exists("pred_ok_list")) pred_ok_list[[as.character(i)]] <- yhat_ok
  if (exists("pred_sd_list")) pred_sd_list[[as.character(i)]] <- yhat_sd
  
  # ----------------------------
  # RMSEs (2017)
  # ----------------------------
  rmse_value    <- sqrt(mean((f.ev - yhat_lk)^2, na.rm = TRUE))
  rmse_value_ok <- sqrt(mean((f.ev - yhat_ok)^2, na.rm = TRUE))
  
  rmse_list    <- c(rmse_list, rmse_value)
  rmse_list_ok <- c(rmse_list_ok, rmse_value_ok)
  sigma2_list  <- c(sigma2_list, sigma2)
  
  print(rmse_list)
  print(rmse_list_ok)
  
  # ----------------------------
  # Total timing
  # ----------------------------
  t_total1 <- Sys.time()
  runtimes_total <- c(runtimes_total, as.numeric(difftime(t_total1, t_total0, units = "secs")))
}



nlpd_list    <- c()
nlpd_list_18 <- c()

for (i in testset) {
  
  # ---------- 2017 ----------
  test_data_path <- paste0(data_dir, "/Turbine", i, "_2017.csv")
  test_data <- read.csv(test_data_path)
  
  f.ev <- test_data$power
  
  yhat    <- pred_list[[as.character(i)]]
  yhat_sd <- pred_sd_list[[as.character(i)]]
  
  # safety
  stopifnot(length(yhat) == length(f.ev), length(yhat_sd) == length(f.ev))
  
  yhat_sd <- pmax(as.numeric(yhat_sd), 1e-12)
  
  nlpds <- 0.5 * log(2*pi*yhat_sd^2) + 0.5 * ((f.ev - as.numeric(yhat))^2) / (yhat_sd^2)
  nlpd  <- mean(nlpds, na.rm = TRUE)
  
  nlpd_list <- c(nlpd_list, nlpd)
  

  if (exists("testset1") && (i %in% testset1)) {
    
    # if you didn't store preds for this turbine, skip safely
    if (!is.null(pred_list_18[[as.character(i)]]) && !is.null(pred_sd_list_18[[as.character(i)]])) {
      
      test_data_path1 <- paste0(data_dir, "/Turbine", i, "_2018.csv")
      test_data1 <- read.csv(test_data_path1)
      
      f.ev_18 <- test_data1$power
      
      yhat_18    <- pred_list_18[[as.character(i)]]
      yhat_sd_18 <- pred_sd_list_18[[as.character(i)]]
      
      stopifnot(length(yhat_18) == length(f.ev_18), length(yhat_sd_18) == length(f.ev_18))
      
      yhat_sd_18 <- pmax(as.numeric(yhat_sd_18), 1e-12)
      
      nlpds18 <- 0.5 * log(2*pi*yhat_sd_18^2) + 0.5 * ((f.ev_18 - as.numeric(yhat_18))^2) / (yhat_sd_18^2)
      nlpd_18 <- mean(nlpds18, na.rm = TRUE)
      
      nlpd_list_18 <- c(nlpd_list_18, nlpd_18)
    }
  }
}

cat("Avg NLPD 2017:", mean(nlpd_list, na.rm = TRUE), "\n")
if (length(nlpd_list_18) > 0) cat("Avg NLPD 2018:", mean(nlpd_list_18, na.rm = TRUE), "\n")
# ============================================================
# WRITE RESULTS
# ============================================================
results_df <- data.frame(
  Turbine_ID    = testset,
  RMSE          = rmse_list,
  RMSE_OK_2017  = rmse_list_ok,
  Sigma2        = sigma2_list,
#  Tau2          = tau2_list,
  NLPD          = nlpd_list,
  Runtime_total = runtimes_total,
  Runtime_fit   = runtimes_fit,
  Runtime_pred17= runtimes_pred_17
)

write.csv(results_df, paste0(output_folder, "/table2-table3.csv"), row.names=FALSE)

if (length(testset1) > 0) {
  results_df_18 <- data.frame(
    Turbine_ID     = testset1,
    RMSE_2018      = rmse_list_18,
    Sigma2         = sigma2_list_18,
  #  Tau2           = tau2_list_18,
    NLPD_2018      = nlpd_list_18,
    Runtime_pred18 = runtimes_pred_18
  )
  write.csv(results_df_18, paste0(output_folder, "/table2_table3.csv"), row.names=FALSE)
}




# ============================================================
# UPDATE final results.csv : Table 2 and Table 3, STGP row
# ============================================================
source(file.path("code", "update_final_results.R"))

# ---- Table 2 (2017) ----
stgp_table2_rmse <- if (nrow(results_df) > 0) mean(results_df$RMSE, na.rm = TRUE) else NA_real_
stgp_table2_nlpd <- if (nrow(results_df) > 0) mean(results_df$NLPD, na.rm = TRUE) else NA_real_
stgp_table2_runtime <- if (nrow(results_df) > 0) mean(results_df$Runtime_total, na.rm = TRUE) else NA_real_

update_final_results(
  method   = "STGP (ours)",
  table_id = "Table 2",
  rmse     = stgp_table2_rmse,
  nlpd     = stgp_table2_nlpd,
  runtime  = stgp_table2_runtime
)

# ---- Table 3 (2018) ----
if (exists("results_df_18") && nrow(results_df_18) > 0) {
  # For 2018, total runtime = fit time from 2017 + pred18 time
  runtime18_total <- results_df_18$Runtime_pred18
  
  if ("Runtime_fit" %in% names(results_df)) {
    idx_match <- match(results_df_18$Turbine_ID, results_df$Turbine_ID)
    runtime18_total <- results_df$Runtime_fit[idx_match] + results_df_18$Runtime_pred18
  }
  
  stgp_table3_rmse <- mean(results_df_18$RMSE_2018, na.rm = TRUE)
  stgp_table3_nlpd <- mean(results_df_18$NLPD_2018, na.rm = TRUE)
  stgp_table3_runtime <- mean(runtime18_total, na.rm = TRUE)
  
  update_final_results(
    method   = "STGP (ours)",
    table_id = "Table 3",
    rmse     = stgp_table3_rmse,
    nlpd     = stgp_table3_nlpd,
    runtime  = stgp_table3_runtime
  )
}