update_final_results <- function(method, table_id, rmse = NA, nlpd = NA, runtime = NA,
                                 version = "x",
                                 final_csv = file.path("results", "final results.csv")) {
  if (!file.exists(final_csv)) {
    stop("final results.csv does not exist. Please create the template first.")
  }
  
  df <- read.csv(final_csv, header = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  df[is.na(df)] <- ""
  
  find_row <- function(table_label, method_name) {
    table_row <- which(df[, 1] == table_label)
    if (length(table_row) == 0) stop(paste("Could not find", table_label))
    
    next_table_row <- which(df[, 1] %in% c("Table 2", "Table 3", "Table 4") & seq_len(nrow(df)) > table_row[1])
    end_row <- if (length(next_table_row) == 0) nrow(df) else next_table_row[1] - 1
    
    method_rows <- which(df[, 1] == method_name)
    method_rows <- method_rows[method_rows > table_row[1] & method_rows <= end_row]
    if (length(method_rows) == 0) stop(paste("Could not find method", method_name, "under", table_label))
    
    method_rows[1]
  }
  
  r <- find_row(table_id, method)
  
  if (table_id %in% c("Table 2", "Table 3")) {
    df[r, 2] <- ifelse(is.na(rmse),    "", sprintf("%.2f", rmse))
    df[r, 3] <- ifelse(is.na(nlpd),    "", sprintf("%.2f", nlpd))
    df[r, 4] <- ifelse(is.na(runtime), "", sprintf("%.2f", runtime))
  } else if (table_id == "Table 4") {
    if (version == "x")       df[r, 2] <- ifelse(is.na(rmse), "", sprintf("%.2f", rmse))
    else if (version == "xs") df[r, 3] <- ifelse(is.na(rmse), "", sprintf("%.2f", rmse))
    else stop("version must be 'x' or 'xs'")
  } else {
    stop("table_id must be one of: Table 2, Table 3, Table 4")
  }
  
  write.table(df, final_csv, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE)
}
