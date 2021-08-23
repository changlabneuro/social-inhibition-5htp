library(R.matlab)
library(FactoMineR)

combine_data_labels <- function(data, labels) {
  all_data <- data
  
  for (name in names(labels)) {
    all_data$name = labels$name
  }
  
  return(all_data)
}

raw_data_to_frame <- function(raw_data) {
  mat_data <- data.frame(raw_data$data)
  names(mat_data) <- unlist(raw_data$data.header)
  return(mat_data)
}

extract_labels <- function(raw_data) {
  categories <- unlist(raw_data$label.categories)
  entries <- unlist(raw_data$label.entries)
  indices <- raw_data$label.indices
  
  label_entries <- matrix(0, nrow(indices), ncol(indices))
  
  for (i in 1:ncol(label_entries)) {
    label_entries[, i] = entries[indices[, i]]
  }
  
  labels <- data.frame(label_entries)
  names(labels) <- categories
  
  return(labels)
}

raw_data <- readMat("/Users/Nick/Desktop/hwwa/pca_data.mat")
labels <- extract_labels(raw_data)
data <- raw_data_to_frame(raw_data)
all_data <- cbind(data, labels)

keep_categories <- c("drug", "trial_type")
keep_columns <- c(names(data), keep_categories)

mfa_data <- all_data[keep_columns]
any_nan = apply(apply(mfa_data[names(data)], 2, is.nan), 1, any)

res.mfa <- MFA(mfa_data[!any_nan,],
               group = c(ncol(data), length(keep_categories)),
               type = c("s", "n"),
               name.group = c("saccade_info", "behavior"),
               num.group.sup = 1,
               graph = FALSE)

quanti_var <- get_mfa_var(res.mfa)
