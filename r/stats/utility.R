read_pair <- function(data_path, label_path) {
  data_csv <- read.csv(data_path)
  label_csv <- read.csv(label_path)
  combined <- cbind(data.frame(data_csv), data.frame(label_csv))
  return(combined)
}

read_pair_dir <- function(dir_path) {
  data_path <- file.path(dir_path, 'data__.csv')
  label_path <- file.path(dir_path, 'labels__.csv')
  return(read_pair(data_path, label_path))
}

require_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

write_csv <- function(data, directory, fname) {
  require_dir(directory)
  write.csv(data, file = file.path(directory, fname))
}

data_root <- function() {
  return('D:\\data\\hwwa')
}

code_root <- function() {
  return('C:\\Users\\nick\\source\\matlab\\hwwa\\r\\stats')
}

export_path <- function(join_with) {
  return(file.path(data_root(), 'export/approach_avoid/behavior', join_with))
}