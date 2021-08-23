sesh_subdir <- '071421'
data_file <- 'data__.csv'
label_file <- 'labels__.csv'
data_type <- 'polar_landing_points'

export_root <- export_path(file.path(sesh_subdir, data_type))
data_path <- file.path(export_root, data_file)
label_path <- file.path(export_root, label_file)
out_path <- file.path(export_root, 'out')

combined <- read_pair(data_path, label_path)

each <- unique(combined$correct)
f_list <- list()
i <- 1

for (corr in each) {
  match <- combined$correct == corr
  subset <- combined[match,]
  stat <- aov.circular(circular(subset$theta), subset$drug)
  f_list[[i]] <- stat
  i = i + 1
}