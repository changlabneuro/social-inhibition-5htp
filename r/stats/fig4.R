sesh_subdir <- '062921'
data_file <- 'data__.csv'
label_file <- 'labels__.csv'
data_type <- 'time_bw'
out_path <- file.path(export_root, 'out')

export_root <- export_path(file.path(sesh_subdir, data_type))
rt_data <- read_pair(file.path(export_root, data_file), file.path(export_root, label_file))

anova_res <- anova_test(rt_data, time_bw ~ drug, effect.size="pes")
write_csv(anova_res, out_path, 'anova.csv')