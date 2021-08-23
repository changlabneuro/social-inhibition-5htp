sesh_subdir <- '062921'
data_file <- 'data__.csv'
label_file <- 'labels__.csv'
data_type <- 'pcorr'

export_root <- export_path(file.path(sesh_subdir, data_type))
data_path <- file.path(export_root, data_file)
label_path <- file.path(export_root, label_file)
out_path <- file.path(export_root, 'out')

combined <- read_pair(data_path, label_path)

pcorr_res <- anova_test(combined, pcorr ~ scrambled_type*trial_type, effect.size="pes")
write_csv(pcorr_res, out_path, 'anova.csv')