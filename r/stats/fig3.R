sesh_subdir <- '062921'
data_file <- 'data__.csv'
label_file <- 'labels__.csv'
data_type <- 'trial_sequence_pcorr'
out_path <- file.path(export_root, 'out')

export_root <- export_path(file.path(sesh_subdir, data_type))
pcorr_data <- read_pair(file.path(export_root, data_file), file.path(export_root, label_file))

anova_res <- anova_test(pcorr_data, pcorr ~ prev_correct_sequence*trial_type, effect.size="pes")
write_csv(anova_res, out_path, 'anova.csv')