library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\trial_sequence_pcorr_nn')
aov_mdl <- aov(pcorr ~ trial_type*prev_trial_type, ps)
tukey_res <- tukey_hsd(aov_mdl)
e2 <- eta_squared(aov_mdl)