library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\pcorr_jj')
aov_mdl <- aov(pcorr ~ target_image_category*trial_type, ps)
tukey_res <- tukey_hsd(aov_mdl)
e2 <- eta_squared(aov_mdl)