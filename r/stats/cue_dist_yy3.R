library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\cue_dist_yy3')
aov_mdl <- aov(cue_dist ~ correct * scrambled_type * target_image_category * monkey * drug, ps)
tukey_res <- tukey_hsd(aov_mdl)
e2 <- eta_squared(aov_mdl)