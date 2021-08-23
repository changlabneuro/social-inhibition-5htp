library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\norm_pcorr')
aov_mdl <- aov(pcorr ~ scrambled_type*target_image_category, ps)
tukey_res <- tukey_hsd(aov_mdl)
e2 <- eta_squared(aov_mdl)