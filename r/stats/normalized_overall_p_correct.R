library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081621\\norm_pcorr')
#mdl <- anova_test(ps, pcorr ~ drug*trial_type*scrambled_type)
aov_mdl <- aov(pcorr ~ trial_type*scrambled_type, ps)
tukey_res <- tukey_hsd(aov_mdl)
e2 <- eta_squared(aov_mdl)