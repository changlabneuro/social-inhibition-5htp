library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\prop_fixated_yy')
mdl <- wilcox_test(ps, prop_fixated ~ drug)
ef <- wilcox_effsize(ps, prop_fixated ~ drug)