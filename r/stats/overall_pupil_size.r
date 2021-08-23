library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081621\\pupil_size')
mdl <- wilcox_test(ps, pupil_size ~ drug)
ef <- wilcox_effsize(ps, pupil_size ~ drug)