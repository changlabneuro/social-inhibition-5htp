library(rstatix)

ps <- read_pair_dir('D:\\data\\hwwa\\export\\approach_avoid\\behavior\\081721\\fix_pattern_freq_3_go_trial')
mdl <- wilcox_test(ps, fix_pattern_freq ~ drug)
ef <- wilcox_effsize(ps, fix_pattern_freq ~ drug)

# nogo 1: ef = 0.115
# nogo 2: ef = 0.116
# nogo 3: ef = 0.152
#   go 1: ef = 0.28
#   go 2: ef = 0.26
#   go 3: ef = 0.16