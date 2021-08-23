function [labels, curr_inds, prev_inds] = label_gonogo_switchiness(labels)

[labels, curr_inds, prev_inds] = hwwa.label_switchiness( labels, 'trial_type' );

end