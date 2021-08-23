function [labels, curr_inds, prev_inds] = label_image_category_switchiness(labels)

[labels, curr_inds, prev_inds] = hwwa.label_switchiness( labels, 'target_image_category' );

end