function labels = label_prev_correct(labels)

[curr_inds, prev_inds] = hwwa.find_nback( labels, 1, 'date' );

correct_cat = 'prev_correct';
incorrect_lab = sprintf( '%s_false', correct_cat );
correct_lab = sprintf( '%s_true', correct_cat );

prev_is_correct = find( labels, 'correct_true', prev_inds );
prev_is_incorrect = find( labels, 'correct_false', prev_inds );

% If the previous trial is correct, current will be labeled
% prev_correct_true, etc.
curr_is_correct = prev_is_correct + 1;
curr_is_incorrect = prev_is_incorrect + 1;

% BUT, ignore subsequent trials that "spill over" into the next date.
curr_is_correct = intersect( curr_is_correct, curr_inds );
curr_is_incorrect = intersect( curr_is_incorrect, curr_inds );

addcat( labels, correct_cat );
setcat( labels, correct_cat, correct_lab, curr_is_correct );
setcat( labels, correct_cat, incorrect_lab, curr_is_incorrect );

end