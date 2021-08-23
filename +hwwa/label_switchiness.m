function [labels, curr_inds, prev_inds] = label_switchiness(labels, orig_cat)

rep_cat = sprintf( 'switch_%s', orig_cat );

[repeat_inds, switch_inds, curr_inds, prev_inds] = ...
  hwwa.find_repeats( labels, orig_cat, 'date' );

addcat( labels, rep_cat );

repeat_values = cellstr( labels, orig_cat, repeat_inds );
repeat_values = cellfun( @(x) sprintf('repeated_%s', x), repeat_values, 'un', 0 );

setcat( labels, rep_cat, repeat_values, repeat_inds );

prev_non_repeat_values = cellstr( labels, orig_cat, switch_inds-1 );
switch_from = cellfun( @(x) sprintf('previous_%s', x), prev_non_repeat_values, 'un', 0 );

setcat( labels, rep_cat, switch_from, switch_inds );

end