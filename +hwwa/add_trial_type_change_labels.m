function labs = add_trial_type_change_labels( labs )

if ( isempty(labs) ), return; end

trial_types = fullcat( labs, 'trial_type' );

changed = cell( size(trial_types) );

changed_lab = 'type_changed';
changed_cat = 'type_changed';

changed_false_lab = sprintf( '%s_false', changed_lab );
changed_true_lab = sprintf( '%s_true', changed_lab );

prev_lab = trial_types{1};
changed{1} = changed_false_lab;

stp = 2;
N = numel( trial_types );

while ( stp <= N )
  curr_lab = trial_types{stp};
  
  did_change = ~strcmp( curr_lab, prev_lab );
  
  if ( did_change )
    changed{stp} = changed_true_lab;
  else
    changed{stp} = changed_false_lab;
  end
  
  stp = stp + 1;
  
  prev_lab = curr_lab;
end

addcat( labs, changed_cat );
setcat( labs, changed_cat, changed );

end