function labels = hwwa_n_go_trial_sequence(labels, each_I, precede_trial_type, precede_correct_type, n, allow_repeated)

seq_label = sprintf( 'prev_%s_%s_%d', precede_correct_type, precede_trial_type, n );
seq_cat = 'prev_correct_sequence';

tt_label = sprintf( 'prev_%s', precede_trial_type );
tt_cat = 'prev_trial_type';

n_label = sprintf( 'prev_%d', n );
n_cat = 'prev_sequence_length';

corr_label = sprintf( 'prev_%s', precede_correct_type );
corr_cat = 'prev_correct';

tt_corr_cat = 'prev_trial_type_sequence_length';
tt_corr_label = sprintf( 'prev_%s_%d', precede_trial_type, n );

switch_tt_cat = 'switch_trial_type';
same_tt_label = 'same_trial_type';
switch_tt_label = 'switch_trial_type';

addsetcat( labels, seq_cat, 'no_prev_correct_sequence' );
addcat( labels, tt_cat );
addcat( labels, n_cat );
addcat( labels, corr_cat );
addcat( labels, tt_corr_cat );
addcat( labels, switch_tt_cat );

for i = 1:numel(each_I)
  seq = 0;
  stp = 1;
  each_ind = each_I{i};
  
  while ( stp <= numel(each_ind) )
    ind = each_ind(stp);
    match = strcmp( ...
        cellstr(labels, {'correct', 'trial_type'}, ind) ...
      , {precede_correct_type, precede_trial_type} ...
    );
    
    if ( all(match) )
      seq = seq + 1;
    else
      seq = 0;
    end
    
    if ( seq == n && ind < size(labels, 1) )
      curr_tt = cellstr( labels, 'trial_type', ind+1 );
      
      if ( strcmp(curr_tt, precede_trial_type) )
        setcat( labels, switch_tt_cat, same_tt_label, ind+1 );
      else
        setcat( labels, switch_tt_cat, switch_tt_label, ind+1 );
      end
      
      setcat( labels, seq_cat, seq_label, ind+1 );
      setcat( labels, tt_cat, tt_label, ind+1 );
      setcat( labels, n_cat, n_label, ind+1 );
      setcat( labels, corr_cat, corr_label, ind+1 );
      setcat( labels, tt_corr_cat, tt_corr_label, ind+1 );
      
      if ( allow_repeated )
        seq = 0;
      end
    end
    
    stp = stp + 1;
  end
end

end