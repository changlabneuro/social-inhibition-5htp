function labs = add_prev_trial_correct_labels( labs )

addcat( labs, 'prevcorrect' );

N = size( labs, 1 );

if ( N == 0 )
  return;
end

current_correct = fullcat( labs, 'correct' );
prev_correct = cell( size(current_correct) );

prev_correct{1} = 'prevcorrect_NaN';

for i = 2:N
  
  lab = current_correct{i-1};
  was_correct = strcmp( lab, 'correct_true' );
  was_incorrect = strcmp( lab, 'correct_false' );
  was_nan = ~was_correct && ~was_incorrect;
  
  if ( was_nan )
    prev_correct{i} = 'prevcorrect_NaN';
  elseif ( was_correct )
    prev_correct{i} = 'prevcorrect_true';
  else
    prev_correct{i} = 'prevcorrect_false';
  end
end

setcat( labs, 'prevcorrect', prev_correct );

end