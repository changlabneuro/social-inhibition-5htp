function tbl = to_hBayesDM_table(labels)

monks = cellstr( labels, 'monkey' );
[~, ~, subj_id] = unique( monks );

go_trial = hwwa.find_go( labels );
nogo_trial = hwwa.find_nogo( labels );

trial_types = nan( size(labels, 1), 1 );
trial_types(go_trial) = 1;
trial_types(nogo_trial) = 2;

go_choice = find( labels, 'go_choice' );
nogo_choice = find( labels, 'nogo_choice' );

key_pressed = nan( size(labels, 1), 1 );
key_pressed(go_choice) = 1;
key_pressed(nogo_choice) = 0;

outcome = nan( size(labels, 1), 1 );
outcome(find(labels, 'correct_true')) = 1;
outcome(find(labels, 'correct_false')) = -1;

m = [ subj_id, trial_types, key_pressed, outcome ];
tbl = array2table( m, 'VariableNames', {'subjID', 'cue', 'keyPressed', 'outcome'} );

end