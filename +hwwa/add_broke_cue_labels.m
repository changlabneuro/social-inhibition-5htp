function labs = add_broke_cue_labels(labs)

did_break_ind = find( labs, 'broke_cue_fixation' );
did_not_break = setdiff( 1:size(labs, 1), did_break_ind );

addcat( labs, 'broke_cue' );
setcat( labs, 'broke_cue', 'broke_true', did_break_ind );
setcat( labs, 'broke_cue', 'broke_false', did_not_break );

prune( labs );

end