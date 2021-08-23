function dirs = collapse_left_saccades(dirs, labels)

assert_ispair( dirs, labels );
left_trials = find( labels, 'center-left' );
dirs(left_trials, 1) = -dirs(left_trials, 1);

end