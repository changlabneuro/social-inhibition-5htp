function labs = add_correct_labels(labs)

import shared_utils.assertions.*;

assert__isa( labs, 'fcat', 'labels' );

% initiated trials are those without error, or those where the error is
% a look away from the 
correct_ind = find( labs, {'no_errors'} );
incorrect_ind = setdiff( 1:size(labs, 1), correct_ind );

addcat( labs, 'correct' );
setcat( labs, 'correct', 'correct_true', correct_ind );
setcat( labs, 'correct', 'correct_false', incorrect_ind );

end