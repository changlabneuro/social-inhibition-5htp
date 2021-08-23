function labs = add_initiated_labels(labs)

import shared_utils.assertions.*;

assert__isa( labs, 'fcat', 'labels' );

% initiated trials are those without error, or those where the error is
% a look away from the 
did_initiate_ind = find( labs, {'wrong_go_nogo', 'no_errors'} );
no_initiated_ind = setdiff( 1:size(labs, 1), did_initiate_ind );

addcat( labs, 'initiated' );
setcat( labs, 'initiated', 'initiated_true', did_initiate_ind );
setcat( labs, 'initiated', 'initiated_false', no_initiated_ind );

end