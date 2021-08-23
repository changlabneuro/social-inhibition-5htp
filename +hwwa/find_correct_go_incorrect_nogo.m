function ind = find_correct_go_incorrect_nogo(labels, varargin)

go = hwwa.find_correct_go( labels, varargin{:} );
nogo = find( labels, {'nogo_trial', 'correct_false'}, varargin{:} );

ind = union( go, nogo );

end