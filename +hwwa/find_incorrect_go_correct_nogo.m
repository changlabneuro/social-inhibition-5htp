function ind = find_incorrect_go_correct_nogo(labels, varargin)

go = hwwa.find_incorrect_go( labels, varargin{:} );
nogo = find( labels, {'nogo_trial', 'correct_true'}, varargin{:} );
ind = union( go, nogo );

end