function mask = find_correct(labels, varargin)

mask = find( labels, 'correct_true', varargin{:} );

end