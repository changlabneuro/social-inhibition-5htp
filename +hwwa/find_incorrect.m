function mask = find_incorrect(labels, varargin)

mask = find( labels, 'correct_false', varargin{:} );

end