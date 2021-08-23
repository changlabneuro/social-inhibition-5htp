function mask = find_incorrect_go(labels, varargin)

mask = hwwa.find_go( labels, find(labels, 'correct_false', varargin{:}) );

end