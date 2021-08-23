function mask = find_correct_go(labels, varargin)

mask = hwwa.find_go( labels, find(labels, 'correct_true', varargin{:}) );

end