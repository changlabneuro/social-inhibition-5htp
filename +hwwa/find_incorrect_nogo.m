function ind = find_incorrect_nogo(labels, varargin)

ind = hwwa.find_nogo( labels, find(labels, 'correct_false', varargin{:}) );

end