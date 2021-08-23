function ind = find_nogo(labels, varargin)

ind = find( labels, 'nogo_trial', varargin{:} );

end