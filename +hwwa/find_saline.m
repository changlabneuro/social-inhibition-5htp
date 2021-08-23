function mask = find_saline(labels, varargin)

mask = find( labels, 'saline', varargin{:} );

end