function result = single_prop_initiated(labels, varargin)

if ( isempty(varargin) )
  denom = rows( labels );
else
  denom = numel( varargin{1} );
end

result = numel( find(labels, 'initiated_true', varargin{:}) ) / denom;

end