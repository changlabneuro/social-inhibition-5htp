function [curr_inds, prev_inds] = find_nback(labels, n_back, specificity, mask)

if ( nargin < 4 )
  mask = rowmask( labels );
end

if ( nargin < 3 )
  specificity = {};
end

validateattributes( n_back, {'numeric'}, {'integer'}, mfilename, 'n_back' );

I = findall_or_one( labels, specificity, mask );

curr_inds = cell( numel(I), 1 );
prev_inds = cell( size(curr_inds) );

for i = 1:numel(I)
  if ( numel(I{i}) <= n_back )
    continue;
  end
  
  curr_inds{i} = I{i}(1+n_back:end);
  prev_inds{i} = I{i}(1:end-n_back);
end

curr_inds = vertcat( curr_inds{:} );
prev_inds = vertcat( prev_inds{:} );

end