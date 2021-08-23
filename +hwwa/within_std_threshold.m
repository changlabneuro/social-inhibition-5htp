function [in_bounds, I] = within_std_threshold(data, labels, n_devs, each, mask)

if ( nargin < 5 )
  mask = rowmask( labels );
end

assert_ispair( data, labels );
validateattributes( data, {'double'}, {'vector'}, mfilename, 'data' );

I = findall( labels, each, mask );
in_bounds = false( rows(data), 1 );

for i = 1:numel(I)
  subset_data = data(I{i});
  mu = nanmean( subset_data );
  sigma = nanstd( subset_data );
  
  lb = mu - sigma * n_devs;
  ub = mu + sigma * n_devs;
  
  in_bounds(I{i}) = subset_data >= lb & subset_data <= ub;
end

end