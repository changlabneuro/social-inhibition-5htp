function [out_data, out_labels, bin_inds] = sliding_window_summary(data, labels, each, window_size, step_size, func, mask)

assert_ispair( data, labels );

if ( nargin < 7 )
  mask = rowmask( labels );
end

I = findall( labels, each, mask );

out_data = cell( size(I) );
out_labels = fcat();
bin_inds = cell( size(I) );

for i = 1:numel(I)
  ind = I{i};
  window_inds = shared_utils.vector.slidebin( 1:numel(ind), window_size, step_size, true );
  
  if ( nargout > 2 )
    tmp_inds = cellfun( @(x) ind(x), window_inds(:), 'un', 0 );
  else
    tmp_inds = {};
  end
  
  tmp_data = cell( numel(window_inds), 1 );
  
  for j = 1:numel(window_inds)
    windowed_ind = ind(window_inds{j});
    result = func( data, labels, windowed_ind );
    
    assert( rows(result) == 1, 'Output of summary function must produce 1 row.' );
    
    tmp_data{j} = result;
    append1( out_labels, labels, windowed_ind );
  end
  
  out_data{i} = vertcat( tmp_data{:} );
  bin_inds{i} = tmp_inds;
end

out_data = vertcat( out_data{:} );
bin_inds = vertcat( bin_inds{:} );

assert_ispair( out_data, out_labels );

end