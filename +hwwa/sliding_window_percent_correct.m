function [props, out_labels, bin_inds] = sliding_window_percent_correct(labels, each, window_size, step_size, mask)

if ( nargin < 5 )
  mask = rowmask( labels );
end

I = findall( labels, each, mask );

props = cell( size(I) );
out_labels = fcat();
bin_inds = cell( size(I) );

for i = 1:numel(I)
  ind = I{i};
  window_inds = shared_utils.vector.slidebin( 1:numel(ind), window_size, step_size, true );
  pcorr_dat = zeros( numel(window_inds), 1 );
  
  if ( nargout > 2 )
    tmp_inds = cellfun( @(x) ind(x), window_inds(:), 'un', 0 );
  else
    tmp_inds = {};
  end
  
  for j = 1:numel(window_inds)
    windowed_ind = ind(window_inds{j});
    
    n_corr = numel( find(labels, 'correct_true', windowed_ind) );
    n_incorr = numel( find(labels, 'correct_false', windowed_ind) );

    pcorr_dat(j) = n_corr / (n_corr + n_incorr);
    append1( out_labels, labels, windowed_ind );
  end
  
  props{i} = pcorr_dat;
  bin_inds{i} = tmp_inds;
end

props = vertcat( props{:} );
bin_inds = vertcat( bin_inds{:} );

assert_ispair( props, out_labels );

end