function show_slope_permutation_test_performance(ps, each_I, ids)

displayed = false( size(ps) );

for i = 1:numel(ids)
  ind = ids(i).index;
  matches = cellfun( @(x) any(ismember(x, ind)), each_I );
  assert( nnz(matches) == 1, 'Expected 1 match; got %d.', nnz(matches) );
  
  if ( ~displayed(matches) )
    p_str = sprintf( 'p-slope-comparison = %0.2f', ps(matches) );
    ax = ids(i).axes;
    text( ax, min(get(ax, 'xlim')), max(get(ax, 'ylim')), p_str );
  end
end

end