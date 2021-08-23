function [samples, min_edges, out_labels, out_inds] = ...
time_slide_bin_apply(data, time, labels, each, window_size, step_size, func, mask)

assert_ispair( data, labels );
assert_ispair( time, labels );
validateattributes( time, {'double'}, {'vector'}, mfilename, 'time' );

if ( nargin < 8 )
  mask = rowmask( data );
end

[out_labels, each_I] = keepeach( labels', each, mask );

min_edges = 0:step_size:ceil(max(time));

if ( max(min_edges) < max(time) )
  min_edges(end+1) = min_edges(end) + step_size;
end

max_edges = min_edges + window_size;
samples = cell( numel(each_I), numel(min_edges) );

keep_out_inds = nargout > 3;

if ( keep_out_inds )
  out_inds = cell( size(samples) );
end

for i = 1:numel(each_I)
  subset_times = time(each_I{i});

  for j = 1:numel(min_edges)
    within_edge = subset_times >= min_edges(j) & subset_times < max_edges(j);

    subset_ind = each_I{i}(within_edge);
    result = func( data(subset_ind, :), labels, subset_ind );
    assert( rows(result) == 1, 'Output of bin function must produce one row.' );
    samples{i, j} = result;

    if ( keep_out_inds )    
      out_inds{i, j} = subset_ind;
    end
  end
end

assert_ispair( samples, out_labels );
assert( numel(min_edges) == size(samples, 2) );

end