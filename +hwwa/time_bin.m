function [samples, edges, labels] = time_bin(data, time, labels, each, window_size, mask)

assert_ispair( data, labels );
assert_ispair( time, labels );
validateattributes( time, {'double'}, {'vector'}, mfilename, 'time' );

if ( nargin < 6 )
  mask = rowmask( data );
end

[labels, each_I] = keepeach( labels, each, mask );

num_bins = ceil( max(time) / window_size );
edges = linspace( 0, num_bins * window_size, num_bins+1 );

min_edges = edges(1:end-1);
max_edges = edges(2:end);

samples = cell( numel(each_I), numel(edges)-1 );

for i = 1:numel(each_I)
  subset_times = time(each_I{i});
  
  nearest_edges = arrayfun( @(x) find(x >= min_edges & x < max_edges), subset_times );
  
  for j = 1:numel(nearest_edges)
    samples{i, nearest_edges(j)}(end+1, :) = data(each_I{i}(j), :);
  end
end

end