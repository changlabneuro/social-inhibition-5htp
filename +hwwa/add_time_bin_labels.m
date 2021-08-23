function labels = add_time_bin_labels(labels, time_stamps, bin_size)

assert_ispair( time_stamps, labels );

adjusted_time = time_stamps - min( time_stamps );
max_time = max( adjusted_time );

n_bins = ceil( max_time / bin_size );

edges = (0:n_bins) * bin_size;
assert( edges(end) >= max_time );

addcat( labels, 'trial_time_bin' );

for i = 1:n_bins
  is_within_time = adjusted_time >= edges(i) & adjusted_time < edges(i+1);
  
  bin_label = sprintf( 'time_%d_%d', edges(i), edges(i+1) );
  
  setcat( labels, 'trial_time_bin', bin_label, find(is_within_time) );
end

end