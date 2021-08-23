function bin_inds = bin_indices(mask, bin_size, step_size)

max_rows = numel( mask );
start = 1;
stop = min( max_rows, bin_size );
bin_idx = 1;

bin_inds = {};

while ( stop <= max_rows )
  bin_inds{bin_idx} = mask(start:stop);
  
  if ( stop == max_rows )
    break;
  end
  
  start = start + step_size;
  stop = min( start + bin_size - 1, max_rows );
  
  bin_idx = bin_idx + 1;
end

end