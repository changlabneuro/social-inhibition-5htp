function labs = add_trial_bin_labels(labs, bin_size, step_size, mask)

if ( nargin < 4 )
  mask = rowmask( labs );
end

is_sliding = bin_size ~= step_size;

max_rows = numel( mask );
start = 1;
stop = min( max_rows, bin_size );
bin_idx = 1;

addcat( labs, 'trial_bin' );

if ( is_sliding )
  newlabs = fcat();
end

while ( stop <= max_rows )
  ind = mask( start:stop );
  
  trial_bin_label = sprintf( 'trial_bin__%d', bin_idx );
  
  if ( is_sliding )
    subset = prune( keep(copy(labs), ind) );
    setcat( subset, 'trial_bin', trial_bin_label );
    append( newlabs, subset );
  else
    setcat( labs, 'trial_bin', trial_bin_label, ind );
  end
  
  if ( stop == max_rows )
    break;
  end
  
  start = start + step_size;
  stop = min( start + bin_size - 1, max_rows );
  
  bin_idx = bin_idx + 1;
end

if ( is_sliding )
  resize( labs, rows(newlabs) );
  assign( labs, newlabs, 1:rows(labs) );
end

end