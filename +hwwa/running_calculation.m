function [pcorr_dat, pcorr_labs] = running_calculation(data, labels, specificity, func, varargin)

defaults = struct();
defaults.mask = rowmask( labels );
defaults.trial_bin_size = 1;
defaults.trial_step_size = 1;
defaults.use_cumulative = false;

validateattributes( func, {'function_handle'}, {'scalar'}, mfilename, 'func' );
assert_ispair( data, labels );

params = hwwa.parsestruct( defaults, varargin );

mask = params.mask;
trial_bin_size = params.trial_bin_size;
trial_step_size = params.trial_step_size;
use_cumulative = params.use_cumulative;

I = findall( labels, specificity, mask );
pcorr_dat = nan( numel(I), 1e4 );
max_ind = 0;
pcorr_labs = fcat();

for i = 1:numel(I)  
  trial_ind = I{i};
  
  max_rows = numel( trial_ind );
  start = 1;
  stop = min( max_rows, trial_bin_size );
  bin_idx = 1;
  inds = {};

  if ( trial_bin_size == 1 && trial_step_size == 1 )
    inds = arrayfun( @(x) x, 1:numel(trial_ind), 'un', 0 );
  else
    while ( stop <= max_rows )
      inds{bin_idx} = start:stop;

      if ( stop == max_rows ), break; end

      start = start + trial_step_size;
      stop = min( start + trial_bin_size - 1, max_rows );

      bin_idx = bin_idx + 1;
    end
  end
  
  assert( numel(inds) <= size(pcorr_dat, 2), 'P-corr dimension is too small.' );
  
  for j = 1:numel(inds)
    current_I = inds{j};
    use_I = trial_ind(current_I);
    
    if ( use_cumulative )
      if ( j > 1 )
        prev = pcorr_dat(i, 1:j-1);
      else
        prev = [];
      end
      
      pcorr_dat(i, j) = func( prev, data, use_I );
    else
      pcorr_dat(i, j) = func( data, use_I );
    end
  end
  
  max_ind = max( j, max_ind );
  
  append1( pcorr_labs, labels, trial_ind );
end

pcorr_dat = pcorr_dat(:, 1:max_ind);

end