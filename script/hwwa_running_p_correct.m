function [pcorr_dat, pcorr_labs] = hwwa_running_p_correct(labels, specificity, varargin)

defaults = struct();
defaults.mask = rowmask( labels );
defaults.trial_bin_size = 1;
defaults.trial_step_size = 1;
defaults.use_cumulative = true;

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
  
  n_corr_tot = 0;
  n_incorr_tot = 0;
  
  for j = 1:numel(inds)
    current_I = inds{j};
    use_I = trial_ind(current_I);
    
    n_corr = numel( find(labels, 'correct_true', use_I) );
    n_incorr = numel( find(labels, 'correct_false', use_I) );
    
    pcorr = n_corr / (n_corr + n_incorr);
    
    total_n = n_incorr + n_incorr_tot + n_corr_tot + n_corr;
    
    if ( use_cumulative )
      pcorr_dat(i, j) = (n_corr + n_corr_tot) / total_n;
    else
      pcorr_dat(i, j) = pcorr;
    end
    
    n_corr_tot = n_corr_tot + n_corr;
    n_incorr_tot = n_incorr_tot + n_incorr;
  end
  
  max_ind = max( j, max_ind );
  
  append1( pcorr_labs, labels, trial_ind );
end

pcorr_dat = pcorr_dat(:, 1:max_ind);

end