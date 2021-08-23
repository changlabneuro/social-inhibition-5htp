function hwwa_run_plot_drug_learning(outs)

conf = hwwa.config.load();

use_files = [ hwwa.get_image_5htp_days(), hwwa.get_image_saline_days() ];
use_files = cellstr( hwwa.to_date(use_files) );

if ( nargin < 1 || isempty(outs) )
  outs = hwwa_load_basic_behav_measures( ...
      'config', conf ...
    , 'files_containing', use_files ...
    , 'trial_bin_size', 50 ...
    , 'trial_step_size', 1 ...
    , 'is_parallel', true ...
  );
end

cumulative_percent_correct( outs );
cumulative_rt( outs );

% average_level_percent_correct( outs );
% average_level_rt( outs );

end

function cumulative_percent_correct(outs)

% kinds = { 'social_minus_scrambled', 'social_vs_scrambled', 'threat_vs_appetitive' };
kinds = { 'threat_vs_appetitive' };
per_monks = [ true ];
% is_normed = [ true, false ];
is_normed = false;

C = dsp3.numel_combvec( kinds, per_monks, is_normed );

for i = 1:size(C, 2)
  shared_utils.general.progress( i, size(C, 2) );
  
  kind = kinds{C(1, i)};
  is_per_monkey = per_monks(C(2, i));
  is_sal_normalized = is_normed(C(3, i));
  
  rt = outs.rt;
  labels = outs.labels';
  norm_func = ternary( is_sal_normalized, @hwwa.saline_normalize, 'no_norm' );
  
  mask = get_base_mask( labels );
  rt = indexpair( rt, labels, mask );
  
  monk_dir = ternary( is_per_monkey, 'per_monk', 'collapsed_monk' );
  norm_dir = ternary( is_sal_normalized, 'norm', 'no-norm' );
  base_subdir = sprintf( '%s_%s_%s', kind, monk_dir, norm_dir );
  
  hwwa_plot_learning_running_p_correct( rt, labels ...
    , 'base_subdir', base_subdir ...
    , 'do_save', true ...
    , 'trial_bin_size', 25 ...
    , 'trial_step_size', 1 ...
    , 'colored_lines_are', kind ...
    , 'combine_days', false ...
    , 'is_per_monkey', is_per_monkey ...
    , 'is_per_drug', true ...
    , 'is_per_day', false ...
    , 'is_rt', false ...
    , 'norm_func', norm_func ...
  );
end
end

function cumulative_rt(outs)

% kinds = { 'social_minus_scrambled', 'social_vs_scrambled', 'threat_vs_appetitive' };
kinds = { 'threat_vs_appetitive' };
per_monks = [ true ];
% is_normed = [ true, false ];
is_normed = [ false ];

C = dsp3.numel_combvec( kinds, per_monks, is_normed );

for i = 1:size(C, 2)
  shared_utils.general.progress( i, size(C, 2) );
  
  kind = kinds{C(1, i)};
  is_per_monkey = per_monks(C(2, i));
  is_sal_normalized = is_normed(C(3, i));
  
  if ( strcmp(kind, 'scrambled_minus_social') )
    lims = [-0.5, 0.5];
  else
    lims = [0, 0.5];
  end
  
  norm_func = ternary( is_sal_normalized, @hwwa.saline_normalize, 'no_norm' );
  
  if ( is_sal_normalized )
    lims = [];
  end
  
  rt = outs.rt;
  labels = outs.labels';
  
  mask = get_base_mask( labels );
  rt = indexpair( rt, labels, mask );
  
  monk_dir = ternary( is_per_monkey, 'per_monk', 'collapsed_monk' );
  norm_dir = ternary( is_sal_normalized, 'norm', 'no-norm' );
  base_subdir = sprintf( '%s_%s_%s', kind, monk_dir, norm_dir );

  hwwa_plot_learning_running_p_correct( rt, labels' ...
    , 'base_subdir', base_subdir ...
    , 'do_save', true ...
    , 'trial_bin_size', 25 ...
    , 'trial_step_size', 1 ...
    , 'colored_lines_are', kind ...
    , 'combine_days', false ...
    , 'is_per_monkey', is_per_monkey ...
    , 'is_rt', true ...
    , 'is_per_drug', true ...
    , 'is_per_day', false ...
    , 'line_ylims', lims ...
    , 'line_xlims', [0, 100] ...
    , 'norm_func', norm_func ...
  );
end
end

function mask = get_base_mask(labels)

mask = findnone( labels, {'021819', 'ephron'} );
mask = union( mask, get_ephron_mask(labels) );

end

function mask = get_ephron_mask(labels)

days = {'041019', '041219', '041619', '041819', '042319', '042519' ...
  , '040819', '041119', '041519', '041719', '042219', '042419' ...
  , '042919', '043019', '050219', '050319' ...
};

mask = fcat.mask( labels ...
  , @find, days ...
  , @find, 'ephron' ...
);

end

%%  average, percent correct

function average_level_percent_correct(outs)

kinds = { 'social_vs_scrambled', 'threat_vs_appetitive' };

for i = 1:numel(kinds)
  kind = kinds{i};

  hwwa_plot_average_learning( outs.rt, outs.labels' ...
    , 'base_subdir', sprintf('collapsed_%s', kind) ...
    , 'do_save', false ...
    , 'colored_lines_are', kind ...
    , 'is_per_monkey', false ...
    , 'is_per_drug', true ...
  );
end

end

%%  average, n trials rt

function average_level_rt(outs)

kinds = { 'social_vs_scrambled', 'threat_vs_appetitive' };

for i = 1:numel(kinds)
  kind = kinds{i};

  hwwa_plot_average_learning( outs.rt, outs.labels' ...
    , 'base_subdir', sprintf('collapsed_%s', kind) ...
    , 'do_save', false ...
    , 'colored_lines_are', kind ...
    , 'is_per_monkey', false ...
    , 'is_rt', true ...
    , 'is_per_drug', true ...
  );
end

end
