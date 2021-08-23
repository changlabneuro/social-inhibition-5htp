function hwwa_run_plot_nondrug_cumulative_learning(outs)

conf = hwwa.config.load();

if ( nargin < 1 || isempty(outs) )
  outs = hwwa_load_basic_behav_measures( ...
      'config', conf ...
    , 'files_containing', cellstr(hwwa.to_date(hwwa.get_learning_days())) ...
    , 'trial_bin_size', 50 ...
    , 'trial_step_size', 50 ...
    , 'is_parallel', true ...
  );
end

cumulative_percent_correct( outs );
% cumulative_rt( outs );

% average_level_percent_correct( outs );
% average_level_rt( outs );

end

function cumulative_percent_correct(outs)

% kinds = { 'social_minus_scrambled', 'social_vs_scrambled', 'threat_vs_appetitive' };
kinds = { 'threat_vs_appetitive' };
per_monks = [ true, false ];
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
    , 'combine_days', true ...
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

mask = [];

if ( count(labels, 'ephron') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'ephron' ...
    , @findor, hwwa.get_ephron_learning_days() ...
    ));
end

if ( count(labels, 'hitch') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'hitch' ...
    , @findor, hwwa.get_hitch_learning_days() ...
    ));
end

if ( count(labels, 'tar') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'tar' ...
    , @findor, hwwa.get_tarantino_learning_days() ...
    ));
end

% Only initiated trials.
mask = add_initiated_mask( labels, mask );

end

function mask = add_initiated_mask(labels, mask)

mask = fcat.mask( labels, mask ...
  , @find, 'initiated_true' ...
);

gonogo_error = find( labels, 'wrong_go_nogo' );
ok_trial = find( labels, 'correct_true' );
error_or_ok = union( gonogo_error, ok_trial );

mask = intersect( mask, error_or_ok );

end