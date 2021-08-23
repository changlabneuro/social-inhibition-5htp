function hwwa_run_plot_face_perinea_learning(outs)

conf = hwwa.config.load();

% use_files = { '18-F', '25-F', '27-F', '04-M', '24-F', '26-F', '28-F', '05-M' };
use_files = { '16-J', '17-J', '18-J', '21-J', '22-J', '23-J', '18-M', '19-M', '20-M', '21-M', '22-M' };

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

average_level_percent_correct( outs );
average_level_rt( outs );

end

function cumulative_percent_correct(outs)

kinds = { 'scrambled_minus_social', 'social_vs_scrambled', 'threat_vs_appetitive' };
per_monks = [ true, false ];

C = dsp3.numel_combvec( kinds, per_monks );

for i = 1:size(C, 2)
  shared_utils.general.progress( i, size(C, 2) );
  
  kind = kinds{C(1, i)};
  is_per_monkey = per_monks(C(2, i));
  
  hwwa_plot_learning_running_p_correct( outs.rt, outs.labels' ...
    , 'base_subdir', sprintf('%s%s', kind, ternary(is_per_monkey, 'per_monk', 'collapsed_monk')) ...
    , 'do_save', true ...
    , 'trial_bin_size', 1 ...
    , 'trial_step_size', 1 ...
    , 'colored_lines_are', kind ...
    , 'combine_days', true ...
    , 'is_per_monkey', is_per_monkey ...
    , 'is_per_drug', false ...
    , 'is_rt', false ...
  );
end
end

function cumulative_rt(outs)

kinds = { 'scrambled_minus_social', 'social_vs_scrambled', 'threat_vs_appetitive' };
per_monks = [ true, false ];

C = dsp3.numel_combvec( kinds, per_monks );

for i = 1:size(C, 2)
  shared_utils.general.progress( i, size(C, 2) );
  
  kind = kinds{C(1, i)};
  is_per_monkey = per_monks(C(2, i));
  
  if ( strcmp(kind, 'scrambled_minus_social') )
    lims = [-0.5, 0.5];
  else
    lims = [0, 0.5];
  end

  hwwa_plot_learning_running_p_correct( outs.rt, outs.labels' ...
    , 'base_subdir', sprintf('%s%s', kind, ternary(is_per_monkey, 'per_monk', 'collapsed_monk')) ...
    , 'do_save', true ...
    , 'trial_bin_size', 1 ...
    , 'trial_step_size', 1 ...
    , 'colored_lines_are', kind ...
    , 'combine_days', true ...
    , 'is_per_monkey', is_per_monkey ...
    , 'is_rt', true ...
    , 'is_per_drug', false ...
    , 'line_ylims', lims ...
    , 'line_xlims', [0, 100] ...
  );
end
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
    , 'is_per_drug', false ...
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
    , 'is_per_drug', false ...
  );
end

end