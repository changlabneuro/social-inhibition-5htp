conf = hwwa.config.load();

% use_files = { '18-F', '25-F', '27-F', '04-M', '24-F', '26-F', '28-F', '05-M' };
use_files = { '16-J', '17-J', '18-J', '21-J', '22-J', '23-J', '18-M', '19-M', '20-M', '21-M', '22-M' };

outs = hwwa_load_basic_behav_measures( ...
    'config', conf ...
  , 'files_containing', use_files ...
  , 'trial_bin_size', 50 ...
  , 'trial_step_size', 1 ...
  , 'is_parallel', true ...
);

date_dir = datestr( now, 'mmddyy' );
plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', date_dir );

shared_utils.io.require_dir( plot_p );

%%

labels = outs.labels';
rt = outs.rt;

addcat( labels, 'image_type' );

%%  rt

base_plot_subdir = '';

do_save = false;

is_per_monkey = true;
is_n_as_trial = true;
is_per_image_type = true;

c_is_across_days = [true];
c_is_per_gender = [false];
c_is_per_image_category = [true];
c_is_per_trial_time_bin = [false];
c_is_per_drug = [true];
c_is_per_trial_bin = [false];

C = dsp3.numel_combvec( c_is_across_days, c_is_per_gender ...
  , c_is_per_image_category, c_is_per_trial_time_bin, c_is_per_drug, c_is_per_trial_bin );

for idx = 1:size(C, 2)
  
plot_subdir = base_plot_subdir;

is_across_days = c_is_across_days(C(1, idx));
is_per_gender = c_is_per_gender(C(2, idx));
is_per_image_category = c_is_per_image_category(C(3, idx));
is_per_trial_time_bin = c_is_per_trial_time_bin(C(4, idx));
is_per_drug = c_is_per_drug(C(5, idx));
is_per_trial_bin = c_is_per_trial_bin(C(6, idx));

if ( ~is_per_image_category && ~is_per_gender )
  continue;
end

if ( is_across_days )
  base_prefix = 'across_days';
else
  base_prefix = '';
end

base_prefix = sprintf( '%s%s', base_prefix, ternary(is_n_as_trial, 'trial_n', '') );

uselabs = labels';
usedat = rt;

assert_ispair( usedat, uselabs );

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true', 'go_trial', 'correct_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pl = plotlabeled.make_common();

pl.x_order = { 'appetitive', 'threat', 'scrambled' };

pl.add_points = false;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.main_line_width = 3;
% pl.y_lims = [0.2, 0.34];
pl.y_lims = [0.2, 0.52];

fcats = {};
xcats = {};
pcats = { 'day' };

if ( is_per_image_type )
  gcats = { 'image_type' };
else
  gcats = { 'cue_delay' };
end

if ( is_per_monkey )
  fcats{end+1} = 'monkey';
  pcats{end+1} = 'monkey';
end

if ( is_per_drug )
  gcats{end+1} = 'drug';
end

if ( is_per_image_category )
  xcats{end+1} = 'target_image_category';
  plot_subdir = sprintf( '%s_per_image_category', plot_subdir );
end

if ( is_per_gender )
  xcats{end+1} = 'gender';
  plot_subdir = sprintf( '%s_per_gender', plot_subdir );
end

if ( is_per_trial_time_bin )
  gcats{end+1} = 'trial_time_bin';
  
  trial_time_bins = combs( uselabs, 'trial_time_bin', mask );
  trial_time = cellfun( @(x) str2double(x(strfind(x, 'time_')+5:max(strfind(x, '_'))-1)) ...
    , trial_time_bins );
  [~, sorted_ind] = sort( trial_time );
  
  pl.group_order = trial_time_bins(sorted_ind);
  
  plot_subdir = sprintf( '%s_per_trial_time_bin', plot_subdir );
end

spec = unique( cshorzcat(xcats, gcats, pcats, {'day'}) );

if ( is_n_as_trial )
  rtlabs = uselabs(mask);
  rtdat = usedat(mask);
else
  [rtlabs, I] = keepeach( uselabs', spec, mask );
  rtdat = rowop( usedat, I, @nanmedian );
end

if ( is_across_days )
  collapsecat( rtlabs, {'day'} );
else
  fcats{end+1} = 'day';
end

[figs, axs, I] = pl.figures( @errorbar, rtdat, rtlabs, fcats, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(figs)
    dsp3.req_savefig( figs(i), fullfile(plot_p, plot_subdir) ...
      , prune(rtlabs(I{i})), unique(cshorzcat(fcats, pcats, xcats, gcats)) ...
      , sprintf('rt_%s', base_prefix) );
  end
end

end

%%  rt over time

base_plot_subdir = 'per_day_panels';

do_save = true;
is_n_as_trial = true;
is_per_day = true;

plot_subdir = base_plot_subdir;

uselabs = labels';
usedat = rt;

assert_ispair( usedat, uselabs );

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true', 'go_trial', 'correct_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

if ( is_n_as_trial )
  rtlabs = uselabs(mask);
  rtdat = usedat(mask);
else
  error( 'Not yet implemented' );
%   [rtlabs, I] = keepeach( uselabs', spec, mask );
%   rtdat = rowop( usedat, I, @nanmedian );
end

if ( ~is_per_day )
  collapsecat( rtlabs, 'day' );
end

bins = combs( rtlabs, 'trial_bin' );
bin_ns = cellfun( @(x) x(numel('trial_bin__')+1:end), bins, 'un', 0 );
bin_ns = str2double( bin_ns );
[~, sorted_I] = sort( bin_ns );

pl = plotlabeled.make_common();
pl.x_order = bins(sorted_I);

fcats = { 'day', 'monkey' };
xcats = { 'trial_bin' };
pcats = { 'monkey', 'day', 'target_image_category' };
gcats = {};

[figs, axs, I] = pl.figures( @errorbar, rtdat, rtlabs, fcats, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(I)
    f = figs(i);
    
    dsp3.req_savefig( f, fullfile(plot_p, plot_subdir) ...
      , prune(rtlabs(I{i})), unique(cshorzcat(fcats, pcats, gcats)) ...
      , sprintf('rt_over_time%s', base_prefix) );
  end
end

%%  p correct, social

base_plot_subdir = '';

do_save = true;

is_per_monkey = true;
is_per_image_type = true;
is_per_run = false;

c_is_across_days = [false];
c_is_per_image_category = [true, false];
c_is_per_gender = [true, false];
c_is_per_trial_time_bin = [false];
c_is_per_drug = [true];

C = dsp3.numel_combvec( c_is_across_days, c_is_per_image_category ...
  , c_is_per_gender, c_is_per_trial_time_bin, c_is_per_drug );

for idx = 1:size(C, 2)
  
plot_subdir = base_plot_subdir; 

is_across_days = c_is_across_days(C(1, idx));
is_per_image_category = c_is_per_image_category(C(2, idx));
is_per_gender = c_is_per_gender(C(3, idx));
is_per_trial_time_bin = c_is_per_trial_time_bin(C(4, idx));
is_per_drug = c_is_per_drug(C(5, idx));

if ( ~is_per_gender && ~is_per_image_category )
  continue;
end

if ( is_across_days )
  base_prefix = 'across_days';
else
  base_prefix = '';
end

uselabs = labels';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pcorr_spec = { 'trial_type', 'cue_delay' };

if ( is_per_drug )
  pcorr_spec{end+1} = 'drug';
end

if ( is_per_run ), pcorr_spec{end+1} = 'date'; end
if ( ~is_across_days ), pcorr_spec{end+1} = 'day'; end
if ( is_per_monkey ), pcorr_spec{end+1} = 'monkey'; end

if ( is_per_image_category )
  pcorr_spec{end+1} = 'target_image_category';
  plot_subdir = sprintf( '%s_per_image_category', plot_subdir );
end

if ( is_per_gender )
  pcorr_spec{end+1} = 'gender';
  plot_subdir = sprintf( '%s_per_gender', plot_subdir );
end

if ( is_per_trial_time_bin )
  pcorr_spec{end+1} = 'trial_time_bin';
  plot_subdir = sprintf( '%s_per_trial_time_bin', plot_subdir );
end

[pcorr_labs, I] = keepeach( uselabs', pcorr_spec, mask );
pcorr_dat = rowzeros( numel(I) );

for i = 1:numel(I)
  n_corr = numel( find(uselabs, 'correct_true', I{i}) );
  n_incorr = numel( find(uselabs, 'correct_false', I{i}) );

  pcorr_dat(i) = n_corr / (n_corr + n_incorr);
end

pl = plotlabeled.make_common();

pl.x_order = { 'appetitive', 'threat', 'scrambled' };
pl.y_lims = [0., 1];
pl.summary_func = @(x) nanmedian(x, 1);

fcats = {};
xcats = {};

pl.color_func = @hot;

if ( is_across_days )
  xcats = {};
  gcats = { 'day' };
  pcats = { 'trial_type' };
  
  if ( is_per_monkey )
    fcats{end+1} = 'monkey';
    pcats{end+1} = 'monkey';
  end
  
  collapsecat( pcorr_labs, 'day' );
  
  func = @errorbar;
else
  xcats = {};
  gcats = { 'day' };
  pcats = { 'trial_type', 'monkey' };
  
  func = @errorbar;
end

if ( is_per_image_type )
  pcats{end+1} = 'image_type';
  fcats{end+1} = 'monkey';
end

if ( is_per_image_category )
  xcats{end+1} = 'target_image_category';
end

if ( is_per_gender )
  xcats{end+1} = 'gender';
end

if ( is_per_drug )
  gcats{end+1} = 'drug';
  
  if ( ~is_across_days )
    pcats{end+1} = 'drug';
  end
end

if ( is_per_trial_time_bin )
  gcats{end+1} = 'trial_time_bin';
  
  trial_time_bins = combs( pcorr_labs, 'trial_time_bin' );
  trial_time = cellfun( @(x) str2double(x(strfind(x, 'time_')+5:max(strfind(x, '_'))-1)) ...
    , trial_time_bins );
  [~, sorted_ind] = sort( trial_time );
  
  pl.group_order = trial_time_bins(sorted_ind);
  
  fcats{end+1} = 'day';
end

[fs, axs, I] = pl.figures( func, pcorr_dat, pcorr_labs, fcats, xcats, gcats, pcats );

shared_utils.plot.match_ylims( axs );
shared_utils.plot.fullscreen( fs );

if ( ~is_across_days )
  for i = 1:numel(axs)

    handles = findobj( axs(i), 'type', 'ErrorBar' );
    colors = flipud( cool( numel(handles) ) );

    for j = 1:numel(handles)
      handles(j).Color = colors(j, :);
      handles(j).LineWidth = 2;
    end
  end
end

% pl.bar( pcorr_dat, pcorr_labs, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(fs)
    dsp3.req_savefig( fs(i), fullfile(plot_p, plot_subdir), pcorr_labs(I{i}) ...
      , unique(cshorzcat(xcats, gcats, pcats, fcats)) ...
      , sprintf('percent_correct_%s', base_prefix) );
  end
end

end

%%

trial_bin_size = 1;
trial_step_size = 1;
use_cumulative = true;

uselabs = labels';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true', 'no_errors', 'wrong_go_nogo'} ...
  , @findnot, {'tar', '011619'} ...
);

% pcorr_spec = { 'trial_type', 'date' };
pcorr_spec = { 'trial_type', 'date', 'target_image_category' };

I = findall( uselabs, pcorr_spec, mask );
pcorr_dat = nan( numel(I), 1e3 );
max_ind = 0;
pcorr_labs = fcat();

for i = 1:numel(I)
  shared_utils.general.progress( i, numel(I) );
  
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
    
    n_corr = numel( find(uselabs, 'correct_true', use_I) );
    n_incorr = numel( find(uselabs, 'correct_false', use_I) );
    
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
  
  append1( pcorr_labs, uselabs, trial_ind );
end

pcorr_dat = pcorr_dat(:, 1:max_ind);

%%

base_plot_subdir = 'per_day_panels_cumulative_per_image_category';
base_prefix = '';

plot_subdir = base_plot_subdir; 

do_save = true;
is_per_day = true;

if ( ~is_per_day )
  collapsecat( pcorr_labs, 'day' );
end

fcats = { 'monkey', 'day' };
gcats = { 'trial_type' };
pcats = { 'monkey', 'day', 'target_image_category' };

pl = plotlabeled.make_common();

pl.x_order = bins(sorted_I);
pl.y_lims = [0, 1];
pl.summary_func = @(x) nanmedian(x, 1);
pl.add_x_tick_labels = false;
pl.add_errors = false;

[figs, axs, I] = pl.figures( @lines, pcorr_dat, pcorr_labs, fcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(figs)
    dsp3.req_savefig( figs(i), fullfile(plot_p, plot_subdir), pcorr_labs(I{i}) ...
      , unique(cshorzcat(pcats, fcats)) ...
      , sprintf('percent_correct_%s', base_prefix) );
  end
end


%%  p correct, social, over time

base_plot_subdir = 'per_day_panels_slide_window';
base_prefix = '';

do_save = true;
is_per_day = true;

plot_subdir = base_plot_subdir; 

uselabs = labels';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pcorr_spec = { 'trial_type', 'date', 'trial_bin', 'target_image_category' };

[pcorr_labs, I] = keepeach( uselabs', pcorr_spec, mask );
pcorr_dat = rowzeros( numel(I) );

for i = 1:numel(I)
  n_corr = numel( find(uselabs, 'correct_true', I{i}) );
  n_incorr = numel( find(uselabs, 'correct_false', I{i}) );

  pcorr_dat(i) = n_corr / (n_corr + n_incorr);
end

if ( ~is_per_day )
  collapsecat( pcorr_labs, 'day' );
end

fcats = { 'monkey', 'day' };
xcats = { 'trial_bin' };
gcats = { 'trial_type' };
pcats = { 'monkey', 'day', 'target_image_category' };

bins = combs( pcorr_labs, 'trial_bin' );
bin_ns = cellfun( @(x) x(numel('trial_bin__')+1:end), bins, 'un', 0 );
bin_ns = str2double( bin_ns );
[~, sorted_I] = sort( bin_ns );

pl = plotlabeled.make_common();

% pl.x_order = bins(sorted_I);
pl.y_lims = [0, 1];
pl.summary_func = @(x) nanmedian(x, 1);
pl.add_x_tick_labels = false;
pl.add_errors = false;

[figs, axs, I] = pl.figures( @errorbar, pcorr_dat, pcorr_labs, fcats, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(figs)
    dsp3.req_savefig( figs(i), fullfile(plot_p, plot_subdir), pcorr_labs(I{i}) ...
      , unique(cshorzcat(pcats, fcats)) ...
      , sprintf('percent_correct_%s', base_prefix) );
  end
end

%%  p correct, social drug

do_save = true;

base_plot_subdir = '';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pcorr_spec = { 'date', 'trial_type', 'cue_delay' };
plot_subdir = base_plot_subdir; 

is_per_drug = true;
is_per_day = false;
is_per_image_category = true;
is_per_monkey = true;

fcats = {};
gcats = {};
pcats = { 'trial_type' };
xcats = {};

if ( is_per_drug )
  if ( is_per_day )
    pcats{end+1} = 'drug'; 
  else
    gcats{end+1} = 'drug';
  end
end
if ( is_per_day ), gcats{end+1} = 'day'; end
if ( is_per_monkey )
  fcats{end+1} = 'monkey';
  pcats{end+1} = 'monkey'; 
end

if ( is_per_image_category )
  pcorr_spec{end+1} = 'target_image_category';
  xcats{end+1} = 'target_image_category';
  
  plot_subdir = sprintf( '%s_per_image_category', plot_subdir );
end

[pcorr_labs, I] = keepeach( uselabs', pcorr_spec, mask );
pcorr_dat = rowzeros( numel(I) );

for i = 1:numel(I)
  n_corr = numel( find(uselabs, 'correct_true', I{i}) );
  n_incorr = numel( find(uselabs, 'correct_false', I{i}) );

  pcorr_dat(i) = n_corr / (n_corr + n_incorr);
end

pl = plotlabeled.make_common();

[fs, axs, I] = pl.figures( @errorbar, pcorr_dat, pcorr_labs, fcats, xcats, gcats, pcats );

shared_utils.plot.match_ylims( axs );
shared_utils.plot.fullscreen( fs );

if ( is_per_day )
  for i = 1:numel(axs)

    handles = findobj( axs(i), 'type', 'ErrorBar' );
    colors = flipud( cool( numel(handles) ) );

    for j = 1:numel(handles)
      handles(j).Color = colors(j, :);
      handles(j).LineWidth = 2;
    end
  end
end

if ( do_save )
  for i = 1:numel(fs)
    dsp3.req_savefig( fs(i), fullfile(plot_p, plot_subdir), pcorr_labs(I{i}) ...
      , unique(cshorzcat(xcats, gcats, pcats, fcats)) ...
      , sprintf('percent_correct_%s', base_prefix) );
  end
end