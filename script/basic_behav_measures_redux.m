conf = hwwa.config.load();

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
labels_p = hwwa.get_intermediate_dir( 'labels', conf );
events_p = hwwa.get_intermediate_dir( 'events', conf );
% use_files = { '16-J', '17-J', '18-J', '21-J', '22-J', '23-J', '25-J' };
% use_files = { '25-J' };
% use_files = { '29-J', '30-J', '31-J', '01-F' };
% use_files = { '05-F', '06-F', '07-F', '08-F', '11-F', '13-F', '14-F' };
% use_files = { '05-F', '06-F', '07-F', '08-F', '11-F' };
% use_files = { '18-F', '19-F', '21-F', '24-F', '25-F', '26-F', '27-F' };
use_files = { '18-F', '15-F', '27-F', '04-M', '24-F', '26-F', '28-F', '05-M' };

label_mats = hwwa.require_intermediate_mats( [], labels_p, use_files );

summary = labeled();

delays = 0.1:0.01:0.5;
n_delays = 3;
[g_starts, g_stops] = hwwa.bin_delays( delays, n_delays );

g_starts = [ 0.1, 0.2, 0.3, 0.4 ];
g_stops = [0.19, 0.29, 0.39, 0.5 ];

timed_trial_bin_size = 15*60; % 15 minutes

for i = 1:numel(label_mats)
  hwwa.progress( i, numel(label_mats), mfilename );
  
  labels = shared_utils.io.fload( label_mats{i} );
  unified = shared_utils.io.fload( fullfile(unified_p, labels.unified_filename) );
  events_file = shared_utils.io.fload( fullfile(events_p, labels.unified_filename) );
  
  rt = [ unified.DATA(:).reaction_time ];
  rt = rt(:);
  
  trial_start_times = events_file.event_times(:, events_file.event_key('new_trial'));
  
  labs = labels.labels;
  
  hwwa.add_day_labels( labs );
  hwwa.add_group_delay_labels( labs, g_starts, g_stops );
  hwwa.add_data_set_labels( labs );
  hwwa.add_drug_labels( labs );
  hwwa.fix_image_category_labels( labs );
  hwwa.add_time_bin_labels( labs, trial_start_times, timed_trial_bin_size );
  hwwa.split_gender_expression( labs );
  
  addcat( labs, 'monkey' );
  setcat( labs, 'monkey', lower(unified.opts.META.monkey) );
  
  addcat( summary, getcats(labs) );
  addcat( labs, getcats(summary) );
  
  prune( labs );
  
  lab = labeled( rt, labs );
  
  append( summary, lab );
end

date_dir = datestr( now, 'mmddyy' );

stats_p = fullfile( conf.PATHS.data_root, 'analyses', 'behavior', date_dir );
plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', date_dir );

shared_utils.io.require_dir( plot_p );
shared_utils.io.require_dir( stats_p );

%%

interleaved_days = { '100118', '100318' };
across_blocks_days = { '100218', '100418' };

summary_labs = getlabels( summary );

addcat( summary_labs, 'delay_manipulation' );
setcat( summary_labs, 'delay_manipulation', 'interleaved', find(summary_labs, interleaved_days) );
setcat( summary_labs, 'delay_manipulation', 'across_blocks', find(summary_labs, across_blocks_days) );

summary_dat = getdata( summary );

assert_ispair( summary_dat, summary_labs );

prune( summary_labs );

%%

bw_days = { '020119', '013019' };
lum_days = { '012919', '013119' };

addcat( summary_labs, 'image_type' );
setcat( summary_labs, 'image_type', 'b+w', find(summary_labs, bw_days) );
setcat( summary_labs, 'image_type', 'luminance-matched', find(summary_labs, lum_days) );

%%  rt

base_plot_subdir = '';

do_save = true;

is_per_monkey = true;
is_n_as_trial = true;
is_per_image_type = true;

c_is_across_days = [false];
c_is_per_gender = [true, false];
c_is_per_image_category = [true, false];
c_is_per_trial_time_bin = [false];

C = dsp3.numel_combvec( c_is_across_days, c_is_per_gender ...
  , c_is_per_image_category, c_is_per_trial_time_bin );

for idx = 1:size(C, 2)
  
plot_subdir = base_plot_subdir;

is_across_days = c_is_across_days(C(1, idx));
is_per_gender = c_is_per_gender(C(2, idx));
is_per_image_category = c_is_per_image_category(C(3, idx));
is_per_trial_time_bin = c_is_per_trial_time_bin(C(4, idx));

if ( ~is_per_image_category && ~is_per_gender )
  continue;
end

if ( is_across_days )
  base_prefix = 'across_days';
else
  base_prefix = '';
end

base_prefix = sprintf( '%s%s', base_prefix, ternary(is_n_as_trial, 'trial_n', '') );

uselabs = summary_labs';
usedat = summary_dat;

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

%%  rt, lines each day

do_save = true;
is_across_days = false;
is_per_monkey = true;
is_n_as_trial = true;
is_per_image_type = true;
is_per_gender = true;
is_per_image_category = false;

base_plot_subdir = 'first_5';

if ( is_across_days )
  base_prefix = 'across_days';
else
  base_prefix = '';
end

plot_subdir = base_plot_subdir;

base_prefix = sprintf( '%s%s', base_prefix, ternary(is_n_as_trial, 'trial_n', '') );

usedat = summary_dat;
uselabs = summary_labs';

assert_ispair( usedat, uselabs );

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true', 'go_trial', 'correct_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pl = plotlabeled.make_common();

pl.x_order = { 'appetitive', 'threat', 'scrambled' };
pl.y_lims = [0.2, 0.4];
pl.summary_func = @(x) nanmedian(x, 1);

fcats = {};
xcats = {};

pl.color_func = @hot;

if ( is_across_days )
  gcats = { 'day' };
  pcats = { 'trial_type' };
  
  if ( is_per_monkey )
    fcats{end+1} = 'monkey';
    pcats{end+1} = 'monkey';
  end
  
  func = @errorbar;
else
  gcats = { 'day' };
  pcats = { 'trial_type', 'monkey' };
  
  func = @errorbar;
end

if ( is_per_image_category )
  xcats{end+1} = 'target_image_category';
  plot_subdir = sprintf( '%s_per_image_category', plot_subdir );
end

if ( is_per_gender )
  xcats{end+1} = 'gender';
  plot_subdir = sprintf( '%s_per_gender', plot_subdir );
end

if ( is_per_image_type )
  pcats{end+1} = 'image_type';
end

spec = unique( cshorzcat(xcats, gcats, pcats, {'day'}) );

if ( is_n_as_trial )
  rtdat = usedat(mask);
  rtlabs = uselabs(mask);
else
  [rtlabs, I] = keepeach( uselabs', spec, mask );
  rtdat = rowop( usedat, I, @nanmedian );
end

[fs, axs, I] = pl.figures( func, rtdat, rtlabs, fcats, xcats, gcats, pcats );

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
    dsp3.req_savefig( fs(i), fullfile(plot_p, plot_subdir), rtlabs(I{i}), cshorzcat(pcats, fcats) ...
      , sprintf('rt_day_lines_%s', base_prefix) );
  end
end

%%  n initiated

do_save = false;

uselabs = summary_labs';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true'} ...
);

[init_labs, I] = keepeach( uselabs', {'date'}, mask );
init_dat = rowzeros( numel(I) );

for i = 1:numel(I)
  init_dat(i) = count( uselabs, 'initiated_true', I{i} );
end

pl = plotlabeled.make_common();

pl.add_points = true;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.panel_order = { 'saline' };

fcats = {};
xcats = { 'drug' };
gcats = { 'correct', 'cue_delay' };
pcats = { 'delay_manipulation' };

[fs, axs, I] = pl.figures( @bar, init_dat, init_labs, fcats, xcats, gcats, pcats );
shared_utils.plot.match_ylims( axs );
shared_utils.plot.fullscreen( fs );

% pl.bar( pcorr_dat, pcorr_labs, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(fs)
    dsp3.req_savefig( fs(i), plot_p, pcorr_labs(I{i}), cshorzcat(pcats, fcats), 'n_initiated' );
  end
end

%%  p correct, social

base_plot_subdir = '';

do_save = true;

is_per_monkey = true;
is_per_image_type = true;
is_per_run = false;

c_is_across_days = [true, false];
c_is_per_image_category = [true, false];
c_is_per_gender = [true, false];
c_is_per_trial_time_bin = [false];

C = dsp3.numel_combvec( c_is_across_days, c_is_per_image_category ...
  , c_is_per_gender, c_is_per_trial_time_bin );

for idx = 1:size(C, 2)
  
plot_subdir = base_plot_subdir; 

is_across_days = c_is_across_days(C(1, idx));
is_per_image_category = c_is_per_image_category(C(2, idx));
is_per_gender = c_is_per_gender(C(3, idx));
is_per_trial_time_bin = c_is_per_trial_time_bin(C(4, idx));

if ( ~is_per_gender && ~is_per_image_category )
  continue;
end

if ( is_across_days )
  base_prefix = 'across_days';
else
  base_prefix = '';
end

uselabs = summary_labs';

mask = fcat.mask( uselabs ...
  , @find, {'pilot_social', 'initiated_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
);

pcorr_spec = { 'trial_type', 'cue_delay' };

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

%%  p correct
%
% by day, across days

do_save = false;

uselabs = summary_labs';

% mask = fcat.mask( uselabs, @find ...
%   , {'tarantino_delay_check', 'initiated_true', 'interleaved', 'across_blocks'} );

mask = fcat.mask( uselabs ...
  , @find, {'tarantino_drug_check', 'initiated_true'} ...
  , @find, hwwa_get_first_days() ...
);

[pcorr_labs, I] = keepeach( uselabs', {'date', 'trial_type', 'cue_delay'}, mask );
pcorr_dat = rowzeros( numel(I) );

for i = 1:numel(I)
  n_corr = numel( find(uselabs, 'correct_true', I{i}) );
  n_incorr = numel( find(uselabs, 'correct_false', I{i}) );

  pcorr_dat(i) = n_corr / (n_corr + n_incorr);
end

pl = plotlabeled.make_common();
pl.fig = figure(2);

pl.add_points = true;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.group_order = { 'saline' };

fcats = {};
% xcats = { 'trial_type' };
% gcats = { 'correct', 'cue_delay' };
% pcats = { 'delay_manipulation', 'drug' };

xcats = { 'cue_delay' };
gcats = { 'drug' };
pcats = { 'trial_type', 'correct' };

[fs, axs, I] = pl.figures( @bar, pcorr_dat, pcorr_labs, fcats, xcats, gcats, pcats );
shared_utils.plot.match_ylims( axs );
shared_utils.plot.fullscreen( fs );

% pl.bar( pcorr_dat, pcorr_labs, xcats, gcats, pcats );

if ( do_save )
  for i = 1:numel(fs)
    dsp3.req_savefig( fs(i), plot_p, pcorr_labs(I{i}), cshorzcat(pcats, fcats), 'percent_correct' );
  end
end

%%  broken cues 
%
% by day, across days

do_save = true;
prefix = 'tarantino';

% uselabs = summary_labs';

% mask = fcat.mask( uselabs, @find ...
%   , {'tarantino_delay_check', 'interleaved', 'across_blocks'} );

mask = fcat.mask( uselabs ...
  , @find, hwwa_get_first_days() ...
);

[broke_labs, I] = keepeach( uselabs', {'date', 'trial_type', 'cue_delay'}, mask );

broke_ind = rowzeros( size(broke_labs, 1), 1 );

for i = 1:numel(I)
  broke_cue_fix_ind = find( uselabs, 'broke_cue_fixation', I{i} );
  did_not_break_ind = setdiff( I{i}, find(uselabs, 'no_fixation') );
  
  broke_ind(i) = numel(broke_cue_fix_ind) / numel(did_not_break_ind);
end

pl = plotlabeled.make_common();

pl.add_points = true;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.group_order = { 'saline' };
pl.y_lims = [0, 1];

% xcats = { 'trial_type' };
% gcats = { 'correct', 'cue_delay' };
% pcats = { 'drug' };

xcats = { 'cue_delay' };
gcats = { 'drug' };
pcats = { 'trial_type', 'correct' };

pl.bar( broke_ind, broke_labs, xcats, gcats, pcats );

shared_utils.plot.fullscreen( gcf );

if ( do_save )
  dsp3.req_savefig( gcf, plot_p, broke_labs', cshorzcat(pcats), sprintf('%s_abort_rate', prefix) );
end

%%  d' / criterion

delay_group = 'cue_delay';

specificity = { 'date', 'trial_type', delay_group };

uselabs = summary_labs';
usedat = summary_dat;

assert_ispair( usedat, uselabs );

% mask = fcat.mask( uselabs, @find ...
%   , {'tarantino_delay_check', 'interleaved', 'across_blocks'} );

mask = fcat.mask( uselabs ...
  , @find, {'tarantino_drug_check'} ...
  , @find, hwwa_get_first_days() ...
);

[pcorr_labs, I] = keepeach( uselabs', specificity, mask );
data = zeros( numel(I), 1 );

for i = 1:numel(I)  
  n_correct = numel( find(uselabs, 'no_errors', I{i}) );
  n_init = numel( find(uselabs, 'initiated_true', I{i}) );
  
  data(i) = n_correct / n_init;
end

I = findall( pcorr_labs, {'drug'} );

d_prime_labs = fcat.like( pcorr_labs );
d_prime_data = zeros( size(data, 1)/2, 1 );
C_data = zeros( size(d_prime_data) );
B_data = zeros( size(C_data) );
C_prime = zeros( size(C_data) );
stp = 1;

for i = 1:numel(I)
  go_ind = find( pcorr_labs, 'go_trial', I{i} );
  nogo_ind = find( pcorr_labs, 'nogo_trial', I{i} );
  
  go_percent = data(go_ind);
  nogo_percent = (1 - data(nogo_ind));
  
  z_hit = norminv( go_percent, mean(go_percent), std(go_percent) );
  z_fa = norminv( nogo_percent, mean(nogo_percent), std(nogo_percent) );
  
  d_prime = z_hit - z_fa;
  
  c_dat = -0.5 .* (z_hit + z_fa);
  
  b = exp( c_dat .* d_prime );
  c_prime = c_dat ./ d_prime;
  
  assign_seq = stp:stp+numel(d_prime)-1;
  
  d_prime_data(assign_seq) = d_prime;
  C_data(assign_seq) = c_dat;
  B_data(assign_seq) = b;
  C_prime(assign_seq) = c_prime;
  
  stp = stp + numel(d_prime);
  
  append( d_prime_labs, pcorr_labs, go_ind );
end

%%

kind = 'd_prime';

do_save = true;

switch ( kind )
  case 'd_prime'
    pltdat = d_prime_data;
  case 'criterion'
    pltdat = C_data;
  otherwise
    error( 'unrecognized kind "%s".', kind );
end

pltlabs = d_prime_labs';

pl = plotlabeled.make_common();

pl.add_points = true;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.panel_order = { 'saline' };

xcats = { 'trial_type' };
gcats = { 'correct' };
pcats = { 'delay_manipulation', 'drug' };

pl.bar( pltdat, pltlabs, xcats, gcats, pcats );

shared_utils.plot.fullscreen( gcf );

if ( do_save )
  dsp3.req_savefig( gcf, plot_p, pltlabs', cshorzcat(pcats, gcats), kind );
end


