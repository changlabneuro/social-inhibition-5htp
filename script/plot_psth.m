conf = hwwa.config.load();

label_p = hwwa.get_intermediate_dir( 'labels' );
psth_p = hwwa.get_intermediate_dir( 'psth' );

psth_mats = hwwa.require_intermediate_mats( psth_p );

evt = 'go_target_onset';
% evt = 'go_target_acquired';
% evt = 'go_nogo_cue_onset';
% evt = 'reward_onset';
% evt = 'go_target_offset';

do_norm = false;
% baseline_evt = 'go_nogo_cue_onset';
baseline_evt = 'fixation_on';
pre_baseline = -0.3;
post_baseline = 0;
norm_func = @minus;
% norm_func = @rdivide;

delays = 0.1:0.01:0.5;
n_delays = 3;
[g_starts, g_stops] = hwwa.bin_delays( delays, n_delays );

do_z = false;
z_within = { 'id' };

do_thresh = false;
thresh_within = { 'id' };
ndevs = 3;

psth_data = [];
psth_labs = fcat();

p_removed = [];
p_removed_labs = fcat();

for i = 1:numel(psth_mats)
  hwwa.progress( i, numel(psth_mats) );
  
  psth_file = shared_utils.io.fload( psth_mats{i} );
  un_filename = psth_file.unified_filename;
  labs_file = shared_utils.io.fload( fullfile(label_p, un_filename) );
  
  psth = psth_file.psth(evt);
  psth_d = psth.data;
  
  if ( do_z )
    psth_d = hwwa.zscore( psth_d, psth.labels, z_within );
  end
  
  if ( do_thresh )
    [psth_d, inds, ps] = hwwa.stdthresh( psth_d, psth.labels, thresh_within, ndevs );
    ps = cellfun( @(x) perc(x), ps );
  end
  
  if ( do_norm )
    baseline_psth = psth_file.psth(baseline_evt);
    baseline_t = baseline_psth.time > pre_baseline & baseline_psth.time <= post_baseline;
    baseline_d = nanmean( baseline_psth.data(:, baseline_t), 2 );
    
    psth_d = norm_func( psth_d, baseline_d );
  end
  
  hwwa.merge_unit_labs( labs_file.labels, psth.labels );
  hwwa.add_day_labels( labs_file.labels );
  hwwa.add_data_set_labels( labs_file.labels );
  hwwa.add_broke_cue_labels( labs_file.labels );
  hwwa.add_group_delay_labels( labs_file.labels, g_starts, g_stops );
  
  psth_data = [ psth_data; psth_d ];
  append( psth_labs, labs_file.labels );
  
  t = psth.time;
  psth_t = t;
end

hwwa.add_region_labels( psth_labs, 'WB' );

prune( psth_labs );

save_p = fullfile( conf.PATHS.data_root, 'plots', 'spike', datestr(now, 'mmddyy'), evt );

%%

pfclabs = prune( only(psth_labs', {'ro1', 'dlpfc'}) );

keepeach( pfclabs, 'id' );

[y, I] = keepeach( pfclabs', {'drug', 'region'} );

counts = Container( cellfun(@numel, I), SparseLabels.from_fcat(y) );

table( counts, 'drug', 'region' )

[c, C] = countrows( pfclabs, {'drug', 'region'} );

%%

[y, I] = keepeach( psth_labs', 'id' );
grand_means = nanmean( rownanmean(psth_data, I), 2 );

figure(1); clf();

plt = Container.from( labeled(grand_means, y) );

plt = only( plt, {'dlpfc', 'ro1'} );

pl = ContainerPlotter();
pl.x_label = 'Grand mean sp/s';

pl.hist( plt, 20, [], 'drug' );

%%

do_save = false;

pl = plotlabeled();
pl.error_func = @plotlabeled.nansem;
pl.summary_func = @plotlabeled.nanmean;
pl.smooth_func = @(x) smooth(x, 10);
pl.add_smoothing = true;
pl.add_errors = true;
pl.one_legend = true;
pl.x = t;
pl.fig = figure(2);
% pl.group_order = { 'grouped_delay__0.1-0.22', 'grouped_delay__0.23-0.35' };
pl.group_order = { '5-htp' };
pl.panel_order = { 'go', 'nogo', 'correct_false' };

% panels_are = { 'id', 'trial_type' };
% lines_are = { 'trial_outcome' };
% lines_are = { 'broke_cue' };
% panels_are = { 'id' };

% panels_are = { 'id', 'trial_type', 'drug', 'region' };
% lines_are = { 'correct' };

panels_are = { 'id', 'trial_type', 'region' };
lines_are = { 'drug', 'date', 'correct' };

fnames_are = { 'id', 'drug' };
figs_are = { 'id' };

specificity = unique( [lines_are, panels_are, 'id'] );

psth_labeled = labeled( psth_data, psth_labs );

prune( only(psth_labeled, {'dlpfc', 'initiated_true', 'ro1'}) );

only( psth_labeled, '0514-WB16-MS08-1' ); % classic dlpfc
% only( psth_labeled, '0510-WB13-MS05-7' ); % go action
% only( psth_labeled, '0508-WB13-MS05-1' ); % cue response
% only( psth_labeled, '0510-WB11-MS03-1' ); % poss. nogo
% only( psth_labeled, '0508-WB14-MS06-2' ); % poss. nogo
% only( psth_labeled, '0509-WB13-MS05-1' ); % both outcomes
% only( psth_labeled, '0508-WB11-MS03-2' );

% only( psth_labeled, '0511-WB09-2' );  % classic dlpfc (closest)
% only( psth_labeled, '0509-WB13-6' );  % decrease due to go; nogo correct higher than go incorrect
% only( psth_labeled, '0509-WB15-8' );  % decrease due to go; nogo correct higher than go incorrect
% only( psth_labeled, '0508-WB15-8' ); % planning for go action
% only( psth_labeled, '0508-WB15-8' ); % planning for go action
% only( psth_labeled, '0508-WB11-4' ); % correct go (pre-response)



% ind = find( psth_labeled, '051018' );
% to_keep = setdiff( 1:size(psth_labeled, 1), ind );
% keep( psth_labeled, to_keep );

% good_site_510 = only( psth_labeled', '0510-WB13-6' );
% keep( psth_labeled, to_keep );
% append( psth_labeled, good_site_510 );

% meaned = Container.from( psth_labeled );
% meaned = each1d( meaned, union(specificity, {'date', 'correct'}), @rowops.nanmean );
% meaned = collapse( meaned, {'error', 'trial_outcome'} );
% meaned = meaned({'correct_true'}) - meaned({'correct_false'});
% psth_labeled = labeled.from( meaned );

% collapsecat( psth_labeled, {'id'} );

I = findall( psth_labeled, figs_are );

for i = 1:numel(I)

plt = prune( psth_labeled(I{i}) );

ind = find( plt, 'no_choice' );

keep( plt, setdiff(1:size(plt, 1), ind) );

axs = pl.lines( plt, lines_are, panels_are );

arrayfun( @(x) set(x, 'nextplot', 'add'), axs );
arrayfun( @(x) plot(x, [0; 0], get(x, 'ylim'), 'k--'), axs );
arrayfun( @(x) xlabel(x, sprintf('Time (s) from %s', strrep(evt, '_', ' '))), axs );

fname = joincat( plt, fnames_are );
fname = sprintf( 'psth_%s', fname );
full_fname = fullfile( save_p, fname );

if ( do_save )
  shared_utils.io.require_dir( save_p );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

end

xlim( axs, [-0.5, 1] );

%%  plot means across time

do_save = false;

t_start = -0.1;
t_stop = 0.3;
t_ind = psth_t >= t_start & psth_t <= t_stop;

pl = plotlabeled();
pl.error_func = @plotlabeled.nansem;
pl.summary_func = @plotlabeled.nanmean;
pl.add_errors = true;
pl.one_legend = true;
pl.x_tick_rotation = 0;
% pl.group_order = { 'grouped_delay__0.1-0.22', 'grouped_delay__0.23-0.35' };
pl.group_order = { '5-htp', 'saline' };

specificity = unique( [lines_are, panels_are, 'id'] );

plt_data = nanmean( psth_data(:, t_ind), 2 );

plt = labeled( plt_data, psth_labs );
only( plt, 'ro1' );
only( plt, {'initiated_true', 'dlpfc'} );

collapsecat( plt, 'id' );

prune( plt );

figs_are = { 'id' };
x_are = { 'correct' };
% groups_are = { 'gcue_delay' };
groups_are = { 'drug' };
panels_are = { 'trial_type' };

I = findall( plt, figs_are );

for i = 1:numel(I)
  subset_plt = plt(I{i});
  
  bar( pl, subset_plt, x_are, groups_are, panels_are );
end


%%  plot targ onset

pl = ContainerPlotter();
pl.error_function = @(x, y) plotlabeled.nansem(x);
pl.summary_function = @(x, y) nanmean( x, 1 );
pl.smooth_function = @(x) smooth(x, 7);
pl.add_smoothing = true;
pl.x = t;
pl.marker_size = 8;
pl.add_ribbon = true;
pl.compare_series = true;
pl.p_correct_type = 'fdr';

panels_are = { 'id' };
lines_are = { 'correct' };

specificity = unique( [lines_are, panels_are, 'id'] );

psth_labeled = labeled( psth_data, psth_labs );

I = findall( psth_labeled, 'id' );

for i = 1:numel(I)

plt = prune( psth_labeled(I{i}) );

ind = find( plt, 'no_choice' );

keep( plt, setdiff(1:size(plt, 1), ind) );

axs = pl.plot( Container.from(plt), lines_are, panels_are );

arrayfun( @(x) set(x, 'nextplot', 'add'), axs );
arrayfun( @(x) plot(x, [0; 0], get(x, 'ylim'), 'k--'), axs );
arrayfun( @(x) xlabel(x, sprintf('Time (s) from %s', strrep(evt, '_', ' '))), axs );

fname = strjoin( incat(plt, specificity), '_' );
fname = sprintf( 'psth_%s', fname );
full_fname = fullfile( save_p, fname );

shared_utils.io.require_dir( save_p );
shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );

end

%%  cue onset | target onset

do_save = true;

pl = ContainerPlotter();
pl.error_function = @(x, y) plotlabeled.nansem(x);
pl.summary_function = @(x, y) nanmean( x, 1 );
pl.smooth_function = @(x) smooth(x, 7);
pl.add_smoothing = true;
pl.x = t;
pl.marker_size = 8;
pl.add_ribbon = true;
pl.compare_series = true;
pl.p_correct_type = 'fdr';

panels_are = { 'id', 'trial_type' };
lines_are = { 'trial_outcome' };

specificity = unique( [lines_are, panels_are, 'id'] );

psth_labeled = labeled( psth_data, psth_labs );

I = findall( psth_labeled, 'id' );

for i = 1:numel(I)

plt = prune( psth_labeled(I{i}) );

ind = find( plt, 'no_choice' );

keep( plt, setdiff(1:size(plt, 1), ind) );

axs = pl.plot( Container.from(plt), lines_are, panels_are );

arrayfun( @(x) set(x, 'nextplot', 'add'), axs );
arrayfun( @(x) plot(x, [0; 0], get(x, 'ylim'), 'k--'), axs );
arrayfun( @(x) xlabel(x, sprintf('Time (s) from %s', strrep(evt, '_', ' '))), axs );

fname = strjoin( incat(plt, specificity), '_' );
fname = sprintf( 'psth_%s', fname );
full_fname = fullfile( save_p, fname );

if ( do_save )
  shared_utils.io.require_dir( save_p );
  shared_utils.plot.save_fig( gcf, full_fname, {'epsc', 'png', 'fig'}, true );
end

end