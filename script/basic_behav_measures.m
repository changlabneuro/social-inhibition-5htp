conf = hwwa.config.load();

conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

unified_p = hwwa.get_intermediate_dir( 'unified', conf );
labels_p = hwwa.get_intermediate_dir( 'labels', conf );

label_mats = hwwa.require_intermediate_mats( [], labels_p, [] );

summary = labeled();

delays = 0.1:0.01:0.5;
n_delays = 3;
[g_starts, g_stops] = hwwa.bin_delays( delays, n_delays );

% g_starts = delay_starts;
% g_stops = delay_stops;

g_starts = [ 0.1, 0.2, 0.3, 0.4 ];
g_stops = [0.19, 0.29, 0.39, 0.5 ];

for i = 1:numel(label_mats)
  hwwa.progress( i, numel(label_mats), mfilename );
  
  labels = shared_utils.io.fload( label_mats{i} );
  unified = shared_utils.io.fload( fullfile(unified_p, labels.unified_filename) );
  
  rt = [ unified.DATA(:).reaction_time ];
  rt = rt(:);
  
  labs = labels.labels;
  
  hwwa.add_day_labels( labs );
  hwwa.add_group_delay_labels( labs, g_starts, g_stops );
  hwwa.add_data_set_labels( labs );
  
  prune( labs );
  
  lab = labeled( rt, labs );
  
  append( summary, lab );
end

date_dir = datestr( now, 'mmddyy' );

stats_p = fullfile( conf.PATHS.data_root, 'analyses', 'behavior', date_dir );
save_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', date_dir );
shared_utils.io.require_dir( save_p );
shared_utils.io.require_dir( stats_p );

%%

interleaved_days = { '100118', '100318' };
across_blocks_days = { '100218', '100418' };

summary_labs = getlabels( summary );

addcat( summary_labs, 'delay_manipulation' );
setcat( summary_labs, 'delay_manipulation', 'interleaved', find(summary_labs, interleaved_days) );
setcat( summary_labs, 'delay_manipulation', 'across_blocks', find(summary_labs, across_blocks_days) );

%%  p correct, saline and 5-htp

do_save = false;

pcorr = prune( only(summary', 'tarantino_delay_check') );

labs = getlabels( pcorr );

I = findall( labs, 'date' );

dat = [];
plabs = fcat();

for i = 1:numel(I)
  
  subset = labs(I{i});
  
  only( subset, 'initiated_true' );
  
  [y, corr_i] = keepeach( subset', {'trial_type'} );
  
  for j = 1:numel(corr_i)
    
    ind = corr_i{j};
    
    n_corr = numel( intersect(ind, find(subset, 'correct_true')) );
    n_incorr = numel( intersect(ind, find(subset, 'correct_false')) );
    
    n_total = n_corr + n_incorr;
    
    p_corr = n_corr / n_total;
    
    dat = [ dat; p_corr ];
  end
  
  append( plabs, y );
end

colors = containers.Map();
colors('050818') = 'r';
colors('050918') = 'c';
colors('051018') = 'b';
colors('051118') = 'g';
colors('051418') = 'k';

pl = plotlabeled();
to_plt = labeled( dat, plabs );
pl.add_points = true;
pl.error_func = @plotlabeled.sem;
pl.points_are = { 'day' };
pl.marker_size = 18;
pl.marker_type = '*';
pl.points_color_map = colors;

bar( pl, to_plt, 'trial_type', 'drug', 'correct' );

fname = joincat( prune(to_plt'), {'trial_type', 'drug', 'correct'} );

shared_utils.plot.save_fig( gcf(), fullfile(save_p, fname), {'epsc', 'png', 'fig'} );

%% p corr grouped delays

do_save = true;

labs = hwwa.add_day_labels( getlabels(summary) );

only( labs, {'5-htp', 'saline'} );
only( labs, {'050918', '051418', '051118', '050818', '051018'} );

prune( labs );

g_starts = [0.1, 0.2, 0.3, 0.4];
g_stops = [0.19, 0.29, 0.39, 0.5];

hwwa.add_group_delay_labels( labs, g_starts, g_stops );

I = findall( labs, {'date', 'gcue_delay'} );

dat = [];
plabs = fcat();

for i = 1:numel(I)
  subset = labs(I{i});
  
  only( subset, 'initiated_true' );
  
  [y, corr_i] = keepeach( subset', {'trial_type'} );
  
  for j = 1:numel(corr_i)
    
    ind = corr_i{j};
    
    n_corr = numel( intersect(ind, find(subset, 'correct_true')) );
    n_incorr = numel( intersect(ind, find(subset, 'correct_false')) );
    
    n_total = n_corr + n_incorr;
    
    p_corr = n_corr / n_total;
    
    dat = [ dat; p_corr ];
  end
  
  append( plabs, y );
end

delay_strs = combs( plabs, 'gcue_delay' );
stop_inds = cellfun( @(x) strfind(x, '-'), delay_strs );

delays = arrayfun( @(x, y) str2double(x{1}(numel('grouped_delay__')+1:y-1)), delay_strs, stop_inds );
[~, sorted_ind] = sort( delays );

pl = plotlabeled();
to_plt = labeled( dat, plabs );

pl.group_order = delay_strs( sorted_ind );
pl.error_func = @plotlabeled.sem;
pl.one_legend = true;

bar( pl, to_plt, 'trial_type', 'gcue_delay', 'drug' );

fname = joincat( prune(to_plt'), {'trial_type', 'drug', 'correct'} );


%%  rt, per cue delay

ind = intersect( find(summary.data > 0), find(summary, {'no_errors', 'wrong_go_nogo'}) );

rt = prune( summary(ind) );

collapsecat( rt, 'cue_delay' );
only( rt, 'ro1' );

rt = eachindex( rt', {'date', 'gcue_delay', 'trial_outcome'}, @rownanmean );

% delay_strs = rt('cue_delay');
% delays = shared_utils.container.cat_parse_double( 'delay__', delay_strs );
% [~, sort_i] = sort( delays );

delay_strs = combs( rt, 'gcue_delay' );
stop_inds = cellfun( @(x) strfind(x, '-'), delay_strs );

delays = arrayfun( @(x, y) str2double(x{1}(numel('grouped_delay__')+1:y-1)), delay_strs, stop_inds );
[~, sorted_ind] = sort( delays );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;
pl.group_order = delay_strs(sorted_ind);
pl.bar( rt, 'drug', 'gcue_delay', 'trial_outcome' );

fname = strjoin( incat(rt, {'date', 'drug'}), '_' );
fname = sprintf( 'rt_%s', fname );

shared_utils.plot.save_fig( gcf, fullfile(save_p, fname), {'fig', 'epsc', 'png'}, true );

%%  rt anova

full_stats_p = fullfile( stats_p, 'rt' );
shared_utils.io.require_dir( full_stats_p );

grps = { 'cue_delay', 'drug' };

hwwa.write_anova_tables( rt, grps, full_stats_p, '' );

%%  means + devs

grps = { 'cue_delay' };
hwwa.write_anova_summary_tables( rt, grps, full_stats_p );

%%  rt, across cue delays

ind = intersect( find(summary.data > 0), find(summary, {'no_errors', 'wrong_go_nogo'}) );

rt = summary( ind );
collapsecat( rt, 'cue_delay' );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;
pl.group_order = { 'delay__0.01' };
pl.bar( rt, 'drug', 'trial_outcome', 'cue_delay' );

%%  n initiated

specificity = { 'date' };

[labs, I, C] = keepeach( getlabels(summary), specificity );
data = zeros( numel(I), 1 );
initiated_trials = find( summary, {'no_errors', 'wrong_go_nogo'} );
for i = 1:numel(I)
  data(i) = numel( intersect(I{i}, initiated_trials) );
end

n_initiated = labeled( data, labs );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;

pl.bar( n_initiated, 'drug', 'trial_type', 'trial_outcome' );

shared_utils.plot.save_fig( gcf, fullfile(save_p, 'n_initiated_all_trials') ...
  , {'fig', 'epsc', 'png'}, true );

saline_data = data( find(n_initiated, 'saline') );
serotonin_data = data( find(n_initiated, '5-htp') );

sal_means = mean( saline_data );
sal_devs = std( saline_data );
htp_means = mean( serotonin_data );
htp_devs = std( serotonin_data );

[~, p] = ttest( saline_data, serotonin_data );

%%  broken cues 

do_save = false;

labels = getlabels( summary );

% prune( only(labels, 'ro1') );
prune( only(labels, hwwa_get_cron_days()) );

[y, I] = keepeach( labels', {'date', 'trial_type', 'gcue_delay'} );
data = zeros( size(y, 1), 1 );

for i = 1:numel(I)
  initiated_ind = intersect( I{i}, find(labels, 'initiated_true') );
  broke_cue_fix_ind = intersect( I{i}, find(labels, 'broke_cue_fixation') );
  
  did_not_break_ind = setdiff( I{i}, find(labels, 'no_fixation') );
  
%   data(i) = numel(broke_cue_fix_ind) / numel(initiated_ind);
  data(i) = numel(broke_cue_fix_ind) / numel(did_not_break_ind);
end

p_broke_cue = labeled( data, y );

delay_strs = combs( p_broke_cue, 'gcue_delay' );
% stop_inds = cellfun( @(x) strfind(x, '-'), delay_strs );

if ( false )
  delays = arrayfun( @(x, y) str2double(x{1}(numel('grouped_delay__')+1:y-1)), delay_strs, stop_inds );
  [~, sorted_ind] = sort( delays );
end

pl = plotlabeled.make_common();
pl.error_func = @plotlabeled.sem;
% pl.panel_order = delay_strs(sorted_ind);

pl.bar( p_broke_cue, 'cue_delay', 'drug', 'trial_type' );

fname = strjoin( incat(p_broke_cue, {'drug', 'trial_type', 'cue_delay'}), '_' );
fname = sprintf( 'p_broke_cue_fix_%s', fname );

if ( do_save )
  shared_utils.plot.save_fig( gcf, fullfile(save_p, fname), {'fig', 'epsc', 'png'}, true );
end

%%  broke cue anova

full_stats_p = fullfile( stats_p, 'broke_cue' );
shared_utils.io.require_dir( full_stats_p );

grps = { 'drug', 'cue_delay' };

hwwa.write_anova_tables( p_broke_cue, grps, full_stats_p, 'no_cue_delay' );

%%  means + devs

grps = { 'cue_delay' };
hwwa.write_anova_summary_tables( p_broke_cue, grps, full_stats_p );

%%  percent correct

specificity = { 'date', 'trial_type', 'cue_delay' };

[labs, I, C] = keepeach( getlabels(summary), specificity );
data = zeros( numel(I), 1 );
correct_trials = find( summary, 'no_errors' );
initiated_trials = find( summary, {'no_errors', 'wrong_go_nogo'} );
for i = 1:numel(I)
  n_correct = numel( intersect(I{i}, correct_trials) );
  n_init = numel( intersect(I{i}, initiated_trials) );
  
  data(i) = n_correct / n_init;
end

p_correct = labeled( data, labs );

% only( p_correct, {'delay__0.1', 'delay__0.5'} );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;
pl.one_legend = true;
pl.group_order = { 'delay__0.01' };

pl.bar( p_correct, 'drug', 'cue_delay', 'trial_type' );

shared_utils.plot.save_fig( gcf, fullfile(save_p, 'p_correct') ...
  , {'fig', 'epsc', 'png'}, true );

%%

full_stats_p = fullfile( stats_p, 'p_correct' );
shared_utils.io.require_dir( full_stats_p );

grps = { 'trial_type' };

hwwa.write_anova_tables( p_correct, grps, full_stats_p, '' );

%%  means + devs

grps = { 'trial_type' };
hwwa.write_anova_summary_tables( p_correct, grps, full_stats_p );

%% number correct

specificity = { 'date', 'trial_type', 'cue_delay' };

[labs, I, C] = keepeach( getlabels(summary), specificity );
data = zeros( numel(I), 1 );
correct_trials = find( summary, 'no_errors' );
for i = 1:numel(I)
  n_correct = numel( intersect(I{i}, correct_trials) );
  data(i) = n_correct;
end

n_correct = labeled( data, labs );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;
pl.one_legend = true;
pl.group_order = { 'delay__0.01' };

pl.bar( n_correct, 'drug', 'cue_delay', 'trial_type' );

fname = strjoin( incat(n_correct, {'drug', 'trial_type'}), '_' );
fname = sprintf( 'ncorrect_%s', fname );

shared_utils.plot.save_fig( gcf, fullfile(save_p, fname) ...
  , {'fig', 'epsc', 'png'}, true );

%%  d'

do_save = true;

delay_group = 'cue_delay';

specificity = { 'date', 'trial_type', delay_group };

to_dprime = summary';

to_keep = find( ~trueat(to_dprime, find(to_dprime, 'no drug')) );
keep( to_dprime, to_keep );

% prune( only(to_dprime, 'ro1') );

ind = find( ~trueat(to_dprime, find(to_dprime, 'ro1')) );
keep( to_dprime, ind );

% collapsecat( to_dprime, 'drug' );

all_labs = getlabels( to_dprime );
[labs, I, C] = keepeach( getlabels(to_dprime), specificity );
data = zeros( numel(I), 1 );

correct_trials = find( to_dprime, 'no_errors' );
initiated_trials = find( to_dprime, 'initiated_true' );

for i = 1:numel(I)
  n_correct = numel( intersect(I{i}, correct_trials) );
  n_init = numel( intersect(I{i}, initiated_trials) );
  
  data(i) = n_correct / n_init;
end

p_correct = labeled( data, labs );

z_each = { 'drug' };

[~, I] = keepeach( getlabels(p_correct), z_each );

d_prime_labs = fcat.like( labs );
d_prime_data = zeros( size(data, 1)/2, 1 );
C_data = zeros( size(d_prime_data) );
B_data = zeros( size(C_data) );
C_prime = zeros( size(C_data) );
stp = 1;

for i = 1:numel(I)
  go_ind = intersect( I{i}, find(p_correct, 'go_trial') );
  nogo_ind = intersect( I{i}, find(p_correct, 'nogo_trial') );
  
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
  
  append( d_prime_labs, labs(go_ind) );
end

%%

plot_type = 'c';

switch ( plot_type )
  case 'd_prime'
    plt_data = d_prime_data;
  case 'c'
    plt_data = C_data;
  case 'c_prime'
    plt_data = C_prime;
  case 'b'
    plt_data = B_data;
  otherwise
    error( 'Unrecognized plot type "%s".', plot_type );
end

plt_sdt = labeled( plt_data, d_prime_labs );

pl = plotlabeled();
pl.error_func = @plotlabeled.sem;
pl.one_legend = true;
% pl.group_order = { 'grouped_delay__0.1-0.22', 'grouped_delay__0.23-0.35' };
pl.group_order = { 'delay__0.01', 'delay__0.1' };
pl.x_order = pl.group_order;
% pl.y_lims = [0.55, 0.82];

% pl.bar( plt_sdt, 'drug', delay_group, 'trial_type' );
pl.bar( plt_sdt, delay_group, 'drug', 'trial_type' );

means = eachindex( plt_sdt', 'drug', @rowmean );
devs = each( plt_sdt', 'drug', @(x) std(x, [], 1) );

m_t = table( Container.from(means), 'drug' );
d_t = table( Container.from(devs), 'drug' );

saline_data = plt_data( find(d_prime_labs, 'saline') );
serotonin_data = plt_data( find(d_prime_labs, '5-htp') );

sal_mean = mean( saline_data );
htp_mean = mean( serotonin_data );
sal_dev = std( saline_data );
htp_dev = std( serotonin_data );
%%
mean_each = 'drug';
I = findall( d_prime_labs, mean_each );

for i = 1:numel(I)
  subset_one_drug = plt_data(I{i});
  subset_labs = d_prime_labs(I{i});
  
  [y, mean_i] = keepeach( subset_labs, delay_group );
  means = rowop( subset_one_drug, mean_i, @(x) mean(x, 1) );
  devs = rowop( subset_one_drug, mean_i, @(x) std(x, [], 1) );
  
  delays = combs( subset_labs, delay_group );
  [~, sort_ind] = sort( shared_utils.container.cat_parse_double('delay__', delays) );
  
  m = means(sort_ind);
  n = (1:numel(sort_ind))';
  eval_n = 1:0.1:numel(sort_ind);
  
  p = polyfit( n, m(:), 1 );
  
  hold on;
  plot( eval_n, polyval(p, eval_n) );
end
%%
try
  [~, p, ~, stats] = ttest( saline_data, serotonin_data );
  stats = struct2table( stats );
  stats(:, 'p') = { p };
catch err
  warning( err.message );
end

if ( do_save )

  shared_utils.plot.save_fig( gcf, fullfile(save_p, plot_type), {'fig', 'epsc', 'png'}, true );

  writetable( m_t, fullfile(stats_p, ['means_', plot_type, '.csv']), 'WriteRowNames', true );
  writetable( d_t, fullfile(stats_p, ['devs_', plot_type, '.csv']), 'WriteRowNames', true );
  writetable( stats, fullfile(stats_p, ['stats_', plot_type, '.csv']), 'WriteRowNames', true );

end



