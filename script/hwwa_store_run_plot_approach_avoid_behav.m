conf = hwwa.config.load();
behav_outs = hwwa_load_approach_avoid_behavioral_data( conf );
saccade_outs = hwwa_load_approach_avoid_saccade_info( conf );

go_targ_aligned = shared_utils.io.fload( ...
  fullfile(hwwa.gid('processed', conf), 'behavior', 'go_target_aligned.mat') );

%%  Or reload

roi_outs = hwwa_load_approach_avoid_behavior();
assert( rmcat(behav_outs.labels', 'is_outlier') == roi_outs.labels );

%%  event info

event_info = hwwa_load_events( 'files_containing', hwwa.approach_avoid_files() );

%%

cue_on_pup = behav_outs.cue_on_aligned.pupil;
cue_on_t = behav_outs.cue_on_aligned.time;

cue_on_norm_pup = hwwa.time_normalize( cue_on_pup, cue_on_t, [-150, 0] );
cue_on_norm_pup(cue_on_norm_pup == 0) = nan;
cue_on_slice = cue_on_norm_pup(:, cue_on_t >= 0 & cue_on_t <= 800 );
cue_on_max_diff = min( cue_on_slice, [], 2 );


%%

per_rt_quantile = false;
session_level_rt_quantile = true;
num_rt_quantiles = 2;

find_non_outliers = @(labels, varargin) find( labels, 'is_outlier__false', varargin{:} );

use_labs = behav_outs.labels';
hwwa.decompose_social_scrambled( use_labs );
time = behav_outs.run_relative_start_times;

rt = saccade_outs.rt.saccade_based_rt;

slide_each = { 'unified_filename', 'target_image_category', 'trial_type', 'scrambled_type' };

prop_mask = find_non_outliers( use_labs, hwwa.get_approach_avoid_mask(use_labs) );
rt_mask = hwwa.find_correct_go_incorrect_nogo( use_labs, prop_mask );

quant_cat = 'rt-quantile';
make_quant_labs = @(quants) arrayfun(@(x) sprintf('%s__%d', quant_cat, x), quants, 'un', 0);

if ( per_rt_quantile && ~session_level_rt_quantile )
  quants_of = 'unified_filename';
  quants_each = setdiff( slide_each, quants_of );
  [quants, quant_I] = dsp3.quantiles_each( rt, use_labs', num_rt_quantiles ...
    , quants_each, quants_of, rt_mask );
  addsetcat( use_labs, quant_cat, make_quant_labs(quants) );
  slide_each = union( slide_each, {quant_cat} );
end

% Average percent correct + rt per session.
[pcorr, pcorr_labels] = hwwa.percent_correct( use_labs', slide_each, prop_mask );
[mean_rt_labels, mean_rt_I] = keepeach( use_labs', slide_each, rt_mask );
mean_rt = bfw.row_nanmean( rt, mean_rt_I );

if ( per_rt_quantile && session_level_rt_quantile )
  quants_of = 'unified_filename';
  [quants, quant_I] = dsp3.quantiles_each( mean_rt, mean_rt_labels', num_rt_quantiles ...
    , setdiff(slide_each, quants_of), quants_of );
  addsetcat( mean_rt_labels, quant_cat, make_quant_labs(quants) );
  slide_each = union( slide_each, {quant_cat} );
end

%%

base_pup = nanmean( cue_on_pup(:, cue_on_t >= -150 & cue_on_t <= 0), 2 );
inds = findall( use_labs, 'unified_filename' );
quants = bfw.row_nanmedian( base_pup, inds );

set_a = {};
set_b = {};

for i = 1:numel(inds)
  set_a{end+1, 1} = inds{i}(base_pup(inds{i}) < quants(i));
  set_b{end+1, 1} = inds{i}(base_pup(inds{i}) >= quants(i));
end

pre = vertcat( set_a{:} );
post = vertcat( set_b{:} );

ax = subplot( 1, 2, 1 );
hist( ax, cue_on_max_diff(pre), 1e3 );

ax2 = subplot( 1, 2, 2 );
hist( ax2, cue_on_max_diff(post), 1e3 );

%%  Proportion of trials remained fixated on cue / distance from cue center.

start_event_name = 'go_nogo_cue_onset';
stop_event_name = 'go_target_offset';

start_event_cols = arrayfun( ...
  @(x) x.event_times(:, x.event_key(start_event_name)), roi_outs.events, 'un', 0 );
stop_event_cols = arrayfun( ...
  @(x) x.event_times(:, x.event_key(stop_event_name)), roi_outs.events, 'un', 0 );

offs = cellfun( @(start, stop) stop - start, start_event_cols, stop_event_cols, 'un', 0 );
edf_offs = cat_expanded( 1, cellfun(@(x) floor(x * 1000), offs, 'un', 0) );

assert( numel(edf_offs) == size(behav_outs.labels, 1) );

cue_on_x = behav_outs.cue_on_aligned.x;
cue_on_y = behav_outs.cue_on_aligned.y;
test_rects = cat_expanded( 1, arrayfun(@(x) x.nogo_cue, roi_outs.rois, 'un', 0) );

remained_on_cue = false( size(cue_on_x, 1), 1 );
dist_from_cue_center = zeros( size(remained_on_cue) );

for i = 1:size(cue_on_x, 1)
  start = 1;
  stop = min( edf_offs(i), size(cue_on_x, 2) );
  cue_on_p = [ cue_on_x(i, start:stop)', cue_on_y(i, start:stop)' ];
  test_rect = test_rects(i, :);
  ib = shared_utils.rect.inside( test_rect, cue_on_p );
  remained_on_cue(i) = sum( ib ) == numel( ib );
  
  rect_center = shared_utils.rect.center( test_rect );
  from_center = cue_on_p - rect_center;
  dist_from_cue_center(i) = nanmean( vecnorm(from_center, 2, 2) );
end

%%  Stats prop trials remained fixed on cue

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_nogo(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m, true) ...
  , @(m) hwwa.find_correct(l, m) ...
);

prop_each = { 'unified_filename', 'scrambled_type', 'target_image_category' };
[props, prop_labs] = hwwa.maybe_apply_rowop( double(remained_on_cue), plt_labs ...
  , prop_each, @mean, hwwa.make_mask(plt_labs, mask_func) ...
);

rs_outs = dsp3.ranksum( props, prop_labs', {}, '5-htp', 'saline' );
anova_outs = dsp3.anovan( props, prop_labs' ...
  , {}, {'drug', 'scrambled_type', 'target_image_category'} ); 

mu = mean( props );
sig = plotlabeled.sem( props );

if ( true )
  data_tbl = array2table( props, 'VariableNames', {'prop_fixated'} );
  label_tbl = array2table( cellstr(prop_labs), 'VariableNames', getcats(prop_labs) );
  save_table_pair( data_tbl, label_tbl, 'prop_fixated_yy2', '' );
end

%%  Plot proportion of trials remained fixated on cue.

do_save = false;

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_nogo(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m, true) ...
);

prop_each = { 'monkey', 'unified_filename', 'trial_type', 'scrambled_type', 'target_image_category', 'correct' };
[props, prop_labs] = hwwa.maybe_apply_rowop( double(remained_on_cue), plt_labs ...
  , prop_each, @mean, hwwa.make_mask(plt_labs, mask_func) ...
);

[props, prop_labs] = hwwa.saline_normalize( props, prop_labs, setdiff(prop_each, 'unified_filename') );

anova_factors = { 'correct', 'scrambled_type', 'target_image_category' };
anova_outs = dsp3.anovan( props, prop_labs, {'trial_type'}, anova_factors ...
  , 'include_significant_factor_descriptives', true ...
);

hwwa_basic_bar( double(remained_on_cue), plt_labs, 'prop_trials_fixated_on_cue' ...
  , 'mask_func', mask_func ...
  , 'collapse_each', {'unified_filename', 'trial_type', 'scrambled_type', 'target_image_category', 'correct'} ...
  , 'xcats', {'trial_type'} ...
  , 'gcats', {'drug'} ...
  , 'pcats', {} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'y_lims', [0, 1] ...
  , 'do_save', do_save ...
  , 'run_anova', false ...
  , 'anova_each', {} ...
  , 'anova_factors', {'monkey', 'drug', 'correct', 'scrambled_type', 'target_image_category'} ...
);

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'prop_trials_fixated_on_cue' );
  dsp3.save_anova_outputs( anova_outs, save_p, anova_factors );
end

%%  Plot distance from cue center on cue-broken trials

do_save = false;

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) intersect(m, find(~remained_on_cue)) ...
  , @(m) hwwa.find_incorrect_go_correct_nogo(l, m) ...
);

[dist, dist_labs] = hwwa.maybe_apply_rowop( dist_from_cue_center, plt_labs ...
  , {'unified_filename', 'trial_type'}, @mean, hwwa.make_mask(plt_labs, mask_func) ...
);

t_outs = dsp3.ttest2( dist, dist_labs, {'trial_type'}, 'saline', '5-htp' );

hwwa_basic_bar( dist_from_cue_center, plt_labs, 'distance_from_cue_center' ...
  , 'mask_func', mask_func ...
  , 'collapse_each', {'unified_filename', 'trial_type'} ...
  , 'xcats', {'trial_type'} ...
  , 'gcats', {'drug'} ...
  , 'pcats', {} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'do_save', do_save ...
);

%%  descr stats distance from cue center

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m, true) ...
  , @(m) hwwa.find_nogo(l, m) ...
  , @(m) hwwa.find_saline(l, m) ...
);

mean_each = { 'unified_filename' };
[dists, dist_labs] = hwwa.maybe_apply_rowop( dist_from_cue_center, plt_labs ...
  , mean_each, @mean, hwwa.make_mask(plt_labs, mask_func) ...
);

mu = mean( dists );
sig = plotlabeled.sem( dists );

%% anova stats distance from cue center

do_save = false;
do_norm = false;

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m, true) ...
  , @(m) hwwa.find_nogo(l, m) ...
);

mean_each = { 'monkey', 'unified_filename', 'trial_type', 'scrambled_type', 'target_image_category', 'correct' };
[dists, dist_labs] = hwwa.maybe_apply_rowop( dist_from_cue_center, plt_labs ...
  , mean_each, @mean, hwwa.make_mask(plt_labs, mask_func) ...
);

if ( do_norm )
  [dists, dist_labs] = hwwa.saline_normalize( ...
    dists, dist_labs, setdiff(mean_each, 'unified_filename') );
  anova_factors = { 'scrambled_type', 'target_image_category', 'correct' };
  anova_each = { 'trial_type' };
else
  anova_factors = { 'correct', 'scrambled_type', 'target_image_category', 'monkey', 'drug' };
  anova_each = {'trial_type'};
end

if ( true )
  data_tbl = array2table( dists, 'VariableNames', {'cue_dist'} );
  label_tbl = array2table( cellstr(dist_labs), 'VariableNames', getcats(prop_labs) );
  save_table_pair( data_tbl, label_tbl, 'cue_dist_yy3', '' );
end

anova_outs = dsp3.anovan( dists, dist_labs, anova_each, anova_factors ...
  , 'include_significant_factor_descriptives', true ...
);

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'distance_from_cue_center' );
  dsp3.save_anova_outputs( anova_outs, save_p, union(anova_factors, anova_each) );
end

hwwa_basic_bar( dists, dist_labs, 'distance_from_cue_center' ...
  , 'mask_func', @hwwa.default_mask_func ...
  , 'collapse_each', {} ...
  , 'xcats', {'trial_type'} ...
  , 'gcats', {'drug', 'correct'} ...
  , 'pcats', {'scrambled_type', 'target_image_category'} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'do_save', do_save ...
);

%%  microsaccades + regular saccades

x = behav_outs.cue_on_aligned.x;
y = behav_outs.cue_on_aligned.y;
t = behav_outs.cue_on_aligned.time;
microsaccade_start_stops = cell( size(x, 1), 1 );

parfor i = 1:numel(microsaccade_start_stops)
  [deg_x, deg_y] = hwwa.run_px2deg( x(i, :), y(i, :) );
  microsaccade_start_stops(i) = hwwa.find_microsaccades( deg_x, deg_y ...
    , 'max_dur', 60 ...
    , 'amp_thresh', 30 ...
  );
end

[x_deg, y_deg] = hwwa.run_px2deg( x, y );
regular_saccades = hwwa.run_find_saccades( x_deg, y_deg );
num_microsaccades = cellfun( @rows, microsaccade_start_stops );

%%  num saccades per trial

restricted_saccades = regular_saccades;
restricted_tlim = [0, 750];
for i = 1:numel(restricted_saccades)
  t_starts = t(restricted_saccades{i}(:, 1));
  t_stops = t(restricted_saccades{i}(:, 2));
  
  keep_start = t_starts >= restricted_tlim(1) & t_starts < restricted_tlim(2);
  keep_stop = t_stops >= restricted_tlim(1) & t_stops < restricted_tlim(2);
  keep_sacc = keep_start & keep_stop;
  
  restricted_saccades{i}(~keep_sacc, :) = [];
end

num_regular_saccades = cellfun( @rows, restricted_saccades );
regular_saccade_durs = cellfun( @(x) x(:, 2) - x(:, 1), restricted_saccades, 'un', 0 );

%%  debug plot saccades

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_incorrect_go(l, m) ...
  , @(m) intersect(m, find(num_regular_saccades >= 4)) ...
);

trial_mask = mask_func( use_labs, rowmask(use_labs) );
trial = trial_mask(1);

eg_x = x(trial, :);
eg_y = y(trial, :);
eg_starts = regular_saccades{trial}(:, 1);
eg_stops = regular_saccades{trial}(:, 2);

axs = hwwa.plot_saccade( t(1, :), eg_x, eg_y, eg_starts, eg_stops );

%%  plot prop regular saccades on incorrect trials

do_save = false;

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_incorrect_go(l, m) ...
);

max_num_saccs = 3;
use_num_regular_saccades = min( num_regular_saccades, max_num_saccs );

num_saccs = unique( use_num_regular_saccades );

% plot_each = {'trial_type', 'correct', 'scrambled_type', 'target_image_category', 'drug'};
plot_each = {'trial_type', 'drug'};

mask = mask_func( plt_labs, rowmask(plt_labs) );
[plt_labs, plt_I] = keepeach( plt_labs, plot_each, mask );
addcat( plt_labs, 'num_saccades' );

props = [];
prop_labs = fcat();
for i = 1:numel(plt_I)
  p_ind = plt_I{i};
  for j = 1:numel(num_saccs)
    ct = sum( use_num_regular_saccades(p_ind) == num_saccs(j) );
    append( prop_labs, plt_labs, i );
    if ( num_saccs(j) > max_num_saccs-1 )
      ct_label = sprintf( '> %d saccades', max_num_saccs-1 );
    else
      ct_label = sprintf( '%d-saccades', num_saccs(j) );
    end
    setcat( prop_labs, 'num_saccades', ct_label, rows(prop_labs) );
    props(end+1, 1) = ct / numel( p_ind );
  end
end

hwwa_basic_bar( props, prop_labs, 'prop_regular_saccades' ...
  , 'mask_func', @hwwa.default_mask_func ...
  , 'xcats', {'drug'} ...
  , 'gcats', {'num_saccades'} ...
  , 'pcats', {'trial_type', 'correct'} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'y_label', 'Proportion of trials' ...
  , 'do_save', do_save ...
  , 'run_anova', false ...
  , 'anova_each', {} ...
  , 'anova_factors', {'drug', 'trial_type', 'correct', 'scrambled_type', 'target_image_category'} ...
  , 'stacked', true ...
);

%%  compare number of regular saccades

do_save = true;

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_go(l, m) ...
);

hwwa_basic_bar( num_regular_saccades, plt_labs, 'num_regular_saccades' ...
  , 'mask_func', mask_func ...
  , 'collapse_each', {'unified_filename', 'trial_type', 'correct', 'scrambled_type', 'target_image_category'} ...
  , 'collapse_op', @median ...
  , 'xcats', {'drug'} ...
  , 'gcats', {'target_image_category', 'scrambled_type'} ...
  , 'pcats', {'trial_type', 'correct'} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'y_label', 'Number of saccades' ...
  , 'do_save', do_save ...
  , 'run_anova', true ...
  , 'anova_each', {'correct', 'trial_type'} ...
  , 'anova_factors', {'drug', 'scrambled_type', 'target_image_category'} ...
  , 'stacked', false ...
);

%%  incorrect saccade landing point polar plot

num_bins = 32;
first_only = true;
do_save = true;
collapse_left_right = true;
add_smoothing = false;
smooth_func = @(x) conv(x, gausswin(num_bins, 0.25)/sum(gausswin(num_bins, 0.25)), 'same');

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_go(l, m) ...
);

bin_each = { 'drug', 'trial_type', 'correct', 'unified_filename' };
collapse_each = setdiff( bin_each, 'unified_filename' );
bin_mask = mask_func( use_labs, rowmask(use_labs) );
[bin_labs, bin_I, bin_C] = keepeach( use_labs', bin_each, bin_mask );

theta_bins = linspace( 0, 2*pi, num_bins ) - pi;
bin_width = 2 * pi / (num_bins-1);

bin_values = zeros( numel(bin_I), numel(theta_bins) );
bin_counts = zeros( size(bin_values) );
all_thetas = cell( rows(use_labs), 1 );
all_dists = cell( size(all_thetas) );
max_values = zeros( rows(bin_values), 1 );
mean_vs = zeros( numel(bin_I), 2 );

for idx = 1:numel(bin_I)
  bin_ind = bin_I{idx};
  curr_vs = [];
  
  for i = 1:numel(bin_ind)
    ind = bin_ind(i);
    targ_loc = cellstr( use_labs, 'target_placement', ind );
    is_left = strcmp( targ_loc, 'center-left' );
    
    stops = restricted_saccades{ind}(:, 2);
    if ( first_only && ~isempty(stops) )
      stops = stops(1);
    end
    
    targ_rect = roi_outs.rois(ind).go_target;
    stop_x = behav_outs.cue_on_aligned.x(ind, stops);
    stop_y = behav_outs.cue_on_aligned.y(ind, stops);
    targ_center = shared_utils.rect.center( targ_rect );
    to_stop = [ stop_x(:), stop_y(:) ] - targ_center(:)';
    assert( rows(to_stop) == numel(stops) );
    
    if ( collapse_left_right && is_left )
      to_stop(:, 1) = -to_stop(:, 1);
    end

    sacc_dists = vecnorm( to_stop, 2, 2 );
    sacc_dirs = to_stop ./ sacc_dists;
    thetas = atan2( sacc_dirs(:, 2), sacc_dirs(:, 1) );
    non_nan = ~isnan( thetas );
    thetas = thetas(non_nan);
    sacc_dists = sacc_dists(non_nan);
    curr_vs = [ curr_vs; to_stop(non_nan, :) ];
    
    [~, bin] = histc( thetas, theta_bins );
    
    for j = 1:numel(bin)
      bin_values(idx, bin(j)) = bin_values(idx, bin(j)) + sacc_dists(j);
      bin_counts(idx, bin(j)) = bin_counts(idx, bin(j)) + 1;
    end
    
    all_thetas{ind} = thetas;
    all_dists{ind} = sacc_dists;
  end
  
  mean_vs(idx, :) = nanmean( curr_vs, 1 );
  
%   bin_values(idx, :) = bin_values(idx, :) ./ max_values(idx);
  bin_values(idx, :) = bin_values(idx, :) ./ bin_counts(idx, :);
  bin_values(idx, bin_counts(idx, :) == 0) = 0;
  max_values(idx) = max( bin_values(idx, :) );
end

if ( ~isempty(collapse_each) )
  [bin_I, bin_C] = findall( bin_labs, collapse_each );
  bin_values = zeros( numel(bin_I), numel(theta_bins) );
  bin_counts = zeros( size(bin_values) );
  
  for i = 1:numel(bin_I)
    thetas = atan2( mean_vs(bin_I{i}, 2), mean_vs(bin_I{i}, 1) );
    non_nan = ~isnan( thetas );
    
    sacc_dists = vecnorm( mean_vs(bin_I{i}(non_nan), :), 2, 2 );
    [~, bin] = histc( thetas(non_nan), theta_bins );
    
    for j = 1:numel(bin)
      bin_values(i, bin(j)) = bin_values(i, bin(j)) + sacc_dists(j);
      bin_counts(i, bin(j)) = bin_counts(i, bin(j)) + 1;
    end
  end
  
  mean_vs = bfw.row_nanmean( mean_vs, bin_I );
  max_values = repmat( max(mean_vs(:)), numel(bin_I), 1 );
end

if ( add_smoothing )
  for i = 1:rows(bin_values)
    bin_values(i, :) = smooth_func( bin_values(i, :) );
  end
end

if ( exist('axs', 'var') )
  delete( axs );
end

sub_shape = plotlabeled.get_subplot_shape( rows(bin_values) );
axs = gobjects( rows(bin_values) );
for i = 1:rows(bin_values)
  ax = subplot( sub_shape(1), sub_shape(2), i, polaraxes );
  cla( ax );
  axs(i) = ax;
  hold( ax, 'on' );
  
  if ( isempty(collapse_each) )
    ths = vertcat( all_thetas{bin_I{i}} ); 
    dists = vertcat( all_dists{bin_I{i}} );
    vs = [ cos(ths), sin(ths) ] .* dists;
    mean_v = mean( vs, 1 );
    len_v = norm( mean_v );
    th = atan2( mean_v(2), mean_v(1) );
  else
    mean_v = mean_vs(i, :);
    len_v = norm( mean_v );
    th = atan2( mean_v(2), mean_v(1) );
    
    bin_values(i, :) = bin_values(i, :) ./ bin_counts(i, :);
    bin_values(i, bin_counts(i, :) == 0) = 0;
  end
  
  line_hs = gobjects( 0 );
  for j = 1:numel(theta_bins)
    if ( j == numel(theta_bins) )
      j1 = 1;
    else
      j1 = j + 1;
    end
    
    th0 = theta_bins(j) + bin_width * 0.5;
    th1 = theta_bins(j1) + bin_width * 0.5;
    
    bin0 = bin_values(i, j) / max_values(i);
    bin1 = bin_values(i, j1) / max_values(i);
    line_hs(end+1, 1) = polarplot( ax, [th0, th1], [bin0, bin1] );
  end
  set( line_hs, 'color', 'b' );
  
%   polarplot( ax, theta_bins, bin_values(i, :)/max_values(i) );
  h = polarplot( ax, [th, th], [0, len_v/max_values(i)] );
  set( h, 'LineWidth', 5 );
  set( h, 'color', 'r' );
  title( ax, strrep(fcat.strjoin(bin_C(:, i), ' | '), '_', ' ') );
end

test_each = { 'trial_type', 'correct' };
test_a = 'saline';
test_b = '5-htp';
[test_labs, test_I] = keepeach( use_labs', test_each, bin_mask );
ups = zeros( numel(test_I), 2 );
mstats = cell( numel(test_I), 1 );

data_tbl_data = {};
data_tbl_inds = {};

data_tbl_labels = fcat();
data_tbl_mean_data = [];
session_average = true;

for i = 1:numel(test_I)
  ind_a = find( use_labs, test_a, test_I{i} );
  ind_b = find( use_labs, test_b, test_I{i} );
  tot_a = vertcat( all_thetas{ind_a} );
  tot_b = vertcat( all_thetas{ind_b} );
  [up, u2] = watsons_U2_approx_p( tot_a, tot_b );
  ups(i, :) = [up, u2];
  
  store_ind_a = cat_expanded( 1, arrayfun(@(x) repmat(x, size(all_thetas{x}, 1), 1), ind_a, 'un', 0) );
  store_ind_b = cat_expanded( 1, arrayfun(@(x) repmat(x, size(all_thetas{x}, 1), 1), ind_b, 'un', 0) );
  
  dists_a = vertcat( all_dists{ind_a} );
  dists_b = vertcat( all_dists{ind_b} );
  vs_a = [ cos(tot_a), sin(tot_a) ] .* dists_a;
  vs_b = [ cos(tot_b), sin(tot_b) ] .* dists_b;
  
  grp_labels = [ ones(size(dists_a)); zeros(size(dists_b)) ];
  vs = [ vs_a; vs_b ];
  [~, mp, mstat] = manova1( vs, grp_labels );
  mstat.p = mp;
  mstats{i} = mstat;
  
  if ( session_average )
    sesh_I = findall( use_labs, 'unified_filename', test_I{i} );
    curr_a = [];
    curr_b = [];
    
    for j = 1:numel(sesh_I)
      sesh_a = find( use_labs, test_a, sesh_I{j} );
      sesh_b = find( use_labs, test_b, sesh_I{j} );
      collapse_a = mean( vertcat(all_thetas{sesh_a}) );
      collapse_b = mean( vertcat(all_thetas{sesh_b}) );
      append1( data_tbl_labels, use_labs, sesh_I{j} );
      setcat( data_tbl_labels, whichcat(use_labs, test_a), test_a, rows(data_tbl_labels) );
      append1( data_tbl_labels, use_labs, sesh_I{j} );
      setcat( data_tbl_labels, whichcat(use_labs, test_b), test_b, rows(data_tbl_labels) );
      data_tbl_mean_data = [ data_tbl_mean_data; collapse_a; collapse_b ];
      curr_a = [ curr_a; collapse_a ];
      curr_b = [ curr_b; collapse_b ];
    end
    
    [up, u2] = watsons_U2_approx_p( curr_a, curr_b );
    ups(i, :) = [up, u2];
  end
  
  data_tbl_inds{end+1, 1} = [ store_ind_a; store_ind_b ];
  data_tbl_data{end+1, 1} = [ tot_a; tot_b ];
  
  assert( numel(store_ind_a) == numel(tot_a) && numel(store_ind_b) == numel(tot_b) );
end

if ( do_save )  
  tbl_cats = { 'drug', 'correct', 'trial_type' };
  if ( session_average )
    data_tbl = array2table( data_tbl_mean_data, 'VariableNames', {'theta'} );
    label_tbl = array2table( cellstr(data_tbl_labels, tbl_cats), 'VariableNames', tbl_cats );
  else
    data_tbl = array2table( vertcat(data_tbl_data{:}), 'VariableNames', {'theta'} );
    data_inds = vertcat( data_tbl_inds{:} );
    label_tbl = array2table( cellstr(use_labs, tbl_cats, data_inds), 'VariableNames', tbl_cats );
  end
  save_table_pair( data_tbl, label_tbl, 'polar_landing_points', '' );
end

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  plt_params.do_save = true;
  save_p = hwwa.approach_avoid_save_fig( gcf, axs, bin_labs, 'drug' ...
    , 'polar_incorrect_saccade_landing', plt_params );
end

%%  incorrect saccade landing point scatter

first_only = true;
collapse_left_right = true;

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_go(l, m) ...
);

bin_mask = mask_func( use_labs, rowmask(use_labs) );

all_thetas = cell( rows(use_labs), 1 );
all_dists = cell( size(all_thetas) );
  
for i = 1:numel(bin_mask)
  ind = bin_mask(i);
  targ_loc = cellstr( use_labs, 'target_placement', ind );
  is_left = strcmp( targ_loc, 'center-left' );

  stops = restricted_saccades{ind}(:, 2);
  if ( first_only && ~isempty(stops) )
    stops = stops(1);
  end

  targ_rect = roi_outs.rois(ind).go_target;
  stop_x = behav_outs.cue_on_aligned.x(ind, stops);
  stop_y = behav_outs.cue_on_aligned.y(ind, stops);
  targ_center = shared_utils.rect.center( targ_rect );
  to_stop = [ stop_x(:), stop_y(:) ] - targ_center(:)';
  assert( rows(to_stop) == numel(stops) );

  if ( collapse_left_right && is_left )
    to_stop(:, 1) = -to_stop(:, 1);
  end

  sacc_dists = vecnorm( to_stop, 2, 2 );
  sacc_dirs = to_stop ./ sacc_dists;
  thetas = atan2( sacc_dirs(:, 2), sacc_dirs(:, 1) );
  non_nan = ~isnan( thetas );
  thetas = thetas(non_nan);
  sacc_dists = sacc_dists(non_nan);
  all_thetas{ind} = thetas;
  all_dists{ind} = sacc_dists;
end

scatter_each = { 'correct', 'drug', 'trial_type' };
[scatter_I, scatter_C] = findall( use_labs, scatter_each, bin_mask );

if ( exist('axs', 'var') )
  delete( axs );
end

sub_shape = plotlabeled.get_subplot_shape( numel(scatter_I) );
axs = gobjects( numel(scatter_I), 1 );  

for i = 1:numel(scatter_I)
  ax = subplot( sub_shape(1), sub_shape(2), i, polaraxes );
  cla( ax );
  axs(i) = ax;
  hold( ax, 'on' );
  
  si = scatter_I{i};
  ths = vertcat( all_thetas{si} );
  dists = vertcat( all_dists{si} );
  polarscatter( ax, ths, dists, 0.5 );
  title( ax, strrep(fcat.strjoin(scatter_C(:, i), ' | '), '_', ' ') );
end

rlims = [0, 2e3];
arrayfun( @(x) set(x, 'rlim', rlims), axs, 'un', 0 );

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  plt_params.do_save = true;
  hwwa.approach_avoid_save_fig( gcf, axs, use_labs, 'drug' ...
    , 'polar_saccade_landing_scatter', plt_params );
end

%%  distribution of saccade latencies

first_only = true;

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_correct_go_incorrect_nogo(l, m) ...
);

bin_mask = mask_func( use_labs, rowmask(use_labs) );

latencies = [];
latency_labs = fcat();
  
for i = 1:numel(bin_mask)
  ind = bin_mask(i);
  targ_loc = cellstr( use_labs, 'target_placement', ind );
  is_left = strcmp( targ_loc, 'center-left' );
  
  starts = restricted_saccades{ind}(:, 1);
  if ( first_only && ~isempty(starts) )
    starts = starts(1);
  end
  
  start_ts = behav_outs.cue_on_aligned.time(starts);
  if ( isempty(start_ts) )
    start_ts = [];
  else
    start_ts = start_ts(:);
  end
  
  append1( latency_labs, use_labs, ind, numel(starts) );  
  latencies = [ latencies; start_ts ];
  assert_ispair( latencies, latency_labs );
end

%%  saccade latency stats

collapse_each = { 'unified_filename', 'trial_type' };

[lats, lat_labs] = hwwa.maybe_apply_rowop( ...
  latencies, latency_labs', collapse_each, @mean );
ks_outs = dsp3.kstest2( lats, lat_labs', {'trial_type'}, '5-htp', 'saline' );

%%  plot

plt_latencies = latencies;
plt_labs = latency_labs';

% plt_latencies = saccade_outs.saccade_start_stop(:, 1);
% plt_labs = saccade_outs.labels';

collapse_each = {'unified_filename', 'trial_type', 'correct'};

[plt_latencies, plt_labs] = hwwa.maybe_apply_rowop( ...
  plt_latencies, plt_labs', collapse_each, @mean );
% [plt_latencies, plt_labs] = hwwa.saline_normalize( ...
%   plt_latencies, plt_labs, setdiff(collapse_each, 'unified_filename') );

% mask_func = @(l, m) pipe(m ...
%   , @(m) hwwa.find_go(l, m) ...
% );

mask_func = @hwwa.default_mask_func;

hwwa_basic_hist( plt_latencies, plt_labs, 'time_to_initiate' ...
  , 'mask_func', mask_func ...
  , 'collapse_each', {} ...
  , 'collapse_op', @mean ...
  , 'pcats', {'correct', 'drug', 'trial_type'} ...
  , 'y_label', 'Time to initiate saccade' ...
  , 'do_save', true ...
  , 'run_anova', true ...
  , 'anova_each', {'trial_type'} ...
  , 'anova_factors', {'correct', 'drug'} ...
  , 'run_ranksum', true ...
  , 'ranksum_each', {'correct'} ...
  , 'ranksum_a', 'saline' ...
  , 'ranksum_b', '5-htp' ...
  , 'hist_size', 30 ...
);

%% check side-bias during fixations

fix_info = cell( size(microsaccade_start_stops) );
all_fix_info = cell( size(fix_info) );

fix_dur_thresh = 50;
t_lims = [-150, 0];
fix_x_width_px = 150;
fix_t = t;

for i = 1:numel(fix_info)
  curr_rois = roi_outs.rois(i);
  screen_center = shared_utils.rect.center( curr_rois.nogo_cue );
  
  all_sacc_info = [ microsaccade_start_stops{i}; regular_saccades{i} ];
  curr_fix = find_fixations( all_sacc_info, size(x, 2), fix_dur_thresh );
  curr_fix_p = cat_expanded( 1 ...
    , arrayfun(@(start, stop) [nanmean(x(i, start:stop)), nanmean(y(i, start:stop))] ...
    , curr_fix(:, 1), curr_fix(:, 2), 'un', 0) );
  all_fix_info{i} = [ curr_fix_p, curr_fix(:, 1) ];  
  
  fix_start_t = t(curr_fix(:, 1));
  within_time_bounds = fix_start_t >= t_lims(1) & fix_start_t <= t_lims(2);
  % express relative to center x
  curr_fix_p(:, 1) = curr_fix_p(:, 1) - screen_center(1);
  within_x_bounds = abs( curr_fix_p(:, 1) ) <= fix_x_width_px;
  meets_criteria = within_time_bounds(:) & within_x_bounds(:);
  % keep last fixation if there are multiple.
  meets_criteria = find( meets_criteria, 1, 'last' );
  
  tot_fix_info = [ curr_fix, curr_fix_p ];
  tot_fix_info = tot_fix_info(meets_criteria, :);
  fix_info{i} = tot_fix_info;
end

%%  unpack fixation x and rt

plt_labs = use_labs';
mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_correct_go_incorrect_nogo(l, m) ...
);

fix_x = cellfun( @(x) x(:, 4), fix_info, 'un', 0 );
fix_mask = mask_func( plt_labs, rowmask(plt_labs) );

unpacked_rt = [];
unpacked_fix_x = [];
unpacked_labels = fcat();

for i = 1:numel(fix_mask)
  ind = fix_mask(i);
  tmp_fix_x = fix_x{ind};
  
  unpacked_rt = [ unpacked_rt; repmat(rt(ind), rows(tmp_fix_x), 1) ];
  unpacked_fix_x = [ unpacked_fix_x; tmp_fix_x ];
  append1( unpacked_labels, plt_labs, ind, rows(tmp_fix_x) );
end

keep_non_nan = find( ~isnan(unpacked_fix_x) & ~isnan(unpacked_rt) );
unpacked_rt = unpacked_rt(keep_non_nan);
unpacked_fix_x = unpacked_fix_x(keep_non_nan);
unpacked_labels = unpacked_labels(keep_non_nan);

%%  plot 

do_save = true;

collapse_left = true;
use_fix_x = unpacked_fix_x;
use_fix_rt = unpacked_rt;
use_fix_labels = unpacked_labels';

if ( collapse_left )
  left_ind = find( use_fix_labels, 'center-left' );
  use_fix_x(left_ind) = -use_fix_x(left_ind);
  setcat( use_fix_labels, whichcat(use_fix_labels, 'center-left'), 'center-right' );
end

hwwa_basic_scatter( use_fix_x, use_fix_rt, use_fix_labels', 'fix_x_vs_rt' ...
  , 'mask_func', mask_func ...
  , 'collapse_each', {} ...
  , 'collapse_op', @mean ...
  , 'gcats', {} ...
  , 'pcats', {'target_placement', 'drug'} ...
  , 'marker_size', 3 ...
  , 'x_label', 'Center-relative x fixation position (px)' ...
  , 'y_label', 'RT (s)' ...
  , 'x_summary_func', @nanmedian ...
  , 'x_test_func', [] ...
  , 'y_summary_func', @nanmedian ...
  , 'do_save', do_save ...
  , 'permute', true ...
  , 'permute_each', {'target_placement'} ...
  , 'permute_compare', {'drug'} ...
  , 'perm_iters', 1e3 ...
);

%%  microsacc stats

plt_labs = use_labs';
plt_microsaccades = num_microsaccades;

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m) ...
  , @(m) hwwa.find_nogo(l, m) ...
);

microsacc_mask = hwwa.make_mask( plt_labs, mask_func );

collapse_each = { 'unified_filename', 'correct', 'drug', 'scrambled_type', 'target_image_category', 'monkey' };

[plt_microsaccades, plt_labs] = hwwa.maybe_apply_rowop( ...
    plt_microsaccades, plt_labs, collapse_each, @mean, microsacc_mask );
  
mu = mean( plt_microsaccades );
sig = plotlabeled.sem( plt_microsaccades );

if ( true )
  data_tbl = array2table( plt_microsaccades, 'VariableNames', {'num_ms'} );
  label_tbl = array2table( cellstr(plt_labs), 'VariableNames', getcats(plt_labs) );
  save_table_pair( data_tbl, label_tbl, 'num_ms_yy4', '' );
end

%%  plot microsaccades

do_save_plot = false;
do_save_anova = false;

plt_labs = use_labs';
plt_microsaccades = num_microsaccades;

% mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
%   , @(m) intersect(m, find(~remained_on_cue)) ...
%   , @(m) hwwa.find_incorrect_go_correct_nogo(l, m) ...
% );

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.get_approach_avoid_mask(l, m) ...
  , @(m) hwwa.find_nogo(l, m) ...
);

microsacc_mask = hwwa.make_mask( plt_labs, mask_func );

collapse_each = { 'unified_filename', 'monkey', 'trial_type', 'correct' ...
  , 'scrambled_type', 'target_image_category'};
do_norm = false;

[plt_microsaccades, plt_labs] = hwwa.maybe_apply_rowop( ...
    plt_microsaccades, plt_labs, collapse_each, @mean, microsacc_mask );
  
if ( do_norm )
  [plt_microsaccades, plt_labs] = hwwa.saline_normalize( ...
    plt_microsaccades, plt_labs, setdiff(collapse_each, 'unified_filename') );
  
  anova_factors = { 'correct', 'scrambled_type', 'target_image_category' };
else
  anova_factors = { 'drug', 'monkey', 'correct', 'scrambled_type', 'target_image_category' };
end

anova_each = { 'trial_type' };
anova_outs = dsp3.anovan( plt_microsaccades, plt_labs, anova_each, anova_factors ...
  , 'include_significant_factor_descriptives', true ...
);

if ( do_save_anova )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'num_microsaccades' );
  dsp3.save_anova_outputs( anova_outs, save_p, anova_factors );
end

hwwa_basic_bar( plt_microsaccades, plt_labs, 'num_microsaccades' ...
  , 'mask_func', @hwwa.default_mask_func ...
  , 'collapse_each', {} ...
  , 'collapse_op', @mean ...
  , 'xcats', {'trial_type'} ...
  , 'gcats', {'drug'} ...
  , 'pcats', {'correct', 'scrambled_type', 'target_image_category'} ...
  , 'points_are', 'monkey' ...
  , 'marker_size', 3 ...
  , 'y_label', 'Mean # microsaccades per trial' ...
  , 'do_save', do_save_plot ...
  , 'run_anova', false ...
  , 'anova_each', {'trial_type'} ...
  , 'anova_factors', {'correct', 'scrambled_type', 'target_image_category'} ...
);

%%  Compare pupil size between NPP and current dataset

do_save_stats = false;

npp_root = fullfile( fileparts(hwwa.dataroot()), 'hannah-npp' );
pup_file = fullfile( npp_root, 'pupil_during_fixation_plus_lookdur.mat' );
npp_pup_data = load( pup_file );
npp_pup = npp_pup_data.dat.pupil_size;

pup_labs = use_labs';
addsetcat( pup_labs, 'dose', 'low' );

curr_pup_size = behav_outs.pupil_size;
pup_mask_func = @(l, m) pipe(find_non_outliers(l, m));
pup_mask = pup_mask_func( pup_labs, rowmask(pup_labs) );
curr_pup_size = curr_pup_size(pup_mask);
pup_labs = prune( pup_labs(pup_mask) );

npp_labs = fcat.from( npp_pup.labels );
sal_ind = find( npp_labs, 'saline' );
ser_ind = findnot( npp_labs, 'saline' );

prune( replace(npp_labs, 'saline', 'low') );
prune( replace(npp_labs, 'tarantino', 'tar') );
addcat( npp_labs, 'drug' );
setcat( npp_labs, 'drug', 'saline', sal_ind );
setcat( npp_labs, 'drug', '5-htp', ser_ind );

npp_pup = npp_pup.data;
% npp_pup = cellfun( @nanmean, npp_pup );
npp_pup = cellfun( @nanmean, npp_pup );
assert_ispair( npp_pup, npp_labs );

[curr_norm_pup, curr_labs] = hwwa.saline_normalize( curr_pup_size, pup_labs', {'monkey'} );

% Normalize npp data
npp_norm_pup = nan( size(npp_pup) );
npp_norm_labs = npp_labs';
[npp_norm_I, npp_norm_C] = findall( npp_labs, {'monkeys', 'doses', 'drug'}, find(npp_labs, '5-htp') );
for i = 1:numel(npp_norm_I)
  norm_C = npp_norm_C(:, i);
  search_for = [norm_C(1), {'saline'}];
%   search_for = 'saline';
  sal_ind = find( npp_labs, search_for );
  sal_mean = nanmean( npp_pup(sal_ind) );
  npp_norm_pup(npp_norm_I{i}) = npp_pup(npp_norm_I{i}) ./ sal_mean;
  setcat( npp_norm_labs, 'drug', '5-htp/saline', npp_norm_I{i} );
end

[npp_norm_labs, npp_I] = keepeach( npp_norm_labs', {'days', 'doses', 'monkeys'} ...
  , find(npp_norm_labs, '5-htp/saline') );
npp_norm_pup = bfw.row_nanmean( npp_norm_pup, npp_I );

[curr_norm_labs, curr_I] = keepeach( curr_labs', {'unified_filename'} );
curr_norm_pup = bfw.row_nanmean( curr_norm_pup, curr_I );

keep_curr_labels = cellstr( curr_norm_labs, {'dose', 'monkey'} );
keep_curr_labels(:, end+1) = { 'current-data' };

keep_npp_labels = cellstr( npp_norm_labs, {'doses', 'monkeys'} );
keep_npp_labels(:, end+1) = { 'npp-data' };

new_cats = { 'dose', 'monkey', 'data-set' };
combined_labels = [ fcat.from(keep_curr_labels, new_cats); fcat.from(keep_npp_labels, new_cats) ];
combined_data = [ curr_norm_pup; npp_norm_pup ];
assert_ispair( combined_data, combined_labels );

anova_outs = dsp3.anova1( combined_data, combined_labels, {}, {'data-set'} ...
  , 'mask', find(combined_labels, {'low', 'ephron', 'tar', 'hitch'}) ...
);

t_outs = dsp3.ttest2( combined_data, combined_labels, {}, 'current-data', 'npp-data' ...
  , 'mask', find(combined_labels, {'low', 'ephron', 'tar', 'hitch'}) ...
);

rs_outs = dsp3.ranksum( combined_data, combined_labels, {}, 'current-data', 'npp-data' ...
  , 'mask', find(combined_labels, {'low', 'ephron', 'tar', 'hitch'}) ...
);

if ( do_save_stats )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'pupil_npp/stats' );
  dsp3.save_ttest2_outputs( t_outs, save_p );
  dsp3.save_ranksum_outputs( rs_outs, save_p );
end

%%

do_save = true;

pl = plotlabeled.make_common();
pl.add_points = true;
pl.points_are = { 'monkey' };

mask_func = @(l, m) pipe(m ...
  , @(m) find(l, 'low', m) ...
  , @(m) find(l, {'ephron', 'tar', 'hitch'}, m) ...
);

mask = hwwa.make_mask( combined_labels, mask_func );
pltdat = combined_data(mask);
pltlabs = combined_labels(mask);

% [pltlabs, I] = keepeach( pltlabs', {'data-set', 'dose', 'monkey'} );
% pltdat = bfw.row_nanmean( pltdat, I );

axs = pl.bar( pltdat, pltlabs, 'data-set', {'dose'}, {} );

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  plt_params.do_save = true;
  save_p = hwwa.approach_avoid_save_fig( gcf, axs, pltlabs, 'monkey', 'pupil_npp', plt_params );
end

%%  PCA main measures

do_save = false;

use_rt = rt;
use_pup_size = behav_outs.pupil_size;
use_time_bw = time_bw;

assert_ispair( use_time_bw, time_bw_labels );

rt_mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) find(l, {'go_trial', 'correct_true'}, m) ...
);
pcorr_mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) find(l, {'go_trial'}, m) ...
);
pup_mask_func = @(l, m) pipe(find_non_outliers(l, m));
time_bw_mask_func = @hwwa.default_mask_func;

rt_mask =       rt_mask_func( use_labs, rowmask(use_labs) );
pcorr_mask =    pcorr_mask_func( use_labs, rowmask(use_labs) );
pup_mask =      pup_mask_func( use_labs, rowmask(use_labs) );
time_bw_mask =  time_bw_mask_func( time_bw_labels, rowmask(time_bw_labels) );

norm_each = { 'monkey' };
pcorr_each = { 'unified_filename' };
sesh_each = { 'unified_filename', 'monkey' };

[rt_norm, rt_norm_labs] = hwwa.saline_normalize( use_rt, use_labs', norm_each, rt_mask );
pup_norm = hwwa.saline_normalize( use_pup_size, use_labs', norm_each, pup_mask );

[pcorr, pcorr_labs] = hwwa.percent_correct( use_labs', pcorr_each, pcorr_mask );
[pcorr_norm, pcorr_labs] = hwwa.saline_normalize( pcorr, pcorr_labs, norm_each );

[time_bw_norm, time_bw_norm_labs] = hwwa.saline_normalize( use_time_bw, time_bw_labels, norm_each, time_bw_mask );

all_sessions = cellstr( unique(categorical([...
    combs(rt_norm_labs, sesh_each) ...
  , combs(pcorr_labs, sesh_each) ...
  , combs(time_bw_norm_labs, sesh_each) ...
])', 'rows') )';

pca_labs = fcat();
pca_dat = [];

for i = 1:size(all_sessions, 2)
  sesh_c = all_sessions(:, i);
  rt_ind = find( rt_norm_labs, sesh_c );
  pup_ind = rt_ind;
  pcorr_ind = find( pcorr_labs, sesh_c );
  time_bw_ind = find( time_bw_norm_labs, sesh_c );
  
  mean_rt = nanmean( rt_norm(rt_ind) );
  mean_pup = nanmean( pup_norm(pup_ind) );
  mean_pcorr = nanmean( pcorr_norm(pcorr_ind) );
  mean_time_bw = nanmean( time_bw_norm(time_bw_ind) );
  
  if ( ~any([isnan(mean_rt), isnan(mean_pup), isnan(mean_pcorr), isnan(mean_time_bw)]) )
    pca_dat(end+1, :) = [ mean_rt, mean_pup, mean_pcorr, mean_time_bw ];
    append( pca_labs, fcat.from(sesh_c', sesh_each) );
  end
end

[coeff, score] = pca( pca_dat );
pcs = score(:, 1:2);
cluster_sex = kmeans( pcs, 2, 'Replicates', 1e3 );
cluster_id = kmeans( pcs, 3, 'Replicates', 1e3 );

addsetcat( pca_labs, 'subject_cluster_id' ...
  , arrayfun(@(x) sprintf('subject_cluster_%d', x), cluster_id, 'un', 0) );
addsetcat( pca_labs, 'sex_cluster_id' ...
  , arrayfun(@(x) sprintf('sex_cluster_%d', x), cluster_sex, 'un', 0) );

addcat( pca_labs, 'sex' );
setcat( pca_labs, 'sex', 'male', find(pca_labs, {'tar', 'hitch'}) );
setcat( pca_labs, 'sex', 'female', find(pca_labs, 'ephron') );

[n_subj, n_subj_labels, ~, subj_denom] = counts_of( pca_labs, 'subject_cluster_id', 'monkey' );
chi2_subj = dsp3.chi2_tabular_frequencies( n_subj, n_subj_labels, {}, 'subject_cluster_id', 'monkey' );

[n_sex, n_sex_labels] = counts_of( pca_labs, 'sex_cluster_id', 'sex' );
chi2_sex = dsp3.chi2_tabular_frequencies( n_sex, n_sex_labels, {}, 'sex_cluster_id', 'sex' );

colors_sex = hsv( 2 );
colors_id = hsv( 3 );
colors_actual_sex = summer( 2 );
colors_actual_id = summer( 3 );

ax1 = subplot( 1, 2, 1 );
cla( ax1 ); hold( ax1, 'on' );
scatter( ax1, score(:, 1), score(:, 2), 40, colors_sex(cluster_sex, :) );
[unique_sexes, ~, sex_ind] = unique( pca_labs(:, 'sex') );
h_sex = gscatter( ax1, score(:, 1), score(:, 2), sex_ind, colors_actual_sex(sex_ind, :), [], 5 );
legend( h_sex, unique_sexes );
xlabel( ax1, 'pc1' ); ylabel( ax1, 'pc2' );
title( ax1, sprintf('sex - chi2 = %0.2f, p = %0.2f', chi2_sex.chi2, chi2_sex.p) );

ax2 = subplot( 1, 2, 2 );
cla( ax2 ); hold( ax2, 'on' );
scatter( ax2, score(:, 1), score(:, 2), 40, colors_id(cluster_id, :) );
[unique_subjects, ~, subj_ind] = unique( pca_labs(:, 'monkey') );
h_subj = gscatter( ax2, score(:, 1), score(:, 2), subj_ind, colors_actual_id(subj_ind, :), [], 5 );
legend( h_subj, unique_subjects );
xlabel( ax2, 'pc1' ); ylabel( ax2, 'pc2' );
title( ax2, sprintf('subject - chi2 = %0.2f, p = %0.2f', chi2_subj.chi2, chi2_subj.p) );

if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  plt_params.do_save = true;
  save_p = hwwa.approach_avoid_save_fig( gcf, [ax1, ax2], pca_labs, 'monkey', 'pca_main_measures', plt_params );
end

%%  export to R

do_save = true;

[t_labs, t_cats] = cellstr( pca_labs );
label_tbl = array2table( t_labs, 'VariableNames', t_cats );
data_tbl = array2table( pca_dat, 'VariableNames', {'rt', 'pupil_size', 'pcorr', 'time_bw'} );

if ( do_save )
  params = hwwa.get_common_make_defaults( hwwa.get_common_plot_defaults() );
  save_p = hwwa.approach_avoid_data_path( params, 'export' );
  dsp3.req_writetable( label_tbl, save_p, pca_labs, 'monkey', 'labels_' );
  dsp3.req_writetable( data_tbl, save_p, pca_labs, 'monkey', 'data_' );
end

%%  broken fix heatmap

monitor_info = hwwa.monitor_constants();
monitor_width_px = monitor_info.horizontal_resolution_px;
monitor_height_px = monitor_info.vertical_resolution_px;

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) intersect(m, find(~remained_on_cue)) ...
  , @(m) hwwa.find_incorrect_go_correct_nogo(l, m) ...
);

mask = mask_func( use_labs, rowmask(use_labs) );

heat_labs = fcat();
xs = [];
ys = [];

for i = 1:numel(mask)
  info = all_fix_info{mask(i)};
  info_t = fix_t(info(:, 3));
  keep_fix = find( info_t >= 0 ); % keep after cue onset
  
  xs = [ xs; info(keep_fix, 1)/monitor_width_px ];
  ys = [ ys; info(keep_fix, 2)/monitor_height_px ];
  append1( heat_labs, use_labs, mask(i), numel(keep_fix) );
end

assert_ispair( xs, heat_labs );
assert_ispair( ys, heat_labs );

x_lims = [ -0.2, 1.2 ];
y_lims = [ -0.2, 1.2 ];

stp = 0.01;
win = 0.01;

heat_map_each = { 'trial_type', 'drug', 'correct', 'unified_filename' };

[heat_maps, heat_map_labs, x_edges, y_edges] = hwwa_make_gaze_heatmap( ...
    xs, ys, heat_labs', heat_map_each, x_lims, y_lims, stp, win ...
  , 'mask', rowmask(heat_labs) ...
);

%%  plot broken fix heatmap

cx_edges = x_edges - 0.5;
cy_edges = y_edges - 0.5;

cx_deg = normalized_edges_to_degrees( cx_edges, monitor_width_px );
cy_deg = normalized_edges_to_degrees( cy_edges, monitor_height_px );
overlay_rect = [0, 0, 1, 1] - 0.5;
overlay_rect([1, 3]) = normalized_edges_to_degrees( overlay_rect([1, 3]), monitor_width_px );
overlay_rect([2, 4]) = normalized_edges_to_degrees( overlay_rect([2, 4]), monitor_height_px );

hwwa_plot_target_onset_heatmap( heat_maps, heat_map_labs', cy_deg, cx_deg ...
  , 'do_save', true ...
  , 'match_c_lims', true ...
  , 'c_lims', [0, 0.05] ...
  , 'custom', true ...
  , 'pcats', {'drug'} ...
  , 'fcats', {} ...
  , 'overlay_rects', {overlay_rect} ...
);

%%  hBayesDM data

do_save_tbl = true;

drug_type = 'saline';
trial_types = { 'nogo_trial' }; 

table_path = fullfile( hwwa.dataroot ...
  , sprintf('public/hBayesDM-%s-%s.txt', drug_type, strjoin(trial_types)) );

tbl = hwwa.to_hBayesDM_table( use_labs );

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) find(l, drug_type, m) ...
  , @(m) find(l, trial_types, m) ...
);

base_mask = hwwa.get_approach_avoid_base_mask( use_labs, mask_func );

tbl = tbl(base_mask, :);
tbl_labels = prune(use_labs(base_mask));

nans = isnan( table2array(varfun(@isnan, tbl)) );
assert( ~any(nans(:)) );

if ( do_save_tbl )
  writetable( tbl, table_path );
end

%%  sliding window

time_window_size = 240;
time_step_size = 20;

time_prop_bin_func = @(data, labels, mask) hwwa.single_percent_correct(labels, mask);

[time_props, time_prop_edges, time_prop_labels, time_bin_inds] = ...
  hwwa.time_slide_bin_apply( time, time, use_labs', slide_each ...
    , time_window_size, time_step_size, time_prop_bin_func, prop_mask );

% Sliding window rt
rt_func = @(rt, labels, inds) nanmedian( rt );
  
[time_rt, time_rt_edges, time_rt_labels, time_rt_bin_inds] = ...
  hwwa.time_slide_bin_apply( rt, time, use_labs', slide_each ...
    , time_window_size, time_step_size, rt_func, rt_mask );
  
% Sliding window num initiated

time_init_mask = find_non_outliers( use_labs );
time_init_func = @(~, labels, mask) hwwa.single_num_initiated(labels, mask);

[time_num_init, ~, time_num_init_labels] = ...
  hwwa.time_slide_bin_apply( time, time, use_labs', slide_each ...
  , time_window_size, time_step_size, time_init_func, time_init_mask );
  
time_props = cellfun( @nanmedian, time_props ); 
time_rt = cellfun( @nanmedian, time_rt );
time_num_inits = cellfun( @nanmedian, time_num_init );
time_num_inits = time_num_inits / time_window_size;

%%  num initiated over time

trial_bin_size = 25;
trial_step_size = 25;
each_I = findall( use_labs, 'unified_filename' );

num_init_func = @hwwa.single_prop_initiated;
[num_init_over_trials, trial_bin_labels] = ...
  trial_bin_apply( use_labs', num_init_func, each_I, trial_bin_size, trial_step_size );

num_init_over_trials(cellfun(@isempty, num_init_over_trials)) = {nan};
num_init_over_trials = cell2mat( num_init_over_trials );

% mask_func = @(l, m) find_non_outliers(l, m);
mask_func = @(l, m) m;

num_init = num_init_over_trials;
t = 1:size(num_init_over_trials, 2);
num_init_labs = trial_bin_labels';

hwwa_plot_num_initiated_over_time( num_init, t, num_init_labs' ...
  , 'mask_func', mask_func ...
  , 'y_label', 'Proportion initiated' ...
  , 'y_lims', [0.5, 1.2] ...
  , 'do_save', true ...
);

%%  time between initiated

[time_bw, time_bw_labels] = ...
  time_between_initiated( event_info.events, event_info.labels', event_info.event_names );

%%  time between initiated

time_bw_outs = hwwa_plot_time_between_initiated( time_bw, time_bw_labels' ...
  , 'do_save', false ...
);

rs_outs = dsp3.ranksum( time_bw_outs.means, time_bw_outs.labels', {}, '5-htp', 'saline' );

if ( false )
  data_tbl = array2table( time_bw_outs.means, 'VariableNames', {'time_bw'} );
  label_tbl = array2table( cellstr(time_bw_outs.labels), 'VariableNames', getcats(time_bw_outs.labels) );
  save_table_pair( data_tbl, label_tbl, 'time_bw', '' );
end

%%  time b/w initiated vs p correct

hwwa_time_between_initiated_vs_p_correct( time_bw, time_bw_labels', use_labs' ...
  , 'pcorr_mask_func', find_non_outliers ...
  , 'each', {'unified_filename', 'monkey', 'trial_type'} ...
  , 'normalize', true ...
  , 'do_save', true ...
  , 'percent_change', true ...
  , 'plot_cats', {{'drug'}, {'trial_type'}} ...
  , 'per_panel_corr', true ...
);

%%  quantiles x percent correct y

for_quant_type = 'max_pup_diff';

switch ( for_quant_type )
  case 'max_pup_diff'
    for_quant = cue_on_max_diff;
  case 'rt'
    for_quant = rt;
  case 'base_pupil_size'
    for_quant = behav_outs.pupil_size;
  otherwise
    error( 'Unhandled for quant type.' );
end

find_non_nan = @(l, m, x) intersect(find_non_outliers(l, m), find(~isnan(x)));
find_non_nan_quant = @(l, m) find_non_nan(l, m, for_quant);

for_quant_labels = use_labs';
% for_quant_each = union( slide_each, {'unified_filename'} );
for_quant_each = { 'unified_filename' };
for_quant_of = {};
for_quant_mask = hwwa.get_approach_avoid_base_mask( for_quant_labels, find_non_nan_quant );
num_quants = 20;

% [quant_pcorr, quant_pcorr_labels, quant_match_inds] = ...
%   hwwa.quantiles_x_percent_correct_y( for_quant, for_quant_labels ...
%   , for_quant_each, for_quant_of, num_quants, for_quant_mask );
% quant_x = bfw.row_nanmean( for_quant, quant_match_inds );

[quants, each_I] = dsp3.quantiles_each( for_quant, for_quant_labels, num_quants ...
  , for_quant_each, for_quant_of, for_quant_mask );
dsp3.add_quantile_labels( for_quant_labels, quants, 'x-quantile' );

[quant_pcorr, quant_pcorr_labels] = ...
  hwwa.percent_correct( for_quant_labels', union(slide_each, 'x-quantile'), for_quant_mask );
quant_x = nan( size(quant_pcorr) );

%%

per_sts = trufls;
per_tts = trufls;
cs = dsp3.numel_combvec( per_sts, per_tts );

for i = 1:size(cs, 2)
  c = cs(:, i);

  mask_func = find_non_outliers;

  hwwa_relate_quantiles_x_percent_correct_y( quant_x, quant_pcorr, quant_pcorr_labels' ...
    , 'mask_func', mask_func ...
    , 'do_save', true ...
    , 'per_scrambled_type', per_sts(c(1)) ...
    , 'per_trial_type', per_tts(c(2)) ...
    , 'x_label', for_quant_type ...
    , 'quantile_index_x', true ...
    , 'x_lims', [0, 21] ...
    , 'permutation_test', true ...
  );
end

%%  binned baseline pupil vs. p-correct

per_sts = false;
per_tts = true;
per_monks = true;

num_pup_quantiles = 20;

mask_func = find_non_outliers;
cs = dsp3.numel_combvec( per_sts, per_tts, per_monks );

for i = 1:size(cs, 2)
  shared_utils.general.progress( i, size(cs, 2) );
  c = cs(:, i);
  
  hwwa_baseline_pupil_quantile_vs_p_correct( behav_outs.pupil_size, use_labs' ...
    , 'mask_func', mask_func ...
    , 'per_scrambled_type', per_sts(c(1)) ...
    , 'per_trial_type', per_tts(c(2)) ...
    , 'per_monkey', per_monks(c(3)) ...
    , 'do_save', true ...
    , 'permutation_test', true ...
    , 'permutation_test_iters', 1e2 ...
    , 'num_quantiles', num_pup_quantiles ...
    , 'y_lims', [0.1, 1] ...
  );
end

%%  pupil vs num initiated

labels = use_labs';

per_sts = trufls;
per_tts = trufls;
per_monks = trufls;
metric_names = { 'baseline', 'max_diff' };

mask_func = find_non_outliers;
cs = dsp3.numel_combvec( per_sts, per_tts, per_monks, metric_names );

for i = 1:size(cs, 2)
  c = cs(:, i);
  
  pupil_metric_name = metric_names{c(4)};

  pupil_metric = ternary( strcmp(pupil_metric_name, 'baseline') ...
    , behav_outs.pupil_size, cue_on_max_diff ...
  );
  
  hwwa_relate_num_initiated_to_pupil( pupil_metric, labels ...
    , 'mask_func', mask_func ...
    , 'per_scrambled_type', per_sts(c(1)) ...
    , 'per_trial_type', per_tts(c(2)) ...
    , 'per_monkey', per_monks(c(3)) ...
    , 'pupil_metric', pupil_metric_name ...
    , 'permutation_test', true ...
    , 'do_save', true ...
  );
end

%%  overall p correct

do_save_pcorr = true;

% mask_func = find_non_outliers;
% mask_func = @hwwa.default_mask_func;

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_5htp(l, m) ...
  , @(m) hwwa.find_nogo(l, m) ...
);

plt_labs = use_labs';
% setcat( plt_labs, 'monkey', 'tar_hitch', find(plt_labs, {'tar', 'hitch'}) );
% setcat( plt_labs, 'monkey', 'ephron_hitch', find(plt_labs, {'ephron', 'hitch'}) );

% apply_social_housing_and_gender_labels( plt_labs );

pcorr_outs = hwwa_overall_p_correct( plt_labs' ...
  , 'mask_func', mask_func ...
  , 'do_save', true ...
  , 'anova_factors', {'drug', 'trial_type'} ...
  , 'anova_each', {} ...
  , 'norm_anova_factors', {'scrambled_type', 'trial_type'} ...
  , 'norm_anova_each', {} ...
  , 'per_drug', true ...
  , 'per_target_image_category', false ...
  , 'include_normalized', true ...
  , 'include_raw', true ...
  , 'pcorr_each', {'unified_filename', 'drug', 'target_image_category', 'scrambled_type'} ...
  , 'errorbar_cats', {{'trial_type'}, {'drug', 'scrambled_type'}, {}} ...
  , 'y_label', 'prop. correct' ...
  , 'points_are', {'monkey', 'scrambled_type'} ...
  , 'marker_size', 4 ...
);

% anova_outs = dsp3.anovan( pcorr_outs.normalized_pcorr, pcorr_outs.normalized_labels' ...
%   , {'target_image_category', 'scrambled_type', 'trial_type', 'drug'} );

% [t_labs, t_cats] = cellstr( pcorr_outs.normalized_labels );
% label_tbl = array2table( t_labs, 'VariableNames', t_cats );
% data_tbl = array2table( pcorr_outs.normalized_pcorr, 'VariableNames', {'pcorr'} );

anova_outs = dsp3.anovan( pcorr_outs.raw_corr, pcorr_outs.raw_labels' ...
  , {'drug'}, {'target_image_category', 'scrambled_type'} );

[t_labs, t_cats] = cellstr( pcorr_outs.raw_labels );
label_tbl = array2table( t_labs, 'VariableNames', t_cats );
data_tbl = array2table( pcorr_outs.raw_corr, 'VariableNames', {'pcorr'} );

if ( do_save_pcorr )
  save_table_pair( data_tbl, label_tbl, 'pcorr_mm', '' );
end

%%  anova normalized overall pcorr

do_save_anova = false;

use_nested = false;

if ( use_nested )
  nesting = zeros( 4 );
  nesting(2, 1) = true;

  anova_factors = {'social-housing', 'gender', 'trial_type', 'scrambled_type'};
  anova_outs = dsp3.anovan( pcorr_outs.normalized_pcorr, pcorr_outs.normalized_labels ...
    , {}, anova_factors ...
    , 'anovan_inputs', {'display', 'off', 'varnames', anova_factors, 'model', 'full', 'nested', nesting} ...
    , 'run_multcompare', false ...
  );
else
  anova_factors = { 'social-housing', 'trial_type', 'scrambled_type'};
  anova_outs = dsp3.anovan( pcorr_outs.normalized_pcorr, pcorr_outs.normalized_labels ...
    , {}, anova_factors ...
  );
end

if ( do_save_anova )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'pcorr' );
  dsp3.save_anova_outputs( anova_outs, save_p, anova_factors );
end

%% histogram

do_save_stats = true;

% plt_labs = pcorr_outs.normalized_labels';
% pcorr = pcorr_outs.normalized_pcorr;

plt_labs = pcorr_outs.raw_labels';
pcorr = pcorr_outs.raw_corr;

test_I = findall( plt_labs, 'drug' );
test_labs = fcat();
test_ps = [];

for i = 1:numel(test_I)
  monks = combs( plt_labs, 'monkey', test_I{i} );
  inds = bfw.pair_combination_indices( numel(monks) );
  
  for j = 1:size(inds, 1)
    monk_a = monks{inds(j, 1)};
    monk_b = monks{inds(j, 2)};
    
    ind_a = find( plt_labs, monk_a, test_I{i} );
    ind_b = find( plt_labs, monk_b, test_I{i} );
    
    [~, ksp] = kstest2( pcorr(ind_a), pcorr(ind_b) );
    append1( test_labs, plt_labs, test_I{i} );
    setcat( test_labs, 'monkey', sprintf('%s_%s', monk_a, monk_b), rows(test_labs) );
    
    test_ps = [ test_ps; ksp ];
  end
end

stat_tbl = array2table( test_ps, 'variablenames', {'p'} );
stat_tbl.Properties.RowNames = fcat.strjoin( test_labs(:, {'monkey', 'drug'})', ' | ' );

if ( do_save_stats )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  save_p = hwwa.approach_avoid_data_path( plt_params, 'plots', 'pcorr/stats/kstest2' );
  dsp3.req_writetable( stat_tbl, save_p, test_labs, {'monkey', 'drug'} );
end

%%  grouped histogram each animal percent correct

plt_labs = pcorr_outs.raw_labels';
pcorr = pcorr_outs.raw_corr;

[figs, all_axs, fig_I] = hwwa.plot_grouped_histogram( pcorr, plt_labs ...
  , {'drug'}, {'drug'}, 'monkey' ...
  , 'hist_inputs', {'binwidth', 0.01} ...
);

shared_utils.plot.match_xlims( vertcat(all_axs{:}) );

do_save = true;
if ( do_save )
  plt_params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
  plt_params.do_save = true;
  for i = 1:numel(figs)
    ind = fig_I{i};
    save_p = hwwa.approach_avoid_save_fig( figs{i}, all_axs{i} ...
      , prune(plt_labs(ind)), {'drug', 'monkey'} ...
      , 'pcorr_hist', plt_params );
  end
end

%%  mean target distance

mask_func = @(l, m) fcat.mask(l, find_non_outliers(l, m) ...
  , @find, 'nogo_trial' ...
  , @find, 'correct_true' ...
);

labs = use_labs';
dist = go_targ_aligned.mean_targ_dist;

hwwa_plot_mean_target_distance( dist, labs ...
  , 'mask_func', mask_func ...
  , 'do_save', true ...
  , 'mean_each', {'unified_filename', 'drug', 'target_image_category', 'trial_type'} ...
  , 'bar_cats', {{'drug'}, {'target_image_category'}, {'trial_type', 'correct'}} ...
);

%% go trial sequence

labs = use_labs';
each_I = findall( labs, 'unified_filename' ...
  , find_non_outliers(labs, hwwa.get_approach_avoid_mask(labs)) );

seq_labs = fcat();
ns = 1:3;
make_incorrect = false( size(ns) );
ns(end+1) = 1;
make_incorrect(end+1) = true;

allow_repeated = false;
allow_duplicates = false;

if ( ~allow_duplicates )
  ns = sort( ns, 'descend' );
  marked = false( rows(labs), 1 );
  miss_label = 'no_prev_correct_sequence';
end

for i = 1:numel(ns)
  if ( make_incorrect(i) )
    corr_lab = 'correct_false';
  else
    corr_lab = 'correct_true';
  end
  
  go_n = hwwa_n_go_trial_sequence( labs', each_I, 'go_trial', corr_lab, ns(i), allow_repeated );
  nogo_n = hwwa_n_go_trial_sequence( labs', each_I, 'nogo_trial', corr_lab, ns(i), allow_repeated ); 
  
  if ( ~allow_duplicates )
    found_seq = union( findnone(go_n, miss_label), findnone(nogo_n, miss_label) );
    keep_seq = intersect( found_seq, find(~marked) );
    rm_seq = setdiff( found_seq, keep_seq );
    
    setcat( go_n, 'prev_correct_sequence', miss_label, rm_seq );
    setcat( nogo_n, 'prev_correct_sequence', miss_label, rm_seq );
    
    assert( ~any(marked(found_seq)) );
    marked(found_seq) = true;
  end
  
  append( seq_labs, go_n );
  append( seq_labs, nogo_n );
end

find_incorrect1 = @(l, m) find(l, {'prev_correct_false_nogo_trial_1', 'prev_correct_false_go_trial_1'}, m);
% setcat( seq_labs, 'prev_correct_sequence', 'prev_correct_false_1', incorr_ind0 );

%%  plot go trial sequence

pg_lab_ = @(kind, n) sprintf( 'prev_%s_go_trial_%d', kind, n );
png_lab_ = @(kind, n) sprintf( 'prev_%s_nogo_trial_%d', kind, n );

pg_lab = @(n) pg_lab_('correct_true', n);
png_lab = @(n) png_lab_('correct_true', n);

make_cs = @(n) { ...
  {pg_lab(n), 'go_trial'}, {pg_lab(n), 'nogo_trial'} ...
  , {png_lab(n), 'go_trial'}, {png_lab(n), 'nogo_trial'} ...
};

cs = cat_expanded( 1, [make_cs(1), make_cs(2), make_cs(3)] );

% mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
%   , @(m) union(find_incorrect1(l, m), hwwa.find_many(l, cs, m)) ...
% );

% mask_func = find_non_outliers;

has_len = findnone( seq_labs, '<prev_sequence_length>' );
seq_len = cellstr( seq_labs, 'prev_sequence_length', has_len );
prev_type = cellstr( seq_labs, 'prev_correct', has_len );
joined_types = cellfun( @(x, y) sprintf('%s_%s', x, y), seq_len, prev_type, 'un', 0 );
addsetcat( seq_labs, 'prev_correct_sequence_length', joined_types, has_len );

mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) findnone(l, '<prev_correct_sequence_length>', m) ...
  , @(m) findnone(l, 'prev_3', m) ...
  , @(m) find(l, 'prev_correct_true', m) ...
  , @(m) hwwa.find_5htp(l, m) ...
  , @(m) hwwa.find_go(l, m) ...
);

%   , @(m) findnone(l, 'prev_1_prev_correct_false', m) ...

pcorr_outs = hwwa_overall_p_correct( seq_labs' ...
  , 'mask_func', mask_func ...
  , 'do_save', false ...
  , 'pcorr_each', {'unified_filename', 'monkey', 'drug' ...
                  , 'prev_correct_sequence_length', 'prev_trial_type', 'trial_type'} ...
  , 'anova_each', {} ...
  , 'anova_factors', {'drug', 'monkey', 'trial_type', 'switch_trial_type', 'prev_correct_sequence_length'} ...
  , 'norm_anova_factors', {'trial_type', 'prev_trial_type'} ...
  , 'norm_anova_each', {'prev_correct_sequence_length'} ...
  , 'per_drug', true ...
  , 'per_target_image_category', false ...
  , 'include_normalized', true ...
  , 'include_raw', true ...
  , 'errorbar_cats', {{'prev_correct_sequence'}, {'drug'}, {'trial_type'}} ...
  , 'norm_t_elements', {'prev_go_trial', 'prev_nogo_trial'} ...
  , 'norm_t_each', {'trial_type', 'prev_correct_sequence_length'} ...
  , 'points_are', {'monkey'} ...
);

anova_outs = dsp3.anovan( pcorr_outs.raw_corr, pcorr_outs.raw_labels' ...
  , {'trial_type'}, {'prev_correct_sequence_length', 'prev_trial_type'} );

if ( true )
  data_tbl = array2table( pcorr_outs.raw_corr, 'VariableNames', {'pcorr'} );
  label_tbl = array2table( cellstr(pcorr_outs.normalized_labels) ...
    , 'VariableNames', getcats(pcorr_outs.raw_labels) );
  save_table_pair( data_tbl, label_tbl, 'trial_sequence_pcorr_pp', '' );
end

%%  p correct vs num initiated

mask_func = find_non_outliers;

hwwa_num_initiated_vs_p_correct( use_labs' ...
  , 'mask_func', mask_func ...
  , 'each', {'drug', 'monkey', 'unified_filename'} ...
  , 'plot_cats',  {{'drug'}, {'monkey'}} ...
  , 'do_save', true ...
  , 'normalize', true ...
  , 'marker_size', 10 ...
  , 'per_panel_corr', true ...
);

%%  pupil vs percent correct

mask_func = find_non_outliers;

hwwa_pupil_vs_percent_correct( behav_outs.pupil_size, use_labs' ...
  , 'mask_func', mask_func ...
  , 'corr_each', {} ...
  , 'pupil_each', {'unified_filename', 'monkey', 'drug', 'trial_type'} ...
  , 'pcorr_each', {'unified_filename', 'monkey', 'drug', 'trial_type'} ...
  , 'corr_match', {'unified_filename', 'monkey', 'drug', 'trial_type'} ...
  , 'scatter_cats', {{'drug'}, {'trial_type'}} ...
  , 'normalize', true ...
  , 'do_save', true ...
  , 'per_panel_corr', true ...
);

%%  rt vs percent correct

% Only correct go.
rt_mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) find(l, {'go_trial', 'correct_true'}, m) ...
);

% Only go.
pcorr_mask_func = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) find(l, {'go_trial'}, m) ...
);

hwwa_rt_vs_p_correct( rt, use_labs' ...
  , 'pcorr_mask_func', pcorr_mask_func ...
  , 'rt_mask_func', rt_mask_func ...
  , 'do_save', true ...
  , 'normalize', true ...
  , 'p_format', '%0.4f' ...
  , 'plot_cats', {{'drug'}, {'trial_type', 'monkey'}} ...
  , 'per_panel_corr', true ...
);

%%  overall pupil size

pupil_outs = hwwa_overall_pupil_size( behav_outs.pupil_size, use_labs' ...
  , 'mask_func', find_non_outliers ...
  , 'do_save', false ...
  , 'per_monkey', false ...
  , 'per_scrambled_type', false ...
  , 'per_trial_type', false ...
  , 'per_correct', false ...
  , 'normalize', false ...
  , 'add_points', true ...
  , 'points_are', 'monkey' ...
  , 'run_level_average', true ...
  , 'signrank_inputs', {1} ... % test against 1.
);

rs_outs = dsp3.ranksum( pupil_outs.pupil_size, pupil_outs.labels', {}, '5-htp', 'saline' );

if ( false )
  data_tbl = array2table( pupil_outs.pupil_size, 'VariableNames', {'pupil_size'} );
  label_tbl = array2table( cellstr(pupil_outs.labels), 'VariableNames', getcats(pupil_outs.labels) );
  save_table_pair( data_tbl, label_tbl, 'pupil_size', '' );
end

%%  overall rt

mask_func = find_non_outliers;

include_nogo = false;

per_cats = true;
per_monks = false;
per_scrambled_types = false;
norms = false;
trial_levels = false;

cs = dsp3.numel_combvec( per_cats, per_monks, per_scrambled_types, norms ...
  , trial_levels );

for i = 1:size(cs, 2)
  
per_cat = per_cats(cs(1, i));
per_monk = per_monks(cs(2, i));
per_scrambled_type = per_scrambled_types(cs(3, i));
do_norm = norms(cs(4, i));
trial_level = trial_levels(cs(5, i));

norm_lims = {[], []};
non_norm_lims = {[0, 0.4], [0.1, 0.3]};

lims = ternary( do_norm, norm_lims, non_norm_lims );

if ( trial_level )
  y_lims = lims{1};
else
  y_lims = lims{2};
end

plt_labs = use_labs'; 
apply_social_housing_and_gender_labels( plt_labs );

rt_outs = hwwa_overall_rt( rt, plt_labs' ...
  , 'mask_func', mask_func ...
  , 'do_save', false ...
  , 'per_monkey', per_monk ...
  , 'per_scrambled_type', per_scrambled_type ...
  , 'per_image_category', per_cat ...
  , 'include_nogo', include_nogo ...
  , 'trial_level', trial_level ...
  , 'normalize', do_norm ...
  , 'add_points', true ...
  , 'points_are', 'monkey' ...
  , 'y_lims', y_lims ...
);

end

%%  anova overall rt

use_nested = false;

if ( use_nested )
  nesting = zeros( 3 );
  nesting(2, 1) = true;

  anova_factors = {'social-housing', 'gender', 'scrambled_type'};
  anova_outs = dsp3.anovan( rt_outs.norm_rt, rt_outs.norm_rt_labels ...
    , {}, anova_factors ...
    , 'anovan_inputs', {'display', 'off', 'varnames', anova_factors, 'model', 'full', 'nested', nesting} ...
    , 'run_multcompare', false ...
  );
else
  anova_factors = { 'gender', 'scrambled_type'};
  anova_outs = dsp3.anovan( rt_outs.norm_rt, rt_outs.norm_rt_labels ...
    , {'social-housing'}, anova_factors ...
  );
end

%%  export overall rt

tbl_labels = array2table( cellstr(rt_outs.rt_labels), 'VariableNames', getcats(rt_outs.rt_labels) );
tbl_data = array2table( rt_outs.rt, 'VariableNames', {'rt'} );

save_table_pair( tbl_data, tbl_labels, 'rt_uu2', '' );

%%  p correct vs. rt over time

hwwa_scatter_rt_p_correct_over_time( ...
  pcorr, pcorr_labels', mean_rt, mean_rt_labels', slide_each ...
  , 'permutation_test', false ...
  , 'permutation_test_iters', 1e2 ...
  , 'do_save', true ...
  , 'per_monkey', false ...
  , 'per_rt_quantile', true ...
  , 'mask_func', @(l, m) findnone(l, 'rt-quantile__NaN') ...
  , 'remove_outliers', false ...
  , 'std_threshold', 2 ...
  , 'order_each_func', @(slide_each) setdiff(slide_each, 'rt-quantile') ...
  , 'x_lims', [0.14, 0.3] ...
  , 'y_lims', [0, 1.2] ...
  , 'anova_each', {} ...
  , 'normalize_pcorr', true ...
);

%%  p correct vs. rt over time

hwwa_scatter_rt_p_correct_over_time( ...
  time_props, time_prop_labels, time_rt, time_rt_labels, slide_each ...
);

%%  p correct or rt over time

use_rt_for_over_time = true;

if ( use_rt_for_over_time )
  use_props = time_rt;
  use_prop_edges = time_rt_edges;
  use_prop_labels = time_rt_labels';
  over_time_prefix = 'rt--';
else
  use_props = time_props;
  use_prop_edges = time_prop_edges;
  use_prop_labels = time_prop_labels';
  over_time_prefix = 'pcorr--';
end

hwwa_plot_time_binned_behavior( use_props, use_prop_edges, use_prop_labels' ...
  , 'time_limits', [-inf, 3.5e3] ...
  , 'do_save', true ...
  , 'per_monkey', false ...
  , 'per_drug', true ...
  , 'prefix', over_time_prefix ...
);

%%  num initiated per session

% per_scrambled_types = trufls();
per_scrambled_types = false;
% per_drugs = trufls();
per_drugs = true;
% per_trial_types = trufls();
per_trial_types = false;
cs = dsp3.numel_combvec( per_scrambled_types, per_drugs, per_trial_types );

for i = 1:size(cs, 2)
  
shared_utils.general.progress( i, size(cs, 2) );
  
c = cs(:, i);

hwwa_plot_num_initiated_per_session( use_labs' ...
  , 'per_scrambled_type', per_scrambled_types(c(1)) ...
  , 'per_drug', per_drugs(c(2)) ...
  , 'per_trial_type', per_trial_types(c(3)) ...
  , 'do_save', true ...
  , 'mask_func', find_non_outliers ...
  , 'permutation_test', false ...
  , 'compare_drug', true ...
);

end

%%  p correct vs num initiated next

per_scrambled_types = trufls();
per_drugs = true;
per_trial_types = trufls();
cs = dsp3.numel_combvec( per_scrambled_types, per_drugs, per_trial_types );

for i = 1:size(cs, 2)
  
shared_utils.general.progress( i, size(cs, 2) );
  
c = cs(:, i);

hwwa_plot_n_plus_one_performance( use_labs' ...
  , 'per_scrambled_type', per_scrambled_types(c(1)) ...
  , 'per_drug', per_drugs(c(2)) ...
  , 'per_trial_type', per_trial_types(c(3)) ...
  , 'do_save', true ...
  , 'mask_func', find_non_outliers ...
);

end

%%  amp-velocity tradeoff

vel = nan( rows(use_labs), 1 );
amp = nan( size(vel) );

to_sacc_ind = saccade_outs.aligned_to_saccade_ind;

vel(to_sacc_ind) = saccade_outs.saccade_peak_velocities;
amp(to_sacc_ind) = saccade_outs.saccade_lengths;

vel = hwwa.convert_from_old_degrees_to_corrected_degrees( vel );
amp = hwwa.convert_from_old_degrees_to_corrected_degrees( amp );

correct_go_mask = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_correct_go(l, m) ...
);

incorrect_nogo_mask = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_incorrect_nogo(l, m) ...
);

both_mask = @(l, m) pipe(find_non_outliers(l, m) ...
  , @(m) hwwa.find_correct_go_incorrect_nogo(l, m) ...
);

hwwa_plot_amp_vel_tradeoff( amp, vel, use_labs' ...
  , 'mask_func', correct_go_mask ...
  , 'gcats', {'drug'} ...
  , 'pcats', {'scrambled_type'} ...
  , 'per_monkey', false ...
  , 'do_save', true ...
  , 'permutation_test', true ...
  , 'permutation_test_iters', 1e3 ...
  , 'base_subdir', 'across_all' ...
);

%%  relate pupil size with rt, saccade, velocity, amplitude

% types = { 'amp', 'vel', 'rt' };
% types = { 'vel' };

first_types = { 'amp' };
% types = { 'pupil_max_diff' };
types = { 'vel' };

per_monks = [false];
per_scrambled_types = [true];
cs = dsp3.numel_combvec( first_types, types, per_monks, per_scrambled_types );

use_pupil_size = behav_outs.pupil_size;
% use_pupil_size = cue_on_max_diff;

for i = 1:size(cs, 2)
  
first_measure_type = first_types{cs(1, i)};
second_measure_type = types{cs(2, i)};
per_monk = per_monks(cs(3, i));
per_scrambled_type = per_scrambled_types(cs(4, i));

switch ( first_measure_type )
  case 'pupil'
    first_measure = use_pupil_size;
  case 'amp'
    firt_measure = amp;
  otherwise
    error( 'Unrecognized first measure type "%s".', first_measure_type );
end

switch ( second_measure_type )
  case 'rt'
    second_measure = rt;
  case 'vel'
    second_measure = vel;
  case 'amp'
    second_measure = amp;
  case 'pupil_max_diff'
    second_measure = cue_on_max_diff;
  otherwise
    error( 'Unrecognized second measure type "%s".', second_measure_type );
end

mask_func = ...
  @(l, m) hwwa.find_correct_go_incorrect_nogo(l, find_non_outliers(l, m));

hwwa_relate_to_pupil_size( first_measure, second_measure, use_labs' ...
  , 'second_measure_name', second_measure_type ...
  , 'mask_func', mask_func ...
  , 'do_save', true ...
  , 'per_monkey', per_monk ...
  , 'per_scrambled_type', per_scrambled_type ...
  , 'permutation_test', true ...
  , 'remove_x_outliers', false ...
  , 'outlier_std_thresh', 2 ...
);

end

%%  Pupil traces

pup_labs = use_labs';

use_quantiles = false;

use_aligned_outs = behav_outs.cue_on_aligned;
use_aligned_outs.labels = pup_labs;

quantile_measure = rt;
quantile_measure(quantile_measure == 0) = nan;

num_tiles = 3;
quant_cat = 'rt-quantile';
quant_each = { 'unified_filename', 'trial_type' };
quant_mask_func = @(l, m) fcat.mask(l, find_non_outliers(l, m) ...
  , @(l, arg, m) intersect(arg, m), find(~isnan(quantile_measure)) ...
);
quant_mask_func = @(l, m) ...
  hwwa.find_correct_go_incorrect_nogo(l, quant_mask_func(l, m));

quant_mask = hwwa.get_approach_avoid_base_mask( pup_labs, quant_mask_func );

quant_inds = dsp3.quantiles_each( quantile_measure, pup_labs' ...
  , num_tiles, quant_each, {}, quant_mask );
quant_labs = arrayfun( @(x) sprintf('%s-%d', quant_cat, x), quant_inds, 'un', 0 );
addsetcat( pup_labs, quant_cat, quant_labs );

thresh = 0;
require_crit = false;

if ( use_quantiles )
  wrap_mask_func = @(l, m) fcat.mask(l, find_non_outliers(l, m) ...
    , @findnone, 'rt-quantile-NaN' ...
  );
  pupil_fcats = quant_cat;
else
  wrap_mask_func = find_non_outliers;
  pupil_fcats = {};
end

prefix = ternary( require_crit, sprintf('crit-%d', thresh), '' );

hwwa_plot_pupil_traces( use_aligned_outs ...
  , 'mask_func', wrap_mask_func ...
  , 'time_limits', [0, 800] ...
  , 'do_save', true ...
  , 'smooth_func', @(x) smoothdata(x, 'SmoothingFactor', 0.25) ...
  , 'fcats', pupil_fcats ...
  , 'prefix', prefix ...
  , 'compare_series', true ...
  , 'formats', {'svg'} ...
);

%%  trial duration

trial_dur_labs = use_labs';
start_ts = behav_outs.run_relative_start_times;
trial_durations = get_trial_duration( start_ts, trial_dur_labs );

mask_func = find_non_outliers;

hwwa_plot_trial_durations( trial_durations, trial_dur_labs ...
  , 'mask_func', mask_func ...
  , 'each', {'unified_filename', 'drug', 'correct'} ...
  , 'do_save', true ...
);

%%

function trial_durations = get_trial_duration(start_ts, labels)

assert_ispair( start_ts, labels );
trial_durations = nan( size(start_ts) );

run_I = findall_or_one( labels, 'unified_filename' );
for i = 1:numel(run_I)
  run_ts = start_ts(run_I{i});
  trial_durs = diff( run_ts );
  trial_durs(end+1) = nan;
  trial_durations(run_I{i}) = trial_durs;
end

end

function [data, out_labels] = trial_bin_apply(labels, func, each_I, bin_size, step_size)

data = {};
out_labels = fcat();

for i = 1:numel(each_I)
  bin_inds = hwwa.bin_indices( each_I{i}, bin_size, step_size );
  append1( out_labels, labels, each_I{i} );
  
  for j = 1:numel(bin_inds)
    data{i, j} = func( labels, bin_inds{j} );
  end
end

end

function [time_between, time_between_labels] = ...
  time_between_initiated(events, labels, event_key)

assert( numel(event_key) == size(events, 2) );
assert_ispair( events, labels );

[~, event_inds] = ismember( {'fixation_on', 'go_nogo_cue_onset'}, event_key );
assert( ~any(event_inds == 0), 'Missing some required events.' );

fix_on_ind = event_inds(1);
gng_ind = event_inds(2);

run_I = findall( labels', 'unified_filename' );

time_between = [];
time_between_labels = fcat();

for i = 1:numel(run_I)
  run_ind = run_I{i};
  j = 1;
  last_initiated_ind = 0;
  
  while ( j < numel(run_ind) )
    trial_ind = run_ind(j);
    did_init = strcmp( cellstr(labels, 'initiated', trial_ind), 'initiated_true' );
    
    if ( did_init )
      if ( last_initiated_ind ~= 0 )
        % Last trial offset is fixation onset of the subsequent trial.
        last_init_run_ind = run_ind(last_initiated_ind+1);
        last_trial_offset = events(last_init_run_ind, fix_on_ind);

        curr_trial_onset = events(trial_ind, gng_ind);

        delta = curr_trial_onset - last_trial_offset;
        assert( delta >= 0 );

        time_between(end+1, 1) = delta;
        append( time_between_labels, labels, trial_ind );
      end
      
      last_initiated_ind = j;
    end
    
    j = j + 1;
  end
end

assert_ispair( time_between, time_between_labels );
prune( time_between_labels );

end

function fix_info = find_fixations(all_sacc_info, len, dur_thresh)

is_sacc = false( 1, len );

for i = 1:rows(all_sacc_info)
  start = all_sacc_info(i, 1);
  stop = all_sacc_info(i, 2);
  is_sacc(start:stop) = true;
end

is_fix = ~is_sacc;
[fix_starts, fix_durs] = shared_utils.logical.find_islands( is_fix );

too_short = fix_durs < dur_thresh;
fix_starts(too_short) = [];
fix_durs(too_short) = [];

fix_info = [ ...
    fix_starts(:) ...
  , fix_starts(:) + fix_durs(:) - 1 ...
  , nan(numel(fix_starts), size(all_sacc_info, 2)-2) ...
];

end

function degs = normalized_edges_to_degrees(edges, max_v)

pxs = edges * max_v;
degs_old = hwwa.run_px2deg( pxs, [] );
degs = hwwa.convert_from_old_degrees_to_corrected_degrees( degs_old );

end

function save_table_pair(data_tbl, label_tbl, subdir, fname)

assert( size(data_tbl, 1) == size(label_tbl, 1) );

params = hwwa.get_common_make_defaults( hwwa.get_common_plot_defaults() );
save_p = fullfile( hwwa.approach_avoid_data_path(params, 'export'), subdir );

shared_utils.io.require_dir( save_p );
dsp3.writetable( data_tbl, fullfile(save_p, sprintf('data__%s.csv', fname)) );
dsp3.writetable( label_tbl, fullfile(save_p, sprintf('labels__%s.csv', fname)) );

end

function plt_labs = apply_social_housing_and_gender_labels(plt_labs)

setcat( plt_labs, 'gender', 'male', find(plt_labs, {'tar', 'hitch'}) );
setcat( plt_labs, 'gender', 'female', find(plt_labs, {'ephron'}) );

addcat( plt_labs, 'social-housing' );
setcat( plt_labs, 'social-housing', 'socially-housed', find(plt_labs, 'tar') );
setcat( plt_labs, 'social-housing', 'not-socially-housed', find(plt_labs, {'ephron', 'hitch'}) );

end