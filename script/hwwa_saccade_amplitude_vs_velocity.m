function hwwa_saccade_amplitude_vs_velocity(amp_vel, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @(l, m) m;
defaults.iters = 1e3;
defaults.test_each = {};
defaults.seed = [];
defaults.std_threshold_num_devs = 2;
defaults.mean_each = {};
defaults.marker_size = [];

assert_ispair( amp_vel, labels );

params = hwwa.parsestruct( defaults, varargin);

base_mask = intersect( get_base_mask(labels, params.mask_func), nonnan_indices(amp_vel) );
base_mask = find_all_non_outliers( amp_vel, labels, params.test_each, base_mask, params.std_threshold_num_devs );

if ( ischar(params.mean_each) || ~isempty(params.mean_each) )
  [amp_vel, labels, base_mask] = mean_each( amp_vel, labels', params.mean_each, base_mask );
end

[ps, slope_diffs, labs, mask] = compare_each_amp_vs_vel_slopes( amp_vel, labels, base_mask, params );

if ( params.do_save )
  make_save_table( ps, slope_diffs, labs, params );
end

plot_amp_vs_vel( amp_vel, labels, mask, params );

end

function [data, labels, mask] = mean_each(data, labels, each, mask)

[labels, mean_I] = keepeach( labels', each, mask );
data = bfw.row_nanmean( data, mean_I );
mask = rowmask( data );

end

function keep_ind = find_all_non_outliers(amp_vel, labels, each, mask, num_devs)

each_I = findall_or_one( labels, each, mask );
keep_I = cell( size(each_I) );

for i = 1:numel(each_I)
  vel_non_outliers = find_non_outliers( amp_vel(:, 1), each_I{i}, num_devs );
  amp_non_outliers = find_non_outliers( amp_vel(:, 2), each_I{i}, num_devs );

  keep_I{i} = intersect( vel_non_outliers, amp_non_outliers );
end

keep_ind = unique( vertcat(keep_I{:}) );

end

function make_save_table(ps, slope_diffs, labs, params)

row_labels = fcat.strjoin( combs(labs, params.test_each), ' | ' );
if ( isempty(row_labels) )
  row_labels = { '<undefined>' };
end

col_labels = { 'p', 'slope-diff-sal-minus-5htp' };

t = fcat.table( [ps, slope_diffs], row_labels, col_labels );

save_p = hwwa.approach_avoid_data_path( params, 'plots', 'saccade-amp-vs-vel', 'stats' );
dsp3.req_writetable( t, save_p, labs, params.test_each );

end

function [ps, slope_diffs, out_labels, mask] = compare_each_amp_vs_vel_slopes(amp_vel, labels, mask, params)

each_I = findall_or_one( labels, params.test_each, mask );

ps = nan( numel(each_I), 1 );
slope_diffs = nan( size(ps) );
out_labels = fcat();
mask = cell( size(each_I) );

for i = 1:numel(each_I)
  if ( ~isempty(params.seed) )
    rng_state = rng();
    rng( params.seed );
  end
  
  [p, slope_diff, labs, tmp_mask] = compare_amp_vs_vel_slopes( amp_vel, labels, each_I{i}, params );

  ps(i) = p;
  slope_diffs(i) = slope_diff;
  append( out_labels, labs );
  mask{i} = tmp_mask;
  
  if ( ~isempty(params.seed) )
    rng( rng_state );
  end
end

mask = vertcat( mask{:} );

end

function ind = find_non_outliers(data, mask, num_devs)

mean_dat = nanmean( data(mask) );
dev_dat = nanstd( data(mask) );

below_thresh = data < mean_dat - dev_dat*num_devs;
above_thresh = data > mean_dat + dev_dat*num_devs;

ind = intersect( find(~below_thresh & ~above_thresh), mask );

end

function [p, slope_diff, out_labels, out_mask] = compare_amp_vs_vel_slopes(amp_vel, labels, mask, params)

saline = 'saline';
drug = '5-htp';

% vel_non_outliers = find_non_outliers( amp_vel(:, 1), mask, params.std_threshold_num_devs );
% amp_non_outliers = find_non_outliers( amp_vel(:, 2), mask, params.std_threshold_num_devs );
% 
% mask = intersect( vel_non_outliers, amp_non_outliers );
out_mask = mask;

sal_ind = find( labels, saline, mask );
drug_ind = find( labels, drug, mask );
all_inds = sort( [sal_ind; drug_ind] );

num_inds = count( labels, {saline, drug}, all_inds );
real_abs_diff = model_beta_difference( amp_vel, sal_ind, drug_ind );
real_signed_diff = model_beta_difference( amp_vel, sal_ind, drug_ind, false );

ps = nan( params.iters, 1 );

for i = 1:params.iters
  ind = all_inds(randperm(numel(all_inds)));
  
  sal_ind = ind(1:num_inds(1));
  drug_ind = ind(num_inds(1)+1:end);
  
  beta_diff = model_beta_difference( amp_vel, sal_ind, drug_ind );
  
  ps(i) = beta_diff > real_abs_diff;
end

p = sum( ps ) / params.iters;
slope_diff = real_signed_diff;
out_labels = append1( fcat, labels, all_inds );

end

function beta_diff = model_beta_difference(amp_vel, saline_ind, drug_ind, take_abs)

if ( nargin < 4 )
  take_abs = true;
end

sal = fit_model( amp_vel, saline_ind );
drug = fit_model( amp_vel, drug_ind );

beta_diff = beta1( sal ) - beta1( drug );

if ( take_abs )
  beta_diff = abs( beta_diff );
end

end

function b = nth_beta(model, num)

b = model.Coefficients.Estimate(num);

end

function b = beta1(model)

b = nth_beta( model, 2 );

end

function mdl = fit_model(amp_vel, mask)

mdl = fitlm( amp_vel(mask, 1), amp_vel(mask, 2) );

end

function inds = nonnan_indices(amp_vel)

inds = find( all(~isnan(amp_vel), 2) );

end

function plot_amp_vs_vel(amp_vel, labels, mask, params)

%%

x = amp_vel(mask, 1);
y = amp_vel(mask, 2);
labs = prune( labels(mask) );

pl = plotlabeled.make_common();

if ( ~isempty(params.marker_size) )
  pl.marker_size = params.marker_size;
end

gcats = { 'drug' };
pcats = params.test_each;

[axs, ids] = pl.scatter( x, y, labs, gcats, pcats );
% plotlabeled.scatter_addcorr( ids, x, y );

for i = 1:numel(ids)
  sx = x(ids(i).index);
  sy = y(ids(i).index);
  mdl = fitlm( sx, sy );
  intercept = nth_beta( mdl, 1 );
  beta1 = nth_beta( mdl, 2 );
  ax = ids(i).axes;
  xs = get( ax, 'xtick' );
  ys = polyval( [beta1, intercept], xs );
  hold( ax, 'on' );
  h = plot( ax, xs, ys );
  set( h, 'linewidth', 4 );
  set( h, 'linestyle', '--' );
  line_col = unique( ids(i).series.CData, 'rows' );
  line_col(3) = 0;
  set( h, 'Color', line_col ); 
end

shared_utils.plot.match_xlims( axs );
shared_utils.plot.match_ylims( axs );
shared_utils.plot.xlabel( axs(1), 'Saccade amplitude (deg)' );
shared_utils.plot.ylabel( axs(1), 'Peak saccade velocity (deg/s)' );

if ( params.do_save )
  shared_utils.plot.fullscreen( gcf );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'saccade-amp-vs-vel' );
  dsp3.req_savefig( gcf, save_p, labs, params.test_each );
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels, {} );
mask = mask_func( labels, mask );

end