function hwwa_plot_approach_avoid_behav(rt, labels, varargin)

assert_ispair( rt, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @(labels, mask) mask;
defaults.normalize = false;
defaults.match_figure_limits = true;
defaults.per_monkey = true;
defaults.per_image_category = true;
defaults.rt_kind = 'reaction_time';
defaults.y_lims = [];

params = hwwa.parsestruct( defaults, varargin );

if ( ~hascat(labels, 'scrambled_type') )
  hwwa.decompose_social_scrambled( labels );
end

[rt, labels] = ensure_correct_subset( rt, labels' );

mask = get_base_mask( labels, params.mask_func, true );
any_initiated_mask = get_base_mask( labels, params.mask_func, false );

% num_initiated_per_session( labels', any_initiated_mask, params );

rt_per_image_cat( rt, labels', mask, params );
rt_variance_per_image_cat( rt, labels', mask, params );

% p_correct_social_minus_scrambled( labels', mask, params );
% p_correct_social_scrambled( labels', mask, params );
% p_correct_per_image_cat( labels', mask, params );
% p_correct_merged_image_cats( labels', mask, params );

% p_correct_per_image_cat_stacked( labels', mask, params );
% p_correct_variance_per_image_cat( labels', mask, params );
% rt_per_image_cat( rt, labels', mask, params );

end

function subdir = make_subdir(subdir, params)

norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );
subdir = sprintf( '%s-%s', subdir, norm_subdir );
subdir = fullfile( subdir, params.rt_kind );

end

function rt_variance_per_image_cat(rt, labels, mask, params)

mask = get_rt_mask( labels, mask );

each = { 'drug', 'trial_type', 'scrambled_type' ...
  , 'unified_filename', 'monkey' };

if ( params.per_image_category )
  each{end+1} = 'target_image_category';
end

fcats = {};
pcats = { 'trial_type', 'drug' };
gcats = { 'scrambled_type', 'target_image_category' };

subdir = make_subdir( 'rt_variance', params );

[var_labels, mean_I] = keepeach( labels', each, mask );
var_rt = rowop( rt, mean_I, @(x) nanvar(x, [], 1) );

if ( params.normalize )
  norm_each = setdiff( each, {'unified_filename', 'drug'} );
  [var_rt, var_labels] = hwwa.saline_normalize( var_rt, var_labels', norm_each );
end

pl = plotlabeled.make_common();
box_plot( pl, var_rt, var_labels', fcats, gcats, pcats, subdir, params );

anova_each = { 'trial_type' };
anova_factors = { 'scrambled_type', 'target_image_category' };

if ( ~params.normalize )
  anova_factors{end+1} = 'drug';
end

anova_factors = intersect( anova_factors, each ); 

try
  anova_stats( var_rt, var_labels', anova_each, anova_factors, rowmask(var_rt), subdir, params );
catch err
  warning( err.message );
end

end

function num_initiated_per_session(labels, any_init_mask, params)

each = { 'unified_filename' };
[init_labels, each_I] = keepeach( labels', each, any_init_mask );
num_init = zeros( numel(each_I), 1 );

for i = 1:numel(each_I)
  num_init(i) = numel( find(labels, 'initiated_true', each_I{i}) );
end

pl = plotlabeled.make_common();
gcats = { 'monkey' };
pcats = { 'drug' };

if ( params.do_save )
  box_plot( pl, num_init, init_labels', {}, gcats, pcats, 'num_initiated', params );
end

end

function rt_per_image_cat(rt, labels, mask, params)

%%

mask = get_rt_mask( labels, mask );

each = { 'drug', 'trial_type', 'scrambled_type', 'unified_filename', 'monkey' };

% fcats = { 'trial_type' };
fcats = {};

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
end
if ( params.per_image_category )
  each{end+1} = 'target_image_category';
end

pcats = { 'trial_type', 'monkey', 'scrambled_type', 'target_image_category' };
gcats = { 'drug' };

subdir = make_subdir( 'rt_boxplot', params );

[mean_labels, mean_I] = keepeach( labels', each, mask );
mean_rt = bfw.row_nanmean( rt, mean_I );

if ( params.normalize )
  norm_each = setdiff( each, {'unified_filename', 'drug'} );
  [mean_rt, mean_labels] = hwwa.saline_normalize( mean_rt, mean_labels', norm_each );
end

if ( ~ismember('monkey', fcats) )
  collapsecat( mean_labels, 'monkey' );
end

gcats = intersect( gcats, each );

pl = plotlabeled.make_common();
pl.prefer_multiple_groups = true;
box_plot( pl, mean_rt, mean_labels', fcats, gcats, pcats, subdir, params );

xcats = { 'scrambled_type' };
pl.prefer_multiple_groups = false;
regular_bar( pl, mean_rt, mean_labels', fcats, xcats, setdiff(gcats, xcats) ...
  , setdiff(pcats, xcats), strrep(subdir, 'boxplot', 'bar'), params );

anova_each = {};
anova_factors = { 'scrambled_type', 'target_image_category', 'trial_type' };

if ( ~params.normalize )
  anova_factors{end+1} = 'drug';
end

if ( ismember('monkey', fcats) )
  anova_factors{end+1} = 'monkey';
end

anova_factors = intersect( anova_factors, each );
anova_stats( mean_rt, mean_labels', anova_each, anova_factors, rowmask(mean_rt), subdir, params );

post_hoc_factors = intersect( anova_factors, {'scrambled_type', 'target_image_category'} );
post_hoc_each = { 'trial_type' };
anova_stats( mean_rt, mean_labels', post_hoc_each, post_hoc_factors, rowmask(mean_rt), subdir, params );

end

function p_correct_social_minus_scrambled(labels, mask, params)

%%

each = { 'drug', 'trial_type', 'scrambled_type', 'unified_filename', 'target_image_category' };
  
fcats = {};
pcats = { 'trial_type', 'scrambled_type', 'target_image_category' };
gcats = { 'drug' };

norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );

subdir = 'percent_correct_social_minus_scrambled_boxplot';
subdir = sprintf( '%s-%s', subdir, norm_subdir );

soc_norm_each = setdiff( each, {'scrambled_type'} );
drug_norm_each = setdiff( each, {'unified_filename'} );

[props, prop_labels] = hwwa.percent_correct( labels, each, mask );
[props, prop_labels] = hwwa.social_divide_scrambled( props, prop_labels', soc_norm_each );
[props, prop_labels] = hwwa.saline_normalize( props, prop_labels', drug_norm_each );

pl = plotlabeled.make_common();
box_plot( pl, props, prop_labels, fcats, gcats, pcats, subdir, params );

% anova_factors = { 'target_image_category', 'drug' };
% anova_each = { 'trial_type' };

% anova_stats( props, prop_labels', anova_each, anova_factors, rowmask(props), subdir, params)

ttest_outs = dsp3.ttest2( props, prop_labels, 'trial_type', 'appetitive', 'threat' );

end

function p_correct_social_scrambled(labels, mask, params)

for i = 1:2
  each = { 'drug', 'trial_type', 'scrambled_type', 'unified_filename' };
  
  fcats = {};
  pcats = { 'trial_type', 'drug' };
  gcats = { 'scrambled_type' };
  
  if ( i == 2 )
    each{end+1} = 'run_time_quantile';
    pcats{end+1} = 'run_time_quantile';
  end

  norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );

  subdir = 'percent_correct_boxplot';
  subdir = sprintf( '%s-%s', subdir, norm_subdir );

  [props, prop_labels] = hwwa.percent_correct( labels, each, mask );

  if ( params.normalize )
    norm_each = setdiff( each, {'unified_filename', 'drug'} );
    [props, prop_labels] = hwwa.saline_normalize( props, prop_labels', norm_each );
  end

  collapsecat( prop_labels, 'monkey' );

  pl = plotlabeled.make_common();

  box_plot( pl, props, prop_labels, fcats, gcats, pcats, subdir, params );
end

end

function p_correct_merged_image_cats(labels, mask, params)

merge_scrambled_type_target_image_category( labels, mask );

each = { 'target_image_category', 'drug', 'trial_type', 'scrambled_type', 'unified_filename' };  
norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );

subdir = fullfile( 'percent_correct_3_categories'  );
subdir = sprintf( '%s-%s', subdir, norm_subdir );

[props, prop_labels] = hwwa.percent_correct( labels, each, mask );

if ( params.normalize )
  norm_each = setdiff( each, {'unified_filename', 'drug'} );
  [props, prop_labels] = hwwa.saline_normalize( props, prop_labels', norm_each );
end
  
fcats = { 'monkey' };
pcats = { 'monkey', 'trial_type' };
  
if ( ~params.per_monkey )
  collapsecat( prop_labels, 'monkey' );
end

xcats = { 'target_image_category' };
gcats = { 'drug' };

pl = plotlabeled.make_common();

anova_each = {};
anova_factors = { 'target_image_category', 'trial_type' };

if ( ~params.per_monkey && ~params.normalize )
  anova_factors{end+1} = 'monkey';
end    
if ( ~params.normalize )
  anova_factors{end+1} = 'drug';
end
  
anova_stats( props, prop_labels, anova_each, anova_factors ...
  , rowmask(prop_labels), subdir, params );

regular_bar( pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params );

end

function p_correct_per_image_cat(labels, mask, params)

combinations = [5];

for i = combinations
  each = { 'target_image_category', 'drug', 'trial_type', 'scrambled_type', 'unified_filename' };
  
  include_run_time_quantile = i == 3;
  
  if ( include_run_time_quantile )
    each{end+1} = 'run_time_quantile';
  end
  
  norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );
  
  subdir = fullfile( 'percent_correct_non_stacked', 'per_image_category' );
  subdir = sprintf( '%s-%s', subdir, norm_subdir );

  [props, prop_labels] = hwwa.percent_correct( labels, each, mask );
  
  if ( params.normalize )
    norm_each = setdiff( each, {'unified_filename', 'drug'} );
    [props, prop_labels] = hwwa.saline_normalize( props, prop_labels', norm_each );
  end
  
  fcats = { 'monkey' };
  pcats = { 'monkey', 'trial_type', 'target_image_category' };
  
  if ( include_run_time_quantile )
    pcats{end+1} = 'run_time_quantile';
    fcats{end+1} = 'run_time_quantile';
  end
  
  if ( i == 2 || i == 3 || i == 4 )
    collapsecat( prop_labels, 'monkey' );
  end

  xcats = { 'scrambled_type' };
  gcats = { 'drug' };

  pl = plotlabeled.make_common();
  
  anova_each = {};

  if ( i ~= 4 && i ~= 5 )
    anova_each = { 'trial_type' };
  end
  
  if ( include_run_time_quantile )
    anova_each{end+1} = 'run_time_quantile';
  end
  
  anova_factors = { 'target_image_category', 'scrambled_type' };

  if ( i == 4 || i == 5 )
    anova_factors{end+1} = 'trial_type';
  end
  if ( i == 5 )
    anova_factors{end+1} = 'monkey';
  end

  if ( ~params.normalize )
    anova_factors{end+1} = 'drug';
  end
  
  anova_stats( props, prop_labels, anova_each, anova_factors ...
    , rowmask(prop_labels), subdir, params );

  if ( i == 3 )
    box_plot( pl, props, prop_labels, fcats, [xcats, gcats], pcats, subdir, params );
  else
    regular_bar( pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params );
  end
end

end

function anova_stats(data, labels, each, factors, mask, subdir, params)

anova_outs = dsp3.anovan( data, labels, each, factors ...
  , 'mask', mask ...
  , 'remove_nonsignificant_comparisons', false ...
  , 'dimension', 1:numel(factors) ...
);

if ( params.do_save )
  save_p = get_analysis_path( params, subdir );
  dsp3.save_anova_outputs( anova_outs, save_p, [each, factors] );
end

end

function p_correct_variance_per_image_cat(labels, mask, params)

for i = 1
  each = { 'target_image_category', 'drug', 'trial_type', 'scrambled_type', 'unified_filename' };
  
  norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );
  
  subdir = fullfile( 'variance_percent_correct', 'per_image_category' );
  subdir = sprintf( '%s-%s', subdir, norm_subdir );

  [props, prop_labels] = trial_outcome_proportions( labels, each, mask );
  
  if ( params.normalize )
    norm_each = setdiff( each, 'unified_filename' );
    [props, prop_labels] = hwwa.saline_normalize( props, prop_labels', norm_each );
  end

  fcats = { 'monkey' };
  xcats = { 'scrambled_type' };
  gcats = { 'drug' };
  pcats = { 'monkey', 'trial_type', 'target_image_category' };

  pl = plotlabeled.make_common();
  pl.summary_func = @(x) nanvar( x, [], 1 );
  pl.add_errors = false;

  regular_bar( pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params );
end

end

function p_correct_per_image_cat_stacked(labels, mask, params)

for i = 1:2
  each = { 'target_image_category', 'drug', 'trial_type', 'scrambled_type' };
  
  if ( i == 2 )
    each{end+1} = 'monkey';
  end
  
  subdir = fullfile( 'percent_correct', 'per_image_category' );
  norm_subdir = ternary( params.normalize, 'norm', 'non-norm' );
  
  subdir = sprintf( '%s-%s', subdir, norm_subdir );

  [props, prop_labels] = trial_outcome_proportions( labels, each, mask );
  
  if ( params.normalize )
    norm_each = union( each, 'error' );
    [props, prop_labels] = hwwa.saline_normalize( props, prop_labels', norm_each );
  end

  fcats = { 'monkey' };
  xcats = { 'drug', 'scrambled_type' };
  gcats = { 'error' };
  pcats = { 'monkey', 'trial_type', 'target_image_category' };

  pl = plotlabeled.make_common();
  pl.x_order = { 'saline', 'scrambled', 'saline', 'not-scrambled' };

  if ( params.normalize )
    regular_bar( pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params );
  else
    stacked_bar( pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params );
  end
end

end

function multi_plot(func, pl, rt, prop_labels, fcats, spec, subdir, params)

assert_ispair( rt, prop_labels );

fig_I = findall_or_one( prop_labels, fcats );

store_labels = cell( size(fig_I) );
figs = cell( size(fig_I) );
store_axs = cell( size(fig_I) );

save_spec = cshorzcat( fcats, spec{:} );

for i = 1:numel(fig_I)
  figs{i} = figure(i);
  
  if ( isempty(pl) )
    pl = plotlabeled.make_common();
  end
  
  pl.x_tick_rotation = 30;
  
  plt_props = rt(fig_I{i});
  plt_labels = prune( prop_labels(fig_I{i}) );

  axs = func( pl, plt_props, plt_labels, spec{:} );
  
  if ( ~isempty(params.y_lims) )
    shared_utils.plot.set_ylims( axs, params.y_lims );
  end

  store_labels{i} = plt_labels;
  store_axs{i} = axs(:);
end

store_axs = vertcat( store_axs{:} );

if ( params.match_figure_limits )
  shared_utils.plot.match_ylims( store_axs );
end

if ( params.do_save )
  for i = 1:numel(fig_I)
    save_p = get_plot_path( params, subdir );
    shared_utils.plot.fullscreen( figs{i} );
    dsp3.req_savefig( figs{i}, save_p, store_labels{i}, save_spec );
  end
end

end

function box_plot(pl, data, labels, fcats, gcats, pcats, subdir, params)

multi_plot( @boxplot, pl, data, labels, fcats, {gcats, pcats}, subdir, params );

end

function regular_bar(pl, rt, prop_labels, fcats, xcats, gcats, pcats, subdir, params)

multi_plot( @bar, pl, rt, prop_labels, fcats, {xcats, gcats, pcats}, subdir, params );

end

function stacked_bar(pl, props, prop_labels, fcats, xcats, gcats, pcats, subdir, params)

assert_ispair( props, prop_labels );

fig_I = findall_or_one( prop_labels, fcats );

for i = 1:numel(fig_I)
  if ( isempty(pl) )
    pl = plotlabeled.make_common();
  end
  
  pl.x_tick_rotation = 30;
  
  plt_props = props(fig_I{i});
  plt_labels = prune( prop_labels(fig_I{i}) );

  axs = pl.stackedbar( plt_props, plt_labels, xcats, gcats, pcats );

  if ( params.do_save )
    save_p = get_plot_path( params, subdir );
    shared_utils.plot.fullscreen( gcf );
    dsp3.req_savefig( gcf, save_p, plt_labels, cshorzcat(fcats, gcats, pcats) );
  end
end

end

function [props, prop_labels] = trial_outcome_proportions(labels, each, mask)

each_I = findall( labels, each, mask );
error_types = setdiff( combs(labels, 'error', mask), 'no_fixation' );

prop_labels = fcat();
props = [];

for i = 1:numel(each_I)
  error_inds = findnot( labels, 'no_fixation', each_I{i} );
  denom = numel( error_inds );
  
  for j = 1:numel(error_types)
    num_of_error = double( count(labels, error_types{j}, each_I{i}) );
    props(end+1, 1) = num_of_error / denom;
    
    append1( prop_labels, labels, each_I{i} );
    setcat( prop_labels, 'error', error_types{j}, rows(prop_labels) );
  end
end

end

function [rt, labels] = ensure_correct_subset(rt, labels)

assert_ispair( rt, labels );

days_5htp = hwwa.get_image_5htp_days();
days_sal = hwwa.get_image_saline_days();

use_days = unique( horzcat(days_5htp, days_sal) );

rt = indexpair( rt, labels, findor(labels, use_days) );
prune( labels );

end

function save_figs(save_p, figs, fig_I, labels, filename_cats)

for i = 1:numel(figs)
  f = figs(i);
  ind = fig_I{i};

  shared_utils.plot.fullscreen( f );

  dsp3.req_savefig( f, save_p, prune(labels(ind)), filename_cats );
end

end

function p = get_data_path(params, data_type, varargin)

p = fullfile( hwwa.dataroot(params.config), data_type, 'approach_avoid' ...
  , 'behavior', hwwa.datedir, params.base_subdir, varargin{:} );

end

function p = get_analysis_path(params, varargin)

p = get_data_path( params, 'analyses', varargin{:} );

end

function p = get_plot_path(params, varargin)

p = get_data_path( params, 'plots', varargin{:} );

end

function mask = get_base_mask(labels, mask_func, require_initiated)

mask = hwwa.get_approach_avoid_mask( labels, {}, require_initiated );
mask = mask_func( labels, mask );

end

function mask = get_rt_mask(labels, mask)

go_mask = fcat.mask( labels, mask, @find, {'correct_true', 'go_trial'} );
nogo_mask = fcat.mask( labels, mask, @find, {'correct_false', 'nogo_trial'} );

mask = union( go_mask, nogo_mask );

% mask = fcat.mask( labels, mask ...
%   , @find, {'correct_true', 'go_trial'} ...
% );

end

function merge_scrambled_type_target_image_category(labels, varargin)

scrambled_ind = find( labels, 'scrambled', varargin{:} );
setcat( labels, 'target_image_category', 'scrambled_image', scrambled_ind );

base_ind = fcat.mask( labels, varargin{:} );
setcat( labels, 'scrambled_type', makecollapsed(labels, 'scrambled_type'), base_ind );

end