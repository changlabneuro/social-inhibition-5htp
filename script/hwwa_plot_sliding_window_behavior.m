function hwwa_plot_sliding_window_behavior(data, time, labels, subdir, varargin)

assert_ispair( data, labels );
assert_ispair( time, labels );
validateattributes( data, {'double'}, {'vector'}, mfilename, 'data' );
validateattributes( time, {'double'}, {'vector'}, mfilename, 'time' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @(labels, mask) mask;
defaults.saline_normalize = true;
defaults.over_time_window_size = 60;
defaults.per_monkey = true;
defaults.line_xlims = [];

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

variance_per_image_category( data, labels', mask, subdir, params );
lines_over_time( data, time, labels', mask, subdir, params );
% correlation_with_time( data, time, labels', mask, params );

end

function lines_over_time(data, time, labels, mask, subdir, params)

window_size = params.over_time_window_size;
each = { 'unified_filename', 'target_image_category', 'monkey', 'scrambled_type', 'trial_type' };
[samples, edges, over_time_labels] = hwwa.time_bin( data, time, labels', each, window_size, mask );

med_samples = cellfun( @nanmean, samples );

if ( params.saline_normalize )
  norm_each = setdiff( each, 'unified_filename' );
  [med_samples, over_time_labels] = hwwa.saline_normalize( med_samples, over_time_labels, norm_each );
end

pl = plotlabeled.make_common();
pl.x = edges(1:end-1);
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.25);

fcats = {};

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
else
  collapsecat( over_time_labels, 'monkey' );
end

gcats = { 'drug' };
pcats = { 'trial_type', 'monkey' };

[figs, axs, fig_inds] = pl.figures( @lines, med_samples, over_time_labels, fcats, gcats, pcats );

shared_utils.plot.match_xlims( axs );
shared_utils.plot.match_ylims( axs );

if ( ~isempty(params.line_xlims) )
  shared_utils.plot.set_xlims( axs, params.line_xlims );
end

if ( params.do_save )
  plot_p = hwwa.approach_avoid_data_path( params, 'plots', subdir );  
  
  for i = 1:numel(figs)
    shared_utils.plot.fullscreen( figs(i) );
    dsp3.req_savefig( figs(i), plot_p, over_time_labels(fig_inds{i}), [gcats, pcats, fcats] );
  end
end

end

function correlation_with_time(data, time, labels, mask, params)

subset_data = data(mask);
subset_time = time(mask);
subset_labels = prune( labels(mask) );

gcats = { 'drug' };
pcats = { 'target_image_category', 'scrambled_type', 'trial_type' };

pl = plotlabeled.make_common();

[axs, ids] = pl.scatter( subset_time, subset_data, subset_labels, gcats, pcats );
plotlabeled.scatter_addcorr( ids, subset_time, subset_data );

end

function variance_per_image_category(data, labels, mask, subdir, params)

pl = plotlabeled.make_common();
pl.x_tick_rotation = 30;

is_sal_normalized = params.saline_normalize;

norm_postfix = ternary( is_sal_normalized, 'norm', 'non-norm' );
subdir = sprintf( '%s-variance', subdir );
subdir = sprintf( '%s-%s', subdir, norm_postfix );

var_each = { 'unified_filename', 'monkey', 'scrambled_type', 'target_image_category', 'trial_type' };
var_each = setdiff( var_each, 'target_image_category' );

[var_labels, run_I] = keepeach( labels', var_each, mask );
vars = rowop( data, run_I, @(x) nanvar(x, [], 1) );

if ( is_sal_normalized )
  norm_each = setdiff( var_each, 'unified_filename' );
  [vars, var_labels] = hwwa.saline_normalize( vars, var_labels, norm_each );
end

xcats = { 'scrambled_type' };
gcats = { 'drug' };
pcats = intersect( {'target_image_category', 'trial_type'}, var_each );

% axs = pl.bar( vars, var_labels, xcats, gcats, pcats );
axs = pl.boxplot( vars, var_labels, [xcats, gcats], pcats );

if ( params.do_save )
  plot_p = hwwa.approach_avoid_data_path( params, 'plots', subdir );  
  dsp3.req_savefig( gcf, plot_p, var_labels, [gcats, pcats] );
end

anova_factors = intersect( {'target_image_category', 'scrambled_type'}, var_each );
variance_anova_stats( vars, var_labels', anova_factors, is_sal_normalized, subdir, params );

end

function variance_anova_stats(vars, var_labels, anova_factors, is_sal_normalized, subdir, params)

anova_each = { 'trial_type' };

if ( ~is_sal_normalized )
  anova_factors{end+1} = 'drug';
end

anova_outs = dsp3.anovan( vars, var_labels', anova_each, anova_factors ...
  , 'remove_nonsignificant_comparisons', false ...
  , 'dimension', 1:2 ...
);

ttest_outs = [];
post_hoc_anova_outs = [];

if ( is_sal_normalized )
  ttest_each = { 'trial_type', 'target_image_category' };
  ttest_outs = dsp3.ttest2( vars, var_labels', ttest_each ...
    , 'not-scrambled', 'scrambled' ...
    , 'mask', hwwa.find_nogo(var_labels) ...
  );
else
  post_hoc_anova_each = { 'target_image_category' };
  post_hoc_anova_factors = { 'drug', 'scrambled_type' };
  post_hoc_anova_outs = dsp3.anovan( vars, var_labels', post_hoc_anova_each, post_hoc_anova_factors ...
    , 'mask', hwwa.find_nogo(var_labels) ...
  );
end

if ( params.do_save )
  analysis_p = hwwa.approach_avoid_data_path( params, 'analyses', subdir );
  anova_spec = [ anova_factors, anova_each ];
  
  main_p = fullfile( analysis_p, 'main_effects' );
  post_hoc_p = fullfile( analysis_p, 'post_hoc' );
  
  if ( ~isempty(anova_outs) )
    dsp3.save_anova_outputs( anova_outs, main_p, anova_spec );    
  end
  
  if ( ~isempty(post_hoc_anova_outs) )
    dsp3.save_anova_outputs( post_hoc_anova_outs, post_hoc_p, anova_spec );
  end
  
  if ( ~isempty(ttest_outs) )
    dsp3.save_ttest2_outputs( ttest_outs, post_hoc_p, anova_spec );
  end
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels );
mask = mask_func( labels, mask );

end