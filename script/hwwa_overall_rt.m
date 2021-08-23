function outs = hwwa_overall_rt(rt, labels, varargin)

assert_ispair( rt, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = true;
defaults.per_image_category = false;
defaults.per_monkey = false;
defaults.include_nogo = false;
defaults.trial_level = false;
defaults.normalize = false;
defaults.add_points = false;
defaults.points_are = {};
defaults.y_lims = [0.1, 0.26];

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func, params.include_nogo );
mask = intersect( mask, find(~isnan(rt)) );

if ( params.trial_level )
  mean_rt = rt;
  labs = labels';
  
else
  mean_each = rt_each( params );
  [labs, I] = keepeach( labels', mean_each, mask );
  mean_rt = bfw.row_nanmean( rt, I );
  mask = rowmask( labs );
end

outs = struct();
outs.rt = mean_rt;
outs.rt_labels = labs';

if ( params.normalize )
  norm_each = setdiff( rt_each(params), {'unified_filename'} );
  [mean_rt, labs] = hwwa.saline_normalize( mean_rt, labs', norm_each, mask );
  mask = rowmask( mean_rt );
end

outs.norm_rt = mean_rt;
outs.norm_rt_labels = labs';

box_plots( mean_rt, labs', mask, '', params );

end

function box_plots(rt, labels, mask, subdir, params)

fcats = {};
gcats = { 'drug' };
pcats = { 'trial_type' };

if ( params.per_image_category )
  pcats{end+1} = 'target_image_category';
end
if ( params.per_scrambled_type )
  pcats{end+1} = 'scrambled_type';
end
if ( params.per_monkey )
  fcats{end+1} = 'monkey';
end

pcats = [ pcats, fcats ];

fig_I = findall_or_one( labels, fcats, mask );

if ( params.trial_level )
  subdir = sprintf( '%s-trial-level', subdir );
else
  subdir = sprintf( '%s-session-level', subdir );
end
if ( params.normalize )
  subdir = sprintf( '%s-norm', subdir );
else
  subdir = sprintf( '%s-non-norm', subdir );
end

for i = 1:numel(fig_I)
  plt_rt = rt(fig_I{i});
  plt_labels = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.y_lims = params.y_lims;
  pl.add_points = params.add_points;
  pl.points_are = params.points_are;
  
%   axs = pl.boxplot( plt_rt, plt_labels', gcats, pcats );
  axs = pl.bar( plt_rt, plt_labels', gcats, pcats, {} );
  
  maybe_save_fig( gcf, plt_labels, [gcats, pcats, fcats], subdir, params );
end

anova_factors = ...
  union( intersect(rt_each(params), {'scrambled_type', 'target_image_category'}), {'drug'} );

run_anovas( rt, labels', {'trial_type', 'correct'}, anova_factors ...
  , mask, subdir, params );

end

function run_anovas(rt, labels, each, factors, mask, subdir, params)

factors(isuncat(labels, factors, mask)) = [];

anova_outs = dsp3.anovan( rt, labels', each, factors ...
  , 'mask', mask ...
  , 'remove_nonsignificant_comparisons', false ...
  , 'dimension', 1:numel(factors) ...
);

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'rt', subdir );
  dsp3.save_anova_outputs( anova_outs, save_p, [factors, each] );
end

end

function maybe_save_fig(fig, labels, each, subdir, params)

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'rt', subdir );
  shared_utils.plot.fullscreen( fig );
  dsp3.req_savefig( fig, save_p, labels, each );
end

end

function each = rt_each(params)

each = { 'unified_filename', 'monkey', 'trial_type' };
if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end
if ( params.per_image_category )
  each{end+1} = 'target_image_category';
end
if ( params.per_monkey )
  each{end+1} = 'monkey';
end

end

function mask = get_base_mask(labels, mask_func, include_nogo)

mask = mask_func( labels, hwwa.get_approach_avoid_mask(labels) );
go_mask = hwwa.find_correct_go( labels, mask );

if ( include_nogo )
  nogo_mask = hwwa.find_incorrect_nogo( labels, mask );
  mask = union( go_mask, nogo_mask );
else
  mask = go_mask;
end

end