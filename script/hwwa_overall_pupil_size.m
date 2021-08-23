function outs = hwwa_overall_pupil_size(ps, labels, varargin)

assert_ispair( ps, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = false;
defaults.per_image_category = false;
defaults.per_monkey = false;
defaults.per_trial_type = true;
defaults.per_correct = true;
defaults.normalize = false;
defaults.add_points = false;
defaults.points_are = {};
defaults.run_level_average = false;
defaults.signrank_inputs = {};

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

if ( params.normalize )
  norm_spec = union( get_base_specificity(params), {'monkey'} );
  [ps, labels] = hwwa.saline_normalize( ps, labels', norm_spec, mask );
  mask = rowmask( ps );
end

if ( params.run_level_average )
  spec = union( get_base_specificity(params), {'unified_filename'} );
  [~, each_I] = keepeach( labels, spec, mask );
  ps = bfw.row_nanmean( ps, each_I );
  mask = rowmask( labels );
end

outs = struct();
outs.pupil_size = ps;
outs.labels = labels';

% hists( ps, labels', mask, 'hists', params );
box_plots( ps, labels', mask, 'boxplots', params );

end

function hists(ps, labels, mask, subdir, params)

fcats = {};
pcats = { 'drug', 'correct' };

if ( params.per_trial_type )
  fcats{end+1} = 'trial_type';
end
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

for i = 1:numel(fig_I)
  plt_ps = ps(fig_I{i});
  plt_labels = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.hist_add_summary_line = true;
  axs = pl.hist( plt_ps, plt_labels', pcats );
  
  maybe_save_fig( gcf, plt_labels, [pcats, fcats], subdir, params );
end

end

function [gcats, pcats, fcats] = get_plot_specificity(gcats, pcats, fcats, params)

if ( params.per_trial_type )
  fcats{end+1} = 'trial_type';
end
if ( params.per_image_category )
  pcats{end+1} = 'target_image_category';
end
if ( params.per_scrambled_type )
  pcats{end+1} = 'scrambled_type';
end
if ( params.per_monkey )
  fcats{end+1} = 'monkey';
end
if ( ~params.per_correct )
  gcats = setdiff( gcats, 'correct' );
end

if ( isempty(gcats) )
  gcats = { 'drug' };
  pcats = setdiff( pcats, {'drug'} );
end

end

function box_plots(ps, labels, mask, subdir, params)

fcats = {};
gcats = { 'correct' };
pcats = { 'drug' };
[gcats, pcats, fcats] = get_plot_specificity( gcats, pcats, fcats, params );

pcats = [ pcats, fcats ];

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  plt_ps = ps(fig_I{i});
  plt_labels = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.add_points = params.add_points;
  pl.points_are = params.points_are;
  axs = pl.boxplot( plt_ps, plt_labels', gcats, pcats );
  
  maybe_save_fig( gcf, plt_labels, [gcats, pcats, fcats], subdir, params );
  
  axs = pl.bar( plt_ps, plt_labels', {}, gcats, pcats );
  p = params;
  p.prefix = 'bar--';
  maybe_save_fig( gcf, plt_labels, [gcats, pcats, fcats], subdir, p );
end

try
  rs_tests( ps, labels, mask, [fcats, pcats, gcats], 'saline', '5-htp', subdir, params );
catch err
  warning( err.message );
end

try
  ttests( ps, labels, mask, [fcats, pcats, gcats], 'saline', '5-htp', subdir, params );
catch err
  warning( err.message );
end

try
  signrank1( ps, labels, mask, [fcats, pcats, gcats], subdir, params );
catch err
  warning( err.message );
end

end

function signrank1(data, labels, mask, each, subdir, params)

sr_outs = dsp3.signrank1( data, labels, each ...
  , 'mask', mask ...
  , 'signrank_inputs', params.signrank_inputs ...
);

if ( params.save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'pupil size', subdir );
  dsp3.save_signrank1_outputs( sr_outs, save_p, each );
end

end

function ttests(data, labels, mask, each, a, b, subdir, params)

targ_cat = union( whichcat(labels, a), whichcat(labels, b) );
each = setdiff( each, targ_cat );

ttest_outs = dsp3.ttest2( data, labels, each, a, b ...
  , 'mask', mask ...
);

if ( params.save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'pupil size', subdir );
  dsp3.save_ttest2_outputs( ttest_outs, save_p, [each, targ_cat] );
end

end

function rs_tests(data, labels, mask, each, a, b, subdir, params)

targ_cat = union( whichcat(labels, a), whichcat(labels, b) );
each = setdiff( each, targ_cat );

rs_outs = dsp3.ranksum( data, labels, each, a, b ...
  , 'mask', mask ...
);

if ( params.save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'pupil size', subdir );
  dsp3.save_ranksum_outputs( rs_outs, save_p, [each, targ_cat] );
end

end

function maybe_save_fig(fig, labels, each, subdir, params)

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pupil size', subdir );
  shared_utils.plot.fullscreen( fig );
  dsp3.req_savefig( fig, save_p, labels, each, params.prefix );
end

end

function mask = get_base_mask(labels, mask_func)

mask = mask_func( labels, hwwa.get_approach_avoid_mask(labels) );

end

function spec = get_base_specificity(params)

spec = {};

if ( params.per_trial_type )
  spec{end+1} = 'trial_type';
end
if ( params.per_image_category )
  spec{end+1} = 'target_image_category';
end
if ( params.per_scrambled_type )
  spec{end+1} = 'scrambled_type';
end
if ( params.per_monkey )
  spec{end+1} = 'monkey';
end
if ( params.per_correct )
  spec{end+1} = 'correct';
end
end