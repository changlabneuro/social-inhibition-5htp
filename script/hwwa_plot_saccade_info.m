function hwwa_plot_saccade_info(saccade_outs, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.pre_plot_func = @(varargin) deal(varargin{1:2});
defaults.scatter_xlims = [];
defaults.scatter_ylims = [];
defaults.scatter_rois = {};
defaults.scatter_points_to_right = true;
defaults.funcs = all_function_names();
defaults.normalize = false;
defaults.per_monkey = true;
defaults.per_image_category = true;
defaults.prefer_boxplot = false;
defaults.per_correct_type = false;
defaults.mask_func = @(labels, mask) mask;
params = hwwa.parsestruct( defaults, varargin );

funcs = cellstr( params.funcs );

for i = 1:numel(funcs)
  feval( funcs{i}, saccade_outs, params );
end

% in_bounds_roi_proportions( saccade_outs, params );

% saccade_landing_point_stddev( saccade_outs, params );
% saccade_landing_point_scatter( saccade_outs, params );
% saccade_angle_variance_bars( saccade_outs, params );
% velocity_bars( saccade_outs, params );
% polar_histograms( saccade_outs, params );

end

function names = all_function_names()

names = { 'in_bounds_roi_proportions', 'saccade_landing_point_stddev' ...
  , 'saccade_landing_point_scatter', 'saccade_angle_variance_bars' ...
  , 'velocity_bars', 'polar_histograms' ...
};

end

function [data, labels] = saline_normalize(varargin)

data = varargin{1};
labels = varargin{2};
spec = varargin{3};

norm_spec = setdiff( spec, {'drug', 'day', 'unified_filename'} );
[data, labels] = hwwa.saline_normalize( data, labels, norm_spec );

end

function [fcats, gcats, pcats] = distribute_bar_categories(fcats, xcats, gcats, pcats)

pcats = csunion( pcats, gcats );
gcats = xcats;

end

function in_bounds_roi_proportions(saccade_outs, params)

pltlabs = saccade_outs.in_bounds_roi_labels';

% mask = find( pltlabs, {'correct_true', 'go_trial'}, get_base_mask(pltlabs) );
% mask = find( pltlabs, {'correct_false', 'nogo_trial'}, get_base_mask(pltlabs) );
mask = hwwa.find_correct_go_incorrect_nogo( pltlabs, get_base_mask(pltlabs) );

mask = params.mask_func( pltlabs, mask );

spec = { 'monkey', 'scrambled_type', 'unified_filename', 'roi', 'trial_type', 'drug' };

if ( params.per_image_category )
  spec{end+1} = 'target_image_category';
end

[pltlabs, prop_I] = keepeach( pltlabs, spec, mask );
pltdat = nan( size(prop_I) );

for i = 1:numel(prop_I)
  pltdat(i) = pnz( saccade_outs.in_bounds_roi(prop_I{i}) ) * 100;
end

fcats = { 'trial_type' };
pcats = { 'drug' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats{end+1} = 'monkey';
end

multiple_figures = true;
subdirs = {'in_bounds_roi_proportions'};

xcats = intersect( {'target_image_category', 'scrambled_type'}, spec );
gcats = {'trial_type', 'roi'};

if ( params.per_correct_type )
  gcats{end+1} = 'correct';
end

if ( params.prefer_boxplot )
  [fcats, gcats, pcats] = distribute_bar_categories( fcats, xcats, gcats, pcats );
  plot_cats = { fcats, gcats, pcats };
  plot_func = @boxplot;
else
  plot_cats = { fcats, xcats, gcats, pcats };
  plot_func = @bar;
end

plot_outs = dsp3.multi_plot( plot_func, pltdat, pltlabs ...
  , plot_cats{:} ...
  , 'match_limits', true ...
  , 'multiple_figures', multiple_figures ...
  , 'consolidate_outputs', true ...
  , 'pre_plot_func', require_function_handle( params.pre_plot_func ) ...
  , 'configure_pl_func', @(pl) set_property(pl, 'prefer_multiple_groups', true) ...
  , 'post_plot_func', @(varargin) post_plot(varargin, multiple_figures, subdirs, params) ...
);

if ( multiple_figures )
  prefix = func2str( plot_func );
  maybe_save_multiple( plot_outs, subdirs, params, prefix );
end

if ( params.normalize )
  norm_each = setdiff( spec, {'unified_filename'} );
  [pltdat, pltlabs] = hwwa.saline_normalize( pltdat, pltlabs, norm_each );
end

anova_each = csunion( fcats, {'trial_type'} );
anova_factors = intersect( {'target_image_category', 'scrambled_type' ...
  , 'roi', 'drug'}, spec );

if ( params.normalize )
  anova_factors = setdiff( anova_factors, 'drug' );
end

in_bounds_anova_stats( pltdat, pltlabs, anova_each, anova_factors ...
  , rowmask(pltdat), subdirs, params );

end

function in_bounds_anova_stats(data, labels, each, factors, mask, subdirs, params)

if ( ischar(factors) || numel(factors) == 1 )
  anova_outs = dsp3.anova1( data, labels', each, factors ...
    , 'mask', mask ...
    , 'remove_nonsignificant_comparisons', false ...
  );
else
  anova_outs = dsp3.anovan( data, labels', each, factors ...
    , 'mask', mask ...
    , 'remove_nonsignificant_comparisons', false ...
    , 'dimension', 1:numel(factors) ...
  );
end

if ( params.do_save )
  save_p = get_save_p( params, 'analyses', subdirs{:} );
  dsp3.save_anova_outputs( anova_outs, save_p, csunion(each, factors) );
end

end

function [dispersion, out_labels] = calculate_dispersion(points, labels, each, mask)

[out_labels, I] = keepeach( labels', each, mask );
dispersion = nan( numel(I), 1 );

for i = 1:numel(I)
  subset_pts = points(I{i}, :);
  centroid = nanmean( subset_pts, 1 );
  
  dist = sqrt( (subset_pts(:, 2) - centroid(2)).^2 + (subset_pts(:, 1) - centroid(1)).^2 );
  dispersion(i) = nanstd( dist, [] );
end

end

function saccade_landing_point_stddev(saccade_outs, params)

pltdat = saccade_outs.saccade_stop_points_to_right;
pltlabs = saccade_outs.labels';

mask = find( pltlabs, {'correct_true', 'go_trial'}, get_base_mask(pltlabs) );

disp_each = { 'monkey', 'target_image_category', 'scrambled_type', 'unified_filename' };
[pltdat, pltlabs] = calculate_dispersion( pltdat, pltlabs', disp_each, mask );

multiple_figures = true;
subdirs = {'saccade_stop_point_stddev'};

fcats = {};
xcats = { 'target_image_category', 'scrambled_type' };
gcats = { 'trial_type', 'correct' };
pcats = { 'drug' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats{end+1} = 'monkey';
end

if ( params.prefer_boxplot )
  [fcats, gcats, pcats] = distribute_bar_categories( fcats, xcats, gcats, pcats );
  plot_cats = { fcats, gcats, pcats };
  plot_func = @boxplot;
else
  plot_cats = { fcats, xcats, gcats, pcats };
  plot_func = @bar;
end

plot_outs = dsp3.multi_plot( plot_func, pltdat, pltlabs ...
  , plot_cats{:} ...
  , 'match_limits', true ...
  , 'multiple_figures', multiple_figures ...
  , 'consolidate_outputs', true ...
  , 'pre_plot_func', require_function_handle( params.pre_plot_func ) ...
  , 'post_plot_func', @(varargin) post_plot(varargin, multiple_figures, subdirs, params) ...
);

if ( multiple_figures )
  prefix = func2str( plot_func );
  maybe_save_multiple( plot_outs, subdirs, params, prefix );
end

end

function func = require_function_handle(func)

if ( ~isa(func, 'function_handle') )
  validateattributes( func, {'char'}, {}, mfilename, 'function string' );
  func = str2func( func );
end

end

function [dist, labels] = centroid_relative_distance(pts, labels, spec, mask)

centroid = nanmean( pts, 1 );
dist = sqrt( (pts(:, 1) - centroid(1)).^2 + (pts(:, 2) - centroid(2)).^2 );

end

function saccade_landing_point_scatter(saccade_outs, params)

if ( params.scatter_points_to_right )
  pltdat = saccade_outs.saccade_stop_points_to_right;
else
  pltdat = saccade_outs.saccade_stop_points;
end

pltlabs = saccade_outs.labels';

mask = get_base_mask( pltlabs );
correct_mask = find( pltlabs, {'correct_true', 'go_trial'}, mask );
mask = correct_mask;

subdirs = {'saccade_stop_point_scatter'};

fcats = { 'monkey' };
gcats = { 'drug' };
pcats = { 'monkey', 'target_image_category', 'scrambled_type', 'trial_type', 'correct' };

fig_I = findall_or_one( pltlabs, fcats, mask );

spec = [ fcats, gcats, pcats ];

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.x_lims = params.scatter_xlims;
  pl.y_lims = params.scatter_ylims;
  
  subset_data = pltdat(fig_I{i}, :);
  subset_labels = prune( pltlabs(fig_I{i}) );
  
  [subset_data, subset_labels] = params.pre_plot_func( subset_data, subset_labels );
  axs = pl.scatter( subset_data(:, 1), subset_data(:, 2), subset_labels, gcats, pcats );
  
  if ( ~isempty(params.scatter_rois) )
    shared_utils.plot.prevent_legend_autoupdate( gcf );
    for j = 1:numel(axs)
      for k = 1:numel(params.scatter_rois)
        bfw.plot_rect_as_lines( axs(j), params.scatter_rois{k} );
      end
    end
  end
  
  inputs = { gcf, [], subset_labels, spec };
  
  post_plot( inputs, false, subdirs, params );
end

end

function saccade_angle_variance_bars(saccade_outs, params)

pltdat = saccade_outs.saccade_angles_to_right;
pltlabs = saccade_outs.labels';

mask = get_base_mask( pltlabs );
correct_mask = find( pltlabs, {'correct_true', 'go_trial'}, mask );
incorrect_mask = find( pltlabs, {'correct_false', 'nogo_trial'}, mask );
mask = union( correct_mask, incorrect_mask );
% mask = correct_mask;

% prune( collapsecat(pltlabs, 'monkey') );

multiple_figures = true;
subdirs = {'saccade_variance_in_angle'};

plot_outs = dsp3.multi_plot( @bar, pltdat, pltlabs ...
  , {'monkey'} ...  % figures
  , {'target_image_category', 'scrambled_type'} ... % x
  , {'trial_type', 'correct'} ... % g
  , {'monkey', 'drug'} ...  % p
  , 'match_limits', true ...
  , 'multiple_figures', multiple_figures ...
  , 'mask', mask ...
  , 'consolidate_outputs', true ...
  , 'configure_pl_func', @(pl) set_property(pl, 'summary_func', @(x) nanvar(x, [], 1)) ...
  , 'pre_plot_func', params.pre_plot_func ...
  , 'post_plot_func', @(varargin) post_plot(varargin, multiple_figures, subdirs, params) ...
);

if ( multiple_figures )
  maybe_save_multiple( plot_outs, subdirs, params );
end


end

function velocity_bars(saccade_outs, params)

pltdat = saccade_outs.saccade_velocities;
pltlabs = saccade_outs.labels';

mask = get_base_mask( pltlabs );
correct_mask = find( pltlabs, {'correct_true', 'go_trial'}, mask );
incorrect_mask = find( pltlabs, {'correct_false', 'nogo_trial'}, mask );
% mask = union( correct_mask, incorrect_mask );
mask = correct_mask;

prune( collapsecat(pltlabs, 'monkey') );

multiple_figures = true;
subdirs = {'saccade_velocity'};

plot_outs = dsp3.multi_plot( @bar, pltdat, pltlabs ...
  , {'monkey'} ...  % figures
  , {'target_image_category', 'scrambled_type'} ... % x
  , {'trial_type', 'correct'} ... % g
  , {'monkey', 'drug'} ...  % p
  , 'match_limits', true ...
  , 'multiple_figures', multiple_figures ...
  , 'mask', mask ...
  , 'consolidate_outputs', true ...
  , 'pre_plot_func', params.pre_plot_func ...
  , 'post_plot_func', @(varargin) post_plot(varargin, multiple_figures, subdirs, params) ...
);

if ( multiple_figures )
  maybe_save_multiple( plot_outs, subdirs, params );
end

end

function polar_histograms(saccade_outs, params)

pltdat = saccade_outs.saccade_angles_to_right;
pltlabs = saccade_outs.labels';

mask = get_base_mask( pltlabs );

multiple_figures = true;
subdirs = {'saccade_direction'};

plot_outs = dsp3.multi_plot( @polarhistogram, pltdat, pltlabs ...
  , {'monkey', 'drug'} ...  % figures
  , {'monkey', 'trial_type', 'correct', 'drug'} ... % panels
  , 'plot_func_inputs', {1e2, 'norm', 'prob'} ...
  , 'match_limits', true ...
  , 'multiple_figures', multiple_figures ...
  , 'mask', mask ...
  , 'consolidate_outputs', true ...
  , 'pre_plot_func', params.pre_plot_func ...
  , 'post_plot_func', @(varargin) post_plot(varargin, multiple_figures, subdirs, params) ...
);

if ( multiple_figures )
  maybe_save_multiple( plot_outs, subdirs, params );
end

end

function post_plot(post_plot_inputs, multiple_figures, subdirs, params)

if ( multiple_figures )
  return
end

[fig, labels, spec] = dsp3.util.post_plot.fig_labels_specificity( post_plot_inputs{:} );

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = get_save_p( params, 'plots', subdirs{:} );
  dsp3.req_savefig( fig, save_p, labels, spec );
end

end

function maybe_save_multiple(plot_outs, subdirs, params, prefix)

if ( nargin < 4 )
  prefix = '';
end

if ( params.save )
  save_p = get_save_p( params, 'plots', subdirs{:}, params.base_subdir );
  shared_utils.plot.fullscreen( plot_outs.figs );
  dsp3.req_savefigs( plot_outs.figs, save_p, plot_outs.labels, plot_outs.full_specificity, prefix );
end

end

function save_p = get_save_p(params, kind, varargin)

save_p = fullfile( hwwa.dataroot(params.config), kind ...
  , 'approach_avoid', 'behavior', dsp3.datedir, varargin{:} ); 

end

function mask = get_base_mask(labels)

mask = fcat.mask( labels, hwwa.get_approach_avoid_mask(labels) ...
  , @find, 'initiated_true' ...
);

end