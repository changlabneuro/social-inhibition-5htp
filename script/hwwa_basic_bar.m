function hwwa_basic_bar(data, labels, plot_type, varargin)

assert_ispair( data, labels );
validateattributes( data, {'double'}, {'2d', 'column'}, mfilename, 'data' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.y_lims = [];
defaults.y_label = '';
defaults.xcats = {};
defaults.gcats = {};
defaults.pcats = {};
defaults.collapse_each = {};
defaults.collapse_op = @(x) nanmean(x, 1);
defaults.mask_func = @hwwa.default_mask_func;
defaults.points_are = {};
defaults.marker_size = 1;
defaults.run_anova = false;
defaults.anova_each = {};
defaults.anova_factors = {};
defaults.stacked = false;
params = hwwa.parsestruct( defaults, varargin );

base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

pl = plotlabeled.make_common();
pl.y_lims = params.y_lims;

if ( ~isempty(params.points_are) )
  pl.add_points = true;
  pl.points_are = params.points_are;
  pl.marker_size = params.marker_size;
end

[data, labels] = hwwa.maybe_apply_rowop( ...
  data, labels, params.collapse_each, params.collapse_op, base_mask );

xcats = params.xcats;
gcats = params.gcats;
pcats = params.pcats;

if ( params.stacked )
  axs = pl.stackedbar( data, labels, xcats, gcats, pcats );
else
  axs = pl.bar( data, labels, xcats, gcats, pcats );
end

if ( ~isempty(params.y_label) && ~isempty(axs) )
  ylabel( axs(1), params.y_label );
end

anova_outs = [];
if ( params.run_anova )
  anova_outs = dsp3.anovan( data, labels, params.anova_each, params.anova_factors );
end

if ( params.do_save )
  shared_utils.plot.fullscreen( gcf );
  plot_spec = unique( [cellstr(xcats), cellstr(gcats), cellstr(pcats)] );
  save_p = hwwa.approach_avoid_save_fig( gcf, axs, labels, plot_spec, plot_type, params );
  
  if ( ~isempty(anova_outs) )
    dsp3.save_anova_outputs( anova_outs, fullfile(save_p, 'stats'), plot_spec );
  end
end

end