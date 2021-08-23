function hwwa_basic_hist(data, labels, plot_type, varargin)

assert_ispair( data, labels );
validateattributes( data, {'double'}, {'2d', 'column'}, mfilename, 'data' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.y_lims = [];
defaults.y_label = '';
defaults.pcats = {};
defaults.collapse_each = {};
defaults.collapse_op = @(x) nanmean(x, 1);
defaults.mask_func = @hwwa.default_mask_func;
defaults.marker_size = 1;
defaults.run_anova = false;
defaults.anova_each = {};
defaults.anova_factors = {};
defaults.hist_size = [];
defaults.add_summary_line = true;
defaults.summary_func = @median;
defaults.run_ranksum = false;
defaults.ranksum_each = {};
defaults.ranksum_a = '';
defaults.ranksum_b = '';
params = hwwa.parsestruct( defaults, varargin );

base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

pl = plotlabeled.make_common();
pl.y_lims = params.y_lims;
pl.hist_add_summary_line = params.add_summary_line;
pl.summary_func = params.summary_func;

[data, labels] = hwwa.maybe_apply_rowop( ...
  data, labels, params.collapse_each, params.collapse_op, base_mask );

pcats = params.pcats;
axs = pl.hist( data, labels, pcats, params.hist_size );

if ( ~isempty(params.y_label) && ~isempty(axs) )
  ylabel( axs(1), params.y_label );
end

anova_outs = [];
if ( params.run_anova )
  anova_outs = dsp3.anovan( data, labels, params.anova_each, params.anova_factors );
end

rs_outs = [];
if ( params.run_ranksum )
  rs_outs = dsp3.ranksum( data, labels, params.ranksum_each, params.ranksum_a, params.ranksum_b );
end

if ( params.do_save )
  shared_utils.plot.fullscreen( gcf );
  plot_spec = unique( cellstr(pcats) );
  save_p = hwwa.approach_avoid_save_fig( gcf, axs, labels, plot_spec, plot_type, params );
  stat_p = fullfile( save_p, 'stats' );
  
  if ( ~isempty(anova_outs) )
    dsp3.save_anova_outputs( anova_outs, stat_p, plot_spec );
  end
  if ( ~isempty(rs_outs) )
    dsp3.save_ranksum_outputs( rs_outs, stat_p, plot_spec );
  end
end

end