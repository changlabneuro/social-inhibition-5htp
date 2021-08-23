function hwwa_basic_scatter(data0, data1, labels, plot_type, varargin)

assert_ispair( data0, labels );
validateattributes( data0, {'double'}, {'2d', 'column'}, mfilename, 'data' );

assert_ispair( data1, labels );
validateattributes( data1, {'double'}, {'2d', 'column'}, mfilename, 'data' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.x_summary_func = [];
defaults.x_test_func = [];
defaults.y_summary_func = [];
defaults.y_test_func = [];
defaults.x_summary_format = 'M = %0.2f';
defaults.y_summary_format = 'M = %0.2f';
defaults.x_test_format = 'P = %0.2f';
defaults.y_test_format = 'P = %0.2f';
defaults.y_lims = [];
defaults.x_label = '';
defaults.y_label = '';
defaults.gcats = {};
defaults.pcats = {};
defaults.collapse_each = {};
defaults.collapse_op = @(x) nanmean(x, 1);
defaults.mask_func = @hwwa.default_mask_func;
defaults.marker_size = 1;
defaults.permute = false;
defaults.permute_each = {};
defaults.permute_compare = {};
defaults.perm_iters = 1e2;
params = hwwa.parsestruct( defaults, varargin );

base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

pl = plotlabeled.make_common();
pl.y_lims = params.y_lims;
pl.marker_size = params.marker_size;

[data0, ~] = hwwa.maybe_apply_rowop( ...
  data0, labels', params.collapse_each, params.collapse_op, base_mask );
[data1, labels] = hwwa.maybe_apply_rowop( ...
  data1, labels', params.collapse_each, params.collapse_op, base_mask );

gcats = params.gcats;
pcats = params.pcats;
[axs, ids] = pl.scatter( data0, data1, labels, gcats, pcats );
corr_hs = plotlabeled.scatter_addcorr( ids, data0, data1 );

if ( ~isempty(params.x_label) && ~isempty(axs) )
  xlabel( axs(1), params.x_label );
end

if ( ~isempty(params.y_label) && ~isempty(axs) )
  ylabel( axs(1), params.y_label );
end

if ( ~isempty(params.x_summary_func) )
  add_summary_lines( ids, data0, params.x_summary_func, params.x_summary_format ...
    , params.x_test_func, params.x_test_format, true );
end

if ( ~isempty(params.y_summary_func) )
  add_summary_lines( ids, data1, params.y_summary_func, params.y_summary_format ...
    , params.y_test_func, params.y_test_format, false );
end

if ( params.permute )
  [ps, p_labels] = hwwa.permute_slope_differences( ...
    data0, data1, labels, params.perm_iters, params.permute_each, params.permute_compare );
  d = 10;
end

if ( params.do_save )
  shared_utils.plot.fullscreen( gcf );
  plot_spec = unique( [cellstr(gcats), cellstr(pcats)] );
  hwwa.approach_avoid_save_fig( gcf, axs, labels, plot_spec, plot_type, params );
end

end

function add_summary_lines(ids, data, summary_func, summary_format ...
  , test_func, test_format, is_x)

for i = 1:numel(ids)
  ax = ids(i).axes;
  summary = summary_func( data(ids(i).index) );
  
  if ( ~isempty(test_func) )
    p = test_func( data(ids(i).index) );
  end
  
  if ( is_x )
    hs = shared_utils.plot.add_vertical_lines( ax, summary );
  else
    hs = shared_utils.plot.add_horizontal_lines( ax, summary );
  end
  
  colors = unique( ids(i).series.CData, 'rows' );
  if ( rows(colors) == 1 )
    set( hs, 'color', colors );
  end
  
  if ( is_x )
    xp = summary;
    yp = min( get(ax, 'ylim') );
  else
    xp = min( get(ax, 'xlim') );
    yp = summary;
  end
  
  use_text = sprintf( summary_format, summary );
  if ( ~isempty(test_func) )
    test_text = sprintf( test_format, p );
    use_text = sprintf( '%s %s', use_text, test_text );
  end
  
  th = text( ax, xp, yp, use_text );
  
  if ( rows(colors) == 1 )
%     set( th, 'color', colors );
  end
end

end