function hwwa_time_between_initiated_vs_p_correct(t_bw, t_bw_labels, pcorr_labels, varargin)

assert_ispair( t_bw, t_bw_labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.pcorr_mask_func = @hwwa.default_mask_func;
defaults.t_bw_mask_func = @hwwa.default_mask_func;
defaults.each = { 'unified_filename', 'monkey', 'drug' };
defaults.plot_cats = { {'monkey'}, {'drug'} };
defaults.normalize = false;
defaults.marker_size = 5;
defaults.percent_change = false;
defaults.per_panel_corr = false;

params = hwwa.parsestruct( defaults, varargin );
pcorr_mask = hwwa.get_approach_avoid_base_mask( pcorr_labels, params.pcorr_mask_func );
t_bw_mask = hwwa.get_approach_avoid_base_mask( t_bw_labels, params.t_bw_mask_func );

[t_bw_labels, t_bw_I] = keepeach( t_bw_labels', params.each, t_bw_mask );
t_bw = bfw.row_nanmean( t_bw, t_bw_I );

[pcorr, pcorr_labs] = hwwa.percent_correct( pcorr_labels', params.each, pcorr_mask );

if ( params.normalize )
  [t_bw, t_bw_labels] = hwwa.saline_normalize( t_bw, t_bw_labels', norm_each(params.each) );
  [pcorr, pcorr_labs] = hwwa.saline_normalize( pcorr, pcorr_labs', norm_each(params.each) ); 
end

%%

matched_t_bw = nan( size(pcorr) );
[each_I, each_C] = findall( pcorr_labs, params.each );
for i = 1:numel(each_I)
  match_ind = find( t_bw_labels, each_C(:, i) );
  assert( numel(match_ind) == numel(each_I{i}) );
  matched_t_bw(each_I{i}) = t_bw(match_ind);
end
t_bw = matched_t_bw;

if ( params.percent_change )
  t_bw = (t_bw - 1) * 1e2;
  pcorr = (pcorr - 1) * 1e2;
end

%%

pl = plotlabeled.make_common();
pl.marker_size = params.marker_size;

gcats = params.plot_cats{1};
pcats = params.plot_cats{2};
[axs, ids] = pl.scatter( t_bw, pcorr, pcorr_labs, gcats, pcats );

if ( params.per_panel_corr )
  ax_set = arrayfun( @(x) {x}, axs );
  id_set = arrayfun( @(x) {x}, ids );
else
  ax_set = { axs };
  id_set = { ids };
end

for i = 1:numel(ax_set)
  if ( params.normalize )
    hwwa.add_normalized_corr( ax_set{i}(1), id_set{i}, t_bw, pcorr );
  else
    hwwa.add_5htp_saline_all_corr( ax_set{i}(1), id_set{i}, t_bw, pcorr );
  end
end

shared_utils.plot.xlabel( axs(1), 'Time between initiated' );
shared_utils.plot.ylabel( axs(1), '% Correct' );

try
  [ps, p_labs] = hwwa.permute_slope_differences( t_bw, pcorr, pcorr_labs', 1e3, gcats, pcats );
catch err
  warning( err.message );
end

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, pcorr_labs, [gcats, pcats], 'time_bw_vs_pcorr', params );
end

end

function each = norm_each(each)

each = setdiff( each, {'unified_filename'} );

end