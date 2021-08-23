function hwwa_rt_vs_p_correct(rt, labels, varargin)

assert_ispair( rt, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.rt_mask_func = @hwwa.default_mask_func;
defaults.pcorr_mask_func = @hwwa.default_mask_func;
defaults.each = { 'unified_filename', 'monkey', 'drug' };
defaults.plot_cats = { {'monkey', 'drug'}, {'trial_type'} };
defaults.normalize = false;
defaults.p_format = '';
defaults.r_format = '';
defaults.per_panel_corr = false;

params = hwwa.parsestruct( defaults, varargin );
rt_mask = hwwa.get_approach_avoid_base_mask( labels, params.rt_mask_func );
pcorr_mask = hwwa.get_approach_avoid_base_mask( labels, params.pcorr_mask_func );

[rt_labs, rt_I] = keepeach( labels', params.each, rt_mask );
mean_rt = bfw.row_nanmean( rt, rt_I );
[pcorr, pcorr_labs] = hwwa.percent_correct( labels', params.each, pcorr_mask );

if ( params.normalize )
  [mean_rt, rt_labs] = hwwa.saline_normalize( mean_rt, rt_labs', norm_each(params.each) );
  [pcorr, pcorr_labs] = hwwa.saline_normalize( pcorr, pcorr_labs', norm_each(params.each) );
end

% Match rt labels to pcorr labels.
[match_I, match_C] = findall( pcorr_labs, params.each );
matched_rt = nan( size(pcorr) );

for i = 1:numel(match_I)
  rt_ind = find( rt_labs, match_C(:, i) );
  assert( numel(rt_ind) == numel(match_I{i}), 'Failed to match RT and p-corr labels.' );
  matched_rt(match_I{i}) = mean_rt(rt_ind);
end

%%

pl = plotlabeled.make_common();

gcats = params.plot_cats{1};
pcats = params.plot_cats{2};
[axs, ids] = pl.scatter( mean_rt, pcorr, pcorr_labs, gcats, pcats );

if ( params.per_panel_corr )
  ax_set = arrayfun( @(x) {x}, axs );
  id_set = arrayfun( @(x) {x}, ids );
else
  ax_set = { axs };
  id_set = { ids };
end

for i = 1:numel(ax_set)  
  if ( params.normalize )
    hwwa.add_normalized_corr( ax_set{i}(1), id_set{i}, mean_rt, pcorr ...
      , 'p_format', params.p_format, 'r_format', params.r_format );
  else
    hwwa.add_5htp_saline_all_corr( ax_set{i}(1), id_set{i}, mean_rt, pcorr );
  end
end

shared_utils.plot.xlabel( axs(1), 'RT' );
shared_utils.plot.ylabel( axs(1), '% Correct' );

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, pcorr_labs, [gcats, pcats], 'rt_vs_pcorr', params );
end

end

function each = norm_each(each)

each = setdiff( each, {'unified_filename'} );

end