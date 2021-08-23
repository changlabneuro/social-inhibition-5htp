function hwwa_num_initiated_vs_p_correct(labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.each = { 'unified_filename', 'monkey', 'drug' };
defaults.plot_cats = { {'monkey'}, {'drug'} };
defaults.normalize = false;
defaults.marker_size = 5;
defaults.per_panel_corr = false;

params = hwwa.parsestruct( defaults, varargin );
mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

[init, init_labs] = hwwa.num_initiated( labels', params.each, mask );
[pcorr, pcorr_labs] = hwwa.percent_correct( labels', params.each, mask );
assert( init_labs == pcorr_labs );

if ( params.normalize )
  [init, init_labs] = hwwa.saline_normalize( init, init_labs', norm_each(params.each) );
  [pcorr, pcorr_labs] = hwwa.saline_normalize( pcorr, pcorr_labs', norm_each(params.each) );
  assert( pcorr_labs == init_labs );
end

%%

pl = plotlabeled.make_common();
pl.marker_size = params.marker_size;

gcats = params.plot_cats{1};
pcats = params.plot_cats{2};
[axs, ids] = pl.scatter( init, pcorr, init_labs, gcats, pcats );

if ( params.per_panel_corr )
  ax_sets = arrayfun( @identity, axs, 'un', 0 );
  id_sets = arrayfun( @identity, ids, 'un', 0 );
else
  ax_sets = { axs(1) };
  id_sets = { ids(1) };
end

for i = 1:numel(ax_sets)
  if ( params.normalize )
    hwwa.add_normalized_corr( ax_sets{i}, id_sets{i}, init, pcorr );
  else
    hwwa.add_5htp_saline_all_corr( ax_sets{i}, id_sets{i}, init, pcorr );
  end
end

shared_utils.plot.xlabel( axs(1), 'Num initiated' );
shared_utils.plot.ylabel( axs(1), '% Correct' );

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, init_labs, [gcats, pcats], 'num_init_vs_pcorr', params );
end

end

function each = norm_each(each)

each = setdiff( each, {'unified_filename'} );

end