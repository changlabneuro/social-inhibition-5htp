function hwwa_plot_num_initiated_over_time(num_init, t, labels, varargin)

assert_ispair( num_init, labels );
assert( numel(t) == size(num_init, 2), 'Time does not match num initiated matrix.' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.each = { 'unified_filename', 'drug' };
defaults.plot_cats = { {'drug'}, {} };
defaults.normalize = false;
defaults.y_label = '';
defaults.y_lims = [];

params = hwwa.parsestruct( defaults, varargin );
mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

[labs, I] = keepeach( labels', params.each, mask );
mean_init = bfw.row_nanmean( num_init, I );

if ( params.normalize )
  [mean_init, labs] = hwwa.saline_normalize( mean_init, labs', norm_each(params.each) );
end

%%

pl = plotlabeled.make_common();
pl.x = t;

gcats = params.plot_cats{1};
pcats = params.plot_cats{2};

axs = pl.lines( mean_init, labs, gcats, pcats );
shared_utils.plot.ylabel( axs(1), params.y_label );

if ( ~isempty(params.y_lims) )
  shared_utils.plot.set_ylims( axs, params.y_lims );
end

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, labs, [gcats, pcats] ...
    , 'num_initiated_over_time', params );
end

end

function each = norm_each(each)

each = setdiff( each, {'unified_filename'} );

end