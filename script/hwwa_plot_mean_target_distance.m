function hwwa_plot_mean_target_distance(dist, labels, varargin)

assert_ispair( dist, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.mean_each = { 'unified_filename', 'drug', 'scrambled_type', 'target_image_category', 'trial_type' };
defaults.bar_cats = { {'drug'}, {'scrambled_type'}, {'target_image_category', 'trial_type', 'correct'} };

params = hwwa.parsestruct( defaults, varargin );
mask = get_base_mask( labels, params.mask_func );

[mean_labs, mean_I] = keepeach( labels', params.mean_each, mask );
mean_dist = bfw.row_nanmean( dist, mean_I );
run_bar( mean_dist, mean_labs, rowmask(mean_dist), params.bar_cats, 'bar', params );

end

function run_bar(data, labels, mask, cats, subdir, params)

%%

xcats = cats{1};
gcats = cats{2};
pcats = cats{3};

pl = plotlabeled.make_common();

d = data(mask);
l = prune( labels(mask) );

axs = pl.bar( d, l, xcats, gcats, pcats );

if ( params.do_save )
  save_plot( gcf, axs, l, [xcats, gcats, pcats], subdir, params );
end

end

function mask = get_base_mask(labels, mask_func)

mask = mask_func( labels, hwwa.get_approach_avoid_mask(labels) );

end

function save_plot(fig, axs, labels, cats, subdir, params)

shared_utils.plot.fullscreen( fig );
shared_utils.plot.match_ylims( axs );
save_p = hwwa.approach_avoid_data_path( params, 'plots', 'target_distance', subdir );
dsp3.req_savefig( fig, save_p, labels, cats, params.prefix );

end