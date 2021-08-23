function outs = hwwa_plot_time_between_initiated(durs, labels, varargin)

assert_ispair( durs, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.each = { 'unified_filename', 'drug' };
defaults.plot_cats = { {'drug'}, {}, {} };
defaults.normalize = false;

params = hwwa.parsestruct( defaults, varargin );
mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

[labs, I] = keepeach( labels', params.each, mask );
mean_durs = bfw.row_nanmean( durs, I );

if ( params.normalize )
  [mean_durs, labs] = hwwa.saline_normalize( mean_durs, labs', norm_each(params.each) );
end

%%

pl = plotlabeled.make_common();
pl.points_are = 'monkey';
pl.add_points = true;

xcats = params.plot_cats{1};
gcats = params.plot_cats{2};
pcats = params.plot_cats{3};

axs = pl.bar( mean_durs, labs, xcats, gcats, pcats );

shared_utils.plot.ylabel( axs(1), 'Time between trials (s)' );

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, labs, [xcats, gcats, pcats], 'time_between_trials', params );
end

outs = struct();
outs.means = mean_durs;
outs.labels = labs';

%%

rs_each = setdiff( horzcat(params.plot_cats{:}), {'drug'} );
rs_outs = dsp3.ranksum( mean_durs, labs', rs_each, 'saline', '5-htp' );

if ( params.do_save )
  stat_p = hwwa.approach_avoid_data_path( params, 'plots', 'time_between_trials/stats' );
  dsp3.save_ranksum_outputs( rs_outs, stat_p );
end

end

function each = norm_each(each)

each = setdiff( each, {'unified_filename'} );

end