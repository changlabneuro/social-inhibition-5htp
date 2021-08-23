function hwwa_relate_to_pupil_size(pup, y, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_monkey = false;
defaults.per_scrambled_type = false;
defaults.second_measure_name = 'rt';
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e3;
defaults.remove_x_outliers = false;
defaults.outlier_std_thresh = 2;
defaults.seed = 0;

params = hwwa.parsestruct( defaults, varargin );

assert_ispair( pup, labels );
assert_ispair( y, labels );

mask = get_base_mask( labels, params.mask_func );

social_v_scrambled( pup, y, labels', mask, params );

end

function social_v_scrambled(amp, vel, labels, mask, params)

fcats = { 'correct' };
% pcats = { 'trial_type', 'scrambled_type' };

pcats = { 'trial_type' };
% gcats = {};
gcats = { 'drug' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  subdir = 'per_monkey';
else
  subdir = 'across_monkeys';
end

if ( params.per_scrambled_type )
  pcats{end+1} = 'scrambled_type';
  subdir = sprintf( '%s-per_scrambled_type', subdir );
end

if ( params.remove_x_outliers )
  outliers_each = union( unique([fcats, gcats, pcats]), {'unified_filename'} );
  non_outlier_inds = non_outlier_indices( amp, labels, params.outlier_std_thresh, outliers_each, mask );
  amp = amp(non_outlier_inds);
  vel = vel(non_outlier_inds);
  labels = prune( labels(non_outlier_inds) );
  mask = rowmask( amp );
end

pcats = columnize( union(pcats, fcats) )';

make_scatters( amp, vel, labels', mask, fcats, gcats, pcats, subdir, params );

end

function ind = non_outlier_indices(data, labels, thresh, each, mask)

each_I = findall( labels, each, mask );
means = bfw.row_nanmean( data, each_I );
devs = rowop( data, each_I, @nanstd );

non_outlier_inds = cell( size(each_I) );

for i = 1:numel(each_I)
  ub = means(i) + devs(i) * thresh;
  lb = means(i) - devs(i) * thresh;
  x = data(each_I{i});
  is_outlier = x < lb | x > ub;
  non_outlier_inds{i} = each_I{i}(~is_outlier);
end

ind = vertcat( non_outlier_inds{:} );

end

function make_scatters(amp, vel, labels, mask, fcats, gcats, pcats, subdir, params)

mask = intersect( mask, find(~isnan(amp) & ~isnan(vel)) );

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.y_lims = [0, 1];
  
  a = amp(fig_I{i});
  v = vel(fig_I{i});
  l = prune( labels(fig_I{i}) );
  
  [axs, ids] = pl.scatter( a, v, l, gcats, pcats );
  [hs, store_stats] = plotlabeled.scatter_addcorr( ids, a, v );
  
  if ( params.permutation_test )
    [ps, labs, each_I] = hwwa.permute_slope_differences( a, v, l ...
      , params.permutation_test_iters, pcats, gcats );
    
    hwwa.show_slope_permutation_test_performance( ps, each_I, ids );
  end
  
  xlabel( axs(1), 'Pupil size' );
  ylabel( axs(1), params.second_measure_name );
  
  subdir_ = fullfile( subdir, params.second_measure_name );
  maybe_save_fig( gcf, l, [fcats, pcats, gcats], subdir_, params );
end

end

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params ...
    , 'plots', 'pupil_size_vs_measure', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec );
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels );
mask = mask_func( labels, mask );

end