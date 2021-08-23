function hwwa_relate_num_initiated_to_pupil(pupil, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = true;
defaults.per_trial_type = true;
defaults.per_monkey = false;
defaults.y_lims = [];
defaults.prefix = '';
defaults.pupil_metric = 'baseline';
defaults.seed = 0;
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e3;
params = hwwa.parsestruct( defaults, varargin );

assert_ispair( pupil, labels );

base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

init_each = num_init_each( params );
[num_init, num_init_labels, each_I] = hwwa.num_initiated( labels', init_each, base_mask );
mean_pup = bfw.row_nanmean( pupil, each_I );

fcats = intersect( init_each, {'monkey'} );
pcats = intersect( init_each, {'scrambled_type', 'trial_type'} );
pcats = union( pcats, fcats );
gcats = 'drug';

subdir = params.pupil_metric;

plot_scatter( mean_pup, num_init, num_init_labels', fcats, pcats, gcats, subdir, params );

end

function each = num_init_each(params)

each = { 'unified_filename' };
if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end
if ( params.per_trial_type )
  each{end+1} = 'trial_type';
end
if ( params.per_monkey )
  each{end+1} = 'monkey';
end

end

function plot_scatter(x, y, labels, fcats, pcats, gcats, subdir, params)

fig_I = findall_or_one( labels, fcats );
for i = 1:numel(fig_I)
  tmp_x = x(fig_I{i});
  tmp_y = y(fig_I{i});
  
  lab = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.y_lims = params.y_lims;
  [axs, ids] = pl.scatter( tmp_x, tmp_y, lab, gcats, pcats );
  plotlabeled.scatter_addcorr( ids, tmp_x, tmp_y );  
  
  if ( params.permutation_test )
    if ( ~isempty(params.seed) )
      rng( params.seed );
    end
    
    [ps, p_labs, each_I] = ...
      hwwa.permute_slope_differences( tmp_x, tmp_y, lab' ...
      , params.permutation_test_iters, pcats, gcats );
    
    hwwa.show_slope_permutation_test_performance( ps, each_I, ids );
  end
  
  hwwa.maybe_save_fig( gcf, lab, [fcats, pcats], plot_data_type(), subdir, params );
end

end

function dt = plot_data_type()
dt = 'pupil_v_num_initiated';
end