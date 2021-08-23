function hwwa_relate_quantiles_x_percent_correct_y(x, y, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.x_label = '';
defaults.y_label = 'P-correct';
defaults.per_scrambled_type = false;
defaults.per_trial_type = false;
defaults.x_lims = [];
defaults.y_lims = [0, 1];
defaults.quantile_index_x = false;
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e2;
defaults.seed = 0;

params = hwwa.parsestruct( defaults, varargin );

assert_ispair( x, labels );
assert_ispair( y, labels );

base_spec = get_base_specificity( params );
base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

plot_info = struct();
plot_info.fcats = intersect( base_spec, {'trial_type'} );
plot_info.pcats = intersect( base_spec, {'trial_type', 'scrambled_type', 'x-quantile'} );
plot_info.gcats = {'drug'};
plot_info.mask = base_mask;
plot_info.data_type = 'quantiles-x-vs-percent-correct';
plot_info.subdir = params.x_label;

if ( params.quantile_index_x )
  [labels, x] = quantile_index_x( labels' );
end

plot_scatter( x, y, labels', plot_info, params );

end

function plot_scatter(x, y, labels, plot_info, params)

fcats = plot_info.fcats;
gcats = plot_info.gcats;
pcats = plot_info.pcats;
mask = plot_info.mask;

pcats = union( pcats, fcats );
tot_spec = unique( [fcats, gcats, pcats] );

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  ind = fig_I{i};
  subx = x(ind);
  suby = y(ind);
  sub_labels = prune( labels(ind) );
  
  pl = plotlabeled.make_common();
  pl.y_lims = params.y_lims;
  pl.x_lims = params.x_lims;
  
  [axs, ids] = pl.scatter( subx, suby, sub_labels, gcats, pcats );
  plotlabeled.scatter_addcorr( ids, subx, suby );
  
  if ( params.permutation_test )
    if ( ~isempty(params.seed) )
      rng( params.seed );
    end
    
    [ps, plabs, p_I] = hwwa.permute_slope_differences( subx, suby, sub_labels' ...
      , params.permutation_test_iters, setdiff(tot_spec, gcats), gcats );
    hwwa.show_slope_permutation_test_performance( ps, p_I, ids );
  end
  
  xlabel( axs(1), params.x_label );
  ylabel( axs(1), params.y_label );
  
  hwwa.maybe_save_fig( gcf, sub_labels, tot_spec, plot_info.data_type ...
    , plot_info.subdir, params );
end

end

function [labels, x] = quantile_index_x(labels)

quant_labs = cellstr( labels, 'x-quantile' );
quant_nums = fcat.parse( quant_labs, 'x-quantile__' );
collapsecat( labels, 'x-quantile' );
x = quant_nums;

end

function spec = get_base_specificity(params)

spec = {'x-quantile'};

if ( params.per_scrambled_type )
  spec{end+1} = 'scrambled_type';
end
if ( params.per_trial_type )
  spec{end+1} = 'trial_type';
end

end