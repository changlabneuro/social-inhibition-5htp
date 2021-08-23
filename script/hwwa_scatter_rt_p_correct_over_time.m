function hwwa_scatter_rt_p_correct_over_time(pcorr, pcorr_labels, rt, rt_labels, slide_each, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e3;
defaults.prefix = '';
defaults.per_monkey = true;
defaults.per_rt_quantile = false;
defaults.seed = 0;
defaults.remove_outliers = false;
defaults.std_threshold = 2;
defaults.outlier_each_func = @(slide_each) setdiff(slide_each, {'unified_filename'});
defaults.order_each_func = @(slide_each) slide_each;
defaults.x_lims = [];
defaults.y_lims = [];
defaults.normalize_pcorr = false;
defaults.anova_each = {};
params = hwwa.parsestruct( defaults, varargin );

assert_ispair( pcorr, pcorr_labels );
assert_ispair( rt, rt_labels );

mask = get_base_mask( pcorr_labels, params.mask_func );

order_each = params.order_each_func( slide_each );
[rt, rt_labels] = order_rt( pcorr_labels, rt, rt_labels', order_each, mask );
pcorr = pcorr(mask);
pcorr_labels = pcorr_labels(mask);

if ( params.remove_outliers )
  outlier_each = params.outlier_each_func( slide_each );
  is_rt_outlier = find_outliers( rt, rt_labels, outlier_each, params.std_threshold );
  is_p_outlier = find_outliers( pcorr, pcorr_labels, outlier_each, params.std_threshold );  
  
  non_outlier = find( ~is_rt_outlier & ~is_p_outlier );
%   non_outlier = find( ~is_p_outlier );
  rt = indexpair( rt, rt_labels, non_outlier );
  pcorr = indexpair( pcorr, pcorr_labels, non_outlier );
end

scatter_saline_v_5htp( pcorr, pcorr_labels', rt, rt_labels', slide_each, params );

end

function [rt, rt_labels, pcorr_labels] = order_rt(pcorr_labels, rt, rt_labels, slide_each, mask)

[pcorr_I, pcorr_C] = findall( pcorr_labels, slide_each, mask );
assert( unique(cellfun(@numel, pcorr_I)) == 1, 'Expected 1 trial per `slide_each` index.' );

rt_inds = nan( size(pcorr_I) );
missing_cats = setdiff( getcats(rt_labels), getcats(pcorr_labels) );

for i = 1:numel(pcorr_I)
  rt_ind = find( rt_labels, pcorr_C(:, i) );
  
  assert( numel(rt_ind) == 1, 'Expected 1 trial to match "%s"; got %d.' ...
    , strjoin(pcorr_C(:, i)), numel(rt_ind) );
  
  for j = 1:numel(missing_cats)
    missing_val = cellstr( rt_labels, missing_cats{j}, rt_ind );
    addsetcat( pcorr_labels, missing_cats{j}, missing_val, pcorr_I{i} );
  end
  
  rt_inds(i) = rt_ind;
end

rt = rt(rt_inds, :);
rt_labels = rt_labels(rt_inds);

end

function scatter_saline_v_5htp(pcorr, pcorr_labels, rt, rt_labels, slide_each, params)

assert( isequal(size(pcorr), size(rt)), 'Sizes must match between rt and pcorr.' );

data_type = 'rt_vs_pcorr';

num_cols = size( pcorr, 2 );
repmat( pcorr_labels, num_cols );
repmat( rt_labels, num_cols );

pcorr_all = pcorr(:);
rt_all = rt(:);

nan_mask = ~isnan( rt_all ) & ~isnan( pcorr_all );

assert_ispair( pcorr_all, pcorr_labels );
assert_ispair( rt_all, rt_labels );

rt_use = rt_all(nan_mask);
pcorr_use = pcorr_all(nan_mask);
pcorr_labels = pcorr_labels(find(nan_mask));

fcats = intersect( add_specificity({}, params), slide_each );
pcats = intersect( {'trial_type'}, slide_each );
gcats = { 'drug' };

pcats = union( pcats, fcats );
total_spec = unique( [fcats, gcats, pcats] );

fig_I = findall_or_one( pcorr_labels, fcats );

for i = 1:numel(fig_I)
  figure(i);
  
  subset_rt = rt_use(fig_I{i});
  subset_pcorr = pcorr_use(fig_I{i});
  subset_labels = prune( pcorr_labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.x_lims = params.x_lims;
  pl.y_lims = params.y_lims;
  [axs, ids] = pl.scatter( subset_rt, subset_pcorr, subset_labels, gcats, pcats );
  plotlabeled.scatter_addcorr( ids, subset_rt, subset_pcorr );
  
  if ( params.permutation_test )
    if ( ~isempty(params.seed) )
      rng( params.seed );
    end
    
    [ps, labs, each_I] = ...
      hwwa.permute_slope_differences( subset_rt, subset_pcorr, subset_labels' ...
      , params.permutation_test_iters, pcats, gcats );
    
    hwwa.show_slope_permutation_test_performance( ps, each_I, ids );
  end
  
  xlabel( axs(1), 'Reaction time' );
  ylabel( axs(1), '% correct' );
  
  hwwa.maybe_save_fig( gcf, subset_labels', total_spec ...
    , data_type, 'scatter', params );
  
  if ( params.normalize_pcorr )
    [subset_pcorr, subset_labels] = ...
      hwwa.saline_normalize( subset_pcorr, subset_labels', total_spec );
  end
  
  axs = pl.boxplot( subset_pcorr, subset_labels, gcats, pcats );
  hwwa.maybe_save_fig( gcf, subset_labels', [fcats, gcats, pcats] ...
    , data_type, 'boxplot', params );
end

if ( params.normalize_pcorr )
  [pcorr_use, pcorr_labels] = hwwa.saline_normalize( pcorr_use, pcorr_labels' ...
    , total_spec );
end

anova_each = params.anova_each;
pcorr_anova_stats( pcorr_use, pcorr_labels, anova_each, total_spec, data_type, '', params );

end

function pcorr_anova_stats(pcorr, labels, each, factors, data_type, subdir, params)

%%

factors = setdiff( factors, each );
factors(isuncat(labels, factors)) = [];

anova_outs = dsp3.anovan( pcorr, labels, each, factors ...
  , 'remove_nonsignificant_comparisons', false ...
);

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', data_type, subdir );
  dsp3.save_anova_outputs( anova_outs, save_p, [factors, each] );
end

end

function is_outlier = find_outliers(data, labels, each, std_threshold)

assert_ispair( data, labels );
each_I = findall( labels, each );
means = bfw.row_nanmean( data, each_I );
devs = rowop( data, each_I, @nanstd );

is_outlier = false( size(data) );
for i = 1:numel(each_I)
  ind = each_I{i};
  subset = data(ind);
  lb = means(i) - devs(i) * std_threshold;
  ub = means(i) + devs(i) * std_threshold;
  is_outlier(ind) = subset < lb | subset > ub;  
end

end

function spec = add_specificity(spec, params)

if ( params.per_monkey )
  spec = union( spec, {'monkey'} );
end
if ( params.per_rt_quantile )
  spec = union( spec, {'rt-quantile'} );
end

end

function mask = get_base_mask(labels, mask_func)

mask = mask_func( labels, rowmask(labels) );

end