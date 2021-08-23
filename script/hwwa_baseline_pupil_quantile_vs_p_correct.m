function hwwa_baseline_pupil_quantile_vs_p_correct(pupil, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = true;
defaults.per_trial_type = true;
defaults.per_monkey = false;
defaults.num_quantiles = 10;
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e3;
defaults.y_lims = [];
defaults.prefix = '';
params = hwwa.parsestruct( defaults, varargin );

base_mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

assert_ispair( pupil, labels );

quantiles = pupil_quantile( pupil, labels, params.num_quantiles, base_mask );
dsp3.add_quantile_labels( labels, quantiles, quant_cat() );

[pcorr, pcorr_labels] = hwwa.percent_correct( labels, pcorr_each(params), base_mask );

fcats = intersect( pcorr_each(params), {'monkey'} );
pcats = intersect( pcorr_each(params), {'scrambled_type', 'trial_type'} );
pcats = union( pcats, fcats );
gcats = 'drug';

subdir = '';

plot_quantile_scatters( pcorr, pcorr_labels, rowmask(pcorr), fcats, pcats, gcats, subdir, params );
% plot_quantile_lines( pcorr, pcorr_labels, rowmask(pcorr), fcats, pcats, gcats, subdir, params );

end

function each = pcorr_each(params)

each = { 'unified_filename', quant_cat() };
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

function plot_quantile_scatters(data, labels, mask, fcats, pcats, gcats, subdir, params)

fig_I = findall_or_one( labels, fcats, mask );
for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  
  y = data(fig_I{i});
  lab = prune( labels(fig_I{i}) );
  x = fcat.parse( cellstr(lab, quant_cat), sprintf('%s__', quant_cat) );

  [axs, ids] = pl.scatter( x, y, lab, gcats, pcats );
  plotlabeled.scatter_addcorr( ids, x, y );
  
  if ( params.permutation_test )  
    [ps, p_labs, each_I] = ...
      hwwa.permute_slope_differences( x, y, lab' ...
      , params.permutation_test_iters, pcats, gcats );
    hwwa.show_slope_permutation_test_performance( ps, each_I, ids );
  end
  
  if ( ~isempty(params.y_lims) )
    shared_utils.plot.set_ylims( axs, params.y_lims );
  end
  
  maybe_save_fig( gcf, lab, [fcats, pcats], fullfile(subdir, 'scatter'), params );
end

end

function plot_quantile_lines(x, labels, mask, fcats, pcats, gcats, subdir, params)

fig_I = findall_or_one( labels, fcats, mask );
for i = 1:numel(fig_I)
  plt = x(fig_I{i});
  lab = prune( labels(fig_I{i}) );
  
  quants = combs( labels, quant_cat() );
  quant_n = fcat.parse( quants, sprintf('%s__', quant_cat) );
  assert(~any(isnan(quant_n)));
  [~, order] = sort( quant_n );
  
  pl = plotlabeled.make_common();
  pl.y_lims = params.y_lims;
  pl.x_order = quants(order);
  axs = pl.errorbar( plt, lab, quant_cat(), gcats, pcats );
  
  maybe_save_fig( gcf, lab, [fcats, pcats], fullfile(subdir, 'lines'), params );
end

end

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pupil_quantile', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec, params.prefix );
end

end

function quantiles = pupil_quantile(pupil, labels, num_quantiles, mask)

quantiles = dsp3.quantiles_each( pupil, labels, num_quantiles, 'unified_filename', {}, mask );

end

function c = quant_cat()

c = 'pupil_quantile';

end