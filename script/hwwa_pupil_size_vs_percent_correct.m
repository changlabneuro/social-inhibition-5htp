function hwwa_pupil_size_vs_percent_correct(pupil_size, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.num_quantiles = 30;
defaults.mask_func = @(labels) rowmask(labels);
defaults.apply_trial_std_threshold = true;
defaults.trial_std_threshold_num_devs = 2;
defaults.prefix = '';
defaults.per_monkey = true;

params = hwwa.parsestruct( defaults, varargin );

assert_ispair( pupil_size, labels );

mask = get_base_mask( labels, params.mask_func );

plot_hists( pupil_size, labels', mask, params );

end

function plot_hists(pupil_size, labels, mask, params)


quant_each = { 'monkey' };
pcorr_each = { 'trial_type', 'drug', 'monkey' };

quant_I = findall( labels, quant_each, mask );

bin_cat = 'pupil_bin';

hist_labs = fcat();
hist_dat = {};
tot_num_selected = {};

for i = 1:numel(quant_I)
  subset_pupil = pupil_size(quant_I{i});
  quants = quantile( subset_pupil, params.num_quantiles-1 );
  edges = [ -inf, quants, inf ];
  
  for j = 1:numel(edges)-1
    within_edge = subset_pupil >= edges(j) & subset_pupil < edges(j+1);
    subset_mask = quant_I{i}(within_edge);
    
    [pcorr_dat, pcorr_labs, num_selected] = hwwa.percent_correct( labels, pcorr_each, subset_mask );
    addcat( pcorr_labs, bin_cat );
    setcat( pcorr_labs, bin_cat, sprintf('%s__%d', bin_cat, j) );
    
    hist_dat{end+1, 1} = pcorr_dat;
    tot_num_selected{end+1, 1} = num_selected;
    append( hist_labs, pcorr_labs );
  end
end

hist_dat = vertcat( hist_dat{:} );
tot_num_selected = vertcat( tot_num_selected{:} );

if ( params.apply_trial_std_threshold )
  n_devs = params.trial_std_threshold_num_devs;

  threshold_I = findall( hist_labs, pcorr_each );
  for i = 1:numel(threshold_I)
    subset_num_selected = tot_num_selected(threshold_I{i});
    dev = std( subset_num_selected );
    med = mean( subset_num_selected );
    oob = subset_num_selected < (med - dev*n_devs) | subset_num_selected > (med + dev*n_devs);

    threshold_I{i}(oob) = [];
  end

  hist_mask = vertcat( threshold_I{:} );
  params.prefix = sprintf( '%s%s%d', params.prefix, 'with-threshold-', n_devs );
else
  hist_mask = rowmask( hist_labs );
  params.prefix = sprintf( '%s%s', params.prefix, 'no-threshold' );
end

%%

% plot_hist_bars( hist_dat, hist_labs, bin_cat, params );
plot_hist_scatters( hist_dat, hist_labs, bin_cat, hist_mask, params );

end

function plot_hist_scatters(hist_dat, hist_labs, bin_cat, mask, params)

%%

pl = plotlabeled.make_common();

fcats = {};

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
end

gcats = { 'drug' };
pcats = [ {'trial_type'}, fcats ];

bin_prefix = sprintf( '%s__', bin_cat );
xs = fcat.parse( cellstr(hist_labs, bin_cat), bin_prefix );
fig_I = findall_or_one( hist_labs, fcats, mask );
figs = gobjects( numel(fig_I), 1 );
axs = cell( size(figs) );

for i = 1:numel(fig_I)
  f = figure(i);
  figs(i) = f;
  
  X = xs(fig_I{i});
  Y = hist_dat(fig_I{i});
  
  [axs{i}, ids] = pl.scatter( X, Y, hist_labs(fig_I{i}), gcats, pcats );
  plotlabeled.scatter_addcorr( ids, X, Y );
end

axs = vertcat( axs{:} );
shared_utils.plot.match_ylims( axs );

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pupil-vs-pcorr' );
  for i = 1:numel(figs)
    shared_utils.plot.fullscreen( figs(i) );
    dsp3.req_savefig( figs(i), save_p, prune(hist_labs(fig_I{i})), [fcats, pcats], params.prefix );
  end
end

end

function plot_hist_bars(hist_dat, hist_labs, bin_cat, params)

bin_prefix = sprintf( '%s__', bin_cat );
bins = combs( hist_labs, bin_cat );
[~, order] = sort( fcat.parse(bins, bin_prefix) );
bin_order = bins(order);

pl = plotlabeled.make_common();
pl.x_order = bin_order;

fcats = { 'monkey' };
xcats = { bin_cat };
gcats = {};
pcats = [ {'trial_type', 'drug'}, fcats ];

[figs, axs, I] = pl.figures( @bar, hist_dat, hist_labs, fcats, xcats, gcats, pcats );

end

function mask = get_base_mask(labels, mask_func)

mask = mask_func( labels );

end