function hwwa_pupil_vs_percent_correct(pupil, labels, varargin)

assert_ispair( pupil, labels );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.pcorr_each = {};
defaults.pupil_each = {};
defaults.corr_each = {};
defaults.corr_match = {};
defaults.scatter_cats = { {'drug'}, {'trial_type'} };
defaults.normalize = true;
defaults.per_panel_corr = false;

params = hwwa.parsestruct( defaults, varargin );

mask = hwwa.get_approach_avoid_base_mask( labels, params.mask_func );

[pcorr, pcorr_labs] = hwwa.percent_correct( labels', params.pcorr_each, mask );
[pup_labs, pup_I] = keepeach( labels', params.pupil_each, mask );
pup = bfw.row_nanmean( pupil, pup_I );

if ( params.normalize )
  [pcorr, pcorr_labs] = ...
    hwwa.saline_normalize( pcorr, pcorr_labs', prune_norm_each(params.pcorr_each) );
  [pup, pup_labs] = ...
    hwwa.saline_normalize( pup, pup_labs', prune_norm_each(params.pupil_each) );
end

[pcorr_x, pup_y, labs] = match( pcorr, pcorr_labs', pup, pup_labs', params.corr_each, params.corr_match );
% [x, y, labs] = match( pup, pup_labs', pcorr, pcorr_labs', params.corr_each, params.corr_match );

%%

x = pup_y;
y = pcorr_x;

scatter( x, y, labs', params.scatter_cats, 'pupil', 'percent-correct', params );

end

function scatter(x, y, labels, cats, xlab, ylab, params)

pl = plotlabeled.make_common();

gcats = cats{1};
pcats = cats{2};

[axs, ids] = pl.scatter( x, y, labels, gcats, pcats );

if ( params.per_panel_corr )
  ax_set = arrayfun( @(x) {x}, axs );
  id_set = arrayfun( @(x) {x}, ids );
else
  ax_set = { axs };
  id_set = { ids };
end

for i = 1:numel(ax_set)
  if ( params.normalize )
    hwwa.add_normalized_corr( ax_set{i}(1), id_set{i}, x, y );
  else
    hwwa.add_5htp_saline_all_corr( ax_set{i}(1), id_set{i}, x, y );
  end
end

shared_utils.plot.xlabel( axs, xlab );
shared_utils.plot.ylabel( axs, ylab );

for j = 1:numel(ids)
  ser = ids(j).series;
  arrayfun( @(x) set(x, 'SizeData', 6), ser, 'un', 0 );
end

[ps, p_labs] = hwwa.permute_slope_differences( x, y, labels', 1e3, gcats, pcats );

if ( params.do_save )
  hwwa.approach_avoid_save_fig( gcf, axs, labels, [gcats, pcats], 'pupil_vs_pcorr', params );
end

end

function [x, y, labels] = match(data1, labels1, data2, labels2, each, of)

each1 = findall_or_one( labels1, each );

x = [];
y = [];
labels = fcat();

for i = 1:numel(each1)
  [of1, C] = findall( labels1, of, each1{i} );
  
  for j = 1:numel(of1)
    ind1 = of1{j};
    assert( numel(ind1) == 1 );
    
    search_for = C(:, j);
    ind2 = find( labels2, search_for );
    
    for k = 1:numel(ind2)
      x = [ x; data1(ind1) ];
      y = [ y; data2(ind2(k)) ];

      l1 = prune( labels1(ind1) );
      l2 = prune( labels2(ind2(k)) );
      l = join( l1', l2 );

      append( labels, l );
    end
  end
end

assert_ispair( x, labels );
assert_ispair( y, labels );

end

function each = prune_norm_each(each)

each = setdiff( each, {'unified_filename', 'drug'} );

end