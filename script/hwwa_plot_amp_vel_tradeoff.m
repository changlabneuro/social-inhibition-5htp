function hwwa_plot_amp_vel_tradeoff(amp, vel, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_monkey = false;
defaults.permutation_test_iters = 1e3;
defaults.permutation_test = false;
defaults.seed = 0;
defaults.gcats = {};
defaults.pcats = {};

params = hwwa.parsestruct( defaults, varargin );

assert_ispair( amp, labels );
assert_ispair( vel, labels );

mask = get_base_mask( labels, params.mask_func );

revisited( amp, vel, labels', mask, params );

% social_v_scrambled( amp, vel, labels', mask, true, params );
% social_v_scrambled( amp, vel, labels', mask, false, params );

end

function revisited(amp, vel, labels, mask, params)

pl = plotlabeled.make_common();
pl.y_lims = [0, 2e3];
pl.x_lims = [0, 100];

mask = intersect( mask, intersect(find(~isnan(vel)), find(~isnan(amp))) );

a = amp(mask);
v = vel(mask);
l = prune( labels(mask) );

[axs, ids] = pl.scatter( a, v, l, params.gcats, params.pcats );
xlabel( axs(1), 'Amplitude (deg)' );
ylabel( axs(1), 'Peak velocity (deg/s)' );

[~, stats] = plotlabeled.scatter_addcorr( ids, a, v );

try
  [ps, p_labels] = hwwa.permute_slope_differences( a, v, l, params.permutation_test_iters, params.pcats, params.gcats );
catch err
  warning( err.message );
end

maybe_save_fig( gcf, l', [params.gcats, params.pcats], '', params );

end

function social_v_scrambled(amp, vel, labels, mask, is_per_image_cat, params)

fcats = {};
% pcats = { 'trial_type', 'scrambled_type' };
% gcats = { 'drug' };

pcats = { 'trial_type', 'drug' };
gcats = { 'scrambled_type' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  subdir = 'per_monkey';
else
  subdir = 'across_monkeys';
end

if ( is_per_image_cat )
  fcats{end+1} = 'target_image_category';
end

pcats = columnize( union(pcats, fcats) )';

make_scatters( amp, vel, labels', mask, fcats, gcats, pcats, subdir, params );

end

function make_scatters(amp, vel, labels, mask, fcats, gcats, pcats, subdir, params)

mask = intersect( mask, find(~isnan(amp) & ~isnan(vel)) );

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.y_lims = [0, 5e3];
  pl.x_lims = [0, 250];
  
  a = amp(fig_I{i});
  v = vel(fig_I{i});
  l = prune( labels(fig_I{i}) );
  
  [axs, ids] = pl.scatter( a, v, l, gcats, pcats );
  [hs, store_stats] = plotlabeled.scatter_addcorr( ids, a, v );
  
  if ( params.permutation_test )
    rng( params.seed );
    perm_iters = params.permutation_test_iters;
    [ps, labs] = permute_slope_differences( a, v, l, perm_iters, pcats, gcats );
    
    for j = 1:numel(ids)
      sel = ids(j).selectors(numel(gcats)+1:end);
      ax = ids(j).axes;
      
      match_ind = find( labs, sel );
      assert( numel(match_ind) == 1 );
      p = ps(match_ind);
      
      hold( ax, 'on' );
      text( ax, max(get(ax, 'xlim')), max(get(ax, 'ylim')), sprintf('P=%0.2f', p) );
    end
  end
  
%   model = fitlm( amps(ind), vels(ind) );
  
  xlabel( axs(1), 'Saccade amplitude' );
  ylabel( axs(1), 'Saccade peak velocity' );
  
  maybe_save_fig( gcf, l, [fcats, pcats, gcats], subdir, params );
end

end

function [ps, labs] = permute_slope_differences(amp, vel, labels, iters, each, compare)

[labs, each_I] = keepeach_or_one( labels', each );
ps = nan( numel(each_I), 1 );

parfor i = 1:numel(each_I)
  ind = each_I{i};
  
  compare_I = findall( labels, compare, ind );
  assert( numel(compare_I) == 2, 'Expected 2 groups to compare; got %d.' ...
    , numel(compare_I) );
  
  num_a = numel( compare_I{1} );
  null_diffs = zeros( iters, 1 );
  
  [lm1, lm2] = fit_models( amp, vel, compare_I{1}, compare_I{2} );
  real_diff = slope_comparison( lm1, lm2 );
  
  for j = 1:iters
    subset_a = sort( ind(randperm(numel(ind), num_a)) );
    subset_b = setdiff( ind, subset_a );
    
    [lm1, lm2] = fit_models( amp, vel, subset_a, subset_b );
    null_diffs(j) = slope_comparison( lm1, lm2 );
  end
  
  ps(i) = sum( null_diffs > real_diff ) / iters;
end

end

function [lm1, lm2] = fit_models(amp, vel, ind1, ind2)

a1 = amp(ind1);
v1 = vel(ind1);

a2 = amp(ind2);
v2 = vel(ind2);

lm1 = fitlm( a1, v1 );
lm2 = fitlm( a2, v2 );

end

function b = slope_comparison(lm1, lm2)

beta1 = lm1.Coefficients.Estimate(2);
beta2 = lm2.Coefficients.Estimate(2);

b = abs( beta1 - beta2 );

end

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params ...
    , 'plots', 'saccade_amp_v_vel', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec );
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels );
mask = mask_func( labels, mask );
mask = hwwa.find_correct_go_incorrect_nogo( labels, mask );

end