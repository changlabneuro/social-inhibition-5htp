function hwwa_plot_num_initiated_per_session(labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = false;
defaults.per_drug = true;
defaults.per_trial_type = true;
defaults.seed = 0;
defaults.permutation_test = false;
defaults.permutation_test_iters = 1e2;
defaults.compare_drug = true;
defaults.prefix = '';

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

init_each = num_init_each( params );
norm_init_each = cssetdiff( init_each, 'unified_filename' );

p_each = pcorr_each( params );

[num_init, init_labels] = hwwa.num_initiated( labels, init_each, mask );
[normed_init, norm_labels] = ...
  hwwa.saline_normalize( num_init, init_labels, norm_init_each );

[pcorr, pcorr_labels] = hwwa.percent_correct( labels, p_each, mask );

% scatter_pcorr( num_init, init_labels', pcorr, pcorr_labels', params );
plot_raw( num_init, init_labels, rowmask(init_labels), params );
% plot_normalized( normed_init, norm_labels, rowmask(norm_labels), params );

% ttest_outs = ttests( num_init, init_labels', {}, 'saline', '5-htp', rowmask(num_init) );
% if ( params.do_save )
%   save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'num_initiated', 'raw' );
%   dsp3.save_ttest2_outputs( ttest_outs, save_p, 'drug' );
% end

end

function outs = ttests(data, labels, each, a, b, mask)

outs = dsp3.ttest2( data, labels', each, a, b, 'mask', mask );

end

function scatter_pcorr(num_init, init_labels, pcorr, pcorr_labels, params)

[match_init, match_pcorr, match_labels] = ...
  match_init_pcorr( num_init, init_labels, pcorr, pcorr_labels );

fcats = {};
gcats = {};

if ( params.compare_drug )
  if ( params.per_scrambled_type )
    fcats{end+1} = 'scrambled_type';
  end

  gcats{end+1} = 'drug';
else
  if ( params.per_scrambled_type )
    gcats{end+1} = 'scrambled_type';
  end

  fcats{end+1} = 'drug';
end

subdir_prefix = get_subdir_prefix( params );

pcats = { 'trial_type' };
pcats = [ pcats, fcats ];

if ( ~params.per_trial_type )
  [fcats, pcats, gcats] = mult_setdiff( 'trial_type', fcats, pcats, gcats );
end
if ( ~params.per_drug )
  [fcats, pcats, gcats] = mult_setdiff( 'drug', fcats, pcats, gcats );
end

fig_I = findall_or_one( match_labels, fcats );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.marker_size = 10;
  pl.y_lims = [0.1, 1.4];
  pl.x_lims = [100, 1500];
  
  plt_init = match_init(fig_I{i});
  plt_pcorr = match_pcorr(fig_I{i});
  plt_labels = prune( match_labels(fig_I{i}) );

  [axs, ids] = pl.scatter( plt_init, plt_pcorr, plt_labels, gcats, pcats );
  plotlabeled.scatter_addcorr( ids, plt_init, plt_pcorr );
  
  if ( params.permutation_test )
    rng( params.seed );
    [ps, labs, each_I] = hwwa.permute_slope_differences( plt_init, plt_pcorr, plt_labels' ...
      , params.permutation_test_iters, pcats, gcats );
    
    show_permutation_test_perf( ps, labs, each_I, axs, ids );
  end

  shared_utils.plot.xlabel( axs(1), 'Num initiated trials' );
  shared_utils.plot.ylabel( axs(1), 'Percent correct' );
  
  subdir = sprintf( 'num_init_vs_pcorr_%s', subdir_prefix );

  maybe_save_fig( gcf(), plt_labels, [gcats, pcats], subdir, params );
end

end

function prefix = get_subdir_prefix(params)

prefix = params.prefix;

if ( params.per_scrambled_type )
  prefix = sprintf( '%s-per_scrambled_type', prefix );
end
if ( params.per_drug )
  prefix = sprintf( '%s-per_drug', prefix );
end
if ( params.per_trial_type )
  prefix = sprintf( '%s-per_trial_type', prefix );
end

% if ( params.compare_drug )
%   prefix = sprintf( '%s-%s', prefix, 'compare_drug' );
% else
%   prefix = sprintf( '%s-%s', prefix, 'compare_social' );
% end

end

function varargout = mult_setdiff(cs, varargin)

varargout = cell( size(varargin) );

for i = 1:numel(varargin)
  varargout{i} = setdiff( varargin{i}, cs );
end

end

function show_permutation_test_perf(ps, labs, each_I, axs, ids)

displayed = false( size(ps) );

for i = 1:numel(ids)
  ind = ids(i).index;
  matches = cellfun( @(x) any(ismember(x, ind)), each_I );
  assert( nnz(matches) == 1, 'Expected 1 match; got %d.', nnz(matches) );
  
  if ( ~displayed(matches) )
    p_str = sprintf( 'p-slope-comparison = %0.2f', ps(matches) );
    ax = ids(i).axes;
    text( ax, min(get(ax, 'xlim')), max(get(ax, 'ylim')), p_str );
  end
end

end

function [match_init, match_pcorr, match_init_labels] = ...
  match_init_pcorr(num_init, init_labels, pcorr, pcorr_labels)

[match_I, match_C] = findall( init_labels, 'unified_filename' );
match_each = { 'trial_type' };

match_init_labels = fcat();
match_init = [];
match_pcorr = [];

for i = 1:numel(match_I)  
  match_ind = find( pcorr_labels, match_C(:, i) );
  match_each_I = findall( pcorr_labels, match_each, match_ind );
  init_ind = match_I{i};
  
  for j = 1:numel(match_each_I)
    match_ind = match_each_I{j};
    assert( numel(match_ind) == numel(init_ind), 'Non-matching trial subsets.' );
    
    append( match_init_labels, pcorr_labels, match_ind );
    match_init = [ match_init; num_init(init_ind) ];
    match_pcorr = [ match_pcorr; pcorr(match_ind) ];
  end
end

end

function plot_raw(data, labels, mask, params)

fcats = {};
pcats = {};
gcats = { 'drug' };
xcats = {};

subdir = 'raw';

bar_plots( data, labels, fcats, xcats, gcats, pcats, mask, subdir, params );
% box_plots( data, labels, fcats, gcats, pcats, mask, subdir, params );

end

function plot_normalized(data, labels, mask, params)

fcats = {};
pcats = {};
gcats = { 'scrambled_type' };
xcats = {};

subdir = 'normalized';

bar_plots( data, labels, fcats, xcats, gcats, pcats, mask, subdir, params );
% box_plots( data, labels, fcats, gcats, pcats, mask, subdir, params );

end

function box_plots(data, labels, fcats, gcats, pcats, mask, subdir, params)

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  fig = gcf();
  clf( fig );
  
  subset = data(fig_I{i});
  subset_labels = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.add_points = true;
  pl.points_are = { 'monkey' };
  pl.marker_size = 2;
  
  stats_each = [ fcats, gcats, pcats ];
  
  [axs, components] = pl.boxplot( subset, subset_labels, gcats, pcats );
  maybe_save_fig( fig, subset_labels, stats_each, subdir, params );
  
  sr_outs = dsp3.signrank1( subset, subset_labels', stats_each ...
    , 'signrank_inputs', {1} ...
  );
end

end

function bar_plots(data, labels, fcats, xcats, gcats, pcats, mask, subdir, params)

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  fig = gcf();
  clf( fig );
  
  subset = data(fig_I{i});
  subset_labels = prune( labels(fig_I{i}) );
  
  pl = plotlabeled.make_common();
  pl.add_points = true;
  pl.points_are = { 'monkey' };
  pl.marker_size = 2;
  
  axs = pl.bar( subset, subset_labels, xcats, gcats, pcats );
  maybe_save_fig( fig, subset_labels, [fcats, gcats, pcats], subdir, params );
end

end

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'num_initiated', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec );
end

end

function mask = get_base_mask(labels, mask_func)

require_initiated = true;

mask = hwwa.get_approach_avoid_mask( labels, {}, require_initiated );
mask = mask_func( labels, mask );

end

function each = pcorr_each(params)

each = { 'unified_filename', 'monkey', 'trial_type' };

if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end
if ( ~params.per_trial_type )
  each = setdiff( each, 'trial_type' );
end

end

function each = num_init_each(params)

each = { 'unified_filename', 'monkey' };

if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end

end