function hwwa_plot_n_plus_one_performance(labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.per_scrambled_type = false;
defaults.per_drug = true;
defaults.per_trial_type = true;
defaults.prefix = '';

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

init_each = num_init_each( params );
norm_init_each = cssetdiff( init_each, 'unified_filename' );
n_plus_one_initiated( labels, {'unified_filename', 'monkey'} );

mask = intersect( mask, findnone(labels, makecollapsed(labels, 'next_initiated')) );
[props, prop_labels] = proportions_of( labels, init_each, 'next_initiated', mask );

[normed_props, norm_labels] = ...
  hwwa.saline_normalize( props, prop_labels, norm_init_each );

%%

pl = plotlabeled.make_common();

fcats = { 'monkey' };
% plt_mask = find( prop_labels, 'ephron' );
plt_mask = rowmask( prop_labels );

replace( prop_labels, 'correct_false', 'incorrect' );
replace( prop_labels, 'correct_true', 'correct' );
replace( prop_labels, 'next-initiated_false', 'next-uninitiated' );
replace( prop_labels, 'next-initiated_true', 'next-initiated' );

fig_I = findall_or_one( prop_labels, fcats, plt_mask );

for i = 1:numel(fig_I)
  pcats = intersect( init_each, {'monkey', 'trial_type', 'scrambled_type'} );
  gcats = { 'drug' };
  xcats = { 'correct', 'next_initiated' };
  
  plt = props(fig_I{i});
  plt_labs = prune( prop_labels(fig_I{i}) );

  axs = pl.bar( plt, plt_labs, xcats, gcats, pcats );
  
  subdir = get_subdir_prefix( params );
  maybe_save_fig( gcf, plt_labs, [fcats, pcats], subdir, params );
end

end

function labels = n_plus_one_initiated(labels, each)

I = findall( labels, each );
next_cat = 'next_initiated';
addcat( labels, next_cat );

for i = 1:numel(I)
  init = cellstr( labels, 'initiated', I{i} );
  next_init = eachcell( @(x) sprintf('next-%s', x), init(2:end, :) );  
  next_init = [ next_init; {makecollapsed(labels, next_cat)} ];
  setcat( labels, next_cat, next_init, I{i} );
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

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params, 'plots' ...
    , 'n_plus_one_num_initiated', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec );
end

end

function mask = get_base_mask(labels, mask_func)

require_initiated = true;

mask = hwwa.get_approach_avoid_mask( labels, {}, require_initiated );
mask = mask_func( labels, mask );

end

function each = num_init_each(params)

each = { 'unified_filename', 'monkey', 'correct' };

if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end

if ( params.per_trial_type )
  each{end+1} = 'trial_type';
end

end