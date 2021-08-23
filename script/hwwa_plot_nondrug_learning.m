function hwwa_plot_nondrug_learning(behav_outputs, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.xlims = [];
defaults.var_bin_size = 10;
defaults.var_step_size = 1;
defaults.per_monkey = true;

params = hwwa.parsestruct( defaults, varargin );

labels = behav_outputs.labels';

hwwa.decompose_social_scrambled( labels );
add_day_number_labels( labels );

specs = get_specificities();
base_mask = get_base_mask( labels );

for i = 1:numel(specs)
%   plot_errorbar_timecourse( labels', specs{i}, base_mask, params );
%   plot_day_level_timecourse( labels', specs{i}, base_mask, params );
  plot_variance_timecourse( labels', specs{i}, base_mask, params );
end

end

function plot_day_level_timecourse(labels, spec, mask, params)

%%

fig_I = findall_or_one( labels, 'monkey', mask );

spec = csunion( spec, get_base_specificity() );
bin_spec = {'day'};

c_is_smooth = [ false ];
combts = dsp3.numel_combvec( c_is_smooth );

%%
for i = 1:numel(fig_I)

  [pcorr, pcorr_labels] = get_p_correct_data( [], labels, csunion(spec, bin_spec), fig_I{i} );
  [pltdat, pltlabs] = to_matrix( pcorr, pcorr_labels, spec, bin_spec );
  
  assert( ~any(isnan(pcorr)) );

  %%
  
  for j = 1:size(combts, 2)    
    is_smoothed = c_is_smooth( combts(1, j) );
    
    if ( is_smoothed )
      addtl_subdirs = 'smoothed';
    else
      addtl_subdirs = 'non-smoothed';
    end

    pl = plotlabeled.make_common();
    pl.add_errors = false;
    pl.smooth_func = @(x) smooth(x, 20);
    pl.add_smoothing = is_smoothed;
    pl.y_lims = [0, 1];

    gcats = { 'scrambled_type' };
    pcats = { 'target_image_category', 'trial_type', 'monkey' };

    axs = pl.lines( pltdat, pltlabs, gcats, pcats );
    
    if ( ~isempty(params.xlims) )
      shared_utils.plot.set_xlims( axs, params.xlims );
    end
    
    if ( params.do_save )
      save_p = get_base_plot_path( params, 'percent_correct', 'day_level' );
      save_p = fullfile( save_p, addtl_subdirs );
      
      dsp3.req_savefig( gcf, save_p, pltlabs, [pcats, gcats] );
    end
  end
end

end

function plot_variance_timecourse(labels, spec, mask, params)

%%

if ( ~params.per_monkey )
  collapsecat( labels, 'monkey' );
end

fig_I = findall_or_one( labels, 'monkey', mask );

spec = csunion( spec, get_base_specificity() );
bin_spec = { 'day_number', 'trial_bin' };

c_is_smooth = [ false, true ];
combts = dsp3.numel_combvec( c_is_smooth );

%%
for i = 1:numel(fig_I)

  [pcorr, pcorr_labels] = get_p_correct_data( [], labels, csunion(spec, bin_spec), fig_I{i} );
  [mat_corr, pltlabs] = to_matrix( pcorr, pcorr_labels, spec, bin_spec );
  
  bin_inds = hwwa.bin_indices( 1:size(mat_corr, 2), params.var_bin_size, params.var_step_size );
  
%   var_func = @(x) nanvar(x, [], 2);
%   var_func = @(x) nanstd(x, [], 2);
  var_func = @(x) nanstd(x, [], 2) ./ nanmean(x, 2);
  
  bin_corr = cellfun( @(x) var_func(mat_corr(:, x)), bin_inds, 'un', 0 );
  bin_corr = horzcat( bin_corr{:} );
  
  assert( ~any(isnan(pcorr)) );

  %%
  
  for j = 1:size(combts, 2)    
    is_smoothed = c_is_smooth( combts(1, j) );
    
    if ( is_smoothed )
      addtl_subdirs = 'smoothed';
    else
      addtl_subdirs = 'non-smoothed';
    end

    pl = plotlabeled.make_common();
    pl.add_errors = false;
    pl.smooth_func = @(x) smooth(x, 20);
    pl.add_smoothing = is_smoothed;
%     pl.y_lims = [0, 0.5];

    gcats = { 'scrambled_type' };
    pcats = { 'trial_type', 'monkey' };

    axs = pl.lines( bin_corr, pltlabs, gcats, pcats );
    
    if ( ~isempty(params.xlims) )
      shared_utils.plot.set_xlims( axs, params.xlims );
    end
    
    if ( params.do_save )
      save_p = get_base_plot_path( params, 'variance_percent_correct' );
      save_p = fullfile( save_p, addtl_subdirs );
      
      dsp3.req_savefig( gcf, save_p, pltlabs, [pcats, gcats] );
    end
  end
end

end

function plot_errorbar_timecourse(labels, spec, mask, params)

%%

collapsecat( labels, 'monkey' );

fig_I = findall_or_one( labels, 'monkey', mask );

spec = csunion( spec, get_base_specificity() );
% bin_spec = {'date', 'trial_bin'};
bin_spec = { 'day_number', 'trial_bin' };

c_is_smooth = [ true, false ];
combts = dsp3.numel_combvec( c_is_smooth );

%%
for i = 1:numel(fig_I)

  [pcorr, pcorr_labels] = get_p_correct_data( [], labels, csunion(spec, bin_spec), fig_I{i} );
  [pltdat, pltlabs] = to_matrix( pcorr, pcorr_labels, spec, bin_spec );
  
  assert( ~any(isnan(pcorr)) );

  %%
  
  for j = 1:size(combts, 2)    
    is_smoothed = c_is_smooth( combts(1, j) );
    
    if ( is_smoothed )
      addtl_subdirs = 'smoothed';
    else
      addtl_subdirs = 'non-smoothed';
    end

    pl = plotlabeled.make_common();
    pl.add_errors = false;
    pl.smooth_func = @(x) smooth(x, 20);
    pl.add_smoothing = is_smoothed;
    pl.y_lims = [0, 1];

    gcats = { 'scrambled_type' };
    pcats = { 'target_image_category', 'trial_type', 'monkey' };

    axs = pl.lines( pltdat, pltlabs, gcats, pcats );
    
    if ( ~isempty(params.xlims) )
      shared_utils.plot.set_xlims( axs, params.xlims );
    end
    
    if ( params.do_save )
      save_p = get_base_plot_path( params, 'percent_correct' );
      save_p = fullfile( save_p, addtl_subdirs );
      
      dsp3.req_savefig( gcf, save_p, pltlabs, [pcats, gcats] );
    end
  end
end

end

function [new_mat, new_labels] = to_matrix(mat, labels, row_spec, col_spec)
col_spec = cellstr( col_spec );

[t, new_labels, col_labels] = tabular( labels, row_spec, col_spec );

bin_names = combs( col_labels, col_spec );
sorted_idx = sort_bin_names( bin_names, col_spec );

keep( col_labels, sorted_idx );
t = t(:, sorted_idx);

new_mat = nan( size(t) );
for i = 1:numel(t)
  if ( isempty(t{i}) ), continue; end
  new_mat(i) = mat(t{i});
end

end

function sorted_idx = sort_bin_names(bin_names, col_spec)

%%

assert( numel(col_spec) == size(bin_names, 1) );

sort_mat = nan( size(bin_names, 2), size(bin_names, 1) ); % transpose

for i = 1:numel(col_spec)
  spec = col_spec{i};
  bins = bin_names(i, :);
  
  switch ( spec )
    case 'trial_bin'  
      nums = trial_bin_numbers( bins );
    case 'date'
      nums = datenum( bins );
    case 'day'
      nums = datenum( hwwa.to_date(bins) );
    case 'day_number'
      nums = day_numbers( bins );
    otherwise
      error( 'Unhandled sort spec: "%s"', spec );
  end
  
  sort_mat(:, i) = nums;
end

[~, sorted_idx] = sortrows( sort_mat, 1:size(sort_mat, 2) );

end

function nums = trial_bin_numbers(trial_bins)

trial_bin_prefix = 'trial_bin__';
nums = cellfun( @(x) str2double(x(numel(trial_bin_prefix)+1:end)), trial_bins );

end

function nums = day_numbers(bins)

nums = parse_bin_numbers( bins, 'day_');

end

function nums = parse_bin_numbers(bins, prefix)

nums = cellfun( @(x) str2double(x(numel(prefix)+1:end)), bins );

end

function [pcorr, pcorr_labels] = get_p_correct_data(data, labels, data_each, mask)

% Calculate percent correct for each `data_each`
[pcorr, pcorr_labels] = hwwa.percent_correct( labels', data_each, mask );

end

function specs = get_specificities()

% specs = { ...
%   {'target_image_category', 'scrambled_type'} ...
% };

specs = { ...
  {'scrambled_type'} ...
};

end

function spec = get_base_specificity()
spec = { 'date', 'monkey', 'trial_type' };
end

function mask = get_base_mask(labels)

mask = [];

if ( count(labels, 'ephron') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'ephron' ...
    , @findor, hwwa.get_ephron_learning_days() ...
    ));
end

if ( count(labels, 'hitch') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'hitch' ...
    , @findor, hwwa.get_hitch_learning_days() ...
    ));
end

if ( count(labels, 'tar') > 0 )
  mask = union( mask, fcat.mask(labels ...
    , @find, 'tar' ...
    , @findor, hwwa.get_tarantino_learning_days() ...
    ));
end

% Only initiated trials.
mask = add_initiated_mask( labels, mask );

end

function mask = add_initiated_mask(labels, mask)

mask = fcat.mask( labels, mask ...
  , @find, 'initiated_true' ...
);

gonogo_error = find( labels, 'wrong_go_nogo' );
ok_trial = find( labels, 'correct_true' );
error_or_ok = union( gonogo_error, ok_trial );

mask = intersect( mask, error_or_ok );

end

function labels = add_day_number_labels(labels)

addcat( labels, 'day_number' );

monk_I = findall( labels, 'monkey' );
for i = 1:numel(monk_I)
  [date_I, date_C] = findall( labels, 'date', monk_I{i} );
  [~, sorted_ind] = sort( datenum(date_C) );
  
  date_I = date_I(sorted_ind);
  
  for j = 1:numel(date_I)
    setcat( labels, 'day_number', sprintf('day_%d', j), date_I{j} );
  end
end

end

function p = get_base_plot_path(params, varargin)

p = fullfile( hwwa.dataroot(params.config), 'plots', 'behavior' ...
  , hwwa.datedir, 'nondrug_learning', params.base_subdir, varargin{:} );

end