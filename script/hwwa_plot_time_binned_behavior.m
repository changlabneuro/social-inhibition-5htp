function hwwa_plot_time_binned_behavior(data, time, labels, varargin)

validateattributes( data, {'double'}, {'2d', 'ncols', size(data, 2)} ...
 , mfilename, 'data' );

assert( numel(time) == size(data, 2), 'Time does not correspond to data 2nd dimension.' );

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @(labels, mask) mask;
defaults.saline_normalize = true;
defaults.time_limits = [-inf, inf];
defaults.per_monkey = false;
defaults.per_drug = true;

params = hwwa.parsestruct( defaults, varargin );

mask = get_base_mask( labels, params.mask_func );

% social_appetitive_v_threat_normalized( data, labels, time, mask, true, params );
% social_appetitive_v_threat_normalized( data, labels, time, mask, false, params );

appetitive_v_threat( data, labels, time, mask, false, params );
% appetitive_v_threat( data, labels, time, mask, true, params );
% appetitive_v_threat( data, labels, time ...
%   , find(labels, 'not-scrambled', mask), false, params );

% social_v_scrambled( data, labels, time, mask, params );
% social_v_scrambled_normalized( data, labels, time, mask, params );

end

function social_appetitive_v_threat_normalized(data, labels, time, mask, is_saline_norm, params)

subdir = 'social_appetitive_v_threat';

[subset_data, subset_labels, time] = make_subsets( data, labels, time, mask, params );

norm_each = { 'drug', 'target_image_category', 'trial_type', 'monkey' };
[subset_data, subset_labels] = hwwa.scrambled_normalize( subset_data, subset_labels', norm_each );

if ( is_saline_norm )
  [subset_data, subset_labels] = hwwa.saline_normalize( subset_data, subset_labels', norm_each );
  subdir = fullfile( subdir, 'saline_normalized' );
else
  subdir = fullfile( subdir, 'non_saline_normalized' );
end

subdir = fullfile( subdir, 'normalized' );
y_lims = [0.5, 1.5];

pl = plotlabeled.make_common();
pl.x = time;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.2);
pl.add_smoothing = true;
pl.y_lims = y_lims;

fcats = {};
gcats = { 'scrambled_type', 'target_image_category' };
pcats = { 'drug', 'trial_type' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats = union( pcats, 'monkey' );
  subdir = sprintf( '%s_per_monkey', subdir );
end

make_line_figures( pl, subset_data, subset_labels, fcats, gcats, pcats, subdir, params );

end

function appetitive_v_threat(data, labels, time, mask, is_norm, params)

subdir = 'appetitive_v_threat';

[subset_data, subset_labels, time] = make_subsets( data, labels, time, mask, params );

if ( is_norm )
  norm_each = { 'scrambled_type', 'target_image_category', 'trial_type', 'monkey' };
  [subset_data, subset_labels] = hwwa.saline_normalize( subset_data, subset_labels', norm_each );
  
  subdir = fullfile( subdir, 'normalized' );
  y_lims = [0.5, 1.5];
else
  subdir = fullfile( subdir, 'non-normalized' );
  y_lims = [0, 1];
end

pl = plotlabeled.make_common();
pl.x = time;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.2);
pl.add_smoothing = true;
pl.y_lims = y_lims;

fcats = {};
% gcats = { 'target_image_category' };
gcats = { 'drug' };
pcats = { 'trial_type', 'scrambled_type' };

if ( params.per_drug )
%   fcats{end+1} = 'drug';
%   pcats{end+1} = 'drug';
  fcats{end+1} = 'target_image_category';
  pcats{end+1} = 'target_image_category';
end

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats = union( pcats, 'monkey' );
  subdir = sprintf( '%s_per_monkey', subdir );
end

make_line_figures( pl, subset_data, subset_labels, fcats, gcats, pcats, subdir, params );

end

function social_v_scrambled_normalized(data, labels, time, mask, params)

[subset_data, subset_labels, time] = make_subsets( data, labels, time, mask, params );

norm_each = { 'scrambled_type', 'trial_type', 'monkey' };
[subset_data, subset_labels] = hwwa.saline_normalize( subset_data, subset_labels', norm_each );

pl = plotlabeled.make_common();
pl.x = time;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.2);
pl.add_smoothing = true;
pl.y_lims = [0.5, 1.5];

subdir = 'sal_vs_5htp/';

fcats = {};
gcats = { 'scrambled_type' };
pcats = { 'drug', 'trial_type' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats = union( pcats, 'monkey' );
  subdir = sprintf( '%s_per_monkey', subdir );
end

make_line_figures( pl, subset_data, subset_labels, fcats, gcats, pcats, subdir, params );

end

function social_v_scrambled(data, labels, time, mask, params)

[subset_data, subset_labels, time] = make_subsets( data, labels, time, mask, params );

pl = plotlabeled.make_common();
pl.x = time;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.05);
pl.add_smoothing = true;
pl.y_lims = [0, 1];

subdir = 'social_v_scrambled/';

fcats = {};
gcats = { 'scrambled_type' };
pcats = { 'drug', 'trial_type' };

if ( params.per_monkey )
  fcats{end+1} = 'monkey';
  pcats = union( pcats, 'monkey' );
  subdir = sprintf( '%s_per_monkey', subdir );
end

make_line_figures( pl, subset_data, subset_labels, fcats, gcats, pcats, subdir, params );

end

function [data, labels, time] = make_subsets(data, labels, time, mask, params)

t_ind = time >= params.time_limits(1) & time <= params.time_limits(2);

data = data(mask, t_ind);
labels = prune( labels(mask) );
time = time(t_ind);

end

function make_line_figures(pl, subset_data, subset_labels, fcats, gcats, pcats, subdir, params)

fig_I = findall_or_one( subset_labels, fcats );

for i = 1:numel(fig_I)
  fig = gcf();
  clf( fig );
  
  fig_data = subset_data(fig_I{i}, :);
  fig_labels = prune( subset_labels(fig_I{i}) );
  
  [axs, ~, inds] = pl.lines( fig_data, fig_labels, gcats, pcats );
  can_compare = all( cellfun(@numel, inds) == 2 );
  
  if ( can_compare )
    [ps, plot_hs] = dsp3.compare_series( axs, inds, fig_data, @signrank ...
      , 'x', pl.x ...
      , 'fig', gcf ...
      , 'p_levels', 0.05 ...
    );
  
    non_empty_plot_hs = cellfun( @(x) x(cellfun(@(y) ~isempty(y), x)), plot_hs, 'un', 0 );
    cellfun( @(x) cellfun(@(x) set(x, 'markersize', 0.5), x), non_empty_plot_hs );
  end

  if ( params.do_save )
    shared_utils.plot.fullscreen( fig );
    save_p = hwwa.approach_avoid_data_path( params, 'plots', 'over_time', subdir );
    dsp3.req_savefig( fig, save_p, fig_labels, [fcats, gcats, pcats], params.prefix );
  end
end

end

function mask = get_base_mask(labels, mask_func)

mask = hwwa.get_approach_avoid_mask( labels );
mask = mask_func( labels, mask );

end