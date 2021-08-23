function hwwa_plot_drug_behav(rt, labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.norm_func = @no_norm;
defaults.collapse_scrambled_image_category = false;

params = hwwa.parsestruct( defaults, varargin );

hwwa.decompose_social_scrambled( labels );

if ( params.collapse_scrambled_image_category )
  collapse_scrambled_image_category( labels );
end

[rt, labels] = ensure_correct_subset( rt, labels' );
handle_label_repeats( labels );
hwwa.label_prev_correct( labels );

% Plot overall bars
run_plot_bar_percent_correct( rt, labels, params );
run_plot_bar_rt( rt, labels, params );

% Plot lines for each day
% run_plot_percent_correct_per_day( rt, labels, params );
% run_plot_rt_per_day( rt, labels, params );

end

function collapse_scrambled_image_category(labels)

scrambled_ind = find( labels, 'scrambled' );
setcat( labels, 'target_image_category', 'scrambled-image', scrambled_ind );
collapsecat( labels, 'scrambled_type' );

end

function run_plot_percent_correct_per_day(data, labels, params)

spec = get_specificities();

for i = 1:numel(spec)
  plot_lines_per_day( data, labels', spec{i} ...
    , 'percent_correct_per_day', @get_p_correct_data, params );
end

end

function run_plot_rt_per_day(data, labels, params)

spec = get_specificities();

for i = 1:numel(spec)
  plot_lines_per_day( data, labels', spec{i} ...
    , 'rt_per_day', @get_rt_data, params );
end

end

function run_plot_bar_percent_correct(data, labels, params)

bar_specificities = get_specificities();

for i = 1:numel(bar_specificities)
  plot_bar_overall( data, labels', bar_specificities{i} ...
    , 'percent_correct', @get_p_correct_data, params );
end

end

function run_plot_bar_rt(data, labels, params)

bar_specificities = get_specificities();

for i = 1:numel(bar_specificities)
  plot_bar_overall( data, labels', bar_specificities{i} ...
    , 'rt', @get_rt_data, params );
end

end

function [rt, rt_labels] = get_rt_data(data, labels, rt_each, mask)

% Only correct go-trials
rt_mask = find( labels, {'go_trial', 'correct_true'}, mask );

% Calculate average rt for each `rt_each`
[rt_labels, I] = keepeach( labels', rt_each, rt_mask );
rt = rowop( data, I, @(x) nanmean(x, 1) );

end

function [pcorr, pcorr_labels] = get_p_correct_data(data, labels, data_each, mask)

% Calculate percent correct for each `p_each`
[pcorr, pcorr_labels] = hwwa.percent_correct( labels', data_each, mask );

end

function plot_lines_per_day(data, labels, specificity, measure_type, data_func, params)

mask = get_base_mask( labels );

handle_specificity( labels, specificity );

data_each = get_data_each( specificity );

[plt_data, plt_labels] = data_func( data, labels, data_each, mask );

pl = plotlabeled.make_common();
pl.add_errors = false;
pl.per_panel_labels = true;
pl.x_order = { 'appetitive', 'threat', 'scrambled_appetitive' };

fig_cats = { 'monkey' };
xcats = { 'target_image_category' };
gcats = { 'day' };
pcats = { 'drug', 'trial_type', 'monkey', 'trial_type' };

if ( ismember('gender', specificity) )
  xcats{end+1} = 'gender';
end

[figs, axs, fig_I] = pl.figures( @errorbar, plt_data, plt_labels, fig_cats, xcats, gcats, pcats );
shared_utils.plot.match_ylims( axs );

if ( params.do_save )
  base_p = fullfile( get_base_plot_path(params), measure_type );
  save_p = fullfile( base_p, make_specificity_subdir(specificity) );
  
  filename_cats = [ fig_cats, pcats ];
  
  save_figs( save_p, figs, fig_I, plt_labels, filename_cats );
end


end

function data_each = get_data_each(specificity)

data_each = { 'date', 'target_image_category', 'trial_type' };

check_additional = { 'gender', 'scrambled_type', 'switch_trial_type', 'switch_target_image_category' };
use_additional = ismember( check_additional, specificity );

data_each = [ data_each, check_additional(use_additional) ];

end

function mask = handle_repeat_label_mask(labels, specificity, mask)

% Only include trials where the preceding trial was correct.

mask = fcat.mask( labels, mask ...
  , @find, 'prev_correct_true' ...
);

end

function [data, labels] = no_norm(data, labels, spec)
end

function [data, labels] = social_minus_scrambled(data, labels, spec)

[data, labels] = dsp3.summary_binary_op( data, labels, spec ...
  , 'scrambled', 'not-scrambled', @minus, @nanmean );

end

function plot_bar_overall(data, labels, specificity, measure_type, data_func, params)
%%

mask = get_base_mask( labels );
mask = handle_repeat_label_mask( labels, specificity, mask );

handle_specificity( labels, specificity );
norm_spec = csunion( specificity, 'trial_type' );

data_each = get_data_each( specificity );

[plt_data, plt_labels] = data_func( data, labels, data_each, mask );
[plt_data, plt_labels] = feval( params.norm_func, plt_data, plt_labels, norm_spec );
% Appetitie - threat
% [plt_data, plt_labels] = social_minus_scrambled( plt_data, plt_labels ...
%   , setdiff(data_each, 'scrambled_type') );

pl = plotlabeled.make_common();
pl.prefer_multiple_groups = true;
pl.prefer_multiple_xs = true;
pl.x_tick_rotation = 0;
pl.panel_order = { 'go_trial', 'nogo_trial' };
pl.y_lims = [0.3, 1.2];

fig_cats = { 'monkey' };
xcats = { 'target_image_category' };
gcats = { 'drug' };
pcats = { 'scrambled_type', 'trial_type', 'monkey' };

% plt_func = @bar;
plt_func = @boxplot;

if ( ismember('gender', specificity) )
  xcats{end+1} = 'gender';
end

if ( ismember(switch_trial_cat, specificity) )  
  if ( ~ismember('target_image_category', specificity) )
    xcats = switch_trial_cat;
    pcats = setdiff( pcats, 'trial_type' );
    
  else
    pcats = setdiff( pcats, {'trial_type', 'scrambled_type'} );
    pcats{end+1} = switch_trial_cat;
    gcats{end+1} = 'scrambled_type';
  end
  
  plt_func = @errorbar;
  
  plt_data = remove_non_labeled_switch( plt_data, plt_labels, switch_trial_cat );
  
elseif ( ismember(switch_image_cat, specificity) )
  xcats = switch_image_cat;
  pcats = setdiff( pcats, 'scrambled_type' );
  
  plt_data = remove_non_labeled_switch( plt_data, plt_labels, switch_image_cat );
  
  plt_func = @errorbar;
end

if ( isequal(plt_func, @bar) || isequal(plt_func, @errorbar) )
  [figs, axs, fig_I] = pl.figures( plt_func, plt_data, plt_labels, fig_cats, xcats, gcats, pcats );
else
  pl.prefer_multiple_xs = false;
  pl.prefer_multiple_groups = true;
  
  gcats = csunion( gcats, xcats );
  
  if ( ismember('scrambled_type', pcats) )
    pcats = setdiff( pcats, 'scrambled_type' );
    gcats{end+1} = 'scrambled_type';
  end
  
  [figs, axs, fig_I] = pl.figures( plt_func, plt_data, plt_labels, fig_cats, gcats, pcats );
end

if ( params.do_save )
  base_p = fullfile( get_base_plot_path(params), measure_type );
  save_p = fullfile( base_p, make_specificity_subdir(specificity) );
  
  filename_cats = [ fig_cats, pcats ];
  
  save_figs( save_p, figs, fig_I, plt_labels, filename_cats );
end

%%

% factors = { 'target_image_category', 'trial_type', 'scrambled_type' };
% anova_outs = dsp3.anovan( plt_data, plt_labels', fig_cats, factors ...
%   , 'include_per_factor_descriptives', true ...
%   , 'remove_nonsignificant_comparisons', false ...
% );
% 
% a = 'not-scrambled';
% b = 'scrambled';
% ranksum_each = { 'trial_type' };
% 
% rs_outs = dsp3.ranksum( plt_data, plt_labels', ranksum_each, a, b );

%%

% factors = { 'target_image_category' };
% anova_outs = dsp3.anovan( plt_data, plt_labels', 'trial_type', factors ...
%   , 'include_per_factor_descriptives', true ...
%   , 'remove_nonsignificant_comparisons', false ...
%   , 'dimension', 1 ...
% );
% 
% a = 'appetitive';
% b = 'threat';
% ranksum_each = { 'trial_type' };
% 
% rs_outs = dsp3.ranksum( plt_data, plt_labels', ranksum_each, a, b );

end

function [data, labels] = remove_non_labeled_switch(data, labels, switch_cat)

missing_lab = makecollapsed( labels, switch_cat );
data = indexpair( data, labels, findnone(labels, missing_lab) );

end

function [rt, labels] = ensure_correct_subset(rt, labels)

assert_ispair( rt, labels );

days_5htp = hwwa.get_image_5htp_days();
days_sal = hwwa.get_image_saline_days();

use_days = unique( horzcat(days_5htp, days_sal) );

rt = indexpair( rt, labels, findor(labels, use_days) );
prune( labels );

end

function handle_label_repeats(labels)

hwwa.label_gonogo_switchiness( labels );
hwwa.label_image_category_switchiness( labels );

end

function handle_specificity(labels, specificity)

if ( ~ismember('target_image_category', specificity) )
  collapsecat( labels, 'target_image_category' );
end

if ( ~ismember('monkey', specificity) )
  collapsecat( labels, 'monkey' );
end

end

function save_figs(save_p, figs, fig_I, labels, filename_cats)

for i = 1:numel(figs)
  f = figs(i);
  ind = fig_I{i};

  shared_utils.plot.fullscreen( f );

  dsp3.req_savefig( f, save_p, prune(labels(ind)), filename_cats );
end

end

function p = get_base_plot_path(params)

p = fullfile( hwwa.dataroot(params.config), 'plots', 'behavior' ...
  , hwwa.datedir, 'basic_behavior', params.base_subdir );

end

function mask = get_ephron_mask(labels)

days = {'041019', '041219', '041619', '041819', '042319', '042519' ...
  , '040819', '041119', '041519', '041719', '042219', '042419' ...
  , '042919', '043019', '050219', '050319' ...
};

mask = fcat.mask( labels ...
  , @find, days ...
  , @find, 'ephron' ...
);

end

function mask = get_base_mask(labels)

mask = fcat.mask( labels ...
  , @find, {'initiated_true'} ...
  , @findnot, {'tar', '011619'} ...
  , @findnone, 'ephron' ...
  , @findnone, '021819' ... % Eliminate first saline day to match N.
);

mask = union( mask, get_ephron_mask(labels) );

end

function subdir = make_specificity_subdir(specificity)

subdir = '';

for i = 1:numel(specificity)
  format_str = ternary( i < numel(specificity), '%sper_%s__', '%sper_%s' );
  subdir = sprintf( format_str, subdir, specificity{i} );
end

if ( isempty(subdir) )
  subdir = 'across_all';
end

end

function specs = get_specificities()

base_spec = get_base_specificities();
switch_trial_type_spec = get_switch_trial_type_specificities();
switch_image_cat_spec = get_switch_image_category_specificities();

% specs = [ base_spec, switch_trial_type, switch_image_cat ];
% specs = [ base_spec, switch_image_cat_spec, switch_trial_type_spec ];
specs = base_spec;

end

function s = switch_trial_cat

s = 'switch_trial_type';

end

function s = switch_image_cat

s = 'switch_target_image_category';

end

function specs = get_switch_trial_type_specificities(switch_cat)

switch_cat = switch_trial_cat;

specs = { ...
    {switch_cat} ...
  , {switch_cat, 'monkey'} ...
  , {switch_cat, 'scrambled_type'} ...
  , {switch_cat, 'target_image_category'} ...
  , {switch_cat, 'target_image_category', 'scrambled_type'} ...
  , {switch_cat, 'target_image_category', 'monkey'} ...
  , {switch_cat, 'target_image_category', 'monkey', 'scrambled_type'} ...
};

end

function specs = get_switch_image_category_specificities()

switch_cat = switch_image_cat;

specs = { ...
    {switch_cat} ...
  , {switch_cat, 'monkey'} ...
  , {switch_cat, 'scrambled_type'} ...
  , {switch_cat, 'scrambled_type', 'monkey'} ...
};

end

function specs = get_base_specificities()

% specs = { ...
%     {'target_image_category', 'monkey', 'gender'} ...
%   , {'target_image_category', 'gender'} ...
%   , {'target_image_category', 'monkey'} ...
%   , {'target_image_category', 'scrambled_type'} ...
%   , {'target_image_category', 'scrambled_type', 'monkey'} ...
%   , {'scrambled_type', 'monkey'} ...
%   , {'scrambled_type'} ...
%   , {'target_image_category'} ...
% };

% specs = { ...
%     {'target_image_category', 'scrambled_type'} ...
%   , {'target_image_category', 'scrambled_type', 'monkey'} ...
%   , {'scrambled_type'} ...
% };

specs = { ...
  {'scrambled_type'} ...
};

end