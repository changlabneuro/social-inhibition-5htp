function outs = hwwa_overall_p_correct(labels, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.pcorr_each = { 'unified_filename', 'monkey', 'drug', 'trial_type', 'scrambled_type' };
defaults.anova_factors = [];
defaults.anova_each = [];
defaults.norm_anova_factors = [];
defaults.norm_anova_each = [];
defaults.norm_t_each = [];
defaults.norm_t_elements = [];
defaults.per_target_image_category = false;
defaults.per_scrambled_type = false;
defaults.per_drug = true;
defaults.include_raw = true;
defaults.include_normalized = true;
defaults.per_trial_type = false;
defaults.errorbar_cats = {{}, {}, {}};
defaults.points_are = {};
defaults.y_label = '';
defaults.marker_size = [];
params = hwwa.parsestruct( defaults, varargin );

outs = struct();

mask = get_base_mask( labels, params.mask_func );

pcorr_each = get_pcorr_each( params );
[pcorr, pcorr_labels] = hwwa.percent_correct( labels', pcorr_each, mask );

norm_each = setdiff( pcorr_each, {'unified_filename'} );
[norm_corr, norm_labels] = hwwa.saline_normalize( pcorr, pcorr_labels', norm_each );

%%  raw

e_cats = params.errorbar_cats;

if ( params.include_raw )
  % anovas_each = { 'trial_type' };
  anovas_each = {};
  anova_factors = setdiff( pcorr_each, {'unified_filename'} );
  anova_factors = setdiff( anova_factors, anovas_each );

  if ( ~isempty(params.anova_factors) && iscell(params.anova_factors) )
    anova_factors = params.anova_factors;
  end
  if ( ~isempty(params.anova_each) && iscell(params.anova_each) )
    anovas_each = params.anova_each;
  end
  
%   e_xcats = { 'prev_correct_sequence' };
%   e_gcats = { 'drug' };
%   e_pcats = { 'trial_type' };
%   
%   run_bar( pcorr, pcorr_labels', e_cats{1}, e_cats{2}, e_cats{3}, 'raw', params );
%   run_errorbar( pcorr, pcorr_labels', e_cats{1}, e_cats{2}, e_cats{3}, 'raw', params );
% 
%   v_gcats = { 'scrambled_type', 'target_image_category' };
%   v_pcats = { 'trial_type' };
%   run_violin( pcorr, pcorr_labels', v_gcats, v_pcats, 'raw', params );

%   run_anova( pcorr, pcorr_labels', anovas_each, anova_factors, 'raw', params );
  
  outs.raw_corr = pcorr;
  outs.raw_labels = pcorr_labels';
end

%%  normalized

if ( params.include_normalized )
  norm_subdir = 'normalized';

  norm_anovas_each = {};
  norm_anova_factors = setdiff( norm_each, {'unified_filename', 'drug', 'monkey'} );
  norm_anova_factors = setdiff( norm_anova_factors, norm_anovas_each );
  
  if ( ~isempty(params.norm_anova_factors) && iscell(params.norm_anova_factors) )
    norm_anova_factors = params.norm_anova_factors;
  end
  if ( ~isempty(params.norm_anova_each) && iscell(params.norm_anova_each) )
    norm_anovas_each = params.norm_anova_each;
  end
  
  run_errorbar( norm_corr, norm_labels', e_cats{1}, e_cats{2}, e_cats{3}, norm_subdir, params );

%   run_anova( norm_corr, norm_labels', norm_anovas_each, norm_anova_factors, norm_subdir, params );

  norm_ts_each = { 'trial_type' };
  norm_ta = { 'scrambled' };
  norm_tb = { 'not-scrambled' };
  
  if ( ~isempty(params.norm_t_each) )
    norm_ts_each = params.norm_t_each;
    norm_ta = params.norm_t_elements{1};
    norm_tb = params.norm_t_elements{2};
  end
  
%   run_ttests( norm_corr, norm_labels', norm_ts_each, norm_ta, norm_tb, norm_subdir, params );

%   norm_sr_each = { 'trial_type', 'scrambled_type' };
%   run_signrank1( norm_corr, norm_labels', norm_sr_each, norm_subdir, params );
  
  outs.normalized_pcorr = norm_corr;
  outs.normalized_labels = norm_labels';
end

end

function run_signrank1(data, labels, each, subdir, params)

sr_outs = dsp3.signrank1( data, labels', each ...
  , 'signrank_inputs', {1} ...
);

maybe_save_sr_outs( sr_outs, each, subdir, params );

end

function run_ttests(data, labels, each, a, b, subdir, params)

ttest_outs = dsp3.ttest2( data, labels', each, a, b );
maybe_save_ttest_outs( ttest_outs, each, subdir, params );

end

function run_anova(data, labels, each, factors, subdir, params)

anova_outs = dsp3.anovan( data, labels', each, factors ...
  , 'dimension', 1:numel(factors) ...
  , 'remove_nonsignificant_comparisons', false ...
  , 'include_significant_factor_descriptives', true ...
);

maybe_save_anova_outs( anova_outs, [each, factors], subdir, params );

end

function run_errorbar(data, labels, xcats, gcats, pcats, subdir, params)

pl = plotlabeled.make_common();
pl.x_order = { 'saline', '5-htp' };

if ( ~isempty(params.points_are) )
  pl.points_are = params.points_are;
  pl.add_points = true;
end
if ( ~isempty(params.marker_size) )
  pl.marker_size = params.marker_size;
end

axs = pl.errorbar( data, labels, xcats, gcats, pcats );

if ( ~isempty(params.y_label) )
  ylabel( axs(1), params.y_label );
end

if ( params.do_save )
  save_plot( gcf, axs, labels, [xcats, gcats, pcats], fullfile(subdir, 'errorbar'), params );
end

end

function run_bar(data, labels, xcats, gcats, pcats, subdir, params)

pl = plotlabeled.make_common();
pl.x_order = { 'saline', '5-htp' };

axs = pl.bar( data, labels, xcats, gcats, pcats );

if ( params.do_save )
  save_plot( gcf, axs, labels, [xcats, gcats, pcats], fullfile(subdir, 'bar'), params );
end

end

function run_violin(data, labels, gcats, pcats, subdir, params)

%%

pl = plotlabeled.make_common();

axs = pl.violinalt( data, labels, gcats, pcats );

if ( params.do_save )  
  shared_utils.plot.fullscreen( gcf );
  shared_utils.plot.match_ylims( axs );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pcorr', subdir, 'violin' );
  dsp3.req_savefig( gcf, save_p, labels, [gcats, pcats], params.prefix );
end

end

function save_plot(fig, axs, labels, cats, subdir, params)

shared_utils.plot.fullscreen( fig );
shared_utils.plot.match_ylims( axs );
save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pcorr', subdir );
dsp3.req_savefig( fig, save_p, labels, cats, params.prefix );

end

function maybe_save_ttest_outs(ttest_outs, spec, subdir, params)

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'analyses', 'pcorr', subdir );
  dsp3.save_ttest2_outputs( ttest_outs, save_p, spec );
end

end

function maybe_save_anova_outs(anova_outs, spec, subdir, params)

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pcorr', fullfile(subdir, 'stats') );
  dsp3.save_anova_outputs( anova_outs, save_p, spec );
end

end

function maybe_save_sr_outs(sr_outs, spec, subdir, params)

if ( params.do_save )
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pcorr', fullfile(subdir, 'stats') );
  dsp3.save_signrank1_outputs( sr_outs, save_p, spec );
end

end

function each = get_pcorr_each(params)

each = params.pcorr_each;

if ( params.per_target_image_category )
  each{end+1} = 'target_image_category';
end
if ( ~params.per_drug )
  each = setdiff( each, {'drug'} );
end
if ( params.per_trial_type )
  each{end+1} = 'trial_type';
end
if ( params.per_scrambled_type )
  each{end+1} = 'scrambled_type';
end

end

function each = percent_correct_each()

each = { 'unified_filename', 'monkey', 'drug', 'trial_type', 'scrambled_type' };

end

function mask = get_base_mask(labels, mask_func)

mask = mask_func( labels, hwwa.get_approach_avoid_mask(labels) );

end