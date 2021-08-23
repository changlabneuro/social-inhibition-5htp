function hwwa_plot_pupil_traces(aligned_outs, varargin)

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
defaults.mask_func = @hwwa.default_mask_func;
defaults.time_normalize = true;
defaults.time_normalize_range = [-150, 0];
defaults.time_limits = [];
defaults.smooth_func = @(x) x;
defaults.compare_series = true;
defaults.fcats = {};
defaults.prefix = '';
defaults.formats = { 'epsc', 'png', 'fig', 'svg' };
params = hwwa.parsestruct( defaults, varargin );

t = aligned_outs.time(1, :);
pupil = aligned_outs.pupil;
labels = aligned_outs.labels';
assert_ispair( pupil, labels );

if ( params.time_normalize )
  pupil = time_normalize( pupil, t, params.time_normalize_range );
end

base_mask = get_base_mask( labels, params.mask_func );
corr_mask = hwwa.find_correct( labels, base_mask );
base_subdir = get_base_subdir( aligned_outs, params );

plot_soc_nonsoc( pupil, t, labels', corr_mask, base_subdir, params );
plot_5htp_saline( pupil, t, labels', corr_mask, base_subdir, params );
plot_corr_incorr( pupil, t, labels', base_mask, base_subdir, params );

end

function base_subdir = get_base_subdir(aligned_outs, params)

if ( shared_utils.struct.is_field(aligned_outs, 'params.start_event_name') )
  base_subdir = aligned_outs.params.start_event_name;
else
  base_subdir = '';
end

end

function pup = time_normalize(pup, t, trange)

baseline = nanmean( pup(:, t >= trange(1) & t <= trange(2)), 2 );
pup = pup ./ baseline;

end

function plot_corr_incorr(pupil, t, labels, mask, subdir, params)

pcats = { 'trial_type', 'drug' };
gcats = { 'correct' };

subdir = fullfile( subdir, 'compare-correct' );

plot_traces( pupil, t, labels, pcats, gcats, mask, subdir, params );

end

function plot_5htp_saline(pupil, t, labels, mask, subdir, params)

pcats = { 'trial_type', 'scrambled_type' };
gcats = { 'drug' };

subdir = fullfile( subdir, 'compare-drug' );

plot_traces( pupil, t, labels, pcats, gcats, mask, subdir, params );

end

function plot_soc_nonsoc(pupil, t, labels, mask, subdir, params)

pcats = { 'trial_type', 'drug' };
gcats = { 'scrambled_type' };

subdir = fullfile( subdir, 'compare-scrambled-type' );

plot_traces( pupil, t, labels, pcats, gcats, mask, subdir, params );

end

function plot_traces(pupil, t, labels, pcats, gcats, mask, subdir, params)

assert_ispair( pupil, labels );

fcats = params.fcats;
pcats = unique( [pcats, fcats] );

fig_I = findall_or_one( labels, fcats, mask );

for i = 1:numel(fig_I)
  pl = plotlabeled.make_common();
  pl.fig = gcf();
  pl.x = t;
  pl.add_smoothing = true;
  pl.smooth_func = params.smooth_func;
  
  set( gcf, 'Renderer', 'painters' );
  
  pup = pupil(fig_I{i}, :);
  labs = prune( labels(fig_I{i}) );
  
  [axs, ~, inds] = pl.lines( pup, labs, gcats, pcats );
  
  if ( ~isempty(params.time_limits) )
    shared_utils.plot.set_xlims( axs, params.time_limits );
  end
  
  if ( params.compare_series )
    test_func = @extract_ttest2_p;
%     test_func = @signrank;
    
    [ps, plot_hs] = dsp3.compare_series( axs, inds, pup, test_func ...
      , 'x', pl.x ...
      , 'fig', gcf ...
      , 'p_levels', 0.05 ...
    );
  end
  
  maybe_save_fig( gcf, labs, [fcats, pcats, gcats], subdir, params );
end

end

function p = extract_ttest2_p(varargin)

[~, p] = ttest2( varargin{:} );

end

function maybe_save_fig(fig, labels, spec, subdir, params)

if ( params.do_save )
  shared_utils.plot.fullscreen( fig );
  save_p = hwwa.approach_avoid_data_path( params, 'plots', 'pupil', subdir );
  dsp3.req_savefig( fig, save_p, labels, spec, params.prefix, params.formats );
end

end

function mask = get_base_mask(labels, mask_func)

require_initiated = true;

mask = hwwa.get_approach_avoid_mask( labels, {}, require_initiated );
mask = mask_func( labels, mask );

end