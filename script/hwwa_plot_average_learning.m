function hwwa_plot_average_learning(rt, labels, varargin)

defaults = hwwa.get_common_make_defaults();
defaults.colored_lines_are = 'social_vs_scrambled';
defaults.do_save = false;
defaults.base_subdir = '';
defaults.base_prefix = '';
defaults.is_per_monkey = false;
defaults.is_per_drug = false;
defaults.is_rt = false;

params = hwwa.parsestruct( defaults, varargin );

conf = params.config;

plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior' ...
  , datestr(now, 'mmddyy'), 'average_learning' );
plot_p = fullfile( plot_p, ternary(params.is_rt, 'rt', 'percent_correct') );

params.plot_p = plot_p;

plot_average( rt, labels', params )

end

function make_numbered_days(labs, I)

for i = 1:numel(I)
  days = combs( labs, 'day', I{i} );
  
  day_nums = datenum( days, 'mmddyy' );
  
  [~, sorted_I] = sort( day_nums );
  sorted_days = days(sorted_I);
  
  for j = 1:numel(sorted_days)
    replace( labs, days{j}, sprintf('day_%d', j) );
  end
  
  prune( labs );
end

end

function plot_average(rt, labels, params)

%%

switch ( params.colored_lines_are )
  case 'social_vs_scrambled'
    gcats = { 'target_image_category' };
    pcats = { 'scrambled_type' };
    
  case 'threat_vs_appetitive'    
    gcats = { 'scrambled_type' };
    pcats = { 'target_image_category' };
    
  otherwise
    error( 'Unrecognized line type: "%s".', params.colored_lines_are );
end

mask = fcat.mask( labels ...
  , @find, {'pilot_social', 'initiated_true', 'no_errors', 'wrong_go_nogo'} ...
  , @findnot, {'tar', '011619'} ...
);

if ( params.is_per_drug )
  pcats{end+1} = 'drug';
end

spec = { 'date', 'trial_type', 'target_image_category' };

if ( params.is_rt )
  pcorr_dat = rt;
  pcorr_labs = labels';
else
  [pcorr_dat, pcorr_labs] = hwwa.percent_correct( labels, spec, mask );
end

hwwa.decompose_social_scrambled( pcorr_labs );

pl = plotlabeled.make_common();

if ( ~params.is_rt )
  pl.y_lims = [0, 1];
end

if ( ~params.is_per_monkey )
  make_numbered_days( pcorr_labs, findall_or_one(pcorr_labs, 'monkey') );
  collapsecat( pcorr_labs, 'monkey' );
end

fcats = { 'monkey' };
xcats = { 'day' };
pcats = cshorzcat( pcats, {'trial_type', 'monkey'} );

[figs, axs, fig_I] = pl.figures( @errorbar, pcorr_dat, pcorr_labs, fcats, xcats, gcats, pcats );

if ( params.do_save )
  plot_p = fullfile( params.plot_p, params.base_subdir, params.colored_lines_are );
  
  for i = 1:numel(figs)
    dsp3.req_savefig( figs(i), plot_p, pcorr_labs(fig_I{i}) ...
      , unique(cshorzcat(pcats, fcats)) ...
      , sprintf('percent_correct_%s', params.base_prefix) );
  end
end

end