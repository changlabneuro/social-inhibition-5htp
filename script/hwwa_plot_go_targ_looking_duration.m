conf = hwwa.config.load();

use_files = {'13-F', '14-F'};

out = hwwa_check_go_targ_looking_duration( ...
    'files_containing', use_files ...
  , 'is_parallel',      false ...
  , 'config',           conf ...
  , 'is_parallel',      false ...
  , 'start_event',      'go_target_onset' ...
  , 'stop_event',       'iti' ...
  , 'use_reward_offset_time', false ...
);

plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', date_dir, 'lookdur' );

%%

lookdur = out.looking_duration;
labels = out.labels';

%%

plot_subdir = '';
base_prefix = '';

do_save = true;

pltdat = lookdur;
pltlabs = labels';

pl = plotlabeled.make_common();
pl.y_lims = [0, 700];
pl.x_tick_rotation = 0;

mask = fcat.mask( pltlabs ...
  , @find, 'no_errors' ...
  , @findnone, 'ephron' ...
  , @find, 'go_trial' ...
);

fcats = { 'monkey' };
xcats = { 'gender' };
gcats = { 'target_image_category' };
pcats = { 'monkey', 'trial_type', 'day' };

pltlabs = pltlabs(mask);
pltdat = pltdat(mask);

[figs, axs, I] = pl.figures( @bar, pltdat, pltlabs, fcats, xcats, gcats, pcats );

shared_utils.plot.ylabel( axs, 'Looking duration (ms)' );

if ( do_save )
  for i = 1:numel(figs)
    dsp3.req_savefig( figs(i), fullfile(plot_p, plot_subdir) ...
      , prune(pltlabs(I{i})), unique(cshorzcat(fcats, pcats, xcats, gcats)) ...
      , sprintf('lookdur_%s', base_prefix) );
  end
end