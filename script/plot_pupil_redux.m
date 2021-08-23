conf = hwwa.config.load();

conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

labels_p = hwwa.get_intermediate_dir( 'labels', conf );
samples_p = hwwa.get_intermediate_dir( 'edf_trials', conf );

sample_mats = hwwa.require_intermediate_mats( [], samples_p );

puplabs = fcat();
pupdat = [];
t = [];

do_norm = false;

for i = 1:numel(sample_mats)
  hwwa.progress( i, numel(sample_mats), mfilename );
  
  samples_file = shared_utils.io.fload( sample_mats{i} );
  
  un_filename = samples_file.unified_filename;
  
  try
    labels_file = shared_utils.io.fload( fullfile(labels_p, un_filename) );
  catch err
    warning( err.message );
    continue;
  end
  
  if ( ~isKey(samples_file.trials, 'go_nogo_cue_onset') )
    continue;
  end
  
  target = samples_file.trials('go_nogo_cue_onset').samples('pupilSize');
  baseline = samples_file.trials('iti').samples('pupilSize');
  
  t = samples_file.trials('iti').time;
  
  if ( do_norm )
    mean_base = nanmean( baseline, 2 );
    target = target ./ mean_base;
  end
  
  labs = labels_file.labels';
  
  hwwa.add_day_labels( labs );
  hwwa.add_data_set_labels( labs );
  hwwa.add_drug_labels( labs );
  
  prune( labs );
  
  append( puplabs, labs );
  
  pupdat = [ pupdat; target ];
end

plot_p = fullfile( conf.PATHS.data_root, 'plots', 'behavior', date_dir );

%%  time course

pltlabs = puplabs';
pltdat = pupdat;

mask = fcat.mask( pltlabs ...
  , @find, {'tarantino_drug_check', 'initiated_true' } ...
  , @find, 'saline' ...
);

gcats = { 'correct', 'trial_type' };
pcats = { 'drug', 'cue_delay' };

pl = plotlabeled.make_common();
pl.fig = figure(2);

pl.lines( pltdat(mask, :), pltlabs(mask), gcats, pcats );

%%  average

do_save = true;

pltlabs = puplabs';
pltdat = pupdat;

pltdat = nanmean( pltdat, 2 );

mask = fcat.mask( pltlabs ...
  , @find, {'tarantino_drug_check', 'initiated_true'} ...
  , @find, hwwa_get_first_days() ...
);

pl = plotlabeled.make_common();

pl.add_points = true;
pl.points_are = 'day';
pl.marker_size = 8;
pl.marker_type = '*';
pl.panel_order = { 'saline' };

xcats = { 'drug' };
gcats = {};
pcats = {};

% xcats = { 'trial_type' };
% gcats = {  };
% pcats = { 'drug' };

spec = unique( cshorzcat(xcats, gcats, pcats, {'date'}) );

[meanlabs, I] = keepeach( pltlabs', spec, mask );
meandat = rownanmean( pltdat, I );

pl.bar( meandat, meanlabs, xcats, gcats, pcats );

shared_utils.plot.fullscreen( gcf );

if ( do_save )
  dsp3.req_savefig( gcf, plot_p, meanlabs', pcats, 'pupil' );
end