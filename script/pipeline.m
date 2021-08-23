%%  pipeline

conf = hwwa.config.load();
conf.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';

files = { 'test_' };

inputs = { 'overwrite', false, 'files_containing', files, 'config', conf };

%%  save unified trial data

hwwa.make_unified( inputs{:} );

%%  edf -> mat

hwwa.make_edfs( inputs{:} );
hwwa.make_alternate_el_events( inputs{:} );

%%  save spikes + lfp as mat files

hwwa.make_sync_times( inputs{:} );
% hwwa.make_spikes( inputs{:} );
hwwa.make_ms_spikes( inputs{:} );

hwwa.make_lfp( inputs{:} ...
  , 'kind', 'wb' ...
);

%%  save mua spikes

hwwa.make_mua_spikes( inputs{:} );

%%  make trial labels

hwwa.make_labels( inputs{:} );

%%  extract events, convert to eyelink and plex time

hwwa.make_events( inputs{:} );
hwwa.make_el_events( inputs{:} );
hwwa.make_plex_events( inputs{:} );

%%  align pupil size, position to trial events

hwwa.make_edf_trials( inputs{:} ...
  , 'event',      {'iti', 'go_nogo_cue_onset'} ...
  , 'look_back',  [0, -300] ...
  , 'look_ahead', [1e3, 0] ...
);

%%  calculate proportion of samples in bounds, during iti

hwwa.make_iti_bounds( inputs{:} ...
  , 'bin', 1 ...
);

%%  calculate psth, raster aligned to trial events

evts = { 'go_nogo_cue_onset', 'go_target_onset' ...
  , 'reward_onset', 'go_target_offset', 'go_target_acquired', 'fixation_on' };

% evts = { 'fixation_on' };

look_back = repmat( -0.5, 1, numel(evts) );
look_ahead = repmat( 1, 1, numel(evts) );

hwwa.make_psth( inputs{:} ...
  , 'kind',       'sua' ...
  , 'overwrite',  true ...
  , 'append',     false ...
  , 'event',      evts ...
  , 'look_back',  look_back ...
  , 'look_ahead', look_ahead ...
);

%%  align lfp

hwwa.make_aligned_lfp( inputs{:} ...
  , 'kind',       'fp' ...
  , 'event',      evts ...
  , 'look_back',  look_back ...
  , 'look_ahead', look_ahead ...
  , 'window_size', 0.5 ...
);

%%  create raw power

hwwa.make_raw_power( inputs{:} ...
  , 'filter', true ...
);