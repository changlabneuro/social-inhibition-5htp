inputs = hwwa.get_common_make_defaults();

% inputs.config.PATHS.data_root = '/Users/Nick/Desktop/test_data_root2';
inputs.config.PATHS.data_root = '/Volumes/My Passport/NICK/Chang Lab 2016/hww_gng/data';
inputs.config.PATHS.raw_subdirectory = 'raw_redux';

% use_files = arrayfun( @(x) sprintf('12%d', x), [ 11:14, 17, 19 ], 'un', 0 );
% use_files = {'08-', '09-' };
% use_files = { '29-J', '30-J', '31-J', '01-F' };
% use_files = { '24-F', '25-F', '26-F', '27-F' };
% use_files = { '18-M', '19-M', '20-M', '21-M', '22-M' };

% file_nums = [8, 10, 11, 12, 15, 16, 17, 18, 22, 23, 24, 25];
% use_files = arrayfun( @(x) sprintf('%02d-A', x), file_nums, 'un', 0 );

% inputs.files_containing = use_files;
inputs.skip_existing = true;

%%

hwwa.make_unified( inputs );
hwwa.make_edfs( inputs );
hwwa.make_events( inputs );
hwwa.make_el_events( inputs );
hwwa.make_alternate_el_events( inputs );
hwwa.make_labels( inputs );
hwwa.make_trial_data( inputs );
hwwa.make_meta( inputs );

%%

res = hwwa.make_edf_trials( inputs ...
  , 'event',        {'go_target_onset'} ...
  , 'event_subdir', 'edf_events' ...
  , 'look_back',    0 ...
  , 'look_ahead',   6000 ...
  , 'overwrite',    true ...
  , 'append',       true ...
  , 'keep_output',  false ...
  , 'skip_existing', false ...
);