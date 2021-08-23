function trial_data_file = trial_data(files)

hwwa.validatefiles( files, 'unified' );

unified_file = shared_utils.general.get( files, 'unified' );

trial_data_file = struct();
trial_data_file.unified_filename = unified_file.unified_filename;
trial_data_file.trial_data = unified_file.DATA;

end