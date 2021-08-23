function hwwa_run_plot_approach_avoid_behav(varargin)

provided_data = false;

if ( nargin >= 1 && ~ischar(varargin{1}) && ~isempty(varargin{1}) )
  behav_outputs = varargin{1};
  varargin(1) = [];
  provided_data = true;
end

defaults = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );

params = hwwa.parsestruct( defaults, varargin );

if ( ~provided_data )
  behav_outputs = hwwa_load_basic_behav_measures( ...
      'config', params.config ...
    , 'files_containing', get_target_files() ...
    , 'trial_bin_size', 50 ...
    , 'trial_step_size', 50 ...
    , 'is_parallel', true ...
  );
end

hwwa_plot_approach_avoid_behav( behav_outputs.rt, behav_outputs.labels', params );

end

function files = get_target_files()

files = [ hwwa.get_image_5htp_days(), hwwa.get_image_saline_days() ];
files = cellstr( hwwa.to_date(files) );

end