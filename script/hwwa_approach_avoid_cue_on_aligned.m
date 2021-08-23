function outs = hwwa_approach_avoid_cue_on_aligned(varargin)

outs = hwwa_load_edf_aligned( ...
    'start_event_name', 'go_nogo_cue_onset' ...
  , 'look_back', 0 ...
  , 'look_ahead', 1e3 ...
  , 'is_parallel', true ...
  , 'files_containing', hwwa.approach_avoid_files() ...
  , 'error_handler', 'error' ...
  , varargin{:} ...
);

end