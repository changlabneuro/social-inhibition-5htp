function [results, params] = make_edf_trials(varargin)

%   MAKE_EDF_TRIALS -- Make multiple edf-trial intermediate files.
%
%     See also hwwa.make.edf_trials, hwwa.get_common_make_defaults

defaults = hwwa.make.defaults.edf_trials();

params = hwwa.parsestruct( defaults, varargin );

event_subdir = params.event_subdir;
event_subdir = validatestring( event_subdir, {'el_events', 'edf_events'} );

params.event_subdir = event_subdir;

inputs = { 'edf', event_subdir };
output = params.output_directory;

[~, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.edf_trials, params );

end