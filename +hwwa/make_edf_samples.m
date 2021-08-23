function [results, params] = make_edf_samples(varargin)

%   MAKE_EDF_SAMPLES -- Make multiple intermediate edf_samples files.
%
%     See also hwwa.make.trial_data, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = 'edf';
output = 'edf_samples';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.edf_samples );

end