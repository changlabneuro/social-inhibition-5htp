function [results, params] = make_trial_data(varargin)

%   MAKE_LABELS -- Make multiple intermediate trial_data files.
%
%     See also hwwa.make.trial_data, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = 'unified';
output = 'trial_data';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.trial_data );

end