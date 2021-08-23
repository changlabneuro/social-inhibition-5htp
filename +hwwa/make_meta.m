function [results, params] = make_meta(varargin)

%   MAKE_META -- Make multiple intermediate meta files.
%
%     See also hwwa.make.trial_data, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = 'unified';
output = 'meta';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.meta );

end