function [results, params] = make_rois(varargin)

%   MAKE_ROIS -- Make multiple intermediate roi files.
%
%     See also hwwa.make.rois, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = 'unified';
output = 'rois';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.rois );

end