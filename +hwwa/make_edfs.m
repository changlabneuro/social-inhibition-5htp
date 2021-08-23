function [results, params] = make_edfs(varargin)

%   MAKE_EDFS -- Make multiple edf intermediate files.
%
%     See also hwwa.make.edfs, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = 'unified';
output = 'edf';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

% Save as -v7.3
loop_runner.save_func = @(pathstr, var) save( pathstr, 'var', '-v7.3' );

results = loop_runner.run( @hwwa.make.edfs, params.config );

end