function [result, params] = make_alternate_el_events(varargin)

%   MAKE_ALTERNATE_EL_EVENTS -- Make multiple alternate eyelink event 
%     intermediate files.
%
%     See also hwwa.make.alternate_el_events, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = { 'events', 'edf' };
output = 'edf_events';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

result = loop_runner.run( @hwwa.make.alternate_el_events );

end