function [results, params] = make_el_events(varargin)

%   MAKE_EL_EVENTS -- Make multiple eyelink-event intermediate files.
%
%     This function is not recommended, since the algorithm used to 
%     interpolate event times between Matlab and Eyelink clocks is less
%     robust than that used in hwwa.make_alternate_el_events.
%
%     See also hwwa.make_alternate_el_events,
%       hwwa.make.el_events, hwwa.get_common_make_defaults

defaults = hwwa.get_common_make_defaults();

inputs = { 'events', 'edf' };
output = 'el_events';

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, output, defaults, varargin );
loop_runner.func_name = mfilename;

results = loop_runner.run( @hwwa.make.el_events );

end