function [x_deg, y_deg] = run_px2deg(x, y, flip_y, monitor_constants)

if ( nargin < 3 || isempty(flip_y) )
  flip_y = false;
end
if ( nargin < 4 || isempty(monitor_constants) )
  monitor_constants = hwwa.monitor_constants();
end

h = monitor_constants.monitor_height_cm;
d = monitor_constants.z_dist_to_subject_cm;
r = monitor_constants.vertical_resolution_px;

if ( flip_y )
  y = r - y;
end

x_deg = hwwa.px2deg( x, h, d, r );
y_deg = hwwa.px2deg( y, h, d, r );

end