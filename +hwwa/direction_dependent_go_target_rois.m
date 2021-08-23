function out_rois = direction_dependent_go_target_rois(rois, labels, x_displacement, y_displacement)

assert_ispair( rois, labels );

if ( nargin < 3 )
  x_displacement = zeros( rows(rois), 1 );
else
  assert_ispair( x_displacement, labels );
end

if ( nargin < 4 )
  y_displacement = zeros( rows(rois), 1 );
else
  assert_ispair( y_displacement, labels );
end

constants = hwwa.monitor_constants();

left_trials = find( labels, 'center-left' );
right_trials = find( labels, 'center-right' );

left_cx = constants.horizontal_resolution_px / 4;
right_cx = constants.horizontal_resolution_px / 2 + constants.horizontal_resolution_px / 4;

cx = nan( rows(labels), 1 );
cx(left_trials) = left_cx;
cx(right_trials) = right_cx;

cy = repmat( constants.vertical_resolution_px / 2, size(cx) );

roi_widths = shared_utils.rect.width( rois );
roi_heights = shared_utils.rect.height( rois );

use_displacement_x = x_displacement;
use_displacement_x(left_trials) = -use_displacement_x(left_trials);

cx = cx + use_displacement_x;
cy = cy + y_displacement;

x0 = cx - roi_widths/2;
x1 = cx + roi_widths/2;
y0 = cy - roi_heights/2;
y1 = cy + roi_heights/2;

out_rois = [ x0, y0, x1, y1 ];


end