function outs = hwwa_saccade_info(aligned_outs, varargin)

defaults = struct();
defaults.only_first_saccade = true;
params = hwwa.parsestruct( defaults, varargin );

max_num = ternary( params.only_first_saccade, 1, inf );

[x_deg, y_deg] = hwwa.run_px2deg( aligned_outs.x, aligned_outs.y, true );
start_stops = hwwa.run_find_saccades( x_deg, y_deg );

saccades = hwwa.saccade_directions( x_deg, y_deg, start_stops, false );
[saccades, labels, aligned_to_saccade_ind, saccade_to_aligned_ind] = ...
  hwwa.flatten_saccades( saccades, aligned_outs.labels, max_num );

start_stops = hwwa.flatten_saccades( start_stops, aligned_outs.labels, max_num );

outs = struct();

saccade_dirs = hwwa.normalize_saccades( saccades );
saccade_lengths = sqrt( sum(saccades .* saccades, 2) );

amount_ms = aligned_outs.params.look_ahead - aligned_outs.params.look_back;
vel_ratio = 1e3 / amount_ms;

collapsed_left_dirs = hwwa.collapse_left_saccades( saccade_dirs, labels );
saccade_angles_to_right = hwwa.saccade_angles( collapsed_left_dirs );

monitor_constants = hwwa.monitor_constants();

saccade_stop_points = stop_points( start_stops, aligned_outs.x, aligned_outs.y, aligned_to_saccade_ind );
saccade_stop_points(:, 2) = monitor_constants.vertical_resolution_px - saccade_stop_points(:, 2);

saccade_stop_points_to_right = collapse_points_to_right( saccade_stop_points );

outs.saccade_dirs = saccade_dirs;
outs.saccade_lengths = saccade_lengths;
outs.saccade_velocities = saccade_lengths * vel_ratio;
outs.saccade_peak_velocities = start_stops(:, 3);
outs.saccade_angles = hwwa.saccade_angles( saccade_dirs );
outs.saccade_angles_to_right = saccade_angles_to_right;
outs.aligned_to_saccade_ind = aligned_to_saccade_ind;
outs.saccade_to_aligned_ind = saccade_to_aligned_ind;
outs.saccade_start_stop = start_stops;
outs.saccade_stop_points = saccade_stop_points;
outs.saccade_stop_points_to_right = saccade_stop_points_to_right;
outs.labels = labels;

end

function pts = collapse_points_to_right(pts)

consts = hwwa.monitor_constants();
horz_res = consts.horizontal_resolution_px;
cx = horz_res / 2;

pts(:, 1) = abs( pts(:, 1) - cx );

end

function pts = stop_points(start_stops, x, y, aligned_to_saccade_ind)

pts = nan( size(start_stops, 1), 1 );

for i = 1:size(start_stops, 1)
  pts(i, 1) = x(aligned_to_saccade_ind(i), start_stops(i, 2));
  pts(i, 2) = y(aligned_to_saccade_ind(i), start_stops(i, 2));
end

end