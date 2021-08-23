function start_stops = find_saccades(x, y, sr, vel_thresh, min_samples, smooth_func)

%   FIND_SACCADES -- Find start and stop indices of saccades.
%
%     indices = hwwa.find_saccades( x, y, sr, vel_thresh, min_samples, smooth_func )
%     returns indices into columns of `x` and `y` denoting the start(s) and
%     stop(s) of saccades.
%
%     `x` and `y` are matrices of the same size whose whose rows are 
%     observations and columns are samples, assumed to be in units of 
%     degrees of visual angle. `sr` is the sampling rate of `x` and y`. 
%     `vel_thresh` is the minimum velocity threshold (in deg/s) for marking
%     the onset of a saccade. `min_samples` is the minimum duration of a
%     saccade in units samples. `smooth_func` is a handle to a function
%     that receives a row-vector (i.e., a row of `x` or `y`) and returns a
%     smoothed version of that vector.
%
%     The output `indices` is a cell array with the same number of rows as
%     `x` and `y`. Each element contains an Mx2 array of indices into the
%     corresponding row of `x` and `y`. Each row of this sub-array
%     represents a saccade; column 1 is the start index, and column 2 is 
%     the stop index.
%
%     See also smoothdata

assert_rowsmatch( x, y );
num_trials = size( x, 1 );

start_stops = cell( num_trials, 1 );

parfor i = 1:num_trials
  x0 = smooth_func( x(i, :) );
  y0 = smooth_func( y(i, :) );
  
  vx = abs( gradient(x0) ) * sr;
  vy = abs( gradient(y0) ) * sr;
  
  above_thresh_x = vx >= vel_thresh;
  above_thresh_y = vy >= vel_thresh;
  
  x_starts = shared_utils.logical.find_islands( above_thresh_x );
  y_starts = shared_utils.logical.find_islands( above_thresh_y );
  
  [x_starts, x_stops] = find_vel_stops( vx, x_starts, vel_thresh );
  [y_starts, y_stops] = find_vel_stops( vy, y_starts, vel_thresh );
  
  [starts, stops] = merge_intervals( x_starts, x_stops, y_starts, y_stops );
  
  durs = stops - starts;
  within_thresh = durs >= min_samples;
  
  starts = starts(within_thresh);
  stops = stops(within_thresh);
  peak_velocities = peak_velocity( vx, vy, starts, stops );
  
  start_stops{i} = [starts(:), stops(:), peak_velocities];
end

end

function peak_vel = peak_velocity(vx, vy, starts, stops)

peak_vel = zeros( numel(starts), 1 );

for i = 1:numel(starts)
  start = starts(i);
  stop = stops(i);
  
  peak_vel(i) = max( max(vx(start:stop)), max(vy(start:stop)) );
end

end

function [starts, stops] = merge_intervals(x_starts, x_stops, y_starts, y_stops)

if ( isempty(x_starts) )
  starts = y_starts;
  stops = y_stops;
  return;
  
elseif ( isempty(y_starts) )
  starts = x_starts;
  stops = x_stops;
  return;
end

starts = [];
stops = [];

x_stp = 1;
y_stp = 1;

num_x = numel( x_starts );
num_y = numel( y_starts );

while ( x_stp <= num_x && y_stp <= num_y )
  x0 = x_starts(x_stp);
  x1 = x_stops(x_stp);
  
  y0 = y_starts(y_stp);
  y1 = y_stops(y_stp);
  
  if ( x1 < y0 )
    starts(end+1) = x0;
    stops(end+1) = x1;

    x_stp = x_stp + 1;
  elseif ( y1 < x0 )
    starts(end+1) = y0;
    stops(end+1) = y1;

    y_stp = y_stp + 1;
  else    
    use0 = min( x0, y0 );
    use1 = max( x1, y1 );
    
    starts(end+1) = use0;
    stops(end+1) = use1;
   
    x_stp = x_stp + 1;
    y_stp = y_stp + 1;
  end
end

end

function [starts, stops] = find_vel_stops(vel, starts, threshold)

stops = nan( size(starts) );
to_remove = false( size(starts) );

for i = 1:numel(starts)
  j = starts(i);
  
  if ( i < numel(starts) )
    end_point = starts(i+1);
  else
    end_point = numel( vel ) + 1;
  end
  
  while ( j < end_point && vel(j) > threshold )
    j = j + 1;
  end
  
  if ( j > starts(i) && j ~= end_point )
    stops(i) = j;
  else
    to_remove(i) = true;
  end
end

stops = stops(~to_remove);
starts = starts(~to_remove);

end

function [source_starts, merged_source, merged_target] = merge_starts(source_starts, target_starts, threshold)

% nearest_target_starts = bfw.find_nearest( source_starts, target_starts );
non_nan_starts = find( ~isnan(source_starts) );
nearest_target_starts = non_nan_starts(bfw.find_nearest(source_starts(non_nan_starts), target_starts));

merged_source = false( size(source_starts) );
merged_target = false( size(target_starts) );
  
for i = 1:numel(nearest_target_starts)
  source_ind = nearest_target_starts(i);
  offset = abs( source_starts(source_ind) - target_starts(i) );

  if ( offset <= threshold )
    source_starts(source_ind) = max( source_starts(source_ind), target_starts(i) );
    merged_source(source_ind) = true;
    merged_target(i) = true;
  end
end

end

% function debug_plot_starts_stops(starts, stops)
% 
% x_ax = gca();
% 
% hold( x_ax, 'off' );
% plot( x_ax, x0, 'r' ); 
% shared_utils.plot.hold( x_ax, 'on' );
% plot( x_ax, y0, 'b' );
% 
% shared_utils.plot.add_vertical_lines( x_ax, starts, 'k' );
% shared_utils.plot.add_vertical_lines( x_ax, stops, 'k--' );
% 
% end