function [heat_map, x0, y0] = hwwa_gaze_heatmap(x, y, x_lims, y_lims, stp, win, mask)

if ( nargin < 7 )
  mask = rowmask( x );
end

assert_rowsmatch( x, y );

x0 = x_lims(1):stp:x_lims(2);
x1 = x0 + win;

y0 = y_lims(1):stp:y_lims(2);
y1 = y0 + win;

heat_map = zeros( numel(y0), numel(x0) );

num_trials = numel( mask );

for i = 1:num_trials
  x_t = x(mask(i), :);
  y_t = y(mask(i), :);
    
  for j = 1:numel(x0)
    for k = 1:numel(y0)
      ib = x_t >= x0(j) & x_t < x1(j) & y_t >= y0(k) & y_t < y1(k);
      heat_map(k, j) = heat_map(k, j) + sum( ib );
    end
  end
end

end