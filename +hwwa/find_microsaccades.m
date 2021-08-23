function start_stops = find_microsaccades(deg_x, deg_y, varargin)

check_inputs( deg_x, deg_y );

defaults = struct();
defaults.fs = 1e3;
defaults.min_dur = 15;
defaults.max_dur = 40;
defaults.check_finite_pad = 5;
defaults.vel_thresh = 100;
defaults.amp_thresh = inf;
defaults.smooth_func = @(x) filter(gausswin(7), 1, x);

params = hwwa.parsestruct( defaults, varargin );

check_finite_pad = params.check_finite_pad;

start_stops = hwwa.find_saccades( ...
    deg_x(:)' ...
  , deg_y(:)' ...
  , params.fs ...
  , params.vel_thresh ...
  , params.min_dur ...
  , params.smooth_func ...
);

% Collapse overlapping intervals.
for i = 1:numel(start_stops)
  start_stops{i} = make_unique_saccades( start_stops{i}, numel(deg_x) );
end

% Remove microsaccades where NaN or Inf values are within some distance
% of the saccade start / stop
for i = 1:numel(start_stops)
  tf = all_finite( deg_x, deg_y ...
    , pad_starts(start_stops{i}(:, 1), check_finite_pad) ...
    , pad_stops(start_stops{i}(:, 2), check_finite_pad, numel(deg_x)) ...
  );
  start_stops{i}(~tf, :) = [];
end

% Remove microsaccades that are too long.
for i = 1:numel(start_stops)
  durs = start_stops{i}(:, 2) - start_stops{i}(:, 1);
  start_stops{i}(durs > params.max_dur, :) = [];
end

% Remove microsaccades that are too large in amplitude.
for i = 1:numel(start_stops)
  amp = saccade_amplitude( ...
    start_stops{i}(:, 1), start_stops{i}(:, 2), deg_x, deg_y );
  start_stops{i}(amp > params.amp_thresh, :) = [];
end

if ( isempty(deg_x) )
  start_stops(1) = [];
end

end

function a = saccade_amplitude(starts, stops, x, y)
  
a = zeros( size(starts) );
for i = 1:numel(starts)
  p0 = [x(starts(i)), y(starts(i))];
  p1 = [x(stops(i)), y(stops(i))];
  a(i) = norm( p1 - p0 );
end
  
end

function starts = pad_starts(starts, amt)

starts = max( 1, starts - amt );

end

function stops = pad_stops(stops, amt, len)

stops = min( len, stops + amt );

end

function tf = all_finite(x, y, starts, stops)

tf = false( size(starts) );
for i = 1:numel(starts)
  tf(i) = all( isfinite(x(starts(i):stops(i))) & isfinite(y(starts(i):stops(i))) );
end

end

function info = make_unique_saccades(info, len)

[starts, stops] = make_unique( info(:, 1), info(:, 2), len );
info = [ starts(:), stops(:), nan(numel(starts), size(info, 2)-2) ];

end

function [starts, stops] = make_unique(starts, stops, len)

v = false( 1, len );
for i = 1:numel(starts)
  v(starts(i):stops(i)) = true;
end

[starts, durs] = shared_utils.logical.find_islands( v );
stops = starts + durs - 1;

end

function check_inputs(deg_x, deg_y)

if ( ~isempty(deg_x) )
  validateattributes( deg_x, {'double'}, {'vector'}, mfilename, 'deg_x' );
  validateattributes( deg_y, {'double'}, {'vector', 'numel', numel(deg_x)} ...
    , mfilename, 'deg_y' );
else
  assert( isempty(deg_y), ['deg_x and deg_y must be vectors with the same' ...
    , ' number of elements, or both empty.'] );
end

end