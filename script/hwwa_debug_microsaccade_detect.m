fs = 1e3;

trial = 12;
x1 = behav_outs.cue_on_aligned.x(trial, :);
y1 = behav_outs.cue_on_aligned.y(trial, :);
t = behav_outs.cue_on_aligned.time;

[x1, y1] = hwwa.run_px2deg( x1, y1 );

vel_thresh = 100;
min_dur = 15;
max_dur = 40;
check_finite_pad = 5;

smooth_func = @(x) filter(gausswin(7), 1, x);
start_stops = hwwa.find_saccades( x1(:)', y1(:)', fs, vel_thresh, min_dur, smooth_func );

for i = 1:numel(start_stops)
  start_stops{i} = make_unique_saccades( start_stops{i}, numel(t) );
end

for i = 1:numel(start_stops)
  tf = all_finite( x1, y1 ...
    , pad_starts(start_stops{i}(:, 1), check_finite_pad) ...
    , pad_stops(start_stops{i}(:, 2), check_finite_pad, numel(t)) ...
  );
  start_stops{i}(~tf, :) = [];
end

for i = 1:numel(start_stops)
  durs = start_stops{i}(:, 2) - start_stops{i}(:, 1);
  start_stops{i}(durs > max_dur, :) = [];
end

plot_saccade( t, x1, y1, start_stops{1}(:, 1), start_stops{1}(:, 2) );

%%

subj_inds = findall( behav_outs.labels, 'monkey' );
trial = subj_inds{1}(18);

x1 = behav_outs.cue_on_aligned.x(trial, :);
y1 = behav_outs.cue_on_aligned.y(trial, :);
[x1, y1] = hwwa.run_px2deg( x1, y1 );

start_stops = hwwa.find_microsaccades( x1, y1, 'max_dur', 60, 'amp_thresh', 30 );

axs = plot_saccade( t, x1, y1, start_stops{1}(:, 1), start_stops{1}(:, 2) );
shared_utils.plot.set_ylims( axs, [-200, 200] );
shared_utils.plot.set_xlims( axs, [min(t), max(t)] );

%%

x = behav_outs.cue_on_aligned.x;
y = behav_outs.cue_on_aligned.y;
start_stops = cell( size(x, 1), 1 );

parfor i = 1:numel(start_stops)
  [deg_x, deg_y] = hwwa.run_px2deg( x(i, :), y(i, :) );
  start_stops{i} = hwwa.find_microsaccades( deg_x, deg_y ...
    , 'max_dur', 60 ...
    , 'amp_thresh', 30 ...
  );
end

%%

detect_thresh = 10;

[x1, y1] = hwwa.run_px2deg( x1, y1 );

win = gausswin( 100 );
x1 = filter( win, 1, x1 );
y1 = filter( win, 1, y1 );

p = [ x1(:), y1(:) ];

v = diff( p, 1, 1 );
if ( ~isempty(p) )
  v = [ [0, 0]; v ];
end

a = diff( v, 1, 1 );
if ( ~isempty(p) )
  a = [ [0, 0]; a ];
end

v = v .* fs;
a = a .* fs;

speed = vecnorm( v, 2, 2 );
accel_len = vecnorm( a, 2, 2 );
ptr = 1;

while ( ptr <= numel(speed) )
  if ( speed(ptr) > detect_thresh )
    sacc = [];
    while ( ptr <= numel(speed) && speed(ptr) > detect_thresh )
      sacc(end+1, 1) = ptr;
      ptr = ptr + 1;
    end
    
    plot_saccade( t, p(:, 1), p(:, 2), sacc(1), sacc(end) );
    
    d = 10;
  end
  ptr = ptr + 1;
end

%%

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
info = [ starts(:), stops(:) ];

end

function [starts, stops] = make_unique(starts, stops, len)

v = false( 1, len );
for i = 1:numel(starts)
  v(starts(i):stops(i)) = true;
end

[starts, durs] = shared_utils.logical.find_islands( v );
stops = starts + durs - 1;

end

function axs = plot_saccade(t, x, y, starts, stops)

ax1 = subplot( 1, 2, 1 );
cla( ax1 );
plot( ax1, t, x );
title( ax1, 'x' );

for i = 1:numel(starts)
  hold( ax1, 'on' );
  shared_utils.plot.add_vertical_lines( ax1, t(starts(i)) );
  shared_utils.plot.add_vertical_lines( ax1, t(stops(i)), 'r--' );
end

ax2 = subplot( 1, 2, 2 );
cla( ax2 );
plot( ax2, t, y );
title( ax2, 'y' );

for i = 1:numel(starts)
  hold( ax2, 'on' );
  shared_utils.plot.add_vertical_lines( ax2, t(starts(i)) );
  shared_utils.plot.add_vertical_lines( ax2, t(stops(i)), 'r--' );
end

axs = [ ax1, ax2 ];

end