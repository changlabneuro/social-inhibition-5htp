function is_fix = fixation_detection(x, y, t, varargin)

default_fix_params = bfw.make.defaults.raw_fixations();
defaults = struct();
defaults.t1 = default_fix_params.t1;
defaults.t2 = default_fix_params.t2;
defaults.min_dur = 50;

params = hwwa.parsestruct( defaults, varargin );

assert( size(x, 2) == numel(t), 'Number of x samples must match number of time points.' );
assert( isequal(size(x), size(y)), 'Number of x samples must match number of y samples.' );

t1 = params.t1;
t2 = params.t2;
min_dur = params.min_dur;

[num_rows, num_cols] = size( x );
is_fix = false( num_rows, size(x, 2) );

parfor i = 1:num_rows
  use_pos = [x(i, :); y(i, :)];
  is_fix_tmp = is_fixation( use_pos, t, t1, t2, min_dur );      
  is_fix(i, :) = is_fix_tmp(1:num_cols);
end

% hold off;
% plot( use_pos(1, :), 'r' );
% hold on;
% plot( use_pos(2, :), 'b' );
% 
% test = is_fix(1, :);
% l = get( gca, 'ylim' );
% test = test * 0.5 * (diff(l)) + min(l);
% plot( test, 'g', 'linewidth', 0.5 );

end