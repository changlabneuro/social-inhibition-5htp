function saccade_dirs = saccade_directions(x, y, saccade_inds, normalize)

if ( nargin < 4 )
  normalize = true;
end

assert_rowsmatch( x, y );
assert_rowsmatch( x, saccade_inds );

saccade_dirs = cell( size(saccade_inds) );

for i = 1:numel(saccade_inds)
  per_trial_saccades = saccade_inds{i};
  per_trial_dirs = nan( size(per_trial_saccades, 1), 2 );
  
  for j = 1:size(per_trial_saccades, 1)
    dx = x(i, per_trial_saccades(j, 2)) - x(i, per_trial_saccades(j, 1));
    dy = y(i, per_trial_saccades(j, 2)) - y(i, per_trial_saccades(j, 1));
    
    dir = [dx, dy];

    if ( normalize )
      dir = dir ./ sqrt( sum(dir .* dir) );
    end
    
    per_trial_dirs(j, :) = dir;
  end
  
  saccade_dirs{i} = per_trial_dirs;
end

end