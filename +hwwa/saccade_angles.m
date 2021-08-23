function thetas = saccade_angles(dirs)

thetas = atan2( dirs(:, 2), dirs(:, 1) );

end