function pup = time_normalize(pup, t, trange)

assert( numel(t) == size(pup, 2) ...
  , 'Time series does not match pupil matrix col dimension.' );
baseline = nanmean( pup(:, t >= trange(1) & t <= trange(2)), 2 );
pup = pup ./ baseline;

end