function psth = get_mua_psth( data, ndevs )

%   GET_MUA_PSTH -- Convert voltage data to logical index of spike.

N = size( data, 2 );

devs = std( data, [], 2 );
means = mean( data, 2 );
thresh1 = means - (devs .* ndevs);
thresh2 = means + (devs .* ndevs);

thresh1 = repmat( thresh1, 1, N );
thresh2 = repmat( thresh2, 1, N );

psth = data < thresh1 | data > thresh2;

end