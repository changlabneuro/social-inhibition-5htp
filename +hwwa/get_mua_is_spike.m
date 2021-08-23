function data = get_mua_is_spike(data, ndevs)

devs = std( data );
means = mean( data );

thresh1 = means - (devs * ndevs);
thresh2 = means + (devs * ndevs);

data = data < thresh1 | data > thresh2;

end