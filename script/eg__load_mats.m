conf = hwwa.config.load();

datap = fullfile( conf.PATHS.data_root, 'raw' );

mats = shared_utils.io.find( datap, '.mat' );
mats = cellfun( @load, mats, 'un', false );
mats = cellfun( @(x) x.data, mats, 'un', false );

%%

