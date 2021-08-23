function out = hwwa_load_approach_avoid_saccade_info(varargin)

load_p = hwwa.gid( 'processed/saccade_info', varargin{:} );
out = shared_utils.io.fload( fullfile(load_p, 'saccade_info.mat') );

end