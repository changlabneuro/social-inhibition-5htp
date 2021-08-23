function out = hwwa_load_approach_avoid_behavioral_data(varargin)

save_p = hwwa.gid( 'processed/behavior', varargin{:} );
out = shared_utils.io.fload( fullfile(save_p, 'behav.mat') );

end