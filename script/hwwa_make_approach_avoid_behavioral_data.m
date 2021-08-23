function hwwa_make_approach_avoid_behavioral_data(varargin)

defaults = hwwa.get_common_make_defaults();
params = hwwa.parsestruct( defaults, varargin );

conf = params.config;

behav_outs = hwwa_load_approach_avoid_behavior( 'config', conf );
aligned_outs = hwwa_approach_avoid_fix_on_aligned( 'config', conf );

cue_on_aligned = hwwa_approach_avoid_cue_on_aligned( ...
    'config', conf ...
  , 'look_back', -150 ...
  , 'look_ahead', 1.5e3 ...
);

assert( cue_on_aligned.labels == aligned_outs.labels ...
  , 'Expected cue aligned labels to equal fix aligned labels.' );

%%

check_label_eq( behav_outs.labels', aligned_outs.labels' );

%%

pupil_size = hwwa_approach_avoid_make_apply_outlier_labels( aligned_outs, behav_outs );
behav_outs.pupil_size = pupil_size;
behav_outs.cue_on_aligned = struct( ...
  'x', cue_on_aligned.x ...
  , 'y', cue_on_aligned.y ...
  , 'pupil', cue_on_aligned.pupil ...
  , 'time', cue_on_aligned.time(1, :) ...
);

save_p = get_save_p( params );
file_path = fullfile( save_p, get_filename(params) );
shared_utils.io.require_dir( save_p );

save( file_path, 'behav_outs', '-v7.3' );

end

function check_label_eq(bl, al)

c = intersect( getcats(bl), getcats(al) );
miss_b = setdiff( getcats(bl), c );
miss_a = setdiff( getcats(al), c );

miss_b = union( miss_b, 'target_image_category' );
miss_a = union( miss_a, 'target_image_category' );

rmcat( bl, miss_b );
rmcat( al, miss_a );

assert( al == bl, 'Expected behav labels to equal aligned labels.' );

end

function fname = get_filename(params)

fname = 'behav.mat';

end

function save_p = get_save_p(params)

save_p = hwwa.gid( 'processed/behavior', params.config );

end