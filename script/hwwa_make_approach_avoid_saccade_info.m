function hwwa_make_approach_avoid_saccade_info(varargin)

defaults = hwwa.get_common_make_defaults();
params = hwwa.parsestruct( defaults, varargin );

conf = params.config;

behav_outs = hwwa_load_approach_avoid_behavior( 'config', conf );
go_targ_aligned_outs = hwwa_approach_avoid_go_target_on_aligned( 'config', conf );
saccade_outs = hwwa_saccade_info( go_targ_aligned_outs );

behav_ids = categorical( behav_outs.labels, {'unified_filename', 'trial_number'} );
aligned_ids = categorical( go_targ_aligned_outs.labels, {'unified_filename', 'trial_number'} );
assert( isequal(behav_ids, aligned_ids) );

rt_outs = make_behav_outs_sized_rt( behav_outs, go_targ_aligned_outs, saccade_outs );
saccade_outs.rt = rt_outs;

save_p = get_save_p( params );
shared_utils.io.require_dir( save_p );
save( fullfile(save_p, get_filename(params)), 'saccade_outs' );

end

function outs = make_behav_outs_sized_rt(behav_outs, go_targ_aligned_outs, saccade_outs)

saccade_based_rt = hwwa.saccade_based_reaction_time( go_targ_aligned_outs, saccade_outs );
target_entry_based_rt = hwwa.target_entry_based_rt( go_targ_aligned_outs );
movement_time = target_entry_based_rt - saccade_based_rt;

nogo_selectors = { 'nogo_trial', 'correct_false' };

nogo_ind = find( behav_outs.labels, nogo_selectors );
go_ind = hwwa.find_correct_go( behav_outs.labels );

src_rt = behav_outs.rt;

saccade_based_rt = make_adjusted_rt( saccade_based_rt, src_rt, go_ind, nogo_ind );
target_entry_based_rt = make_adjusted_rt( target_entry_based_rt, src_rt, go_ind, nogo_ind );
movement_time = make_adjusted_rt( movement_time, src_rt, go_ind, nogo_ind );

outs = struct();
outs.saccade_based_rt = saccade_based_rt;
outs.target_entry_based_rt = target_entry_based_rt;
outs.movement_time = movement_time;

end

function adjusted_rt = make_adjusted_rt(source_rt, target_rt, go_ind, nogo_ind)

adjusted_rt = target_rt;
adjusted_rt(nogo_ind) = source_rt(nogo_ind) / 1e3;
adjusted_rt(go_ind) = source_rt(go_ind) / 1e3;

end

function save_p = get_save_p(params)

save_p = hwwa.gid( 'processed/saccade_info', params.config );

end

function fname = get_filename(params)

fname = 'saccade_info.mat';

end