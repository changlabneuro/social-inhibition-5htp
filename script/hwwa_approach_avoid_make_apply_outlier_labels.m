function [pupil_size, pupil_labels] = hwwa_approach_avoid_make_apply_outlier_labels(aligned_outs, behav_outs)

n_devs = 2;

pupil_labels = aligned_outs.labels';
pupil_size = nanmean( aligned_outs.pupil, 2 );
mask = find( pupil_labels, 'initiated_true' );

norm_each = { 'unified_filename' };

in_bounds = hwwa.within_std_threshold( pupil_size, pupil_labels, n_devs, norm_each, mask );
hwwa.apply_outlier_labels( pupil_labels, in_bounds, 'unified_filename', behav_outs.labels );
hwwa.apply_outlier_labels( pupil_labels, in_bounds, 'unified_filename', pupil_labels );

end