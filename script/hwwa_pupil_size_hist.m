aligned_outs = hwwa_load_edf_aligned( ...
    'start_event_name', 'go_nogo_cue_onset' ...
  , 'look_back', -150 ...
  , 'look_ahead', 0 ...
  , 'is_parallel', true ...
  , 'files_containing', hwwa.approach_avoid_files() ...
  , 'error_handler', 'error' ...
);

%%

n_devs = 2;

labels = aligned_outs.labels';

pl = plotlabeled.make_common();
pl.hist_add_summary_line = true;
pl.summary_func = @plotlabeled.nanmedian;

pupil_size = nanmean( aligned_outs.pupil, 2 );

mask = fcat.mask( labels ...
  , @find, 'initiated_true' ...
  , @find, 'hitch' ...
);

norm_each = { 'unified_filename' };

% pupil_size = hwwa.median_normalize_pupil( pupil_size, labels, norm_each, mask );
all_in_bounds = hwwa.within_std_threshold( pupil_size, labels, n_devs, norm_each, mask );

subset_pupil = pupil_size(mask, :);
subset_labels = labels(mask);

fcats = { 'monkey', 'drug' };
pcats = { 'monkey', 'scrambled_type', 'target_image_category', 'drug' };

figs = pl.figures( @hist, subset_pupil, subset_labels, fcats, pcats );