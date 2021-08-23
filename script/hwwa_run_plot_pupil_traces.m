conf = hwwa.config.load();
cue_on_aligned = hwwa_approach_avoid_cue_on_aligned( ...
    'config', conf ...
  , 'look_back', -150 ...
  , 'look_ahead', 1.5e3 ...
);

behav_outs = hwwa_load_approach_avoid_behavior( 'config', conf );
aligned_outs = hwwa_approach_avoid_fix_on_aligned( 'config', conf );

assert( cue_on_aligned.labels == aligned_outs.labels ...
  , 'Expected cue aligned labels to equal fix aligned labels.' );

[~, use_labs] = hwwa_approach_avoid_make_apply_outlier_labels( aligned_outs, behav_outs );

%%

find_non_outliers = @(labels, varargin) find( labels, 'is_outlier__false', varargin{:} );

%%

use_aligned_outs = cue_on_aligned;
use_aligned_outs.labels = use_labs';
t = use_aligned_outs.time(1, :);

is_fix = hwwa.fixation_detection( use_aligned_outs.x, use_aligned_outs.y, t );

%%

thresh = 300;
c_start = find( t == 0 );
subset_fix = is_fix(:, c_start:c_start+thresh-1);
crit = sum( subset_fix, 2 ) >= thresh;

%%

all_go = hwwa.find_go( use_aligned_outs.labels );
all_nogo = hwwa.find_nogo( use_aligned_outs.labels );
crit_nogo = intersect( all_nogo, find(crit) );
pass_crit = union( all_go, crit_nogo );

p_nogo_crit = numel( crit_nogo ) / numel( all_nogo );

crit_restrict = ...
  @(l, m) find_non_outliers(l, intersect(m, pass_crit));

%%

require_crit = true;
mask_func = ternary( require_crit, crit_restrict, find_non_outliers );
prefix = ternary( require_crit, sprintf('crit-%d', thresh), '' );

hwwa_plot_pupil_traces( use_aligned_outs ...
  , 'mask_func', mask_func ...
  , 'time_limits', [0, 800] ...
  , 'do_save', true ...
  , 'smooth_func', @(x) smoothdata(x, 'SmoothingFactor', 0.25) ...
  , 'prefix', prefix ...
);