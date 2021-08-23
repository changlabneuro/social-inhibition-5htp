conf = hwwa.config.load();
aligned_outs = hwwa_approach_avoid_go_target_on_aligned( ...
    'config', conf ...
  , 'look_ahead', 2e3 ...
);

if ( isempty(aligned_outs.event_key) )
  go_targ_offset_col = [];
else
  go_targ_offset_col = aligned_outs.event_key{1}('go_target_offset');
end

saccade_outs = hwwa_saccade_info( aligned_outs );

%%

labels = aligned_outs.labels';

to_sacc_ind = saccade_outs.aligned_to_saccade_ind;
assert( numel(unique(to_sacc_ind)) == numel(to_sacc_ind), 'More than 1 saccade per trial.' );

go_target_offsets = aligned_outs.sample_relative_events(:, go_targ_offset_col);
saccade_stop_points = nan( rows(go_target_offsets), 2 );
saccade_stop_points(to_sacc_ind, :) = saccade_outs.saccade_stop_points;

start_stops = nan( rows(go_target_offsets), 3 );
start_stops(to_sacc_ind, :) = saccade_outs.saccade_start_stop;

nogo_cue_rois = hwwa.linearize_roi( aligned_outs, 'nogo_cue' );
go_target_rois = hwwa.linearize_roi( aligned_outs, 'go_target_padded' );
go_targ_displacements = hwwa.linearize_roi( aligned_outs, 'go_target_displacement' );

is_left = trueat( labels, find(labels, 'center-left') );
constants = hwwa.monitor_constants();

hres = constants.horizontal_resolution_px;
vres = constants.vertical_resolution_px;

cl_x = hres / 4;
cr_x = cl_x + hres / 2;

use_xs = nan( rows(labels), 1 );
use_xs(is_left) = cl_x;
use_xs(~is_left) = cr_x;

began_look_to_targ_before_dissappeared = false( rows(labels), 1 );
use_go_target_stop_point_bounds = false;

left_screen_roi = [0, 0, hres/2, vres];
right_screen_roi = [hres/2, 0, hres, vres];

for i = 1:rows(labels)
  go_roi = go_target_rois(i, :);
  w = shared_utils.rect.width( go_roi );
  x = use_xs(i);
  
  displacement = go_targ_displacements(i, :);
  if ( ~is_left(i) )
    displacement([1, 3]) = -displacement([1, 3]);
  end
  
  if ( use_go_target_stop_point_bounds )
    roi = [x-w/2, go_roi(2), x+w/2, go_roi(4)] + displacement;
    
  else
    if ( is_left(i) )
      roi = left_screen_roi;
    else
      roi = right_screen_roi;
    end
  end
    
  stop = saccade_stop_points(i, :);
  
  if ( shared_utils.rect.inside(roi, stop) )
    start_ind = start_stops(i, 1);
    stop_ind = start_stops(i, 2);
    targ_off_ind = go_target_offsets(i);
    
    began_look_to_targ_before_dissappeared(i) = ...
      start_ind <= targ_off_ind && stop_ind > targ_off_ind;
  end
end

%%

pltlabs = labels';

saccade_epoch_labels = repmat( {'saccade-unspecified'}, size(began_look_to_targ_before_dissappeared) );
saccade_epoch_labels(began_look_to_targ_before_dissappeared) = {'saccaded-too-late'};
addsetcat( pltlabs, 'saccade_epoch', saccade_epoch_labels );

mask = hwwa.get_approach_avoid_mask( pltlabs );
% mask = find( pltlabs, 'correct_false', mask );
mask = hwwa.find_incorrect_go( pltlabs, mask );
mask = find( pltlabs, 'saccaded-too-late', mask );

pl = plotlabeled.make_common();
pl.summary_func = @sum;

pcats = { 'drug', 'trial_type', 'correct' };
gcats = { 'saccade_epoch' };
xcats = { 'monkey' };

dat = ones(size(mask));
labs = prune( pltlabs(mask) );

axs = pl.bar( dat, labs, xcats, gcats, pcats );

%%

saccades_to_go_target = shared_utils.rect.inside( go_target_rois, saccade_stop_points );

x = aligned_outs.x;
y = aligned_outs.y;
p_inside_nogo_cue = zeros( rows(x), 1 );
num_looks_to_go = zeros( rows(x), 1 );
go_looking_duration = zeros( rows(x), 1 );
first_go_looking_time = zeros( rows(x), 1 );
time_first_go_look_to_offset = zeros( rows(x), 1 );
time_last_go_look_to_offset = zeros( rows(x), 1 );
last_go_looking_duration = zeros( rows(x), 1 );

for i = 1:rows(x)
  offset_ind = go_target_offsets(i);
  
  if ( ~isnan(offset_ind) )
    nogo_roi = nogo_cue_rois(i, :);
    go_roi = go_target_rois(i, :);
    
    x_ = x(i, 1:offset_ind);
    y_ = y(i, 1:offset_ind);
    
    ib_nogo = shared_utils.rect.inside( nogo_roi, x_, y_ );
    ib_go = shared_utils.rect.inside( go_roi, x_, y_ );
    
    [islands, durs] = shared_utils.logical.find_islands( ib_go );
    
    p_inside_nogo_cue(i) = pnz( ib_nogo );
    num_looks_to_go(i) = numel( islands );
    go_looking_duration(i) = sum( ib_go );
    
    if ( ~isempty(islands) )
      first_go_looking_time(i) = islands(1);
      time_first_go_look_to_offset(i) = offset_ind - islands(1) + 1;
      time_last_go_look_to_offset(i) = offset_ind - islands(end) + 1;
      last_go_looking_duration(i) = durs(end);
    end
  end
end

%%

incorrect_reason = repmat( {'unspecified'}, rows(labels), 1 );

wrong_go_ind = hwwa.find_incorrect_go( labels );
wrong_nogo_ind = hwwa.find_incorrect_nogo( labels );
go_mask = trueat( labels, hwwa.find_go(labels) );

nogo_in_bounds_thresh = 0.8;

for i = 1:numel(wrong_go_ind)
  ind = wrong_go_ind(i);
  
  if ( p_inside_nogo_cue(ind) >= nogo_in_bounds_thresh )
    incorrect_reason{ind} = 'stayed-in-nogo-cue-bounds';
    
  elseif ( num_looks_to_go(ind) > 0 )
    incorrect_reason{ind} = 'insufficient-go-target-fixation';
    
  else
    incorrect_reason{ind} = 'looked-around';
  end
end

for i = 1:numel(wrong_nogo_ind)
  ind = wrong_nogo_ind(i);
  
  if ( saccades_to_go_target(ind) || num_looks_to_go(ind) > 0 )
    incorrect_reason{ind} = 'saccade-to-go-target';
    
  else
    incorrect_reason{ind} = 'looked-to-go-target-region';
  end
end

incorrect_reason(began_look_to_targ_before_dissappeared & go_mask) = {'saccaded-too-late'};

addsetcat( labels, 'incorrect_reason', incorrect_reason );
  
%%  Number

do_save = true;
params = hwwa.get_common_plot_defaults( hwwa.get_common_make_defaults() );
save_p = hwwa.approach_avoid_data_path( params, 'plots', 'error_classification' );

mask = hwwa.get_approach_avoid_mask( labels );
mask = find( labels, 'correct_false', mask );

pl = plotlabeled.make_common();
pl.summary_func = @sum;

pcats = { 'drug', 'trial_type', 'correct' };
gcats = { 'incorrect_reason' };
xcats = { 'monkey' };

dat = ones(size(mask));
labs = prune( labels(mask) );

% percent
p_of = 'incorrect_reason';
p_each = setdiff( unique([pcats, gcats, xcats]), p_of );
[ps, plabs] = proportions_of( labs', p_each, p_of );

% axs = pl.stackedbar( dat, labs, xcats, gcats, pcats );
axs = pl.stackedbar( ps, plabs, xcats, gcats, pcats );

if ( do_save )
  fname = dsp3.req_savefig( gcf, save_p, labels, pcats );
end



