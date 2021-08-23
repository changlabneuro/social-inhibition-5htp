function out = hwwa_check_go_targ_looking_duration(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.start_event = 'reward_onset';
defaults.stop_event = 'iti';
defaults.use_reward_offset_time = true;

inputs = { 'edf_trials', 'events', 'unified', 'labels' };

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, '', defaults, varargin );
loop_runner.convert_to_non_saving_with_output();

results = loop_runner.run( @get_looking_duration, params );

results(~[results.success]) = [];

outputs = [results.output];

out = struct();
out.looking_duration = vertcat( outputs.looking_duration );
out.labels = vertcat( fcat(), outputs.labels );

end

function out = get_looking_duration(files, params)

edf_trials_file = shared_utils.general.get( files, 'edf_trials' );
labels_file = shared_utils.general.get( files, 'labels' );
unified_file = shared_utils.general.get( files, 'unified' );
events_file = shared_utils.general.get( files, 'events' );

start_event = params.start_event;
stop_event = params.stop_event;

labels = labels_file.labels;
update_labels( labels, unified_file );

go_targ = unified_file.opts.STIMULI.go_target;
displacement = unified_file.opts.STIMULI.setup.go_target.displacement;
amount_reward = unified_file.opts.TIMINGS.time_in.reward * 1000;  % ms

placements = labels_file.labels(:, 'target_placement');

assert( isempty(setdiff(placements, {'center-left', 'center-right'})) ...
  , 'Extraneous placement(s): "%s".', strjoin(unique(placements), ' | ') );

% Get left bounding rect
put( go_targ, 'center-left' );
shift( go_targ, -displacement(1), displacement(2) );
left_rect = go_targ.vertices;

% Get right bounding rect
put( go_targ, 'center-right' );
shift( go_targ, displacement(1), displacement(2) );
right_rect = go_targ.vertices;

traces = edf_trials_file.trials(start_event);
x = traces.samples('posX');
y = traces.samples('posY');
t = traces.time;

assert( t(1) == 0, 'Expected start time to be 0; was %d.', t(1) );

if ( params.use_reward_offset_time )
  time_offsets = amount_reward;
  
  assert( t(end) >= time_offsets, 'Not enough samples; expected %d; got %d.' ...
    , time_offsets, t(end) );
else
  start_ts = events_file.event_times(:, events_file.event_key(start_event));
  stop_ts = events_file.event_times(:, events_file.event_key(stop_event));
  
  time_offsets = floor( (stop_ts - start_ts)*1e3 );
end

lookdur = rownan( rows(x) );

is_left = strcmp( placements, 'center-left' );
is_right = strcmp( placements, 'center-right' );

assert( all(is_left | is_right), 'Some placements were not accounted for.' );

for i = 1:numel(placements)
  if ( params.use_reward_offset_time )
    subset_ind = t >= 0 & t <= time_offsets;
  else
    subset_ind = t >= 0 & t <= time_offsets(i);
  end
  
  x1 = x(i, subset_ind);
  y1 = y(i, subset_ind);
  
  if ( is_left(i) )
    use_rect = left_rect;
  else
    use_rect = right_rect;
  end
  
  ib = bfw.bounds.rect( x1, y1, use_rect );
  
  lookdur(i) = nnz( ib );
end

out = struct();
out.looking_duration = lookdur;
out.labels = labels;

end

function labs = update_labels(labs, unified_file)

hwwa.add_day_labels( labs );
hwwa.add_data_set_labels( labs );
hwwa.add_drug_labels( labs );
hwwa.fix_image_category_labels( labs );
hwwa.split_gender_expression( labs );

addcat( labs, 'monkey' );
setcat( labs, 'monkey', lower(unified_file.opts.META.monkey) );

end