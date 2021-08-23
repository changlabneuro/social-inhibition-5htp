function outs = hwwa_load_edf_aligned(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.start_event_name = 'go_target_onset';
defaults.stop_event_name = '';
defaults.look_ahead = 500;
defaults.look_back = 0;

inputs = { 'edf_samples', 'edf_events', 'labels', 'meta', 'rois' };

[params, runner] = hwwa.get_params_and_loop_runner( inputs, '', defaults, varargin );
runner.convert_to_non_saving_with_output();

results = runner.run( @main, params );

outputs = [ results([results.success]).output ];

if ( isempty(outputs) )
  outs = struct();
  outs.labels = fcat();
  outs.x = [];
  outs.y = [];
  outs.pupil = [];
  outs.rois = struct();
  outs.roi_indices = [];
  outs.sample_relative_events = [];
  outs.event_key = {};
  outs.time = [];
else
  if ( ~isempty(params.stop_event_name) )
    outputs = reconcile_sample_matrix_sizes( outputs );
  end
  
  outs = shared_utils.struct.soa( outputs ); 
  outs.roi_indices = make_roi_indices( outputs );
end

outs.params = params;

end

function outputs = reconcile_sample_matrix_sizes(outputs)

cols = arrayfun( @(x) size(x.x, 2), outputs );
max_col = max( cols );

match_fields = { 'x', 'y', 'pupil' };

for i = 1:numel(outputs)
  n = cols(i);
  
  if ( n < max_col )
    for j = 1:numel(match_fields)
      src = outputs(i).(match_fields{j});
      new_size = [ size(src, 1), max_col ];
      dest = nan( new_size );
      dest(:, 1:n) = src;
      outputs(i).(match_fields{j}) = dest;
    end
  end
end

end

function roi_inds = make_roi_indices(outputs)

roi_inds = cell( size(outputs) );
  
for i = 1:numel(outputs)
  roi_inds{i} = repmat( i, rows(outputs(i).x), 1 );
end

roi_inds = vertcat( roi_inds{:} );

end

function check_event(key, kind, name)

if ( ~isKey(key, name) )
  error( 'No matching %s event "%s".', kind, name );
end

end

function [start_ts, stop_ts, one_start] = get_start_stop_times(events_file, params)

start_event = params.start_event_name;
check_event( events_file.event_key, 'start', start_event );

start_ts = events_file.event_times(:, events_file.event_key(start_event));

stop_event = params.stop_event_name;

if ( ~isempty(stop_event) )
  check_event( events_file.event_key, 'stop', stop_event );
  stop_ts = events_file.event_times(:, events_file.event_key(stop_event));
  one_start = false;
else
  stop_ts = start_ts;
  one_start = true;
end

start_ts = start_ts + params.look_back;
stop_ts = stop_ts + params.look_ahead;

end

function out = main(files, params)

samples_file = shared_utils.general.get( files, 'edf_samples' );
events_file = shared_utils.general.get( files, 'edf_events' );
labels_file = shared_utils.general.get( files, 'labels' );
meta_file = shared_utils.general.get( files, 'meta' );
roi_file = shared_utils.general.get( files, 'rois' );

[starts, stops, one_start] = get_start_stop_times( events_file, params );
% assert( one_start, 'Stop events-by-name not yet implemented' );

t = double( samples_file.t );

start_inds = double( bfw.find_nearest(t, starts) );
stop_inds = double( bfw.find_nearest(t, stops) );

amts = stop_inds - start_inds + 1;
max_amt = max( amts );

if ( one_start )
  valid_amts = find( amts == max_amt );
else
  valid_amts = find( amts ~= 1 );
end

num_valid = numel( valid_amts );
use_fields = { 'x', 'y', 'pupil' };

out = struct();

for i = 1:numel(use_fields)
  samples = double( samples_file.(use_fields{i}) );
  
  aligned_samples = nan( numel(start_inds), max_amt );
  
  for j = 1:num_valid
    valid_ind = valid_amts(j);
    start = start_inds(valid_ind);
    stop = stop_inds(valid_ind);
    n_assign = stop - start + 1;
    
    aligned_samples(valid_ind, 1:n_assign) = samples(start:stop);
  end
  
  out.(use_fields{i}) = aligned_samples;
end

sample_relative_events = sample_relative_event_indices( events_file.event_times, starts );

out.labels = make_labels( labels_file, meta_file );
out.rois = roi_file;
out.sample_relative_events = sample_relative_events;
out.event_key = { events_file.event_key };
out.time = params.look_back:params.look_ahead;

end

function events = sample_relative_event_indices(events, starts)

for i = 1:size(events, 2)
  events(:, i) = round( events(:, i) - starts ) + 1;
end

end

function labs = make_labels(labels_file, meta_file)

labs = labels_file.labels';

hwwa.add_day_labels( labs );
hwwa.add_data_set_labels( labs );
hwwa.add_drug_labels_by_day( labs );
hwwa.fix_image_category_labels( labs );
hwwa.split_gender_expression( labs );
hwwa.decompose_social_scrambled( labs );

addcat( labs, 'monkey' );

if ( ~isempty(labs) )
  setcat( labs, 'monkey', lower(meta_file.monkey) );
end

prune( labs );

end