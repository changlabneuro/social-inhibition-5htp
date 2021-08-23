function outs = hwwa_load_basic_behav_measures(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.time_trial_bin_size = 15*60; % 15 minutes
defaults.trial_bin_size = 50;
defaults.trial_step_size = 1;
defaults.num_run_time_quantiles = 4;

inputs = { 'trial_data', 'meta', 'labels', 'events', 'rois', 'events' };

[params, loop_runner] = hwwa.get_params_and_loop_runner( inputs, '', defaults, varargin );
loop_runner.convert_to_non_saving_with_output();

results = loop_runner.run( @main, params );
results(~[results.success]) = [];

outputs = [results.output];
outs = struct();

if ( isempty(outputs) )
  outs.rt = [];
  outs.labels = fcat();
  outs.run_relative_start_times = [];
  outs.rois = struct( [] );
  outs.events = struct( [] );
else
  match_categories( outputs );
  
  outs.rt = vertcat( outputs.rt );
  outs.labels = vertcat( fcat(), outputs.labels );
  outs.run_relative_start_times = vertcat( outputs.run_relative_start_times );
  outs.rois = vertcat( outputs.rois );
  outs.events = vertcat( outputs.events );
  outs.event_file_indices = make_event_file_indices( outputs );
end

end

function tot = make_event_file_indices(outputs)

tot = [];
for i = 1:numel(outputs)
  curr = repmat( i, size(outputs(i).rt, 1), 1 );
  tot = [ tot; curr ];
end

end

function out = main(files, params)

trial_data_file = shared_utils.general.get( files, 'trial_data' );
meta_file = shared_utils.general.get( files, 'meta' );
labels_file = shared_utils.general.get( files, 'labels' );
events_file = shared_utils.general.get( files, 'events' );
rois_file = shared_utils.general.get( files, 'rois' );
events_file = shared_utils.general.get( files, 'events' );

rt = reshape( [trial_data_file.trial_data(:).reaction_time], [], 1 );

trial_start_times = events_file.event_times(:, events_file.event_key('new_trial'));
time_trial_bin_size = params.time_trial_bin_size;
trial_bin_size = params.trial_bin_size;
trial_step_size = params.trial_step_size;

labs = labels_file.labels';

trial_bin_mask = get_trial_bin_mask( labs );

hwwa.add_day_labels( labs );
hwwa.add_data_set_labels( labs );
% hwwa.add_drug_labels( labs );
hwwa.add_drug_labels_by_day( labs );
hwwa.fix_image_category_labels( labs );
hwwa.add_time_bin_labels( labs, trial_start_times, time_trial_bin_size );
hwwa.add_trial_bin_labels( labs, trial_bin_size, trial_step_size, trial_bin_mask );
hwwa.split_gender_expression( labs );
add_run_time_quantile_labels( labs, trial_start_times, params.num_run_time_quantiles );

addcat( labs, 'monkey' );
setcat( labs, 'monkey', lower(meta_file.monkey) );

prune( labs );

out = struct();
out.rt = rt;
out.labels = labs;
out.run_relative_start_times = trial_start_times - min( trial_start_times );
out.rois = repmat( rois_file, numel(rt), 1 );
out.events = events_file;

end

function labels = add_run_time_quantile_labels(labels, start_times, num_quantiles)

boundary_points = linspace( min(start_times), max(start_times), num_quantiles+1 );

quant_cat = 'run_time_quantile';
addcat( labels, quant_cat );

for i = 1:num_quantiles
  lower_bound = boundary_points(i);
  upper_bound = boundary_points(i+1);
  
  if ( i < num_quantiles )
    in_bound = start_times >= lower_bound & start_times < upper_bound;
  else
    in_bound = start_times >= lower_bound & start_times <= upper_bound;
  end
  
  in_bound_ind = find( in_bound );
  setcat( labels, quant_cat, sprintf('%s__%d', quant_cat, i), in_bound_ind );
end

prune( labels );

assert( count(labels, makecollapsed(labels, quant_cat)) == 0 );

end

function match_categories(outputs)

cats = containers.Map();

for i = 1:numel(outputs)
  c = getcats( outputs(i).labels );
  
  for j = 1:numel(c)
    cats(c{j}) = 1;
  end
end

cats = keys( cats );

for i = 1:numel(outputs)
  addcat( outputs(i).labels, cats );
end

end

function mask = get_trial_bin_mask(labels)

mask = fcat.mask( labels ...
  , @find, 'initiated_true' ...
);

gonogo_error = find( labels, 'wrong_go_nogo' );
ok_trial = find( labels, 'correct_true' );
error_or_ok = union( gonogo_error, ok_trial );

mask = intersect( mask, error_or_ok );

end