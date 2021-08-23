function outs = hwwa_load_events(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.event_names = { 'fixation_on', 'go_nogo_cue_onset', 'iti' };

inputs = { 'events', 'labels', 'meta' };

[params, runner] = hwwa.get_params_and_loop_runner( inputs, '', defaults, varargin );
runner.convert_to_non_saving_with_output();

results = runner.run( @main, params );

outputs = [ results([results.success]).output ];

if ( isempty(outputs) )
  outs = struct();
  outs.labels = fcat();
  outs.events = [];
else
  outs = shared_utils.struct.soa( outputs ); 
end

outs.event_names = params.event_names;
outs.params = params;

end

function out = main(files, params)

events_file = shared_utils.general.get( files, 'events' );
labels_file = shared_utils.general.get( files, 'labels' );
meta_file = shared_utils.general.get( files, 'meta' );

event_names = params.event_names;
events = cellfun( @(x) events_file.event_times(:, events_file.event_key(x)) ...
  , event_names, 'un', 0 );
events = horzcat( events{:} );

out.labels = make_labels( labels_file, meta_file );
out.events = events;

end

function labs = make_labels(labels_file, meta_file)

labs = labels_file.labels';
hwwa.tidy_up_labels( labs, meta_file.monkey );
prune( labs );

end