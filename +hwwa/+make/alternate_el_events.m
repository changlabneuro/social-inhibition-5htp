function events_file = alternate_el_events(files)

%   ALTERNATE_EL_EVENTS -- Create eyelink-events file.
%
%     file = ... alternate_el_events( files ) creates the eyelink-events 
%     intermediate file from the events and edf intermediate files contained 
%     in `files`. `files` is a containers.Map or struct with an entry for 
%     'events' and 'edf.
%
%     The output file is of the same format as the 'events' intermediate
%     file, but event times are expressed in terms of Eyelink's clock.
%
%     This function, compared to hwwa.make.el_events, uses a more robust 
%     interpolation function to convert event times to Eyelink's clock,
%     and should be preferred.
%
%     IN:
%       - `files` (containers.Map, struct)
%     FILES:
%       - 'events'
%       - 'edf'
%     OUT:
%       - `events_file` (struct)

hwwa.validatefiles( files, {'events', 'edf'} );

events_file = shared_utils.general.get( files, 'events' );
edf_file = shared_utils.general.get( files, 'edf' );

messages = edf_file.Events.Messages;

edf_sync_inds = cellfun( @(x) ~isempty(strfind(x, 'TRIAL__')), messages.info );
edf_sync_times = messages.time(edf_sync_inds);

mat_event_times = events_file.event_times;

trial_start_col = events_file.event_key( 'new_trial' );
mat_sync_times = mat_event_times(:, trial_start_col);

nmat = numel( mat_sync_times );
nedf = numel( edf_sync_times );

assert( nmat == nedf || nedf == nmat + 1, 'Eyelink and matlab sync times do not correspond.' );

edf_sync_times = edf_sync_times(1:nmat);

edf_event_times = nan( size(mat_event_times) );

for i = 1:size(mat_event_times, 2)
  mat_times = mat_event_times(:, i);  
  
  edf_times = shared_utils.sync.cinterp( mat_times, mat_sync_times, edf_sync_times, true );
%   edf_times = round( edf_times );
  
  edf_event_times(:, i) = edf_times;
end

events_file.event_times = edf_event_times;

end