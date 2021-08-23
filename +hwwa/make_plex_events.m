function make_plex_events(varargin)

defaults = hwwa.get_common_make_defaults();

params = hwwa.parsestruct( defaults, varargin );

conf = params.config;

event_p = hwwa.get_intermediate_dir( 'events', conf );
sync_p = hwwa.get_intermediate_dir( 'sync', conf );
output_p = hwwa.get_intermediate_dir( 'plex_events', conf );

mats = hwwa.require_intermediate_mats( params.files, sync_p, params.files_containing );

for i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  sync_file = shared_utils.io.fload( mats{i} );
  unified_filename = sync_file.unified_filename;
  
  events_file = fullfile( event_p, unified_filename );
  
  if ( ~shared_utils.io.fexists(events_file) )
    fprintf( '\n Skipping "%s" because the events file does not exist.', unified_filename );
    continue;
  end
  
  events = shared_utils.io.fload( events_file );
  
  output_filename = fullfile( output_p, sync_file.unified_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  mat_sync_times = sync_file.mat(:);
  pl2_sync_times = sync_file.plex(:);
  
  event_times = events.event_times;
  pl2_event_times = nan( size(event_times) );
  
  for j = 1:size(event_times, 2)
    mat_times = event_times(:, j);
    
    pl2_event_times(:, j) = hwwa.clock_a_to_b( mat_times, mat_sync_times, pl2_sync_times );
  end
  
  events.event_times = pl2_event_times;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, events, 'events' );
end

end