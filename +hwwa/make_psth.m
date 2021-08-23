function make_psth(varargin)

defaults = hwwa.get_common_make_defaults();

defaults.look_back = -0.5;  % s
defaults.look_ahead = 0.5;
defaults.bin_size = 0.01;
defaults.event = '';
defaults.kind = 'sua';

params = hwwa.parsestruct( defaults, varargin );

kind = params.kind;

switch ( kind )
  case 'sua'
    spike_dir = 'spikes';
  case 'mua'
    spike_dir = 'mua_spikes';
  otherwise
    error( 'Unrecognized kind "%s".', kind );
end

event_p = hwwa.get_intermediate_dir( 'plex_events' );
spike_p = hwwa.get_intermediate_dir( spike_dir );
output_p = hwwa.get_intermediate_dir( 'psth' );

if ( isempty(params.event) )
  error( 'Specify an event name.' );
end

event_names = params.event;
look_back = params.look_back;
look_ahead = params.look_ahead;
bin_size = params.bin_size;

if ( ~iscell(event_names) ), event_names = { event_names }; end

assert( numel(event_names) == numel(look_back) && numel(look_back) == numel(look_ahead) ...
  , 'Number of event names must match number of look back and look ahead.' );

mats = hwwa.require_intermediate_mats( params.files, spike_p, params.files_containing );

parfor i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  spikes = shared_utils.io.fload( mats{i} );
  unified_filename = spikes.unified_filename;
  
  event_file = fullfile( event_p, unified_filename );
  
  if ( ~shared_utils.io.fexists(event_file) )
    fprintf( '\n Skipping "%s" because it is missing an event file.', unified_filename );
    continue;
  end
  
  events = shared_utils.io.fload( fullfile(event_p, unified_filename) );
  
  output_filename = fullfile( output_p, unified_filename );
  
  if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
    continue;
  end
  
  %   reuse existing samples
  if ( shared_utils.io.fexists(output_filename) && params.append )
    psth = shared_utils.io.fload( output_filename );
  else
    psth = struct();
    psth.psth = containers.Map();
    psth.unified_filename = unified_filename;
  end
  
  units = spikes.units;
  
  for j = 1:numel(event_names)
    fprintf( '\n\t %d of %d', j, numel(event_names) );
    
    evt_name = event_names{j};
    evt_times = events.event_times(:, events.event_key(evt_name));
    
    for k = 1:numel(units)
      fprintf( '\n\t\t %d of %d', k, numel(units) );
      
      unit = units(k);
      
      [one_psth, psth_t] = one_unit( unit, evt_times, look_back(j), look_ahead(j), bin_size );
      
      unit_labs = hwwa.get_unit_labels( units(k) );
      
      repmat( unit_labs, numel(evt_times) );
      
      if ( k == 1 )
        all_labs = unit_labs;
        all_psth = one_psth;
      else
        append( all_labs, unit_labs );
        all_psth = [ all_psth; one_psth ];
      end
    end
    
    psth.psth(evt_name) = struct( 'data', all_psth, 'labels', all_labs, 'time', psth_t );
  end
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, psth, 'psth' );
end

end

function [out_psth, t] = one_unit(unit, events, min_t, max_t, bin_size)

for i = 1:numel(events)
  [psth, t] = looplessPSTH( unit.times, events(i), min_t, max_t, bin_size );
  
  if ( i == 1 )
    out_psth = nan( numel(events), numel(t) );
  end
  
  out_psth(i, :) = psth;
end

end