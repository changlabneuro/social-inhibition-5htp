function make_aligned_lfp(varargin)

defaults = hwwa.get_common_make_defaults();

defaults.look_back = -0.5;  % s
defaults.look_ahead = 0.5;
defaults.window_size = 0;
defaults.event = '';
defaults.kind = 'fp';

params = hwwa.parsestruct( defaults, varargin );

kind = params.kind;

switch ( kind )
  case 'fp'
    lfp_dir = 'lfp';
    aligned_dir = 'aligned_lfp';
  case 'wb'
    lfp_dir = 'wb';
    aligned_dir = 'aligned_wb';
  otherwise
    error( 'Unrecognized kind "%s".', kind );
end

event_p = hwwa.get_intermediate_dir( 'plex_events' );
lfp_p = hwwa.get_intermediate_dir( lfp_dir );
output_p = hwwa.get_intermediate_dir( aligned_dir );

if ( isempty(params.event) )
  error( 'Specify an event name.' );
end

event_names = params.event;
look_back = params.look_back;
look_ahead = params.look_ahead;
win_size = params.window_size;

if ( ~iscell(event_names) ), event_names = { event_names }; end

assert( numel(event_names) == numel(look_back) && numel(look_back) == numel(look_ahead) ...
  , 'Number of event names must match number of look back and look ahead.' );

mats = hwwa.require_intermediate_mats( params.files, lfp_p, params.files_containing );

for i = 1:numel(mats)
  hwwa.progress( i, numel(mats), mfilename );
  
  lfp = shared_utils.io.fload( mats{i} );
  unified_filename = lfp.unified_filename;
  
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
    lfp_psth = shared_utils.io.fload( output_filename );
  else
    lfp_psth = struct();
    lfp_psth.psth = containers.Map();
    lfp_psth.unified_filename = unified_filename;
  end
  
  id_times = lfp.time;
  values = lfp.lfp;
  sr = lfp.sample_rate;
  
  for j = 1:numel(event_names)
    evt_name = event_names{j};
    evt_times = events.event_times(:, events.event_key(evt_name));
    
    id_time_indices = get_event_indices( evt_times, id_times );
    
    stp = 1;
    
    n_trials = numel( evt_times );
    n_chans = size( values, 1 );
    
    for k = 1:n_chans
      one_chan = values(k, :);
      
      [aligned, psth_t] = align_lfp( one_chan, id_time_indices ...
        , sr, look_back(j), look_ahead(j), win_size );
      
      labs = fcat.create( 'channel', lfp.channel{k}, 'id', lfp.id{k} );
      
      repmat( labs, n_trials );
      
      if ( k == 1 )
        all_psth = nan( n_trials * n_chans, numel(psth_t) );
        all_labs = labs;
      else
        append( all_labs, labs );
      end
      
      all_psth(stp:stp+n_trials-1, :) = aligned;
      
      stp = stp + n_trials;
    end
    
    lfp_psth.psth(evt_name) = struct( 'data', all_psth ...
      , 'labels', all_labs, 'time', psth_t, 'sample_rate', sr );
  end
  
  lfp_psth.params = params;
  
  shared_utils.io.require_dir( output_p );
  
  hwwa.psave( output_filename, lfp_psth, 'lfp_psth' );
end

end

function [out, t] = align_lfp(values, event_indices, sr, look_back, look_ahead, win_size)

back_offset = look_back * sr;
ahead_offset = look_ahead * sr;

%   if window size is non-zero, align to center of window

window_samples = win_size * sr;

starts = event_indices + back_offset - round(window_samples / 2);
stops = event_indices + ahead_offset - round(window_samples / 2) + window_samples;

fs = 1 / sr;

t = look_back:fs:look_ahead + win_size - fs;

n_samples = ahead_offset - back_offset + window_samples;
n_trials = numel( event_indices );

out = nan( n_trials, n_samples );

non_nan = find( ~isnan(event_indices) );
n_non_nan = numel( non_nan );

for i = 1:n_non_nan
  ind = non_nan(i);
  
  start = starts(ind);
  stop = stops(ind) - 1;
  
  out(ind, :) = values(start:stop);
end

end

function id_time_indices = get_event_indices(event_times, id_times)

id_time_indices = nan( numel(event_times), 1 );

for j = 1:numel(event_times)
  current_time = event_times(j);
  
  if ( isnan(current_time) ), continue; end
  
  [~, index] = histc( current_time, id_times );
  out_of_bounds_msg = ['The id_times do not properly correspond to the' ...
    , ' inputted events.'];
  is_in_bounds = index ~= 0;
  assert( is_in_bounds, out_of_bounds_msg );
  check = abs( current_time - id_times(index) ) < abs( current_time - id_times(index+1) );
  if ( ~check ), index = index + 1; end
  id_time_indices(j) = index;
end

end