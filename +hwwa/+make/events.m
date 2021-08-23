function events_file = events(files)

%   EVENTS -- Create events file.
%
%     file = ... events( files ) creates the events intermediate file from 
%     the unified intermediate file contained in `files`. `files` is a
%     containers.Map or struct with an entry for 'unified'.
%
%     The output file contains a matrix of event times for each
%     trial, and a key identifying columns of the matrix. The matrix is of
%     size M-trials by N-event types.
%
%     IN:
%       - `files` (containers.Map, struct)
%     FILES:
%       - 'unified'
%     OUT:
%       - `events_file` (struct)

hwwa.validatefiles( files, 'unified' );

unified_file = shared_utils.general.get( files, 'unified' );
unified_filename = hwwa.try_get_unified_filename( unified_file );

evts = reconcile_events( unified_file );
evt_names = sort( fieldnames(evts) );

event_times = zeros( numel(evts), numel(evt_names) );
event_key = containers.Map( 'keytype', 'char', 'valuetype', 'double' );

for j = 1:numel(evt_names)
  evt_name = evt_names{j};
  event_times(:, j) = [ evts(:).(evt_name) ];
  event_key(evt_name) = j;
end

has_iti = any( strcmp(event_key, 'iti') );

%   older sessions don't have an iti event time, but the preceding state
%   is reward onset -- so just add the reward state time
if ( has_iti )
  reward_dur = unified_file.opts.TIMINGS.time_in.reward;
  error_dur = unified_file.opts.TIMINGS.time_in.error_go_nogo;

  reward_ts = event_times(:, event_key('reward_onset'));
  error_ts = event_times(:, event_key('go_target_offset'));

  use_error_ts = strcmp( {unified_file.DATA(:).error}, 'wrong_go_nogo' );

  iti_ts = reward_ts + reward_dur;
  iti_ts(use_error_ts(:)) = error_ts(use_error_ts(:)) + error_dur;

  current_n = numel( evt_names );
  event_key('iti') = current_n + 1;
  event_times(:, end+1) = iti_ts;
end

events_file = struct();
events_file.unified_filename = unified_filename;
events_file.event_key = event_key;
events_file.event_times = event_times;

end

function evts = reconcile_events(unified)

fnames = arrayfun( @(x) fieldnames(x.events), unified.DATA, 'un', 0 );
ns = cellfun( @numel, fnames );

if ( numel(unique(ns)) == 1 )
  evts = [ unified.DATA(:).events ];
  return
end

[max_n, I] = max( ns );
fields = fnames{I};

lt_max = find( ns < max_n );

base_evts = { unified.DATA(:).events };

for i = 1:numel(lt_max)
  ind = lt_max(i);
  
  missing = setdiff( fields, fnames{ind} );
  
  for j = 1:numel(missing)
    miss = missing{j};
    
    base_evts{ind}.(miss) = nan;
  end
end

evts = horzcat( base_evts{:} );

end