function trials_file = edf_trials(files, varargin)

%   EDF_TRIALS -- Create edf_trials intermediate file.
%
%     file = ... edf_trials( files, 'event', event_names ) creates the 
%     edf_trials intermediate file from the events and edf files contained 
%     in `files.` `files` is a containers.Map or struct with an entry for 
%     'edf' and <events>, where <events> is a char vector giving the name 
%     of the intermediate subfolder that contains event times expressed in 
%     Eyelink's clock (default is 'el_events'). `event_names` is a 
%     cell array of strings or char vector giving the name(s) of events to 
%     which to align samples.
%
%     The output file contains a struct for each event in `event_names`;
%     for each trial of a given event, the (x, y) position and pupil size are
%     given as vectors of data +/- 500ms relative to the start of an event.
%     The result for each data type is an MxN matrix of M-trials by
%     N-samples. A time field of the struct gives the time of each 
%     aligned sample, relative to -500ms.
%
%     file = ... edf_trials( ..., 'look_back', lb, 'look_ahead', la )
%     changes how far back and ahead to look relative to an event's time 0.
%     For example, if `lb` is 0 and `la` is 1000, then each data vector
%     will contain 1 second's worth of samples, beginning at the event's
%     time 0. If `lb` is scalar, then it is used for each event in
%     `event_names`. Otherwise, the number of `lb` elements must match the
%     number of `event_names`, and each `lb(i)` gives the number of ms to
%     look back for the corresponding i-th event. The same applies for
%     `la`.
%
%     file = ... edf_trials( ..., 'append', tf ) indicates whether to
%     attempt to re-use an existing edf_trials file and append data to it.
%     Default is false.
%
%     EXAMPLE //
%
%     edf_file = hwwa.load1( 'edf' );
%     events_file = hwwa.load1( 'el_events', edf_file.unified_filename );
%     % Create the file aggregate.
%     files = struct( 'edf', edf_file, 'el_events', events_file );
%     event_name = 'go_target_onset';
%     look_back = 0;
%     look_ahead = 1000;
%     edf_trials_file = hwwa.make.edf_trials( files ...
%       , 'look_back', look_back, 'look_ahead', look_ahead );
%
%     See also hwwa.make.defaults.edf_trials
%
%     IN:
%       - `files` (containers.Map, struct)
%     FILES:
%       - 'edf'
%       - <events>
%     OUT:
%       - `events_file` (struct)

defaults = hwwa.make.defaults.edf_trials();
defaults.append = false;  % Don't attempt to load existing file.

params = hwwa.parsestruct( defaults, varargin );

assert( ~isempty(params.event), 'Specify an event name.' );

event_subdir = params.event_subdir;

hwwa.validatefiles( files, {'edf', event_subdir} );

edf_file = shared_utils.general.get( files, 'edf' );
events_file = shared_utils.general.get( files, event_subdir );

unified_filename = hwwa.try_get_unified_filename( edf_file );

event_names = cellstr( params.event );
look_back = params.look_back;
look_ahead = params.look_ahead;

validate_event_parameters( event_names, look_back, look_ahead );

if ( isscalar(look_back) )
  look_back = repmat( look_back, 1, numel(event_names) );
end

if ( isscalar(look_ahead) )
  look_ahead = repmat( look_ahead, 1, numel(event_names) );
end

should_create_file = true;

if ( params.append )
  % Check if a trials_file already exists. If so, load it and append the
  % contents to it. Otherwise, create it.
  trials_fullfile = get_trials_fullfile( unified_filename, params );
  
  if ( shared_utils.io.fexists(trials_fullfile) )
    trials_file = shared_utils.io.fload( trials_fullfile );
    should_create_file = false;
  end
end

if ( should_create_file )
  trials_file = struct();
  trials_file.trials = containers.Map();
  trials_file.unified_filename = unified_filename;
end

for i = 1:numel(event_names)
  event_name = event_names{i};
  insert_one_event( trials_file.trials, edf_file, events_file ...
    , event_name, look_back(i), look_ahead(i) );
end

trials_file.params = params;

end

function insert_one_event(trials, edf, events, event_name, look_back, look_ahead)

sample_types = { 'pupilSize', 'posX', 'posY' };

event_times = events.event_times(:, events.event_key(event_name));

all_samples = containers.Map();
time = edf.Samples.time;

for i = 1:numel(sample_types)
  samples = edf.Samples.(sample_types{i});

  [mat, t, fs] = get_aligned_samples( time, samples, event_times, look_back, look_ahead );

  all_samples(sample_types{i}) = mat;
end

trials(event_name) = struct( 'samples', all_samples, 'time', t, 'sample_rate', fs );

end

function validate_event_parameters(event_names, look_back, look_ahead)

template = 'Number of "%s" must match the number of events, or else be scalar.';

n_lb = numel( look_back );
n_la = numel( look_ahead );
n_evts = numel( event_names );

if ( n_lb ~= 1 )
  assert( n_evts == n_lb, template, 'look back' );
end

if ( n_la ~= 1 )
  assert( n_evts == n_la, template, 'look ahead' );
end

end

function fs = get_sample_frequency(time)

fs = unique( diff(time(~isnan(time))) );
assert( numel(fs) == 1 ...
  , 'Expected one unique inter-sample interval, but got %d', numel(fs) );

end

function [mat, t, fs] = get_aligned_samples(time, samples, events, look_back, look_ahead)

fs = get_sample_frequency( time );

start_ts = round( events + (look_back/fs) );
amount = round( (look_ahead - look_back) / fs );

mat = nan( numel(start_ts), amount + 1 );

start_inds = shared_utils.sync.nearest( time, start_ts );
is_valid_start = ~isnan( start_ts );

offsets = time(start_inds(is_valid_start)) - start_ts(is_valid_start);

assert( max(abs(offsets)) <= 1, ['Start time was more than 1 sample away from' ...
  , ' a time sample in the raw .edf data.'] );

others = find( is_valid_start );

t = look_back:fs:look_back+round( amount*fs );

for i = 1:numel(others)
  other = others(i);
  start_ind = start_inds(other);
  stop_ind = start_ind + amount;
  
  assign_start = 1;
  assign_stop = size( mat, 2 );
  
  if ( stop_ind > numel(samples) )
    assign_stop = assign_stop - (stop_ind - numel(samples) );
    stop_ind = numel( samples );
  end
  
  mat(other, assign_start:assign_stop) = samples(start_ind:stop_ind);
end

end

function trials_fullfile = get_trials_fullfile(unified_filename, params)

output_directory = params.output_directory;
conf = params.config;

trials_fullfile = fullfile( hwwa.get_intermediate_dir(output_directory, conf) ...
  , unified_filename );

end