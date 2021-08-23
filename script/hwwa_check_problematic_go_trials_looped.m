function outs = hwwa_check_problematic_go_trials_looped(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.kind = 'stats';
defaults.target_event = 'go_target_acquired';
defaults.base_subdir = '';

inputs = { 'unified', 'edf_trials', 'edf_events', 'labels' };

[params, runner] = hwwa.get_params_and_loop_runner( inputs, '', defaults, varargin );

kind = validatestring( params.kind, {'stats', 'plot_traces'}, mfilename, 'kind' );

runner.func_name = mfilename;
runner.main_error_handler = 'error';
runner.convert_to_non_saving_with_output();

if ( strcmp(kind, 'stats') )
  results = runner.run( @check_main, params );

  ok_results = results( [results.success] );

  mat = arrayfun( @(x) x.output.mat, ok_results, 'un', 0 );
  labs = arrayfun( @(x) x.output.labels, ok_results, 'un', 0 );
  traces = arrayfun( @(x) x.output.traces, ok_results, 'un', 0 );
  
  if ( ~isempty(ok_results) )
    data_key = ok_results(1).output.data_key;
  else
    data_key = containers.Map();
  end

  mat = vertcat( mat{:} );
  traces = vertcat( traces{:} );
  labs = vertcat( fcat(), labs{:} );

  assert_ispair( mat, labs );

  outs.data = mat;
  outs.data_key = data_key;
  outs.traces = traces;
  outs.labels = labs;
else
  outs = runner.run( @plot_incorrect_traces, params );
end

end

function outs = check_main(files, params)

labels_file = shared_utils.general.get( files, 'labels' );
events_file = shared_utils.general.get( files, 'edf_events' );
unified_file = shared_utils.general.get( files, 'unified' );
edf_trials_file = shared_utils.general.get( files, 'edf_trials' );

fixation_time = unified_file.opts.STIMULI.setup.go_target.target_duration;
go_nogo_time = unified_file.opts.TIMINGS.time_in.go_nogo;

edf_fix_time = fixation_time * 1e3;
edf_go_nogo_time = go_nogo_time * 1e3;

t0_event_name = params.target_event;

labs = labels_file.labels;

target_event = edf_trials_file.trials(t0_event_name);

event_times = events_file.event_times;
event_key = events_file.event_key;

targ_onset_times = event_times(:, event_key('go_target_onset'));
% targ_offset_times = event_times(:, event_key('go_target_offset'));

targ_offset_times = targ_onset_times + edf_go_nogo_time - 16; % minus one frame

t0_times = targ_onset_times;
t1_times = targ_offset_times;

t = target_event.time;
t1 = floor( t1_times - t0_times );
max_t = max( t );

t1(t1 > max_t) = max_t;

t0_ind = t == 0;
t1_ind = arrayfun( @(x) find(t == x), t1, 'un', 0 );

assert( nnz(t0_ind) == 1, 'No time 0 found for "%s".', t0_event_name );

t0_ind = find( t0_ind );

x = target_event.samples('posX');
y = target_event.samples('posY');

go_targ = unified_file.opts.STIMULI.go_target;

targ = go_targ.targets{1};

x_off = targ.x_offset;
y_off = targ.y_offset;
offsets = [ x_off, y_off, x_off, y_off ];
targ_bounds = targ.bounds + targ.padding + offsets;

center_x = mean( targ_bounds([1, 3]) );
center_y = mean( targ_bounds([2, 4]) );

p_ib = rownan( rows(x) );
did_meet_fix_crit = rownan( rows(x) );
p_nan = rownan( rows(x) );
c_dist = rownan( rows(x) );
x_dist = rownan( rows(x) );
devs = rownan( rows(x) );
blink_durs = rownan( rows(x) );
traces = rowcell( rows(x) );

for j = 1:rows(x)
  c_t1_ind = t1_ind{j};
  
  if ( isempty(c_t1_ind) )
    continue; 
  end
  
  select_cols = t0_ind:c_t1_ind;
  
  subset_x = x(j, select_cols);
  subset_y = y(j, select_cols);
  subset_t = t(select_cols);
  
  is_ib = bfw.bounds.rect( subset_x, subset_y, targ_bounds );
  
  [starts, durs] = shared_utils.logical.find_all_starts( is_ib );
  
  did_meet_fix_crit(j) = 0;
  
  for k = 1:numel(starts)
    start1 = starts(k);
    stop1 = start1 + durs(k) - 1;
    
    if ( subset_t(stop1) - subset_t(start1) >= edf_fix_time )
      did_meet_fix_crit(j) = 1;
      
%       if ( strcmp(labs(j, 'correct'), 'correct_false') && strcmp(labs(j, 'trial_type'), 'go_trial') )
%         d = 10;
%       end
      
      break;
    end
  end
  
  nan_x = isnan( subset_x );
  nan_y = isnan( subset_y );
  
  is_nan = nan_x | nan_y;
  
  [~, nan_durs] = shared_utils.logical.find_all_starts( is_nan );
  
  p_ib(j) = pnz( is_ib );
  p_nan(j) = pnz( is_nan );
  
  blink_durs(j) = nanmedian( nan_durs );
  
  c_dist(j) = nanmedian( bfw.distance(subset_x, subset_y, center_x, center_y)  );
  x_dist(j) = nanmedian( abs(subset_x - center_x) );
  devs(j) = max( nanstd(subset_x), nanstd(subset_y) );
  traces{j} = [ subset_x(:), subset_y(:) ];  
  
end

data_key = containers.Map();
data_key('in_bounds') = 1;
data_key('nan') = 2;
data_key('distance') = 3;
data_key('x_distance') = 4;
data_key('std_position') = 5;
data_key('met_fixation_criterion') = 6;
data_key('blink_duration') = 7;

outs = struct();
outs.mat = [ p_ib, p_nan, c_dist, x_dist, devs, did_meet_fix_crit, blink_durs ];
outs.data_key = data_key;
outs.labels = labs;
outs.traces = traces;

end

function outs = plot_incorrect_traces(files, params)

conf = params.config;
base_plot_p = fullfile( conf.PATHS.data_root, 'plots', 'debug', datestr(now, 'mmddyy') );
base_subdir = params.base_subdir;

outs = [];

t0_event_name = params.target_event;

labels_file = shared_utils.general.get( files, 'labels' );
events_file = shared_utils.general.get( files, 'edf_events' );
unified_file = shared_utils.general.get( files, 'unified' );
edf_trials_file = shared_utils.general.get( files, 'edf_trials' );

fixation_time = unified_file.opts.STIMULI.setup.go_target.target_duration;

labs = prune( hwwa.add_day_labels(labels_file.labels) );

is_go_incorrect = find( labs, {'go_trial', 'nogo_choice'} );
use_trials = is_go_incorrect;

target_event = edf_trials_file.trials(t0_event_name);

event_times = events_file.event_times;
event_key = events_file.event_key;

t0 = round( event_times(:, event_key(t0_event_name)) );

event_names = { 'go_target_onset', 'go_target_offset', 'go_target_acquired' };

all_times = nan( rows(labs), numel(event_names) );

for i = 1:numel(event_names)
  all_times(:, i) = round( event_times(:, event_key(event_names{i})) ) - t0;
end

t = target_event.time;
x = target_event.samples('posX');
y = target_event.samples('posY');

onset_time = all_times(:, strcmp(event_names, 'go_target_onset'));
offset_time = all_times(:, strcmp(event_names, 'go_target_offset'));

t0_ind = t == 0;
assert( nnz(t0_ind) == 1, 'No time 0 found for "%s" event.', t0_event_name );

go_targ = unified_file.opts.STIMULI.go_target;

targ = go_targ.targets{1};

x_off = targ.x_offset;
y_off = targ.y_offset;
offsets = [ x_off, y_off, x_off, y_off ];
bounds = targ.bounds + targ.padding + offsets;

is_problematic = is_problematic_go_trial( t, onset_time, offset_time ...
  , x, y, labs, bounds, round(fixation_time*1e3) );

f = figure(1);

title_cats = { 'unified_filename', 'trial_type', 'trial_outcome' };

subdir = char( fcat.strjoin(combs(labs, {'day', 'unified_filename'}), '_') );
subdir = strrep( subdir , '.mat', '' );
full_save_p = fullfile( base_plot_p, base_subdir, subdir );

shared_utils.io.require_dir( full_save_p );

for i = 1:numel(use_trials)
  shared_utils.general.progress( i, numel(use_trials) );
  
  clf( f );
  shared_utils.plot.fullscreen( f );
  
  trial_ind = use_trials(i);
  
  colors = { 'r', 'b' };
  y_axis_lab = { 'x', 'y' };
  
  title_str = combs( labs, title_cats, trial_ind );
  title_str{end+1} = sprintf( 'trial__%d', trial_ind );
  
  title_str = char( strrep(fcat.strjoin(title_str, ' | '), '_', ' ') );
  
  if ( is_problematic(trial_ind) )
    title_str = sprintf( '**%s', title_str );
  end
  
  for k = 1:2
    
    ax = subplot( 1, 2, k );
  
    trial_x = x(trial_ind, :);
    trial_y = y(trial_ind, :);
    
    is_x = strcmp( y_axis_lab{k}, 'x' );
    
    plot_vec = ternary( is_x, trial_x, trial_y );
    
    plot( ax, t, plot_vec, colors{k} );
    
    hold on;

    ax = gca;

    for j = 1:numel(event_names)
      ind = t == all_times(trial_ind, j);

      if ( nnz(ind) == 0 ), continue; end

      shared_utils.plot.add_vertical_lines( ax, t(ind), 'k--' );
    end
    
    bounds_plot_spec = sprintf( '%s--', colors{k} );
    plot_bounds = ternary( is_x, bounds([1, 3]), bounds([2, 4]) );
    shared_utils.plot.add_horizontal_lines( ax, plot_bounds, bounds_plot_spec );
    
    nans = find( isnan(plot_vec) );
    lims = get( ax, 'ylim' );
    
%     for j = 1:numel(nans)
%       plot( ax, t(nans), lims(2), sprintf('%s*', colors{k}) );
%     end
    
    ylabel( ax, y_axis_lab{k} );
    
    title( ax, title_str );
  end
  
  filename = strrep( char(fcat.trim(title_str)), '|', '_' );
  filename = strrep( filename, '**', 'problem__' );
  file_path = fullfile( full_save_p, filename );
  
  shared_utils.plot.save_fig( f, file_path, {'png'}, true );
end


end

function tf = is_problematic_go_trial(t, t0, t1, x, y, labs, bounds, fix_time)

max_t = max( t );
t1(t1 > max_t) = max_t;

is_go_incorrect = find( labs, {'go_trial', 'nogo_choice'} );

tf = rowzeros( rows(labs), 'logical' );

for i = 1:numel(is_go_incorrect)
  trial_ind = is_go_incorrect(i);
  
  start_ind = find( t == t0(trial_ind) );
  stop_ind = find( t == t1(trial_ind) );
  
  subset_x = x(trial_ind, start_ind:stop_ind);
  subset_y = y(trial_ind, start_ind:stop_ind);
  subset_t = t(start_ind:stop_ind);
  
  ib = bfw.bounds.rect( subset_x, subset_y, bounds );
  
  [starts, durs] = shared_utils.logical.find_all_starts( ib );
  
  for j = 1:numel(starts)
    start1 = starts(j);
    stop1 = start1 + durs(j) - 1;
    
    if ( subset_t(stop1) - subset_t(start1) >= fix_time )
      tf(trial_ind) = true;
      break;
    end
  end
end

end