function make_first_saccade_info(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.t1 = 20;
defaults.t2 = 10;
defaults.min_duration = 0.03;
defaults.debug = false;
defaults.pad_beginning = nan;
defaults.events = 'all';

params = shared_utils.general.parsestruct( defaults, varargin );

conf = params.config;

edf_trial_p = hwwa.get_intermediate_dir( 'edf_trials', conf );
output_p = hwwa.get_intermediate_dir( 'first_saccade', conf );

mats = hwwa.require_intermediate_mats( params.files, edf_trial_p, params.files_containing );

parfor i = 1:numel(mats)
  shared_utils.general.progress( i, numel(mats) );
  
  edf_trials_file = shared_utils.io.fload( mats{i} );
  
  unified_filename = edf_trials_file.unified_filename;
  
  try
    saccade_main( edf_trials_file, output_p, params );
  catch err
    hwwa.print_fail_warn( unified_filename, output_p, err.message );
    continue;
  end
end


end

function saccade_main(edf_trials_file, output_p, params)

events = keys( edf_trials_file.trials );

if ( ~strcmp(params.events, 'all') )
  events = csintersect( events, params.events );
end

for i = 1:numel(events)
  try
    full_output_p = fullfile( output_p, events{i} );
    
    unified_filename = edf_trials_file.unified_filename;
    output_filename = fullfile( full_output_p, unified_filename );
  
    if ( hwwa.conditional_skip_file(output_filename, params.overwrite) )
      continue;
    end
    
    subset = edf_trials_file.trials(events{i});
    
    saccade_file = per_event_type( subset, params );
    saccade_file.unified_filename = unified_filename;
    saccade_file.params = params;
    
    if ( params.save )
      shared_utils.io.require_dir( full_output_p );
      shared_utils.io.psave( output_filename, saccade_file, 'saccade_file' );
    end
  catch err
    hwwa.print_fail_warn( edf_trials_file.unified_filename, err.message );
    continue;
  end
end

end

function outs = per_event_type(subset, params)

x = subset.samples('posX');
y = subset.samples('posY');
t = subset.time;

N = size( x, 1 );

dirs = rownan( N );
x_diffs = rownan( N );
start_times = rownan( N );
start_indices = rownan( N );
destinations = nan( N, 2 );

for i = 1:N
  tn_x = x(i, :);
  tn_y = y(i, :);
  
  is_fix = eye_mmv_is_fixation( tn_x, tn_y, t/1e3, params );
  
  if ( ~isnan(params.pad_beginning) )
    is_fix(1:params.pad_beginning) = true;
  end
  
  if ( params.debug )
    debug_check_plot( t, tn_x, tn_y, is_fix );
  end
  
  [starts, durs] = shared_utils.logical.find_starts( ~is_fix, 1, 1 );
  
  if ( numel(starts) < 1 )
    continue;
  end
  
  stops = starts + durs - 1;
  
  start_sacc = starts(1);
  stop_sacc = stops(1);
  
  diff_x = tn_x(stop_sacc) - tn_x(start_sacc);
  diff_y = tn_y(stop_sacc) - tn_y(start_sacc);
  
  if ( isnan(diff_x) )
    continue;
  end
  
  x_diffs(i) = diff_x;
  dirs(i) = ternary( sign(diff_x) > 0, 1, 0 );
  
  start_times(i) = t(start_sacc);
  start_indices(i) = start_sacc;
  destinations(i, :) = [tn_x(stop_sacc), tn_y(stop_sacc)];
end

outs.directions = dirs;
outs.start_times = start_times;
outs.start_indices = start_indices;

end

function debug_check_plot(t, tn_x, tn_y, is_fix)

figure(1);
clf();
plot( t, tn_x ); 
hold on;
plot( t, tn_y, 'r' );

legend( {'x', 'y'} );

fix_samples = find( is_fix );
lims = get( gca, 'ylim' );

scatter( t(fix_samples), repmat(lims(2), size(fix_samples)) );

first_non_fix = find( ~is_fix, 1 );
first_non_fix_t = t(first_non_fix);

plot( [first_non_fix_t; first_non_fix_t], get(gca, 'ylim') );

end

function is_fix = eye_mmv_is_fixation(x, y, time, params)

pos = [ x(:)'; y(:)' ];

t1 = params.t1;
t2 = params.t2;
min_duration = params.min_duration;

%   repositories/eyelink/eye_mmv
is_fix = is_fixation( pos, time(:)', t1, t2, min_duration );
is_fix = logical( is_fix );
is_fix = is_fix(1:numel(time))';

end
