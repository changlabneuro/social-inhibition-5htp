function hwwa_mark_is_saccade(varargin)

defaults = hwwa.get_common_make_defaults();
defaults.t1 = 20;
defaults.t2 = 10;
defaults.min_duration = 0.03;
defaults.debug = false;
defaults.pad_beginning = nan;

params = shared_utils.general.parsestruct( defaults, varargin );

conf = params.config;

edf_trial_p = hwwa.get_intermediate_dir( 'edf_trials', conf );
mats = hwwa.require_intermediate_mats( params.files, edf_trial_p, params.files_containing );

for i = 1:numel(mats)
  shared_utils.general.progress( i, numel(mats) );
  
  edf_trials_file = shared_utils.io.fload( mats{i} );
  
  try
    outs = saccade_main( edf_trials_file, params );
  catch err
    hwwa.print_fail_warn( edf_trials_file.unified_filename, err.message );
    continue;
  end
  
  d = 10;
end


end

function outs = saccade_main(edf_trials_file, params)

post_cue_display = edf_trials_file.trials('go_nogo_cue_onset');

x = post_cue_display.samples('posX');
y = post_cue_display.samples('posY');
t = post_cue_display.time;

N = size( x, 1 );

dirs = rownan( N );
x_diffs = rownan( N );
start_times = rownan( N );

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
end

outs.directions = dirs;
outs.start_times = start_times;

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
